import time

from sage.all import *
import sage.matrix.matrix_integer_dense_hnf as stein_hnf
from sage.rings.integer import GCD_list
from src.utils import prod_all_elems
from src.determ_rng import *
from src.utils import *
import src.timers as timers
from src.stage_linalg_internal import *
from src.stage_filter import assert_rows_multiply_to_correct_ideal

class LatticeManager:
  '''
  Encapsulates logic for:
   -storing elements of the class group lattice and unit group lattice
   -getting properties of them like the class number
   -adding new elements to either lattice
  The caller's job is to feed new relations into the lattice manager and to have all the index-calculus logic .
  '''
  def __init__(self, hnf_relns_mat, relns_mat, rat_relns_mat, cg_elems_mat, ug_elems_mats, nf, mats_for_vollmer, unit_rank, max_dets, num_chars,
               smoothness_bound, true_reg_for_testing, reg_cols_to_drop):
    '''
    :param relns_mat:
    :param cg_elems_mat:
    :param ug_elems_mats: a list of elements for computing the unit group
    :param nf:
    :param mats_for_vollmer: a list of matrices of relations
    :param unit_rank:
    :param max_dets:
    :param num_chars:
    :param smoothness_bound:
    :param true_reg_for_testing:
    '''
    self.num_ug_mats = len(ug_elems_mats)
    if len(reg_cols_to_drop) != len(ug_elems_mats):
      raise ValueError("number of cols to drop must be")
    if self.num_ug_mats != len(mats_for_vollmer):
      raise ValueError("number of matrices in mats_for_vollmer must equal the number in ug_elems_mats")
    if true_reg_for_testing != None:
      print("Non-production warning: lattice manager was given the true regulator for testing purposes")

    self.true_reg_for_testing = true_reg_for_testing
    self.cg_into_hnf = hnf_relns_mat != None
    if asserts_on():
      assert stein_hnf.is_in_hnf_form(hnf_relns_mat, range(relns_mat.ncols()))
    if self.cg_into_hnf:
      print("Lattice Manager: relns_mat is in HNF so its assumed you want to have relns_mat hold a basis of the class group all the way through")

    if not relns_mat.is_square():
      raise ValueError("relns mat should be square")
    if ug_elems_mats[0].nrows() != relns_mat.nrows():
      raise ValueError("number of rows in the elems mat should be the same as number of rows in relations mat")

    # UG
    self.ug_compact_rep_mats = [Matrix(ZZ, 0, mats_for_vollmer[0].nrows())] * self.num_ug_mats
    self.ug_combined_compact_rep_mat = None
    self.ug_elems_mats = ug_elems_mats
    self.ug_combined_elems_mat = None
    self.mats_for_vollmer = mats_for_vollmer
    self.character_mats_dict = [{}]*self.num_ug_mats
    self.inert_prs_for_chars = [{}]*self.num_ug_mats
    self.regulator_lower_bound=0.2
    self.num_initial_units = None
    self.ug_reg_cols_to_drop = reg_cols_to_drop

    # CG
    self.elems_of_V_as_Zcombination_of_U_elems = [] # where U and V are the two mats for Vollmer
    self.cg_rat_relns_mat = rat_relns_mat
    self.hnf_relns_mat = hnf_relns_mat
    self.cg_relns_mat = relns_mat
    self.cg_compact_rep_mat = identity_matrix(ZZ, self.cg_relns_mat.nrows())
    self.cg_elems_mat = cg_elems_mat
    self.cg_character_mats_dict = {}
    self.cg_inert_prs_for_chars = {}

    # Other
    self.nf = nf
    self.max_dets = max_dets
    self.num_chars = num_chars
    self.smoothness_bound = smoothness_bound
    self.unit_rank = unit_rank

  def rank(self, lattice):
    if lattice=='cg' and self.cg_into_hnf:
      return self.hnf_relns_mat.ncols()

  def det(self, lattice, real_field, ug_matrix_index):
    if lattice=='cg' and self.cg_into_hnf:
      return prod_all_elems(self.hnf_relns_mat.diagonal())
    elif lattice=='ug':
      return self.regulator(real_field, ug_matrix_index)

  def ug_largest_elem_div(self, matrix_index):
    return max(self.ug_compact_rep_mats[matrix_index].elementary_divisors())

  def regulator(self, real_field, matrix_index):
    '''
    A naive implementation of Cohen 6.5.7 p357
    I dont think Cohen defines what this algorithm will return if the subgroup these units generate does not have full rank.
    I wish the algorithm would return 0 in this case but its possible this algorithm could return a nonzero number in this case and that its value could depend on ind_drop_col.
    However, Cohen 6.5.9 does not say anything about how to handle this case.
    Without handling this case correctly the condition hR >= z*sqrt(2) is not valid.
    :param real_field:
    :param matrix_index: -1 if this is the regulator of the combined units mat, o/w the index of particular units mat that we want to select
    :return:
    '''
    start_regulator_time = time.process_time()
    ind_drop_col = self.ug_reg_cols_to_drop[matrix_index]

    # step: construct matrix of element embeddings
    if matrix_index >= 0:
      elems = list(self.ug_elems_mats[matrix_index].columns()[0])
    elif matrix_index == -1:
      if self.ug_combined_elems_mat == None:
        raise ValueError("ug_combined_elems_mat is not set")
      elems = list(self.ug_combined_elems_mat.columns()[0])

    if self.ug_compact_rep_mats[matrix_index].nrows()==0:
      return real_field(0)
    #TODO make this return a matrix with abs and log already applied
    embedds = abs_std_embedding(elems, real_field)
    embd_of_elems_mat = Matrix(real_field, embedds, sparse=False).transpose()
    count_zeros=0
    zero_positions=[]
    for r in range(embd_of_elems_mat.nrows()):
      for c in range(embd_of_elems_mat.ncols()):
        # if embd_of_elems_mat[r,c].abs()<0.0000001:
        if embd_of_elems_mat[r,c].abs() != 0:
          count_zeros+=1
          zero_positions.append([r,c])

    def log(x):
      if x > 0:
        return x.log()
      raise ValueError("input x={} <= 0 so we cant take its log.".format(x))

    try:
      embd_of_elems_mat = embd_of_elems_mat.apply_map(log, real_field, sparse=False)
    except AttributeError as e:
      print(embd_of_elems_mat)
      raise e

    if matrix_index >= 0:
      embd_of_units_mat = self.ug_compact_rep_mats[matrix_index] * embd_of_elems_mat
    elif matrix_index == -1:
      embd_of_units_mat = self.ug_combined_compact_rep_mat * embd_of_elems_mat

    # I did a quick test and found no units were equal to 1. Not even worth checking ?
    nonzero_rows = []
    original_nrows = embd_of_units_mat.nrows()
    for i in range(original_nrows):
      for cc in range(embd_of_units_mat.ncols()):
        if embd_of_units_mat[i,cc].abs() > 0.01:
          nonzero_rows.append(i)
          break
    embd_of_units_mat = embd_of_units_mat.matrix_from_rows(nonzero_rows)

    if asserts_on():
      for rr in range(embd_of_units_mat.nrows()):
        row_sum = real_field(0)
        for _ in embd_of_units_mat[rr]:
          row_sum += _
        max_embd_sum = 0.01
        assert row_sum.abs() < 0.01, "Expected unit to have sum of logs have absolute value less than {} but got:\n{}.\nembd_of_units_mat={}\ncurrent index was {}".format(max_embd_sum, row_sum.abs(), embd_of_units_mat[rr], rr)

    embd_mat_nrows = embd_of_units_mat.nrows()
    embd_mat_ncols = embd_of_units_mat.ncols()
    assert embd_mat_ncols == self.unit_rank + 1
    embd_of_units_mat.swap_columns(0, ind_drop_col)
    embd_with_dropped_col = embd_of_units_mat[:,1:embd_mat_ncols]
    embd_submats = [embd_with_dropped_col[i-self.unit_rank:i] for i in range(self.unit_rank, embd_mat_nrows)]
    embd_submat_dets = [s_det(_) for _ in embd_submats]

    found_nonzero_det = False
    prev_reg = real_field(0)
    cand_regulator = real_field(0)
    err_bound = real_field(0.1)
    all_regs = []
    for _ in range(0, len(embd_submats)):
      if prev_reg > err_bound and cand_regulator == real_field(0):
        raise ValueError("Step {}. The previous regulator was nonzero. cand_rglr=0 but it should be a multiple of the true regulator.".format(_))
      subdet = embd_submat_dets[_]
      if subdet.abs() < (self.regulator_lower_bound + err_bound):
        continue
      if found_nonzero_det == False:
        found_nonzero_det = True
        cand_regulator = subdet
        all_regs.append(cand_regulator)
      else:
        prev_reg = cand_regulator
        cand_regulator = robust_rgcd(cand_regulator, subdet, self.regulator_lower_bound)
        if cand_regulator.abs() > 0 and cand_regulator.abs() < (self.regulator_lower_bound + err_bound):
          raise ArithmeticError('Precision exceeded when computing robust rgcd. Subdet index was i={} of {}'.format(_, len(embd_submats)))
        all_regs.append(cand_regulator)
    if cand_regulator.abs() == 0:
      raise ValueError("regulator is equal to 0 exactly. Possibilities:\na) unit submodule does not have full rank (nunits={})\nb) Precision too small (prec={})\nc) a different col needs to be dropped (Matrix_index={}, ind_drop_col={})\nd) All the units are ROUs. This is most likely if the compact rep mat has small height (log2(height)={})".format(
        self.ug_compact_rep_mats[matrix_index].nrows(),
        real_field.precision(),
        matrix_index,
        ind_drop_col,
        RR(self.ug_compact_rep_mats[matrix_index].height()).log(2).n(10)
      ))
    if cand_regulator.abs() > 0 and cand_regulator.abs() < (self.regulator_lower_bound + err_bound):
      raise ArithmeticError('Precision exceeded when computing rgcd')

    timers.linalg_regulator += time.process_time() - start_regulator_time
    return cand_regulator

  def saturate_ug(self, pr_to_saturate, unit_rank, matrix_index):
    chars = self.character_mats_dict[matrix_index][pr_to_saturate]
    chars_mat = (self.ug_compact_rep_mats[matrix_index] * chars).change_ring(GF(pr_to_saturate))
    num_units = self.ug_compact_rep_mats[matrix_index].nrows()

    if chars_mat.rank() == num_units:
      raise ValueError("kernel is empty so either there are no dependencies in the given set of {} units (get more) OR units lattice must already be saturated for pr={}. Unit rank={}".format(num_units, pr_to_saturate, unit_rank))
    kern_mat = chars_mat.left_kernel().basis_matrix()
    mat_of_solns = kern_mat.change_ring(ZZ)

    new_units_compact_rep = mat_of_solns * self.ug_compact_rep_mats[matrix_index]
    sublattice_of_sqr_units = (new_units_compact_rep.stack(pr_to_saturate*identity_matrix(new_units_compact_rep.ncols()))).change_ring(ZZ)
    pr_pow_units = sublattice_of_sqr_units.left_kernel(algorithm='flint', basis='LLL').basis_matrix().change_ring(ZZ)
    pr_root_units = pr_pow_units[:, :new_units_compact_rep.nrows()] * new_units_compact_rep
    pr_root_units /= ZZ(pr_to_saturate)
    pr_root_units = pr_root_units.change_ring(ZZ)

    self.ug_compact_rep_mats[matrix_index] = self.ug_compact_rep_mats[matrix_index].stack(pr_root_units[-1])

  def extend_relns_and_compute_basis(self, new_relns, new_elems, candidate_class_num=None):
    '''
    :param relns_basis_hnf: takes a pair [mat in HNF, transformation mat] and a matrix with the same number of cols as the hnf
    :param new_relns:
    :param sqr_cg_relns_mat: used to solve the system of equations. Not modified
    :return:
    '''
    nn = new_relns.nrows()
    if nn > 50:
      print("WARNING: Stacking more than 50 rows. Trying to compute HNF of a matrix whose nonzero rows are not close to being a square matrix will probably be slow")

    if self.cg_into_hnf:
      pivots = self.hnf_relns_mat.pivots()
      for i in range(new_relns.nrows()):
        self.hnf_relns_mat, pivots = stein_hnf.add_row(self.hnf_relns_mat, new_relns[i, :], pivots, False)

    self.cg_relns_mat = self.cg_relns_mat.stack(new_relns)
    self.cg_elems_mat = self.cg_elems_mat.stack(new_elems)
    self.cg_compact_rep_mat = identity_matrix(self.cg_elems_mat.nrows())

  def add_rows_to_units_lattice_with_vollmer(self, new_relns, new_elems, forced_sat_pr, matrix_index):
    '''
    :param new_relns:
    :param new_elems:
    :param forced_sat_pr: In production this should be None. For saturation testing this can be used to force a particular prime to divide the lattice's det
    :param matrix_index: either 0 or 1
    :return:
    '''
    new_relns_nr = new_relns.nrows()
    if new_relns_nr != new_elems.nrows():
      raise ValueError("number of new elems must equal number of new relns")
    if not self.mats_for_vollmer[matrix_index].is_square():
      raise ValueError("sqr mat is not sqr")
    nc = new_relns.ncols()
    self.num_initial_units = new_relns.nrows()

    # step 1: set the initial matrix
    mat_of_solns, dnm = self.mats_for_vollmer[matrix_index]._solve_flint(new_relns.augment(zero_matrix(new_relns_nr, self.mats_for_vollmer[matrix_index].ncols() - nc)), right=False)

    # step 2a: divide out denominators
    dnms = [-dnm] * new_relns_nr
    for i in range(0, mat_of_solns.nrows()):
      row = mat_of_solns.row(i)
      row_gcd = GCD_list(row)
      if row_gcd > 1:
        # dividing out row gcd is essential. This can probably be sped up a lot
        dnms[i] = dnms[i] / row_gcd
        for cc in range(mat_of_solns.ncols()):
          mat_of_solns[i, cc] = mat_of_solns[i, cc] / row_gcd

    if forced_sat_pr == None:
      forced_sat_pr = 1
    mat_of_solns *= forced_sat_pr
    dnms = [forced_sat_pr * _ for _ in dnms]

    # step 3: build the bottom of the new transformation mat
    mat_of_solns_extended = mat_of_solns.augment(diagonal_matrix(dnms))

    # step 4: build the new transformation mat
    new_ug_compact_rep_mat = self.ug_compact_rep_mats[matrix_index]
    zero_mat = zero_matrix(self.ug_compact_rep_mats[matrix_index].nrows(), new_relns_nr)
    new_ug_compact_rep_mat = new_ug_compact_rep_mat.augment(zero_mat)
    new_ug_compact_rep_mat = new_ug_compact_rep_mat.stack(mat_of_solns_extended)

    if new_ug_compact_rep_mat.ncols() != (self.ug_compact_rep_mats[matrix_index].ncols() + new_relns_nr):
      raise ValueError('wrong num cols in new transformation')

    self.mats_for_vollmer[matrix_index] = self.mats_for_vollmer[matrix_index].augment(zero_matrix(self.mats_for_vollmer[matrix_index].nrows(), new_relns_nr))
    self.mats_for_vollmer[matrix_index] = self.mats_for_vollmer[matrix_index].stack(mat_of_solns_extended)

    self.ug_compact_rep_mats[matrix_index] = new_ug_compact_rep_mat
    self.ug_elems_mats[matrix_index] = self.ug_elems_mats[matrix_index].stack(new_elems)

  def insersection_of_unit_submodules(self):
    if self.ug_compact_rep_mats[0].nrows() != self.ug_compact_rep_mats[1].nrows():
      raise ValueError("bug")
    if self.elems_of_V_as_Zcombination_of_U_elems == []:
      raise ValueError("this function cant be called until elems_of_V_as_Zcombination_of_U_elems is set")
    self.ug_combined_elems_mat = self.ug_elems_mats[0].stack(self.ug_elems_mats[1])

    nr = self.ug_compact_rep_mats[0].nrows()
    if nr == 0:
      raise ValueError("ug_compact_rep_mats should not have 0 rows")
    nc = self.ug_compact_rep_mats[1].ncols()
    zero_mat = zero_matrix(nr, nc)
    ug_combined_compact_rep_mat = self.ug_compact_rep_mats[0].augment(zero_mat)
    ug_combined_compact_rep_mat = ug_combined_compact_rep_mat.stack(
      zero_matrix(self.ug_compact_rep_mats[1].nrows(), self.ug_compact_rep_mats[0].ncols()).augment(self.ug_compact_rep_mats[1]))

    mat_to_stack = Matrix(ZZ,self.elems_of_V_as_Zcombination_of_U_elems)
    mat_to_stack = mat_to_stack.augment(zero_matrix(mat_to_stack.nrows(), 2*self.num_initial_units))
    ug_combined_compact_rep_mat = ug_combined_compact_rep_mat.stack(mat_to_stack)
    self.ug_combined_compact_rep_mat = ug_combined_compact_rep_mat

  def project_one_vollmer_mat_onto_the_other(self):
    if self.mats_for_vollmer[0].ncols() != self.hnf_relns_mat.ncols():
      raise ValueError("This function must be called before the initial units are constructed")
    project_submodule_start = time.process_time()
    for i in range(len(self.ug_compact_rep_mats)):
      if self.ug_compact_rep_mats[i] == None:
        raise ValueError("ug_units_compact_rep_mats must all be set")
    id_mat = identity_matrix(self.mats_for_vollmer[0].nrows())
    ug_trans_V_onto_U_flint_pair = self.mats_for_vollmer[0]._solve_flint(self.mats_for_vollmer[1], right=False)
    ug_trans_V_onto_U = ug_trans_V_onto_U_flint_pair[0] / ug_trans_V_onto_U_flint_pair[1]
    num_inverse_denom1 = 0
    for i in range(ug_trans_V_onto_U.nrows()):
      row = ug_trans_V_onto_U[i]
      if row.denominator() == 1:
        self.elems_of_V_as_Zcombination_of_U_elems.append(list(row) + list(-id_mat[i]))
        num_inverse_denom1 += 1
    timers.linalg_project_submodule = time.process_time() - project_submodule_start

### Non-class methods

def add_rows_with_separate_transformations(hnf_data, kernel_mat, mat_to_stack, original_mat=None):
  '''
  :param hnf_data: assumed to be contain a full-rank hnf
  :param kernel_mat:
  :param mat_to_stack:
  :return:
  '''
  [hnf, U0] = hnf_data
  if original_mat != None:
    if original_mat.ncols() != hnf.ncols():
      raise ValueError("number of columns is wrong")
    assert U0*original_mat == hnf
  ncols = hnf.ncols()
  if ncols != mat_to_stack.ncols():
    raise ValueError("mat_to_stack is the wrong size")
  if asserts_on() and hnf.rank() < ncols:
    raise ValueError("hnf should have full rank")

  transformation = U0.stack(kernel_mat)
  full_hnf = hnf.stack(zero_matrix(transformation.nrows() - hnf.nrows(), hnf.ncols()))
