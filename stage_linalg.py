import time

from sage.all import *

from src.stage_filter import assert_rows_multiply_to_correct_ideal, assert_rat_pr_mat_rows_multiply_to_correct_ideal
from src.io_utils import write_int_mat, read_int_mat, read_nf_elem_mat, read_fb_inds
from src.stage_linalg_internal import *
import src.constants
import src.timers as timers
from src.lat_manager import *

def stage_linalg(elems_mat, nf, relns_mat, compute_reg):
  linalg_start_time = time.process_time()
  use_HNF = False
  compact_rep_mat, elem_divs = magma_hnf(relns_mat, use_HNF, compute_reg)
  # lll_compact_rep_mat, transformation_mat, elem_divs = magma_hnf(relns_mat, use_HNF) # LLL reduced transformation matrix is currently not working
  real_field = RealField(1000, sci_not=True)
  sig = nf.signature()
  unit_rank = sig[0]+sig[1]-1
  if compute_reg:
    reg = cohen_regulator(real_field, elems_mat.columns()[0], compact_rep_mat, unit_rank)
    # reg = cohen_regulator(real_field, elems_mat.columns()[0], lll_compact_rep_mat, unit_rank, transformation_mat)
    sys.stderr.write("reg={} for dimension {}x{} compact rep mat".format(reg.n(10), lll_compact_rep_mat.nrows(),
                                                                         lll_compact_rep_mat.ncols())+os.linesep)
  timers.linalg_total_time = time.process_time()-linalg_start_time
  timers.linalg_total_time += timers.linalg_magma_convert
  timers.linalg_total_time += timers.linalg_magma_hnf_time

def magma_hnf(relns_mat, use_hnf, compute_reg):
  mag_convert_time_start = magma.cputime()
  mag_sp_reln_mat = magma(relns_mat).Matrix()  # Matrix will be Sparse type at first
  timers.linalg_magma_convert = magma.cputime(mag_convert_time_start)
  mag_hnf_time_start = magma.cputime()
  # mag_sp_reln_mat.HermiteForm(nvals=2)
  if use_hnf:
    magma.function_call('HermiteForm', mag_sp_reln_mat, params={'Optimize':False}, nvals=2)
  else:
    elem_divs = mag_sp_reln_mat.ElementaryDivisors()
    elem_divs = [_ for _ in elem_divs if _ != 1]
    print("elem divs {}".format(elem_divs))

    dense_time_start=time.process_time()
    relns_dense=relns_mat.dense_matrix()
    dense_time=time.process_time()-dense_time_start

    sys.stderr.write("height orig {}. Convert time to get this data {:.1f}".format(relns_dense.height(), dense_time) + os.linesep)
    compact_rep_mat_magma = mag_sp_reln_mat.KernelMatrix()
    if compute_reg:
      lll_compact_rep_mat, transformation_mat, rank = magma.function_call('BasisReduction', mag_sp_reln_mat, nvals=3)
      lll_compact_rep_mat = lll_compact_rep_mat.sage()
      transformation_mat = transformation_mat.sage()
      sys.stderr.write("height orig {}\nheight lll reduced {}".format(relns_dense.height(), lll_compact_rep_mat.dense_matrix().height()) + os.linesep)
  timers.linalg_magma_hnf_time = magma.cputime(mag_hnf_time_start)
  compact_rep_mat=compact_rep_mat_magma.sage()
  return compact_rep_mat, elem_divs
  # return lll_compact_rep_mat, transformation_mat, elem_divs

def cohen_regulator(real_field, elems, compact_rep_mat, unit_rank, transformation=None):
  '''
  Cohen 6.5.7 p357
  I dont think Cohen defines what this algorithm will return if the subgroup these units generate does not have full rank.
  I wish the algorithm would return 0 in this case but its possible this algorithm could return a nonzero number in this case and that its value could depend on ind_drop_col.
  However, Cohen 6.5.9 does not say anything about how to handle this case.
  Without handling this case correctly the condition hR >= z*sqrt(2) is not valid.
  :param real_field:
  :return:
  '''
  start_regulator_time = time.process_time()
  ind_drop_col = 0
  mat_height=compact_rep_mat.height()
  log_mat_height=RR(mat_height).log(2).n(10)
  max_entry_str='max entry of compact rep mat'
  sys.stderr.write("{} {}\nlog {} {}".format(max_entry_str, mat_height, max_entry_str, log_mat_height))

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

  #TODO move
  def log(x):
    if x > 0:
      return x.log()
    raise ValueError("input x={} <= 0 so we cant take its log.".format(x))

  try:
    embd_of_elems_mat = embd_of_elems_mat.apply_map(log, real_field, sparse=False)
  except AttributeError as e:
    print(embd_of_elems_mat)
    raise e

  # if asserts_on():
  #   # if any sparse units have been found check these. Even only checking sparse this code is slow so I usually disable
  #   num_units_tested=0
  #   count_units_eq_to_one =0
  #   for rr in range(unit_rank):
  #     print("starting unit {}".format(rr))
  #     len_kern_vec = len([compact_rep_mat[rr,c] for c in range(compact_rep_mat.ncols()) if compact_rep_mat[rr,c]!=0])
  #     if len_kern_vec < 5:
  #       unit = prod_all_elems([elems_mat[c,0]**compact_rep_mat[rr,c] for c in range(kern.ncols()) if compact_rep_mat[rr,c]!=0])
  #       if unit==1:
  #         count_units_eq_to_one+=1
  #         print(count_units_eq_to_one)
  #       assert unit.norm().abs() == 1, "expected the elem to be a unit but its norm was {}".format(unit.norm().abs())
  #       num_units_tested += 1
  #     else:
  #       print("length of kernel vector {}".format(len_kern_vec))
  #   print("Asserted that {} units had norm equal to 1".format(num_units_tested))

  # print("kern vecs:")
  # print(relns_mat.rows())

  #TODO to allow compact reps of elems I would do something like this, then multiply by ug_units_compact_rep_mat
  # embd_of_elems_including_compact_rep_mat = self.ug_elems_compact_rep_mat * embd_of_elems_mat
  if transformation == None:
    embd_of_units_mat = compact_rep_mat * embd_of_elems_mat
  else:
    # embd_of_elems_mat = transformation * embd_of_elems_mat
    embd_of_units_mat = compact_rep_mat.transpose() * embd_of_elems_mat

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
  assert embd_mat_ncols == unit_rank + 1
  embd_of_units_mat.swap_columns(0, ind_drop_col)
  embd_with_dropped_col = embd_of_units_mat[:,1:embd_mat_ncols]
  embd_submats = [embd_with_dropped_col[i-unit_rank:i] for i in range(unit_rank, embd_mat_nrows)]
  embd_submat_dets = [s_det(_) for _ in embd_submats] # using s_det vs determinant() does not make a significant different

  found_nonzero_det = False
  prev_reg = real_field(0)
  cand_regulator = real_field(0)
  err_bound = real_field(0.01)
  all_regs = []
  regulator_lower_bound=0.2
  for _ in range(0, len(embd_submats)):
    if prev_reg > err_bound and cand_regulator == real_field(0):
      raise ValueError("Step {}. The previous regulator was nonzero. cand_rglr=0 but it should be a multiple of the true regulator.".format(_))
    subdet = embd_submat_dets[_]
    if subdet.abs() < (regulator_lower_bound + err_bound):
      continue
    if found_nonzero_det == False:
      found_nonzero_det = True
      cand_regulator = subdet
      all_regs.append(cand_regulator)
    else:
      prev_reg = cand_regulator
      cand_regulator = rgcd(cand_regulator, subdet, regulator_lower_bound)
      if cand_regulator.abs() > 0 and cand_regulator.abs() < (regulator_lower_bound + err_bound):
        raise ArithmeticError('Precision exceeded when computing rgcd. Subdet index was i={} of {}. Cand regulator prior to problem {}'.format(_, len(embd_submats), cand_regulator))
      all_regs.append(cand_regulator)
  if cand_regulator.abs() == 0:
    #   print("Warning: this set of units does not have full rank")
    raise ValueError("regulator is equal to 0 exactly. Possibilities:\na) unit submodule does not have full rank (nunits={})\nb) Precision too small (prec={})\nc) a different col needs to be dropped (Matrix_index={}, ind_drop_col={})\nd) All the units are ROUs. This is more likely if the compact rep mat has small height (log2(height)={})".format(
      compact_rep_mat.nrows(),
      real_field.precision(),
      0,
      ind_drop_col,
      RR(compact_rep_mat.height()).log(2).n(10)
    ))
  if cand_regulator.abs() > 0 and cand_regulator.abs() < (regulator_lower_bound + err_bound):
    raise ArithmeticError('Precision exceeded when computing rgcd')

  timers.linalg_regulator += time.process_time() - start_regulator_time
  return cand_regulator

def stage_linalg_system_solving_and_ug_saturation(Delta, alg_prime_str_dict, fb, log_level, nf, params_filename, rat_fb, sieve_filtered_mat_path):
  linalg_relns_to_load_ratio = get_linalg_relns_to_load_ratio(path_to_params_dir, params_filename, False)
  submat_write_rel_dir = rel_dirname
  sig = nf.signature()
  unit_rank = sig[0]+sig[1]-1
  print("unit rank {}".format(unit_rank))
  largest_det_multiple = ZZ(2)**80  # we factor this multiple so it can't be too big

  #TODO these constants are both important and very hard to tune
  if unit_rank == 1:
    num_ug_additional_rows = 10
  if unit_rank == 2:
    num_ug_additional_rows = 12
  if unit_rank == 3:
    num_ug_additional_rows = 36
  cg_num_initial_relns = 3
  cg_max_additional_relns = 200
  col_inds_to_drop = [0, 0]
  prec = 4000
  [_, _, fb_for_cn, _, alg_prime_str_dict_for_cn, _, _] = fac_base_data(nf, round(nf.bach_bound()))
  cn_formula_approx = class_num_formula_approx(fb_for_cn, alg_prime_str_dict_for_cn, prec)
  # takes too long to compute these invariants b/c its in Pari
  # [test_data_regulator, test_data_class_group, test_data_class_num] = get_invariants_from_number_field_test_data(test_data)
  test_data_regulator = None
  try:
    [class_num, class_group, regulator] = linalg_by_system_solving_and_ug_saturation(log_level, sieve_filtered_mat_path, cn_formula_approx,
                                                                                     linalg_relns_to_load_ratio, submat_write_rel_dir, unit_rank,
                                                                                     pol, fb, rat_fb, Delta, smoothness_bound, alg_prime_str_dict,
                                                                                     test_data_regulator, largest_det_multiple,
                                                                                     num_ug_additional_rows, cg_max_additional_relns,
                                                                                     cg_num_initial_relns, col_inds_to_drop, prec)
  except ValueError as e:
    print_profile(0, max_num_primes, full_rank_timing)
    print("FAIL in stage {}: ValueError failed {}".format("Linalg", e))
    raise e
  except AssertionError as e:
    print_profile(0, max_num_primes, full_rank_timing)
    print("FAIL in stage {}: Assertion failed {}".format("Linalg", e))
    raise e
  return class_group, class_num, regulator, test_data_regulator


def linalg_by_system_solving_and_ug_saturation(log_level, partial_mat_path, full_det_approx, ratio_of_relations_to_load,
                                               submat_write_rel_dir, unit_rank, pol, fb, rat_fb, Delta, smoothness_bound, alg_pr_str_dict, true_reg_for_testing, largest_determinant_multiple, ug_num_units_per_matrix, cg_max_additional_relns, cg_num_initial_relns, reg_cols_to_drop, prec):
  '''
  :param log_level:
  :param partial_mat_path:
  :param full_det_approx:
  :param ratio_of_relations_to_load:
  :param submat_write_rel_dir:
  :param unit_rank:
  :param pol:
  :param fb:
  :param rat_fb:
  :param Delta:
  :param smoothness_bound:
  :param alg_pr_str_dict:
  :param true_reg_for_testing:
  :param largest_determinant_multiple:
  :return:
  '''
  #TODO remove fb param so that task.linalg can be run independently of filtering
  print("\nStage: Linalg\n-------------------------------\n")
  from src.stage_linalg_internal import char_mat

  real_field = sage.rings.real_mpfr.RealField(prec)
  print("Starting Precision={}".format(prec))
  start_linalg_time = time.process_time()

  linalg_setup_time_start = time.process_time()
  nf = NumberField(pol, names='Y')

  if ratio_of_relations_to_load < 1:
    raise ValueError("ratio_of_relations_to_load must be >=1 in order to get a matrix of full rank")

  # Experimentally, when sieving is used rows_to_cols_ratio_for_load needs a large value or we won't get the right value for class_number
  rows_to_cols_ratio_for_load = ratio_of_relations_to_load # HNF seems to be computed slightly faster when this is a little larger than the excess
  if rows_to_cols_ratio_for_load < ratio_of_relations_to_load:
    raise ValueError("rows_to_cols_ratio_for_load should be at least as large as excess")

  sp_reln_mat = read_int_mat(log_level, partial_mat_path + "_relns", None, rows_to_cols_ratio_for_load)

  sp_rat_reln_mat_path = partial_mat_path + "_sp_rat_relns"

  nrows = sp_reln_mat.nrows()
  elems_mat= read_nf_elem_mat(log_level, partial_mat_path + "_elems", None, nrows, nf)

  from src.paths import basepath_artifacts
  with open(basepath_artifacts + partial_mat_path + "_fb_inds", "r") as fp:
    lines = fp.readlines()
    if len(lines)>1:
      raise ValueError("wrong number of lines")
    line = lines[0]
    line = line[1:len(line)-1]
    nonzero_inds = [int(_) for _ in line.split(",")]
  fb_filtered = [fb[_] for _ in range(len(fb)) if _ in nonzero_inds]

  if asserts_on():
    # assert assert_rat_pr_mat_rows_multiply_to_correct_ideal(sp_rat_reln_mat, elems_mat, rat_fb)
    assert assert_rows_multiply_to_correct_ideal(sp_reln_mat, elems_mat, nf, fb_filtered)
    assert nrows == elems_mat.nrows(), "number of rows of submat={} is not equal to {} the number of rows of elems_mat".format(nrows, elems_mat.nrows())

  ncols = sp_reln_mat.ncols()
  if nrows < 2*ncols:
    raise ValueError("Number of rows is less than 2x the number of cols {} vs {} which is needed for the approach here of using 2-random square matrices to have a reasonable propability of working. Try increasing the sieve radius.".format(nrows, ncols))
  print("mat nrows = {}, ncols = {}".format(nrows, ncols))

  # Setup and initial HNF method 1
  nc = sp_reln_mat.ncols()
  row_ind_add_rows_begins_at = nc
  vollmer_mat_0=sp_reln_mat[:nc]

  elems_mat_0 = elems_mat.matrix_from_rows(range(nc))
  if asserts_on():
    assert_rows_multiply_to_correct_ideal(vollmer_mat_0, elems_mat_0, nf, fb_filtered)
    for i in range(elems_mat.nrows()):
      if elems_mat[i] == 0:
        raise ValueError("Cant have 0 elems")
  max_dets = 5
  num_chars = 5

  relns_with_sqr_mat_removed_from_top = sp_reln_mat[nc:]
  # first we attempt to construct a relations matrix with no rows shared with the 1st relations matrix
  piv_rows_for_second_sqr_mat = pivot_cols_fastest_available(relns_with_sqr_mat_removed_from_top.transpose())
  vollmer_mat_1 = relns_with_sqr_mat_removed_from_top.matrix_from_rows(piv_rows_for_second_sqr_mat)
  if len(piv_rows_for_second_sqr_mat) != nc:
    print("UG 2nd matrix construction. Could not find a second full rank matrix. Rank {} vs ncols {}. Num rows examined {}. Trying to patch the non-pivot columns by adding rows from the entire relations matrix".format(len(piv_rows_for_second_sqr_mat), nc, relns_with_sqr_mat_removed_from_top.nrows()))
    nonpiv_col_inds,_ =  probable_nonpivot_columns_fastest_available(vollmer_mat_1)
    for i in nonpiv_col_inds:
      nonzero_row_inds = [_ for _ in range(sp_reln_mat.nrows()) if sp_reln_mat[_,i] !=0]
      relns_with_sqr_mat_removed_from_top = relns_with_sqr_mat_removed_from_top.stack(sp_reln_mat[nonzero_row_inds[-1]])
      elems_mat = elems_mat.stack(elems_mat[nonzero_row_inds[-1]])
      # sp_rat_reln_mat = sp_rat_reln_mat.stack(sp_rat_reln_mat[nonzero_row_inds[-1]])

  piv_rows_for_second_sqr_mat = pivot_cols_fastest_available(relns_with_sqr_mat_removed_from_top.transpose())
  nonpivots_for_second_sqr_mat = sorted(list(set(range(relns_with_sqr_mat_removed_from_top.nrows())) - set(piv_rows_for_second_sqr_mat)))
  vollmer_mat_1 = relns_with_sqr_mat_removed_from_top.matrix_from_rows(piv_rows_for_second_sqr_mat)
  if len(piv_rows_for_second_sqr_mat) != nc:
    raise ValueError("Patching could not find a second full rank matrix. Rank {} vs ncols {}. Num rows examined {}. Num rows of full relations matrix".format(len(piv_rows_for_second_sqr_mat), nc, relns_with_sqr_mat_removed_from_top.nrows(), sp_reln_mat.nrows()))
  else:
    print("Successful patch")

  number_of_rows_used_in_ug = len(piv_rows_for_second_sqr_mat) + vollmer_mat_0.nrows() # number_of_rows_used_in_ug is a measure of how many relations I need sieve in order to compute the unit group
  elems_mat_1=elems_mat[nc:].matrix_from_rows(piv_rows_for_second_sqr_mat)

  print("Converting sparse mats to dense mats")
  mats_for_vollmer = [vollmer_mat_0.dense_matrix(), vollmer_mat_1.dense_matrix()]
  ug_elems_mats = [elems_mat_0, elems_mat_1]
  num_ug_mats = len(ug_elems_mats)

  timers.linalg_setup_time = time.process_time() - linalg_setup_time_start

  ## Contains timer
  print("Step 0: CG initial HNF")
  reln_basis_mat = initial_hnf(mats_for_vollmer[0])

  # rat_relns_mat=sp_rat_reln_mat[:nc]
  rat_relns_mat=None
  lat_manager = LatticeManager(reln_basis_mat, mats_for_vollmer[0], rat_relns_mat, elems_mat_0, ug_elems_mats, nf, mats_for_vollmer, unit_rank, max_dets,
                               num_chars, smoothness_bound, true_reg_for_testing, reg_cols_to_drop)

  ### CG Add initial relns
  print("Step 1: CG Add initial relns")
  initial_linalg_start = time.process_time()
  max_num_relns_needed=cg_num_initial_relns + cg_max_additional_relns
  inds_of_distinct_relns = piv_rows_for_second_sqr_mat[:max_num_relns_needed]
  new_relns = relns_with_sqr_mat_removed_from_top.matrix_from_rows(inds_of_distinct_relns)
  new_elems = elems_mat[nc:].matrix_from_rows(inds_of_distinct_relns)
  new_elems_for_cg_initial = new_elems[:cg_num_initial_relns]
  new_relns_for_cg_initial = new_relns[:cg_num_initial_relns]
  lat_manager.extend_relns_and_compute_basis(new_relns_for_cg_initial, new_elems_for_cg_initial)

  new_elems = new_elems[cg_num_initial_relns:]
  new_relns = new_relns[cg_num_initial_relns:]
  initial_cn = lat_manager.det('cg', None, None)
  if initial_cn == 0:
    raise ValueError("Initial class number is 0")
  lower_bound_fraction = 0.1  # chosen arbitrarily
  regulator_lower_bound = max(lower_bound_fraction*(full_det_approx/initial_cn), 0.2)
  lat_manager.regulator_lower_bound = regulator_lower_bound
  num_cg_relns_used = cg_num_initial_relns
  print("Initial CN is {} and regulator lower bound is {}".format(initial_cn, RR(lat_manager.regulator_lower_bound).n(15)))

  print("Initial relns time {}".format(time.process_time() - initial_linalg_start))

  ### UG Non-saturation initial units
  print("Step 2: UG Add initial relns")
  linalg_step2_ug_initial_relns_start=time.process_time()
  lat_manager.project_one_vollmer_mat_onto_the_other() # used when the intersection of the compact rep submodules is computed

  index_of_nonpiv_rows_for_second_sqr_mat_in_full_mat = [_+nc for _ in nonpivots_for_second_sqr_mat]

  index_of_nonpiv_rows_for_second_sqr_mat_in_full_mat = list(reversed(index_of_nonpiv_rows_for_second_sqr_mat_in_full_mat))
  ug_inds = index_of_nonpiv_rows_for_second_sqr_mat_in_full_mat[:ug_num_units_per_matrix]
  forced_sat_pr = 1
  print("starting to add units with volmer")
  ug_add_rows_start = time.process_time()
  lat_manager.add_rows_to_units_lattice_with_vollmer(sp_reln_mat.matrix_from_rows(ug_inds).dense_matrix(),
                                                      elems_mat.matrix_from_rows(ug_inds), forced_sat_pr, 0)
  lat_manager.add_rows_to_units_lattice_with_vollmer(sp_reln_mat.matrix_from_rows(ug_inds).dense_matrix(),
                                                     elems_mat.matrix_from_rows(ug_inds), forced_sat_pr, 1)
  compute_intersection_of_unit_submodules_before_saturation=False # I want to avoid computing insersection_of_unit_submodules twice b/c its slow
  if compute_intersection_of_unit_submodules_before_saturation:
    lat_manager.insersection_of_unit_submodules()
    initial_combined_reg = lat_manager.regulator(real_field, -1)
  timers.linalg_step1_cg_add_relns = time.process_time() - initial_linalg_start

  print("Num bits of height of units mats are: {} and {}".format(RR(lat_manager.ug_compact_rep_mats[0].height()).log(2).n(15), RR(lat_manager.ug_compact_rep_mats[1].height()).log(2).n(15)))
  timers.linalg_ug_system_solving = time.process_time() - ug_add_rows_start

  number_of_rows_used_in_ug += ug_num_units_per_matrix
  for j in range(num_ug_mats):
    max_num_times_double_prec=13
    multiple=2
    try:
      reg_candidate = lat_manager.regulator(real_field, j)
    except ArithmeticError as e:
      for i in range(max_num_times_double_prec):
        prec = floor(multiple*prec)
        try:
          reg_candidate, real_field = raise_precision(prec, lat_manager, j)
          break
        except ArithmeticError as e:
          continue
      if reg_candidate == 0:
        raise ValueError("candidate regulator is 0 after raising precision {} times. #Units={}".format(max_num_times_double_prec, ug_num_units_per_matrix))

  rg0 = lat_manager.regulator(real_field, 0)
  rg1 = lat_manager.regulator(real_field, 1)

  initial_reg_gcd = rgcd(rg0, rg1, lat_manager.regulator_lower_bound)

  print("initial regulator is {} with {} units".format(lat_manager.det('ug', real_field, 0).n(10), ug_num_units_per_matrix))

  full_multiple = get_full_multiplier(lat_manager.det('cg', None, None), initial_reg_gcd, full_det_approx)
  if full_multiple > largest_determinant_multiple:
    raise ValueError("the number multiplying the regulator is too large to factor, {} bits. Num additional rows {}. rg0={}, rg1={}".format(RR(full_multiple).log(2).n(10), ug_num_units_per_matrix, rg0, rg1))
  if full_multiple == 1:
    lat_manager.insersection_of_unit_submodules()
    timers.linalg_step2_ug_initial_relns=time.process_time() - linalg_step2_ug_initial_relns_start
    return finish_up(lat_manager, log_level, real_field, start_linalg_time, submat_write_rel_dir, number_of_rows_used_in_ug, cg_num_initial_relns)
  facs = full_multiple.factor()
  timers.linalg_step2_ug_initial_relns=time.process_time() - linalg_step2_ug_initial_relns_start

  ### Saturate UG
  print("Step 3: UG Saturation")
  linalg_step3_ug_saturate_start=time.process_time()
  # Initial UG Saturation
  for j in range(num_ug_mats):
    reg_prs_that_need_saturation = [_[0] for _ in facs]
    print("Starting getting character mats for mat {}".format(j))
    for pr in reg_prs_that_need_saturation:
      prs_start_time = time.process_time()
      prs_for_chars = primes_of_order_dividing_four_using_arithmetic_progressions(lat_manager.nf, lat_manager.num_chars, lat_manager.smoothness_bound, pr)
      print("init inerts for pr={} in {}".format(pr, time.process_time() - prs_start_time))

      char_mat_start_time = time.process_time()
      lat_manager.character_mats_dict[j][pr] = char_mat(lat_manager.ug_elems_mats[j], lat_manager.nf, prs_for_chars, pr)
      lat_manager.inert_prs_for_chars[j][pr] = prs_for_chars
      timers.linalg_saturate_setup_prs_and_character_mat += time.process_time() - prs_start_time

    print("Starting UG saturation for mat {}".format(j))
    max_num_saturations=4
    linalg_saturate_individual_mats_start=time.process_time()
    if reg_prs_that_need_saturation != []:
      for i in range(max_num_saturations):
        saturating_primes_start = time.process_time()
        new_reg_prs_that_need_saturation = []
        print("Primes that need saturation {}".format(reg_prs_that_need_saturation))
        for pr in reg_prs_that_need_saturation:
          old_reg = lat_manager.det('ug', real_field, j)

          try:
            lat_manager.saturate_ug(pr, unit_rank, j)
            reg_candidate = lat_manager.det('ug', real_field, j)
            new_reg_prs_that_need_saturation.append(pr)
          except ArithmeticError:
            for i in range(max_num_times_double_prec):
              prec = floor(multiple*prec)
              try:
                reg_candidate, real_field = raise_precision(prec, lat_manager)
                #TODO for now assume saturation always reduces the regulator
                if reg_candidate < old_reg:
                  new_reg_prs_that_need_saturation.append(pr)
                break
              except ArithmeticError:
                continue
        print("For ug mat {} saturating pr={} for the {} time took {:.2f} and gave R={}".format(j, pr, i, time.process_time() - saturating_primes_start, reg_candidate.n(10)))
        if new_reg_prs_that_need_saturation == []:
          break
        reg_prs_that_need_saturation = new_reg_prs_that_need_saturation
    timers.linalg_saturate_individual_mats += time.process_time() - linalg_saturate_individual_mats_start

  print("Finished saturating unit group")
  reg_gcd = rgcd(lat_manager.regulator(real_field, 0), lat_manager.regulator(real_field, 1), lat_manager.regulator_lower_bound)
  intersection_reg_start = time.process_time()
  lat_manager.insersection_of_unit_submodules()
  combined_reg = lat_manager.regulator(real_field, -1)
  timers.linalg_saturate_intersection_reg = time.process_time() - intersection_reg_start
  print("Before Saturation: rgcd of the regulators {}".format(initial_reg_gcd.n(15)))
  print("After Saturation: rgcd of regulators {} combined reg {}".format(reg_gcd.n(15), combined_reg.n(15)))

  timers.linalg_step3_ug_saturate = time.process_time() - linalg_step3_ug_saturate_start

  # Complete the class group
  print("Step 4: Complete the class group")
  linalg_step4_cg_complete=time.process_time()
  # TODO update starting_row_ind
  candidate_class_num = lat_manager.det('cg', None, None)
  full_multiple = get_full_multiplier(candidate_class_num, combined_reg, full_det_approx)
  if full_multiple == 1:
    return finish_up(lat_manager, log_level, real_field, start_linalg_time, submat_write_rel_dir, number_of_rows_used_in_ug, cg_num_initial_relns)

  print("Extending with {} new relations. This requires calling add_rows")
  for i in range(new_relns.nrows()):
    if full_multiple == 1:
      break
    # lat_manager.extend_relns_and_compute_basis(new_relns[i,:], new_elems[i,:], new_rat_relns[i,:])
    lat_manager.extend_relns_and_compute_basis(new_relns[i, :], new_elems[i, :])
    candidate_class_num = lat_manager.det('cg', None, None)
    full_multiple = get_full_multiplier(candidate_class_num, combined_reg, full_det_approx)
    print(candidate_class_num, full_multiple)
  else:
    raise ValueError("Failed to complete UG or CG. Try a new run.")
  num_cg_relns_used += i

  timers.linalg_step4_cg_complete = time.process_time() - linalg_step4_cg_complete
  return finish_up(lat_manager, log_level, real_field, start_linalg_time, submat_write_rel_dir, number_of_rows_used_in_ug, num_cg_relns_used)

def finish_up(lat_manager, log_level, real_field, start_linalg_time, submat_write_rel_dir, number_of_rows_used_in_ug, number_of_rows_used_in_cg):
  print("Number of rows used to construct the units submodule is {}".format(number_of_rows_used_in_ug))
  print("Number of rows used to construct the class group submodule is {}".format(number_of_rows_used_in_cg))
  print("Finished adding rows. Starting smith form".format(type(lat_manager.hnf_relns_mat)))
  elem_divs = list(filter(lambda _: _ != 1, lat_manager.hnf_relns_mat.elementary_divisors()))
  print("class structure {}".format(elem_divs))
  # assert_elem_divs(sqr_submat, submat, nrows, elem_divs, inds_of_added_rows)
  # any line after we stop the timer is just for debugging purposes so it doesn't need to be timed
  timers.linalg_total_time = time.process_time() - start_linalg_time
  write_submat(submat_write_rel_dir, lat_manager.hnf_relns_mat, log_level)
  return [lat_manager.det('cg', None, None), elem_divs, lat_manager.det('ug', real_field, -1)]

def initial_hnf(mat):
  linalg_hnf_start = time.process_time()
  print("Linalg - Starting HNF. Mat density {}")

  if src.constants.is_magma_available:
    print("Starting HNF with Magma")
    hnf = non_modular_hnf(0, mat, algo='magma', transformation=False)
  else:
    # when transformation is needed Sagemath's implementation is broken https://github.com/sagemath/sage/issues/33418. Probably Flint 2.9 fixes it but Sagemath wont be upgrading for a long time https://github.com/sagemath/sage/issues/34102
    hnf = non_modular_hnf(0, mat, algo='flint', transformation=False)
  print("Finished initial hnf. Starting add more relations until the unit group has full rank.")
  timers.linalg_step0_cg_initial_relns = time.process_time() - linalg_hnf_start

  return hnf

def raise_precision(new_prec, lat_manager, ug_mat_index):
  real_field = sage.rings.real_mpfr.RealField(new_prec)
  cand_rglr = lat_manager.det('ug', real_field, ug_mat_index)
  # print("Precision increased to {} gave new cand_rglr {}".format(new_prec, cand_rglr))
  print("Precision increased to {}".format(new_prec))
  return cand_rglr, real_field

def write_submat(submat_write_rel_dir, relns_submat, log_level):
  if submat_write_rel_dir != None:
    submat_path = "{}submat_used_for_hnf_{}x{}.bin".format(submat_write_rel_dir, relns_submat.nrows(), relns_submat.ncols())
    print("Writing the submat to {}".format(submat_path))
    write_int_mat(submat_path, relns_submat, log_level)

def assert_elem_divs(sqr_submat, submat, nrows, elem_divs, rows_inds_used):
  if not asserts_on():
    return False
  print("Warning: time consuming function call - {}".format("assert_elem_divs"))
  nrows = submat.nrows()
  available_row_inds = [x for x in list(range(nrows)) if x not in rows_inds_used]
  assert sqr_submat.ncols() == submat.ncols()
  if asserts_on():
    nrows_to_add = floor(0.2*nrows) # usually adding this number of rows is sufficient to recover the full class group
    row_inds_to_add = get_sample_and_increment_seed(available_row_inds, nrows_to_add)
    submat_easy = sqr_submat.stack(submat.matrix_from_rows(row_inds_to_add))
    assert sqr_submat.ncols() == submat_easy.ncols()
    herm_easy = submat_easy.hermite_form(include_zero_rows=False, transformation=False)
    elem_divs_easy = list(filter(lambda _: _ != 1, herm_easy.elementary_divisors()))
    assert elem_divs == elem_divs_easy, "actual elem divs {} vs easy {}. Easy constructed by adding {} rows".format(elem_divs, elem_divs_easy, nrows_to_add)
  return True

def modular_hnf(log_level, mat, modulus, algo):
  '''
  :param log_level:
  :param mat:
  :param modulus:
  :return:
  '''
  if algo != src.constants.hnf_algo_sage and src.constants.hnf_algo_gp:
    raise ValueError("algo must be sage or gp")
  if algo == src.constants.hnf_algo_sage:
    print("Warning: modular HNF in Sage requested but for now we just do a normal HNF in this case")
    return mat.hermite_form(include_zero_rows=False)
  elif algo == src.constants.hnf_algo_gp and src.constants.hnf_algo_gp:
    gp_mat = gp(mat)
    return gp_mat.hnfmatmod(modulus).sage()
  else:
    raise ValueError("algorithm not supported. Check if constants.py has correct settings for the algorithm requested")

def non_modular_hnf(log_level, mat, algo, transformation):
  # https://doc.sagemath.org/html/en/reference/matrices/sage/matrix/matrix_integer_dense_hnf.html
  # all I know is that in some tests using code of the following form was slower
  # import sage.matrix.matrix_integer_dense_hnf as matrix_integer_dense_hnf
  # return matrix_integer_dense_hnf.hnf_with_transformation(mat, proof=False)

  # Sage's hermite_form accepts the following arguments for `algorithm` and `proof`:
  # algorithm – String. The algorithm to use. Valid options are:
  # * 'default' – Let Sage pick an algorithm (default). Up to 75 rows or columns with no transformation matrix, use pari with flag 0; otherwise, use flint.
  # * 'flint' - use flint
  # * 'ntl' - use NTL (only works for square matrices of full rank!)
  # * 'padic' - an asymptotically fast p-adic modular algorithm, If your matrix has large coefficients and is small, you may also want to try this.
  # * 'pari' - use PARI with flag 1
  # * 'pari0' - use PARI with flag 0
  # * 'pari1' - use PARI with flag 1
  # * 'pari4' - use PARI with flag 4 (use heuristic LLL)
  # proof - (default: True); if proof=False certain determinants are computed using a randomized hybrid p-adic multimodular strategy until it stabilizes twice (instead of up to the Hadamard bound). It is incredibly unlikely that one would ever get an incorrect result with proof=False.
  ncols=mat.ncols()
  if algo == src.constants.hnf_algo_magma and src.constants.is_magma_available:
    return magma(mat).HermiteForm().sage()
  else:
    print("HNF starting, no transformation {} {}".format(mat.nrows(), mat.ncols()))
    if transformation:
      hnf = mat.hermite_form(include_zero_rows=True, transformation=True, algorithm=algo, proof=False)
      trans = hnf * mat.inverse()
      #This assert will take forever for matrices of dimension larger that 100
      assert trans == mat.hermite_form(include_zero_rows=True, transformation=True, algorithm=algo, proof=False)[1]
      return hnf, trans
    else:
      return mat.hermite_form(include_zero_rows=True, transformation=False, algorithm=algo, proof=False)

def print_dets(lat_manager, real_field):
  '''
  Prints the 2 determinants to 50 bits of precision for the class group and 10 bits for the unit group.
  The unit group determinant is of the intersection of the 2 units submodules
  :param lat_manager:
  :param real_field:
  :return:
  '''
  print("class det {} reg det {}".format(lat_manager.det('cg', None, None).n(50), lat_manager.det('ug', real_field,
                                                                                                  -1).n(10)))

