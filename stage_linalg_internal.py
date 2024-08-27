from sage.all import *
from random import sample
from random import randint # note that sage has a function that shadows randint (and doesn't work with random.seed). Be careful changing the order of these imports
from src.stage_init_internal import fac_base_data
from sage.groups.generic import discrete_log

from src.utils import *
from sage.matrix.matrix_integer_dense_hnf import det_padic
from sage.matrix.matrix_integer_dense_hnf import probable_pivot_columns
from sage.matrix.matrix_integer_dense_saturation import random_sublist_of_size
import src.timers as timers
from src.determ_rng import *
from src.utils import *
from sage.misc.randstate import random
import time

def left_kernel_via_hnf(submodule_of_intersection):
  hnf, transformation = submodule_of_intersection.hermite_form(algorithm='flint', transformation=True, include_zero_rows=True)
  rnk = hnf.rank() # if rank is full this is not needed. TODO how much does this cost ?
  return transformation[rnk:]

def get_full_multiplier(cand_h, cand_R, full_det_approx):
  full_multiplier_approx = (cand_h * cand_R) / full_det_approx
  if (full_multiplier_approx - round(full_multiplier_approx)).abs() < 0.1:
    return ZZ(round(full_multiplier_approx))
  else:
    response_str = "Could not compute full multiplier. class_candidate {}, reg={}. Distance from nearest integer {}. Full_det_approx {}. Possibilities\na) Not enough precision in Reg or CN formula\nb) Not enough terms computed in the Euler product".format(cand_h, cand_R, full_multiplier_approx - round(full_multiplier_approx), full_det_approx)
    raise ValueError(response_str)

def submat_full_rank_at_columns(sage_mat, magma_mat):
  '''
  Uses the rows of the matrix to construct a submatrix with rank equal to len(col_inds) so that the columns of the submatrix generate the submatrix of the original matrix spanned by the columns with indices col_inds
  :param sage_mat:
  :param magma_mat:
  :return:
  '''
  if is_running_on_remote() and magma_mat==None:
    raise ValueError("On remote so Magma should be used but magma_mat is None")
  linalg_submat_setup_start = time.process_time()

  if is_running_on_remote():
    pivots = pivot_cols_fastest_available(magma_mat.Transpose())
  else:
    pivots = pivot_cols_fastest_available(sage_mat.transpose())
  full_rank_at_cols_mat = sage_mat.matrix_from_rows(pivots)

  timers.filter_submat_setup += time.process_time()-linalg_submat_setup_start
  return full_rank_at_cols_mat, pivots

def char_mat(elems_submat, nf, character_pr_idls, qq):
  '''
  :param elems_submat:
  :param nf:
  :param character_pr_idls:
  :param qq: the modulus of the Diri character
  :return:
  '''
  if not is_pseudoprime(qq):
    raise ValueError("modulus must be a prime")
  qq = ZZ(qq)
  character_pr_fields = [_.residue_field() for _ in character_pr_idls]
  # character_mat = Matrix(ZZ, nrows=elems_submat.nrows(), ncols=len(character_pr_idls))
  character_mat = Matrix(ZZ, nrows=elems_submat.nrows(), ncols=len(character_pr_idls))

  for cc in range(len(character_pr_idls)):
    idl = character_pr_idls[cc]
    pp = idl.norm()
    residue_field = character_pr_fields[cc]
    power = (pp - 1) / qq
    number_of_ties_to_find_nonresidue = 10
    for j in range(number_of_ties_to_find_nonresidue):
      rand_elem = residue_field.random_element()
      if rand_elem**power != 1:
        break
    else:
      raise ValueError("did not find an element in the residue field that is a qq-nonresidue. The probability of this happening is less than 2^-{} I think".format(number_of_ties_to_find_nonresidue))
    resroot = rand_elem ** power
    for rr in range(elems_submat.nrows()):
      elem = nf(elems_submat[rr, 0])
      rf_elem = residue_field(elem)
      if rf_elem == 0:
        raise ValueError("elem of the residue field is 0 at pr of norm={}".format(pp))
      x = rf_elem ** power
      # maybe it would be better to allow prs of degree larger than 1 and use fflog to compute the solution
      dlog_start = time.process_time()
      if residue_field.degree() == 1:
        # pr = residue_field.characteristic()
        # lg = gp.znlog(x, gp.Mod(resroot, pr))
        lg = discrete_log(x, resroot, ord=qq)  # dog slow
      else:
        raise ValueError("Chars on residue fields of degree larger than 1 are not supported")
        # print("Warning: Computing dlog of character ideals of degree larger than 1. This is slow for large pr")
        lg = discrete_log(x, resroot, ord=qq)  # dog slow
      # print("dlog took {}".format(dlog_start - time.process_time()))
      character_mat[rr, cc] = lg

  return character_mat

