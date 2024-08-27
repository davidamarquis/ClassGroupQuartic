import time

from sage.all import *
from src.magma_utils import *
from random import *

def random_permutations(mat, num):
  for i in range(10 * num):
    ri0 = randint(0, num - 1)
    ri1 = randint(0, num - 1)
    ri2 = randint(0, 1)
    sign = (-1) ** ri2
    mat.add_multiple_of_row(ri0, ri1, sign)
  for i in range(10 * num):
    ri0 = randint(0, num - 1)
    ri1 = randint(0, num - 1)
    ri2 = randint(0, 1)
    sign = (-1) ** ri2
    mat.add_multiple_of_column(ri0, ri1, sign)
  return mat

def mat_with_hnf_so_that_many_smooths_on_diagonal(num_rows, num_cols, full_rank_submat, magma_mat_output):
  '''
  We start with a construction of Allan Steel:
  For each n, we start with the n by n zero matrix and then set each diagonal entry to a random integer between 1 and (n/10).
  We then do (10*n) row operations on the matrix by randomly adding or subtracting a random row from another random row and
  (10*n) column operations on the matrix by randomly adding or subtracting a random column from another random column."
  http://magma.maths.usyd.edu.au/users/allan/mat/hermite.html
  Then we take a full rank sqr submatrix of this.
  :param num_rows:
  :param num_cols:
  :param full_rank_submat: bool
  :param magma_mat_output: bool
  :return:
  '''
  # Observations for num in [50, 100]
  # The resulting input matrix tends to have:
  # num rows more than 0.9*num
  # entries less than 6-digits
  # diagonal of the Hermite form has many non-trivial entries all of which factor into primes less than (n/10).
  mat = Matrix(ZZ, num_rows, num_cols, sparse=True)
  seed(0)

  for i in range(num_rows):
    mat[i,i] = randint(1, int(num_rows/10.0))

  random_permutations(mat, num_rows)

  if full_rank_submat:
    pivs = mat.pivots()
    full_rank_mat = mat.matrix_from_columns(pivs)[: len(pivs)]
    if magma_mat_output:
      return magma(full_rank_mat)
    else:
      return full_rank_mat
  else:
    if magma_mat_output:
      return magma(mat)
    else:
      return mat

def rank_of_allan_steel_rand_mat():
  for i in range(100, 103):
    rand_mat = mat_with_hnf_so_that_many_smooths_on_diagonal(i)
    print([i, rand_mat.rank(), [_.factor() for _ in rand_mat.hermite_form().diagonal() if _ != 0]])

def validate_mat_with_hnf_so_that_many_smooths_on_diagonal():

  nrows=100
  mag_mat_convert_start_time = magma.cputime()
  test_mat0_mag = mat_with_hnf_so_that_many_smooths_on_diagonal(nrows, nrows, False, True)
  print(test_mat0_mag.Rank())
  print("magma mat convert{:.2f}".format(magma.cputime(mag_mat_convert_start_time)))

  print("Done making test mat")

  # test_mat0_sage = mat_with_hnf_so_that_many_smooths_on_diagonal(150, False, False)

  start_time = magma.cputime()
  # pivs0 = mag_pivot_col_inds(test_mat0_mag)
  # pivs_mag = magma_probable_pivot_cols(test_mat0_mag)
  pivs_mag_fast = magma_probable_pivot_cols(test_mat0_mag)
  print("Pivots, magma {:.2f}".format(magma.cputime(start_time)))
  # start_time_sage = time.process_time()
  # pivs_sage=test_mat0_sage.pivots()
  # print("Pivots, Sage default {:.2f}".format(time.process_time() - start_time_sage))

validate_mat_with_hnf_so_that_many_smooths_on_diagonal()
