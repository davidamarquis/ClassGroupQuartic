import time
from functools import *
from sage.all import *
import src.timers as timers
#Magma timing: all functions used by the algorithm that call Magma code should be in this file.
#They should update the total cpu time of the algorithm by adding the decorator @add_magma_time_to_total
def GRH_bound_command(pol):
  '''
  :param pol:
  :param Magma_cg_proof_str: Magma input for the class group algorithm usually 'GRH'
  :param Magma_include_units:
  :param Magma_verbose_level:
  :return:
  '''
  str="R<X> := PolynomialRing(Integers());O := MaximalOrder({}); print GRHBound(O);".format(pol)
  return str

def cg_command(pol, Magma_cg_proof_str='"GRH"', Magma_verbose_level=1):
  '''
  :param pol:
  :param Magma_cg_proof_str: Magma input for the class group algorithm usually 'GRH'
  :param Magma_include_units:
  :param Magma_verbose_level:
  :return:
  '''
  if Magma_verbose_level not in [1, 2, 3]:
    raise ValueError("bad value for verbose")
  str='SetVerbose("ClassGroup", {});'.format(Magma_verbose_level) # 1 for minimal, 2 for some data, 3 for verbose
  # str+='SetClassGroupBounds("GRH");'
  str+='tt := Cputime();'
  # next line will be like 'R<X> := PolynomialRing(Integers());O := MaXimalOrder(X**4 ...);C, m := ClassGroup(O);
  str+="R<X> := PolynomialRing(Integers());O := MaximalOrder({});C, m := ClassGroup(O: Proof:={});".format(pol, Magma_cg_proof_str)
  # str+="IndependentUnits(O);" #Don't compare with units
  # str+="Regulator(O: Current:=true);"

  str+="print C;print Cputime(tt);"
  return str

def add_magma_time_to_total(fn):
  # idea for this is from High Performance Python, 2nd Edition Ex2.5
  @wraps(fn)
  def measure_time(*args, **kwargs):
    t1 = magma.cputime()
    result = fn(*args, **kwargs)
    timers.total_time_cpu += magma.cputime(t1)
    return result
  return measure_time

@add_magma_time_to_total
def magma_mat_from_sparse_sage_mat(sp_reln_mat):
  if sp_reln_mat.is_dense():
    raise ValueError("matrix should be sparse but it is dense")
  mag_convert_start = magma.cputime()
  mag_relns = magma(sp_reln_mat)  # Magma type is MtrxSprs
  timers.filter_convert_time += magma.cputime(mag_convert_start)
  return mag_relns

@add_magma_time_to_total
def compute_class_num(mag_relns):
  from src.utils import prod_all_elems
  sage_class_num_start=time.process_time()
  magma_elem_divs_start=magma.cputime()
  elem_divs = mag_relns.ElementaryDivisors().sage()
  cand_class_num_before = prod_all_elems([_ for _ in elem_divs if _ != 1])
  timers.filter_class_num=time.process_time()-sage_class_num_start+magma.cputime(magma_elem_divs_start)
  return cand_class_num_before

def gelin_poly_and_elem_from_poly(poly, gel_c):
  nf=NumberField(poly, names='Z')
  return gelin_poly_and_elem(nf, gel_c)

def gelin_poly_and_elem(nf, gel_c):
  pol = nf.polynomial()
  pol_ring = pol.parent()
  mag_pol = magma(pol)
  # I could not get Magma to return a list so I could not return both the polynomial and the element. So I modified it to only return the element's coordinates on the power basis
  # [gelin_red_pol, coeffs] = mag_pol.PolRed(10,Var=True)
  elem_coordinates = mag_pol.PolRed(gel_c, 2**6)
  elem_coordinates = list(elem_coordinates.sage())
  elem = nf(elem_coordinates)
  pol = elem.minpoly()
  output_pol = pol_ring(pol.coefficients(sparse=False))
  return output_pol, elem

def load_all_mag_files():
  magma.load('/files3/home/david.marquis/pycharm_project_468/src/MagmaUtils.m')

@cache
def is_running_on_remote():
  try:
    load_all_mag_files()
    return True
  except RuntimeError:
    return False

def error_if_remote():
  if is_running_on_remote():
    raise ValueError("must be run on local")

def error_if_local():
  if not is_running_on_remote():
    raise ValueError("must be run on remote")

def probable_pivot_cols_mag_echelon_form(mat):
  '''
  Deprecated. Converting a matrix with a dimension 100 or more from Magma to Sage takes forever
  :param mat:
  :return:
  '''
  pr=ZZ.random_element(10007, 46000).next_prime()
  # if gcd(pr, cand_class_num)>1:
  #   raise ValueError("can't compute pivots b/c pr={} divides mat's det".format(pr))
  nrows = mat.NumberOfRows().sage()
  ncols = mat.NumberOfColumns().sage()
  mod_mats = magma.KMatrixSpace(GF(pr), nrows, ncols)
  dense_mat=mat.Matrix()
  mod_mat = mod_mats(dense_mat)
  ech_form=mod_mat.EchelonForm()

  convert_to_sage_start=time.process_time()
  rowspace_mat = ech_form.sage()
  print(time.process_time()-convert_to_sage_start)
  pivs=[]
  for rr in range(nrows):
    for cc in range(rr,ncols):
      if rowspace_mat[rr,cc] != 0:
        pivs.append(cc)
        break
  return pivs

@add_magma_time_to_total
def magma_probable_pivot_cols(mat):
  '''
  Takes a Magma mat as input and returns its pivot columns
  :param mat:
  :return:
  '''
  pr=ZZ.random_element(10007, 46000).next_prime()
  pivs=mat.ProbablePivotCols(pr) #this is in the file MagmaUtils.m
  return [_-1 for _ in pivs.sage()]

#TODO move to test
def validate_sparse_mat_pivs():
  mat = Matrix(ZZ, [[1,1,0],[0,0,1]])
  mag_mat = magma(mat).SparseMatrix()
  mag_pivs = magma_sparse_mat_pivots(mag_mat)
  assert [1,2] in mag_pivs
  assert [0,0] in mag_pivs
  mag_piv_col_inds = mag_pivot_col_inds(mag_mat)
  assert mag_piv_col_inds == mat.pivots()

  mat2 = random_matrix(ZZ, 8, 8, algorithm='echelon_form', num_pivots=3)
  mag_mat2 = magma(mat2).SparseMatrix()
  mag_pivs2 = mag_pivot_col_inds(mag_mat2)
  assert mag_pivs2 == mat2.pivots()
  assert len(mag_pivs2)==3

is_running_on_remote() # loads all Mag files

