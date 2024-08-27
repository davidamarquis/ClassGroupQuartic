import json
from sage.all import *


def poly_to_str(poly):
  #TODO write str_to_poly
  # process this poly into a string without '[' and ']' characters
  poly_str = str(poly.coefficients(sparse=False))
  return poly_str[1:len(poly_str)-1]
def poly_to_flint_poly_str(poly):
  coeffs = poly.coefficients(sparse=False)
  poly_str = str(coeffs)
  return str(len(coeffs)) + "  " + poly_str[1:len(poly_str)-1].replace(",", "")

def nf_elem_to_str(elem):
  # process this tuple into a string without '[' and ']' characters
  elem_str = str(list(elem))
  return elem_str[1:len(elem_str)-1]

def nf_elem_to_flint_poly_str(elem):
  coeffs = list(elem)
  elem_str = str(coeffs)
  return str(len(coeffs)) + elem_str[1:len(elem_str)-1].replace(",", "")

def str_to_nf_elem(list_of_strs, nf):
  if len(list_of_strs) == None:
    raise ValueError("input must be a list")
  coeff_strs = [QQ(_) for _ in list_of_strs]
  dg = nf.degree()
  return nf(coeff_strs + [0]*(dg - len(coeff_strs)))

def mat_to_str(mat):
  '''
  Takes a Sage mat as input and constructs a string representation. The string represents each matrix row as a comma-separated list of ints and terminates with a newline character.
  :param mat:
  :return:
  '''
  n_rows = mat.nrows()
  mat_str = ''
  for i in range(n_rows):
    if i!=0:
      mat_str += '\n'
    line_list = list(mat[i])

    # process this tuple into a string without '[' and ']' characters
    line_str = str(line_list)
    mat_str += line_str[1:len(line_str)-1]
  return mat_str

def str_to_mat(mat_str):
  str_rows = mat_str.split("\n")
  rows = [None]*len(str_rows)
  for i in range(0, len(rows)):
    rows[i] = [int(_) for _ in str_rows[i].split(",")]
  return Matrix(ZZ, rows)

def mat_to_flint_mat_str(mat):
  #TODO write flint_mat_str_to_mat function
  '''
  Takes a Sage mat as input and writes each row to a file.
  The string is in flint format: number of rows, a space, number of columns, two spaces, then a space-separated list of coefficients, one row after the other.
  :param mat:
  :return:
  '''
  n_rows = mat.nrows()
  n_cols = mat.ncols()
  mat_str = '{} {} '.format(n_rows, n_cols)
  for i in range(n_rows):
    line_list = list(mat[i])
    line_str = str(line_list)
    mat_str += " "
    mat_str += line_str[1:len(line_str)-1].replace(",", "")
  return mat_str


def padded_coeffs(elem, dg):
  coeffs = list(elem)
  return list(coeffs) + [0]*(dg - len(coeffs))

def mat_to_num_denom(mat):
  dnm = mat.denominator()
  num = dnm * mat
  return [num, dnm]

def str_dict_to_ind_calc_data(dict, nf):
  '''
  :param dict: {"basis":"[[2,3],[4,5]]"}
  :param nf:
  :return:
  '''
  if not "basis" in dict.keys():
    return dict
  list_of_coeff_lists = [_.split(",") for _ in dict["basis"]]
  return [str_to_nf_elem(_, nf) for _ in list_of_coeff_lists]

