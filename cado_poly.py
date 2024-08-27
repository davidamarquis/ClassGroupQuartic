from sage.all import *
from src.utils import *
from src.constants import *
from sage.rings.integer import GCD_list
from src.io_utils import *

class CadoQuartic:
  def __init__(self, pol):
    if pol.degree() not in [2,4]:
      raise ValueError("CadoPoly degree must be 2 or 4 (2 is for quadratic subfield)")
    if pol.degree() == 2:
      print("WARNING: Constructing a CadoPoly with deg(f)=2")
    self.poly = pol
    self.nf=None

  def degree(self):
    return self.poly.degree()

  def print_timing_results(self):
    for i in range(len(self.timing_results_inputs)):
      print(self.timing_results_inputs[i])
      print(self.timing_results_outputs[i])

  def skew(self):
    coeffs = self.poly.coefficients(sparse=False)
    dg_inv = 1/float(self.poly.degree())
    return (abs(coeffs[0])/abs(coeffs[-1]))**dg_inv

  def coefficients(self):
    return self.poly.coefficients(sparse=False)

  def get_nf(self):
    if self.nf==None:
      self.nf = NumberField(self.poly, names='WW')
    return self.nf

def read_cado(path_to_params_dir, params_spec_str=None, verbose=False, dg=4):
  '''
  Returns the path to the polynomial file
  :param path_to_params_dir:
  :param params_spec_str:
  :param verbose:
  :return:
  '''
  poly_path = get_poly_file(path_to_params_dir, True, params_spec_str, verbose)
  coeffs = [None]*(dg+1)
  polys = PolynomialRing(QQ, names='Z')

  with open(poly_path, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      if line.startswith("c0:"):
        coeffs[0] = int(line.split("c0:")[1].rstrip())
      if line.startswith("c1:"):
        coeffs[1] = int(line.split("c1:")[1].rstrip())
      if line.startswith("c2:"):
        coeffs[2] = int(line.split("c2:")[1].rstrip())
      if line.startswith("c3:"):
        coeffs[3] = int(line.split("c3:")[1].rstrip())
      if line.startswith("c4:"):
        coeffs[4] = int(line.split("c4:")[1].rstrip())
  pol = polys(coeffs)
  if pol.degree() < dg:
    raise ValueError('not the right degree of the polynomial')
  return pol

def write_cado(pol, path, is_SI_poly, idl_gens, fb_exp_vec, log_output_bitlength=None, disc_div=None):
  '''
  :param pol: CadoPoly
  :param path:
  :param is_SI_poly:
  :param idl_gens:
  :param fb_exp_vec:
  :param log_output_bitlength:
  :param disc_div: a divisor of the polynomial's discriminant.
  :return:
  '''

  poly_name = os.path.split(path)[1]

  if pol.poly.degree() == 2:
    if disc_div == None:
      print("Warning time consuming operation: discriminant of number field")
      NN = pol.get_nf().discriminant().abs()
    else:
      NN = disc_div
  else:
    if disc_div==None:
      print("Warning time consuming operation: discriminant of number field")
      NN = pol.get_nf().discriminant().abs()
    else:
      NN = disc_div

  if gcd(NN, (2*3*5*7*11*13)**4) > 1:
    NN = NN/gcd(NN, (2*3*5*7*11*13)**4)

  coeffs = pol.coefficients()
  com = GCD_list(coeffs)
  str = ""
  str += "n: {}\n".format(ZZ(NN).abs())
  if len(coeffs)==(4+1):
    str += "c4: {}\n".format(coeffs[4]/com)
    str += "c3: {}\n".format(coeffs[3]/com)
    str += "c2: {}\n".format(coeffs[2]/com)
    str += "c1: {}\n".format(coeffs[1]/com)
    str += "c0: {}\n".format(coeffs[0]/com)
    str += "Y4: {}\n".format(coeffs[4]/com)
    str += "Y3: {}\n".format(coeffs[3]/com)
    str += "Y2: {}\n".format(coeffs[2]/com)
    str += "Y1: {}\n".format(coeffs[1]/com)
    str += "Y0: {}\n".format(coeffs[0]/com)
  elif len(coeffs)==(2+1):
    str += "c2: {}\n".format(coeffs[2]/com)
    str += "c1: {}\n".format(coeffs[1]/com)
    str += "c0: {}\n".format(coeffs[0]/com)
    str += "Y2: {}\n".format(coeffs[2]/com)
    str += "Y1: {}\n".format(coeffs[1]/com)
    str += "Y0: {}\n".format(coeffs[0]/com)

  str += "skew: {}\n".format(pol.skew())

  str += "{}: {}\n".format(idl_coeff_gcd_str, com)
  str += "name: {}\n".format(poly_name)
  if is_SI_poly:
    str += "{}: {}\n".format(idl_least_int_str, idl_gens[0])
    str += "{}: {}\n".format(idl_poly_gen_str, idl_gens[1])
    str += "{}: {}\n".format(idl_fb_exp_vec_str, fb_exp_vec)

  with open(path, "w+") as text_file:
    text_file.write(str)