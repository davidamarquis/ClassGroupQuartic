from sage.all import *
from sage.all import RR
from sage.calculus.var import var
from sage.symbolic.integration.integral import definite_integral

from src.utils import *

def index(poly):
  nf=NumberField(poly, names='Z')
  disc=poly.discriminant()
  Delta=nf.discriminant()
  return (disc/Delta).abs()

def log_index(poly, prec):
  return RR(index(poly)).log(2).n(prec)

def rand_poly(max_coeff, degree, monic=False, constant_coeff_bnd=None):
  '''
  Returns either random poly with bounded coeffs or poly with all intermediate coeffs bounded by one bound and
  the constant term bounded by a different bound
  :param max_coeff:
  :param degree:
  :param monic:
  :param constant_coeff_bnd:
  :return:
  '''
  polys = PolynomialRing(QQ, names='X')
  if monic:
    if constant_coeff_bnd == None:
      constant_coeff_bnd = max_coeff
    return polys([randint(-constant_coeff_bnd, constant_coeff_bnd)] + [randint(-max_coeff, max_coeff) for _ in range(1, degree)] + [1])
  if constant_coeff_bnd != None:
    raise ValueError("larger constant term not supported for non-monics")
  return PolynomialRing(QQ, names='X')([randint(-max_coeff, max_coeff) for _ in range(0, degree+1)])

def rand_monic_eisenstein(max_coeff, degree, pr):
  polys = PolynomialRing(QQ, names='X')
  coeff_bound=floor(max_coeff/pr)
  for i in range(100):
    rand_ints=[randint(-coeff_bound, coeff_bound) for _ in range(0, degree)]
    if rand_ints[0] % pr != 0:
      pol=polys([pr*_ for _ in rand_ints] + [1])
      return pol
  raise ValueError("no poly found")

def rand_reduced_defining_pol(target_Delta, dg, seed, trial_div_bound_exp, is_monic, use_larger_constant):
  '''
  The goal of this function is to construct defining polynomials with bounded coefficients that are "sufficiently random" for the problem of testing the performance of index-calculus algorithms for computing the class group of the number field generated by this defining polynomial.
  One way to do this is to generate random integral polynomials $ff$ so that the ith coefficient c_i satisfies |c_i|<=A for some bound A.
  One problem with this approach is that the height of ff is tiny (compared to, say, the height of the polynomials this function computes)
  There is no known method for taking a general defining polynomial and finding a defining polynomial for an isomorphic number field that has a height as small as this construction.
  However, there _is_ an approach to index-calculus (Biasse & Fieker 2013) whose performance directly depends on the height of the defining polynomial.
  This method can compute class groups of number fields of very large discriminant when they have this special form.
  Thus, it is important to construct defining polynomials that are *more random* than this approach.
  This function does this by using such a polynomial $ff$ as a starting point, then computing a random map phi:F->F where $F=Q[X]/ff$ and computing the numerator of the minimal polynomial of phi(X).
  The resulting polynomial $g$ usually has the same discriminant as ff.
  However, its coefficients are usually large, so we apply the Pari function ```polredabs```.
  For this $g(x) = \prod_i(x-\r_i)$ the $T_2$ norm of $g$ is $T_2(g) = \sum_i |\r_i|^2$.
  The coefficients then satisfy $|c_{n-k}| <= bin(n,k)(T_2(g)/n)^(k/2)$ where ```bin``` is the binomial function.

  For the details of how polredabs computes this polynomial see https://pari.math.u-bordeaux.fr/dochtml/html/General_number_fields.html):
  and https://www.lmfdb.org/knowledge/show/nf.polredabs
  :param Delta: we want a number field so that the (absolute value of the) discriminant is not too far from Delta.
  :param dg: degree of the polynomial
  :param seed:
  :param trial_div_bound_exp: an integer or None if you want to always factor the discriminant. 15 is a reasonable default
  '''
  exit_if_discrim_has_too_many_large_primes= trial_div_bound_exp != None
  if is_monic:
    coeff_bnd = floor(target_Delta ** (1.0 / (2 * dg - 2)))
  if use_larger_constant:
    coeff_bnd = floor(target_Delta ** (1.0 / ((dg-1)*(dg))))
    constant_coeff_bnd = floor(target_Delta ** (1.0 / (dg-1)))
  for i in range(10):
    set_random_seed(seed+i)
    initial_def_poly = rand_poly(coeff_bnd, dg, is_monic, constant_coeff_bnd)
    if not initial_def_poly.degree()==4:
      continue
    if not initial_def_poly.is_irreducible():
      continue
    break
  else:
    raise ValueError("no irreducible quartic found with coeff_bnd {}".format(coeff_bnd))
  disc = ZZ(initial_def_poly.discriminant()).abs()
  if disc==1:
    print(initial_def_poly)
    return None

  if exit_if_discrim_has_too_many_large_primes:
    trial_div_bound = 2**trial_div_bound_exp
    trial_div_facs = factor_trial_division(disc, trial_div_bound)
    trial_div_prs = [_[0] for _ in trial_div_facs]
    disc_lp = trial_div_facs[-1][0]
    if not ZZ(disc_lp).is_prime(proof=False):
      return None
    trial_div_prs += [disc_lp]
    nf = NumberField(initial_def_poly, maximize_at_primes=trial_div_prs, names='t')
  else:
    nf = NumberField(initial_def_poly, names='t')

  # this will factor the polynomial's discrim. Pari docs claim there's a way around it using a list but I couldn't get it working
  reduced_initial_pol = initial_def_poly.parent()([_.sage() for _ in list(pari.polredabs(str(initial_def_poly)))])
  DD = nf.discriminant().abs()
  logD = RR(DD).log(2)
  (r,s) = nf.signature()
  nf_urank=r+s-1
  return reduced_initial_pol, logD, RR(disc).log(2), RR(coeff_bnd).log(2), DD, nf_urank, initial_def_poly

def nonzero_roots(pol):
  roots = [_ for _ in pol.roots() if _[0] != 0]
  return roots

def cohen_size(pol):
  '''
  Cohen's definition of size
  :param pol:
  :return:
  '''
  roots = vector([_[0].abs() for _ in pol.change_ring(CC).roots()])
  roots_sqr_sum = roots.dot_product(roots)
  return roots_sqr_sum**0.5

def log_cohen_size(pol):
  return RR(cohen_size(pol)).log(2)