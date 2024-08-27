import time
from functools import *
from sage.all import *
from src.obj_string_converters import *
import src.timers as timers

def naive_valn(pr_idl, alpha, nrm_alpha):
  '''
  :param pr_idl:
  :param alpha:
  :param nrm_alpha: the absolute value of the norm of alpha
  :return:
  '''
  nrm_alpha = nrm_alpha.abs()
  assert nrm_alpha == alpha.norm().abs()

  nf = pr_idl.number_field()
  ff = pr_idl.residue_class_degree()
  pr = pr_idl.smallest_integer()
  assert pr > 1

  # valuation() throws ValueError if pr is not prime
  nrm_valn = nrm_alpha.valuation(pr)

  # The ceil expression is also a term in Belabas Algorithm 5.14 line 2. Actually, he gives a tighter bound but it requires computing the ideal's smallest integer which I don't want to do.
  max_valn = max(1, ceil(nrm_valn / ff))
  valn = 0
  gen_elem = pr_idl.gens_two()[1]
  ee = pr_idl.ramification_index()
  for k in range(max_valn+2):
    # "reduce(f) Return the canonical reduction of the element of modulo the ideal self"
    #TODO this is the line that changes in the non-naive version of the algorithm
    k0 = ceil(QQ(k)/ee)
    assert nf.ideal(pr**k0, gen_elem**k) == pr_idl**k # Belabas 5.4.2
    red = (pr_idl**k).reduce(alpha)
    if red != 0:
      break
    valn += 1
  valn -= 1 # the last power so that the remainder was 0
  if valn==(max_valn+1):
    raise ValueError("Failed to compute valuation of idl with factorization {} at the prime {}. max_valn={}".format(nf.ideal(alpha).factor(), pr_idl.factor(), max_valn))
  return valn

def two_elem_valn(pr_idl, alpha, nrm_alpha):
  '''
  :param pr_idl:
  :param alpha:
  :param nrm_alpha: the absolute value of the norm of alpha
  :return:
  '''
  nrm_alpha = nrm_alpha.abs()
  assert nrm_alpha == alpha.norm().abs()
  ee = pr_idl.ramification_index()
  if ee > 1:
    raise ValueError("Only the unramified case is currently supported")

  nf = pr_idl.number_field()
  ff = pr_idl.residue_class_degree()
  pr = pr_idl.smallest_integer()
  assert pr > 1

  # valuation() throws ValueError if pr is not prime
  nrm_valn = nrm_alpha.valuation(pr)

  # TODO I can see why the following max_valn holds when <elem> is a prime ideal b/c N(<elem>) = pr^ff. I haven't proved the general case
  # The ceil expression is also a term in Belabas Algorithm 5.14 line 2. Actually, he gives a tighter bound but it requires computing the ideal's smallest integer which I don't want to do.
  max_valn = max(1, ceil(nrm_valn / ff))
  valn = 1
  gen_elem = pr_idl.gens_two()[1]
  for k in range(1,max_valn+2):
    # "reduce(f) Return the canonical reduction of the element of modulo the ideal self"
    #TODO this is the line that changes in the non-naive version of the algorithm
    red = (pr_idl**k).reduce(alpha)

    k0 = ceil(QQ(k)/ee)
    assert nf.ideal(pr**k0, gen_elem**k) == pr_idl**k # Belabas 5.4.2
    int_mod = IntegerModRing(pr**k0)
    int_mod_polys = PolynomialRing(int_mod, names='Z')
    coeffs= list(gen_elem**k)
    mod_poly = int_mod_polys(coeffs)
    if mod_poly == 0:
      raise ValueError("Failed to compute valuation. Check ramification.")

    if mod_poly.degree() < 0:
      rem = 0 #TODO this is just a guess
    else:

      _, rem = int_mod_polys(list(alpha)).quo_rem(mod_poly)
      rem_reduc = int_mod_polys(rem.coefficients(sparse=False))
    # rem should match this
    print([rem_reduc, rem, mod_poly, gen_elem**k, pr**k0, alpha])

    if rem_reduc != 0:
      break
    valn += 1
  valn -= 1 # the last power so that the remainder was 0
  if valn==(max_valn+1):
    raise ValueError("Failed to compute valuation of idl with factorization {} at the prime {}. max_valn={}".format(nf.ideal(alpha).factor(), pr_idl.factor(), max_valn))
  return valn

@cache
def get_tau0(uniformizer, nf, pr):
  pr = ZZ(pr)
  # defined in Belabas prop 5.9
  dg = nf.degree()
  Y = nf.gen()
  mat = Matrix(GF(pr), dg, dg)
  mat_no_mul = Matrix(QQ, dg, dg)
  # max_ord = nf.maximal_order()
  int_basis = [nf(_) for _ in nf.pari_zk()]
  precomputed_int_bas_inv = Matrix(QQ, [list(_) for _ in int_basis])**(-1)

  for i in range(dg):
    coords_no_mult = vector(QQ, list(int_basis[i] * uniformizer))
    mat_no_mul[i] = coords_no_mult

    coords2 = vector(QQ, list(int_basis[i] * uniformizer)) * precomputed_int_bas_inv
    coords2 = coords2.change_ring(ZZ)
    mat[i] = coords2

  # kernel means left_kernel
  arbitrary_kern_vec = mat.kernel().basis()[0]
  assert arbitrary_kern_vec * mat == zero_vector(dg)
  tau0_coeffs = list(arbitrary_kern_vec)
  tau0_coeffs += [0]*(dg - len(tau0_coeffs))
  tau0 = nf(0)
  for i in range(dg):
    tau0 += ZZ(tau0_coeffs[i]) * int_basis[i]
  mod_int = IntegerModRing(pr)
  assert(all([mod_int(_)==0 for _ in list(tau0 * uniformizer)]))
  return tau0

@cache
def get_tau0_using_uniformizer(pr, elem_gen):
  '''
  Finds an element tau0 of the form in Belabas Prop 5.9
  :param pr:
  :param elem_gen:
  :return:
  '''
  pr = ZZ(pr)
  nf = elem_gen.parent()
  ww = nf.integral_basis()
  pr_idl = nf.ideal(pr, elem_gen)
  ff = pr_idl.residue_class_degree()
  ee = pr_idl.ramification_index()

  if not elem_gen.valuation(pr_idl)==1:
    # Step: find uniformizer. Belabas 6.1.1
    if ee == 1:
      # for degree 4 fields the expected number of tries is <= 1/(1-1/p)^4. Worst case is a completely splitting pr=2 and the expected number of tries is 16
      max_tries = 1000
      for _ in range(max_tries+1):
        rand_elem = pr_idl.random_element()
        # rand_elem = rand_idl_elem(pr, elem_gen, ww)
        if rand_elem.norm().abs().valuation(pr) < ff + 1:
          break
      fail_str = "failed to find uniformizer for pr={}. Valuation at pr is {}".format(pr, rand_elem.norm().abs().valuation(pr))
      if not (rand_elem.norm().abs().valuation(pr) < ff + 1):
        raise ValueError(fail_str)
      if _ == max_tries+1:
        raise ValueError(fail_str)
      uniformizer = rand_elem
    else:
      gen = pr_idl.gens_two()[1]
      if gen.valuation(pr_idl)==1:
        uniformizer=gen
      else:
        uniformizer=gen+pr
  else:
    uniformizer = elem_gen
  assert nf.ideal(pr, uniformizer) == pr_idl
  assert uniformizer.valuation(pr_idl)==1

  tau0 = get_tau0(uniformizer, nf, pr)

  # Check anti-uniformizer satisfies the Belabas definition
  tau = tau0/pr
  assert tau.valuation(pr_idl) == -uniformizer.valuation(pr_idl)

  return [tau0, uniformizer]

def belabas_valn_no_tau0(pr_idl, alpha, nrm_alpha):
  elem_gen = pr_idl.gens_two()[1]
  pr = pr_idl.smallest_integer()
  [tau0, uniformizer] = get_tau0_using_uniformizer(pr, elem_gen)
  nf = pr_idl.number_field()
  return belabas_valn(pr_idl, alpha, nrm_alpha, tau0)

def belabas_valn(pr_idl, alpha, nrm_alpha, tau0):
  #TODO next step is to turn pr_idl into 2 inputs (pr, and gen) and to change get_tau so that its randomized
  #TODO after that I want to write a specialized version of this function for the Stickel rep. In that rep the degree is given by the degree of the 2nd generator
  '''
  Implements Belabas 5.11
  :param pr_idl:
  :param alpha:
  :param nrm_alpha: the absolute value of the norm of alpha
  :return:
  '''
  pr_idl.number_field() # this will raise an Error if a rational prime is passed in
  if tau0 == None:
    raise ValueError("tau0 must be set")
  nrm_alpha = nrm_alpha.abs()
  assert nrm_alpha == alpha.norm().abs()

  valn_setup_start = time.process_time()
  nf = pr_idl.number_field()
  ff = pr_idl.residue_class_degree()
  pr = pr_idl.smallest_integer() #TODO is there a faster way ?
  Fp = GF(pr)

  nrm_valn = nrm_alpha.valuation(pr)
  if nrm_valn == 0:
    return 0
  timers.filter_valn_setup += time.process_time() - valn_setup_start
  if pr == 1:
    raise ValueError("input must be prime")
  if not pr.is_prime():
    raise ValueError("pr_idl must be prime")

  # TODO I can see why the following max_valn holds when <elem> is a prime ideal b/c N(<elem>) = pr^ff. I haven't proved the general case
  # Source: 5.12 comment 2 in Belabas
  max_valn = max(1, ceil(nrm_valn / ff))
  valn = 0
  yy = alpha

  int_basis = [nf(_) for _ in nf.pari_zk()]
  precomputed_int_bas_inv = Matrix(QQ, [list(_) for _ in int_basis])**(-1)
  for k in range(0,max_valn+1):
    yy_prime = yy * tau0
    # max_ord = nf.maximal_order()
    # coords = max_ord.coordinates(yy_prime)
    coords2 = vector(QQ, list(yy_prime)) * precomputed_int_bas_inv

    for cc in coords2:
      if Fp(cc)!=0:
        break
    else:
      assert all([Fp(_) == 0 for _ in coords2])
      yy = yy_prime/pr
      valn += 1
      continue
    break
  if valn==(max_valn+1):
    raise ValueError("Failed to compute valuation of idl with factorization {} at the prime {}. max_valn={}".format(nf.ideal(alpha).factor(), pr_idl.factor(), max_valn))
  return valn
