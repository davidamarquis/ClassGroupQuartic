from multiprocessing.dummy import Pool
from subprocess import call

from sage.rings.factorint import factor_trial_division
from sage.matrix.matrix_integer_dense_hnf import probable_pivot_columns
from sage.rings.integer import GCD_list
from src.magma_utils import *
import src.timers as timers

my_env = os.environ
my_env["LD_LIBRARY_PATH"] = ":/files3/home/david.marquis/GCC-12.2.0/lib64"

def set_ld_lib_path():
  if is_running_on_remote():
    return my_env
  return os.environ

def limited_concurrent_commands(commands, pool_size):
  env = set_ld_lib_path()
  pool = Pool(pool_size) # two concurrent commands at a time
  for i, returncode in enumerate(pool.imap(partial(call, shell=True, env=env), commands)):
    if returncode != 0:
      print("%d command failed: %d" % (i, returncode))

def all_indexs(lst, value):
  '''
  Finds all indexs i in the list lst so that lst[i]==value. From https://stackoverflow.com/a/6294205
  :param lst:
  :param value:
  :return:
  '''
  return [i for i, x in enumerate(lst) if x == value]

def is_pure_quartic(poly):
  if poly.degree() != 4:
    return False
  coeffs=poly.coefficients(sparse=False)
  return coeffs[1]==0 and coeffs[2]==0 and coeffs[3]==0

def disc_divisor(poly):
  div = None
  if is_pure_quartic(poly):
    coeffs=poly.coefficients(sparse=False)
    div=coeffs[0]
  return div

def nf_elem_to_poly(elem, pol_ring):
  return pol_ring(list(elem))

def pr_idl_length(pr_idl):
  '''
  A definition of idl length so that the length of all prime ideals above p will be distinct
  :param pr_idl: 
  :return: 
  '''
  Y=pr_idl.number_field().gen()
  second_gen=pr_idl.gens_two()[1]
  return T2_elem_real(second_gen, 10)

def reduced_basis(nf, idl, prec):
  """
  Copied from Sage's function in the number_field class.
  TODO I have another version of this function where I used LLL somewhere. Find that and compare
  Return an LLL-reduced basis for the Minkowski-embedding
  of the ideal

  INPUT:

  -  ``prec`` (default: ``None``) - the precision with which to
     compute the Minkowski embedding.

  OUTPUT:

  An LLL-reduced basis for the Minkowski-embedding of the
  maximal order of a number field, given by a sequence of (integral)
  elements from the field.

  .. NOTE::

      In the non-totally-real case, the LLL routine we call is
      currently PARI's :pari:`qflll`, which works with floating point
      approximations, and so the result is only as good as the
      precision promised by PARI. The matrix returned will always
      be integral; however, it may only be only "almost" LLL-reduced
      when the precision is not sufficiently high.

  EXAMPLES::

      sage: F.<t> = NumberField(x^6-7*x^4-x^3+11*x^2+x-1)
      sage: F.maximal_order().basis()
      [1/2*t^5 + 1/2*t^4 + 1/2*t^2 + 1/2, t, t^2, t^3, t^4, t^5]
      sage: F.reduced_basis()
      [-1, -1/2*t^5 + 1/2*t^4 + 3*t^3 - 3/2*t^2 - 4*t - 1/2, t, 1/2*t^5 + 1/2*t^4 - 4*t^3 - 5/2*t^2 + 7*t + 1/2, 1/2*t^5 - 1/2*t^4 - 2*t^3 + 3/2*t^2 - 1/2, 1/2*t^5 - 1/2*t^4 - 3*t^3 + 5/2*t^2 + 4*t - 5/2]
      sage: CyclotomicField(12).reduced_basis()
      [1, zeta12^2, zeta12, zeta12^3]

  TESTS:

  Check that the bug reported at :trac:`10017` is fixed::

      sage: x = polygen(QQ)
      sage: k.<a> = NumberField(x^6 + 2218926655879913714112*x^4 - 32507675650290949030789018433536*x^3 + 4923635504174417014460581055002374467948544*x^2 - 36066074010564497464129951249279114076897746988630016*x + 264187244046129768986806800244258952598300346857154900812365824)
      sage: new_basis = k.reduced_basis(prec=120)
      sage: [c.minpoly() for c in new_basis]
      [x - 1,
       x^2 - x + 1,
       x^6 + 3*x^5 - 102*x^4 - 103*x^3 + 10572*x^2 - 59919*x + 127657,
       x^6 - 3*x^5 - 102*x^4 + 315*x^3 + 10254*x^2 - 80955*x + 198147,
       x^3 - 171*x + 848,
       x^6 + 171*x^4 + 1696*x^3 + 29241*x^2 + 145008*x + 719104]
      sage: R = k.order(new_basis)
      sage: R.discriminant()==k.discriminant()
      True
  """
  d = nf.absolute_degree()

  if idl == None:
    int_bas = nf.integral_basis()
  else:
    int_bas = idl.integral_basis()

  if nf.is_totally_real():
    from sage.matrix.constructor import matrix
    M = matrix(ZZ, d, d, [[(x*y).trace() for x in int_bas] for y in int_bas])
    T = pari(M).qflllgram()
  else:
    M = nf.minkowski_embedding(int_bas, prec=prec)
    # T is a unimodular transformation matrix so that M*T is an LLL reduced basis
    T = pari(M).qflll()
  if len(int_bas) < d:
    raise ValueError('Somehow Zk is not an integral basis of the ideal')
  if T.ncols() < d:
    raise ValueError('Somehow minkowski embedding matrix is not full rank')
  return [sum([ZZ(T[i][j]) * int_bas[j] for j in range(d)]) for i in range(d)]

def d_log2(xx, prec=15):
  '''
  log2 is a function in the global namespace so I use d_log2
  :param xx:
  :param prec:
  :return:
  '''
  return RR(xx).log(2).n(prec)

def lll_reduced_sieve_poly(idl, prec):
  nf = idl.number_field()
  ww = reduced_basis(nf, idl, prec)
  lll_elem=ww[1]

  nrm=lll_elem.norm()
  if nrm.is_square():
    # an elem of a degree 4 field that has square norm is very likely to be in a quadratic subfield
    lll_elem=ww[2]

  idl_nrm = idl.norm()
  mod_gens = nf(ww[0]), nf(lll_elem)

  norm_form_start = time.process_time()
  allow_degree_two=True
  sieving_poly = idl_nrm**(-1)*norm_form(mod_gens, allow_degree_two)
  # timers.filter_target_new_polys_form += time.process_time()-norm_form_start
  com = GCD_list(sieving_poly.coefficients(sparse=False))
  return sieving_poly/com, com, mod_gens, idl_nrm

def crt_coeffs(prs):
  prod_prs=prod_all_elems(prs)
  return [ZZ((prod_prs/pr)*((prod_prs/pr)**(-1) % pr)) for pr in prs], prod_prs

def crt_elems(elems, prs, nf):
  '''
  :param elems: a list of number field elements
  :param prs: a list of rational primes
  :param nf:
  :return:
  '''
  coeffs, prod_prs =  crt_coeffs(prs)
  big_elem = sum(p*q for p,q in zip(elems, coeffs))
  reduced_elem = reduce_coordinates_of_nf_elem(big_elem, prod_prs)
  return reduced_elem, prod_prs

def reduce_coordinates_of_nf_elem(elem, modulus):
  nf=elem.parent()
  mod_int = IntegerModRing(modulus)
  reduced_elem = nf([ZZ(mod_int(_)) for _ in list(elem)])
  return reduced_elem

def map_pr_idls_to_dedekind_kummer_elements(pr_data, def_poly):
  '''
  Warning: Does not test that the prs don't divide the index. If they divide the index then Dede-kummer does not apply
  Strangely, there's no built-in way to represent prime ideals in terms of the Dedekind-Kummer theorem when the prime does not divide the index.
  :param pr_data: a list of pairs [pr, pr_idl]
  :param def_poly:
  :return:
  '''
  nf=pr_data[0][1].number_field()
  theta = nf.gen()
  gen_elems = [None]*len(pr_data)
  for i in range(len(pr_data)):
    pr = pr_data[i][0]
    pr_idl = pr_data[i][1]
    poly = def_poly.change_ring(GF(pr))
    polys=poly.parent()
    elem_with_dnm_cleared = pr_idl.gens_two()[1]
    elem_with_dnm_cleared *= elem_with_dnm_cleared.denominator()
    sec_gen_poly = polys(list(elem_with_dnm_cleared))

    facs_in_ZZ = [_[0] for _ in poly.factor()]
    for fac in facs_in_ZZ:
      if gcd(sec_gen_poly, fac) > 1:
        try:
          gen_elems[i]=fac.change_ring(nf)(theta)
        except TypeError as e:
          print("failed to convert {} to nf elem".format(fac))
          raise e
        break
    if gen_elems[i]==None:
      raise ValueError("did not find elem that is 0 modulo the pr_idl at i={}. Check if this pr divides the index".format(i))
  return gen_elems

def crt_sieve_poly(pr_idls, second_gen_residues):
  nf = pr_idls[0].number_field()
  pr_idl_norms=[_.norm() for _ in pr_idls]

  crt_elem, prod_prs = crt_elems(second_gen_residues, pr_idl_norms, nf)
  mod_gens=[nf(prod_prs), crt_elem]
  # print("crt_sieve_poly: 2 elements of the ideal {}".format(mod_gens))

  # does not work:
  # assert nf.ideal(mod_gens).mod(prod_all_elems(pr_idls))==0, "{} is not a multiple of {}".format(nf.ideal(mod_gens).factor(), prod_all_elems(pr_idls).factor())

  sieving_poly = prod_prs**(-1) * norm_form(mod_gens, False)
  com = GCD_list(sieving_poly.coefficients(sparse=False))
  return sieving_poly/com, com, mod_gens


def crt_elem_from_idealchinese(nf, pr_idl_norms, pr_idls, prs):
  norm = prod_all_elems(pr_idl_norms)
  # you would think idl.gens_two()[1] would be the element defined by Dedekind kummer but it is not. Have to compute it manually
  second_gen_residues = map_pr_idls_to_dedekind_kummer_elements(list(zip(pr_idl_norms, pr_idls)), nf.polynomial())
  crt_elem = nf.idealchinese([nf.ideal(pr) for pr in prs], second_gen_residues)  # wrong, returns 0
  dnm = crt_elem.denominator()
  mod_int = IntegerModRing(norm)
  crt_elem *= dnm*ZZ(mod_int(dnm**(-1)))
  reduced_elem = nf([ZZ(mod_int(_)) for _ in list(crt_elem)])
  return reduced_elem

def norm_form(module_gens, allow_deg_two=True):
  '''
  Returns a univariate polynomial representing the norms of elements in the set {ll*x+alpha: x \in Z}
  where idl=<ll, alpha>=the output of Pari's idealtwoelt function.
  alpha has the same parent as the defining polynomial of idl's number field
  :return:
  '''
  [gen0, gen1] = module_gens
  nf = gen1.parent()
  X = nf.gen()
  def_pol = nf.polynomial()
  univar_polys = PolynomialRing(QQ, names='YY')
  YY = univar_polys.gen()
  bivar_polys = PolynomialRing(univar_polys, names='XX')
  gen0_as_poly = nf_elem_to_poly(gen0, bivar_polys)
  gen1_as_poly = nf_elem_to_poly(gen1, bivar_polys)
  def_pol_as_bivar_poly = bivar_polys(def_pol.coefficients(sparse=False))

  def_poly_parent = def_pol.parent()
  idl = nf.ideal([gen0_as_poly, gen1_as_poly])
  norm_form_in_bivar_polys = (gen0_as_poly*YY + gen1_as_poly).resultant(def_pol_as_bivar_poly)
  norm_form = def_poly_parent(norm_form_in_bivar_polys.coefficients(sparse=False))
  if not norm_form.is_irreducible():
    # if not allow_deg_two:
    #   raise ValueError("poly has the wrong degree. Norm of the ideal is {}. residue class degrees of the sieving idl is {}".format(idl.norm(), [[_[0].smallest_integer(), _[0].residue_class_degree()] for _ in idl.factor()]))
    print(norm_form.factor())
    lead_coeff = norm_form.lc()
    irreduc_fac = norm_form().factor()[0][0]
    print('Warning: norm_form({}) is generally irreducible. This norm form of an ideal thats a product of residue class degrees {} factored as {} over Q[X]. Replacing the norm form with irreducible polynomial {}'.format(module_gens, [_[0].residue_class_degree() for _ in idl.factor()], irreduc_fac, lead_coeff * irreduc_fac))
    norm_form = lead_coeff * irreduc_fac
  ts = module_gens[0] * YY + module_gens[1]

  # assertions are only run on degree 4 norm forms
  if norm_form.degree()==4:
    for k in range(10):
      abs_nrm = nf(ts(k)).norm().abs()
      abs_ff_val = ZZ(norm_form(k)).abs()
      assert abs_nrm == abs_ff_val
  return norm_form

def multivar_norm_form(module_gens):
  '''
  Used for a1 on Feb 16.
  :param module_gens: two generators of a Z-submodule of Ok
  :return:
  '''
  if len(module_gens) != 2:
    raise ValueError('Line sieving uses 2 generators of a Z-submodule of Ok')

  # based on https://ask.sagemath.org/question/52278/is-there-a-way-to-compute-the-norm-form-of-a-number-ring/
  nf = module_gens[0].parent()
  Rng = PolynomialRing(nf, names=["x{}".format(i) for i in range(2)])
  X = PolynomialRing(Rng, name='X').gen()

  # ts_as_poly represents the set to be sieved as a multivariable polynomial
  ts_as_poly=sum([Rng.gen(i) * module_gens[i].lift()(X) for i in range(2)])
  # type = sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular
  ff = ts_as_poly.resultant(nf.polynomial()(X))(x1=1)

  if asserts_on():
    ts = module_gens[0] * X + module_gens[1]
    for k in range(10):
      abs_nrm = nf(ts(k)).norm().abs()
      abs_ff_val = ZZ(ff(x0=k)).abs()
      assert abs_nrm == abs_ff_val

  qq_pols = PolynomialRing(QQ, names='Z')
  qq_pol = ff.polynomial(Rng.gen(0)).map_coefficients(lambda x: QQ(x))
  return qq_pols(list(qq_pol))

def deg4_norm_form_linear(poly_coeffs, aa, bb):
  '''
  Returns the norm for of [aa, YY-bb] in Q[X]/(poly) where poly is constructed from the poly_coeffs
  Code that generated this was:
  PolynomialRing(QQ,'c',6)
  c0,c1,c2,c3,c4,c5 = R.gens()
  S = PolynomialRing(R, names='X,Y')
  X,Y=S.gens()
  Q = X**4 + c3*X**3 + c2*X**2 + c1*X + c0
  (Y*c4 + X - c5).resultant(Q)
  :param poly_coeffs:
  :param aa:
  :param bb:

  '''
  R = PolynomialRing(QQ,'c',6)
  c0,c1,c2,c3,c4,c5 = R.gens()
  Z = PolynomialRing(R, names='Z').gen()
  form = c4**4*Z**4 + (-c3*c4**3 - 4*c4**3*c5)*Z**3 + (3*c3*c4**2*c5 + 6*c4**2*c5**2 + c2*c4**2)*Z**2 + (-3*c3*c4*c5**2 - 4*c4*c5**3 - 2*c2*c4*c5 - c1*c4)*Z + c3*c5**3 + c5**4 + c2*c5**2 + c1*c5 + c0
  multi_form=form(c0=poly_coeffs[0],c1=poly_coeffs[1], c2=poly_coeffs[2], c3=poly_coeffs[3], c4=aa, c5=bb)
  qq_pols = PolynomialRing(QQ, names='Z')
  qq_pol = multi_form.polynomial(Z).map_coefficients(lambda x: QQ(x))
  return qq_pols(list(multi_form))

@cache
def nonzero_positions(vec):
  return tuple([_ for _ in range(len(vec)) if vec[_] != 0])

def fourth_roots(nf):
  polys = PolynomialRing(nf, names='Z')
  Z=polys.gen()
  return [_[0] for _ in (Z**4-1).roots()]

@cache
def idl_pow(idl, pow):
  return idl**pow
def pivot_cols_fastest_available(mat):
  '''
  :param mat: if is_running_on_remote() is True this should be a Magma matrix. Otherwise it should be a dense Sage matrix
  :return:
  '''
  if is_running_on_remote():
    mag_mat=magma(mat)
    pivs = magma_probable_pivot_cols(mag_mat)
  else:
    # very expensive in Sage 10, faster in Sage 9.5 but still not that fast
    pivs = probable_pivot_columns(mat)
  return pivs

def probable_nonpivot_columns_fastest_available(mat):
  '''
  This should be at least 100X faster than the nonpivots() function in integer_matrix.dense. That function calls HNF.
  :param mat:
  :return:
  '''
  pivs = pivot_cols_fastest_available(mat)
  np = []
  for j in range(mat.ncols()):
    if not (j in pivs):
      np.append(j)
  return list(np), list(pivs)

@cache
def s_det(mat):
  '''
  Returns the determinant of a matrix
  :param mat:
  :return:
  '''
  return mat.determinant().abs()

def pr_idls_with_order_of_residue_field_divisible_by_modulus(nf, num_characters, smoothness_bound, modulus, max_dg):
  '''
  :param nf:
  :param num_characters:
  :param smoothness_bound:
  :param modulus:
  :return:
  '''
  from src.stage_init_internal import primes_in_range
  rat_pr_search_bound = smoothness_bound + ZZ(2)**15 #TODO shouldnt I use the inert primes. I know they've been removed from the FB
  rat_primes = primes_in_range(smoothness_bound+1, rat_pr_search_bound+1, lambda pr: pr**4 % modulus == 1)
  pr_idls = []
  for prime in rat_primes:
    pr_idls_with_order_of_residue_field_divisible_by_modulus = nonzero_pr_facs_of_rat_pr_idl(modulus, nf, prime, max_dg)
    pr_idls += pr_idls_with_order_of_residue_field_divisible_by_modulus
    if len(pr_idls) > num_characters:
      return pr_idls[:num_characters]
  raise ValueError("did not find enough primes for characters for modulus={}".format(modulus))

def prime_in_arithmetic_progression(lower_bound, modulus, constant):
  '''
  :param nf:
  :param num_characters:
  :param lower_bound:
  :param modulus:
  :return:
  '''
  max_mult = ZZ(2)**20
  start_value=ceil(lower_bound/modulus)
  for k in range(start_value, max_mult):
    num = ZZ(constant + k * modulus)
    if num.is_prime(proof=True):
      return num, modulus, constant
  else:
    raise ValueError("could not find prime in arithmetic progression")

def prime_with_jacobi_symbol_1_on_prs(lower_bound, moduli):
  '''
  :param lower_bound:
  :param moduli: a list of primes that does not include 2
  :return:
  '''
  modded_values=[]
  for pr in moduli:
    if pr==2:
      raise ValueError("initial list of prs cannot include 2")
    modded_values.append(1)
  moduli= [4]+moduli
  modded_values=[1]+modded_values
  constant=CRT_list(modded_values, moduli)
  return prime_in_arithmetic_progression(lower_bound, prod_all_elems(moduli), constant)

def nonzero_pr_facs_of_rat_pr_idl(modulus, nf, pr, max_dg=1):
  '''
  Returns the prime ideals P in nf above (pr) so that pr^dg=1 (mod modulus) where dg is the residue class degree
  :param modulus:
  :param nf:
  :param prime:
  :return:
  '''
  prime_facs = [k[0] for k in list(nf.ideal(pr).factor())]
  pr_idls_with_order_of_residue_field_divisible_by_modulus = []
  for pr_idl in prime_facs:
    if len(prime_facs) == 1:
      continue
    dg = pr_idl.residue_class_degree()
    if dg > max_dg:
      continue
    else:
      order_of_res_field = pr**dg - 1
    rem = order_of_res_field % modulus
    if rem != 0:
      continue
    pr_idls_with_order_of_residue_field_divisible_by_modulus.append(pr_idl)
  return pr_idls_with_order_of_residue_field_divisible_by_modulus

def primes_of_order_dividing_four_using_arithmetic_progressions(nf, num_characters, bound, modulus):
  '''
  In the case nf is a quartic field is a subset of prime idls P so that the order of the residue class of P is divisible by modulus.
  :param nf:
  :param num_characters:
  :param bound:
  :param modulus: the prime that we are trying to saturate
  :return:
  '''
  if not is_prime(modulus):
    raise ValueError("modulus must be prime")

  Z = PolynomialRing(GF(modulus), names='Z').gen()
  max_mult = ZZ(2)**20
  prs = []
  start_value=ceil(bound / modulus)

  # Commented out stuff would returns prime idls larger than bound that satisfy Z^4-1=(Z^2-1)(Z^2+1) = 0 (mod modulus).
  # but currently I limit the primes to have degree 1
  if modulus % 4 == 1:
    imag_nums = [ZZ(_[0]) for _ in (Z**2 + 1).roots()]
    nf_i = imag_nums[0]
    for k in range(start_value, max_mult):
      # handle the (rare) case that the smallest element in the 4 arithmetic progressions is not large enough
      # if -1 + k * modulus <= bound:
      #   continue

      # if ZZ(-nf_i + k * modulus).is_prime(proof=False):
      #   prs += nonzero_pr_facs_of_rat_pr_idl(modulus, nf, -nf_i + k * modulus)
      # if ZZ(-1 + k * modulus).is_prime(proof=False):
      #   prs += nonzero_pr_facs_of_rat_pr_idl(modulus, nf, -1 + k * modulus)
      if ZZ(1 + k * modulus).is_prime(proof=False):
        prs += nonzero_pr_facs_of_rat_pr_idl(modulus, nf, 1 + k * modulus)
      # if ZZ(nf_i + k * modulus).is_prime(proof=False):
      #   prs += nonzero_pr_facs_of_rat_pr_idl(modulus, nf, nf_i + k * modulus)
      if len(prs) > num_characters:
        break
  elif modulus % 4 == 3:
    for k in range(start_value, max_mult):
      # handle the (rare) case that the smallest element in the 4 arithmetic progressions is not large enough
      # if -1 + k * modulus <= bound:
      #   continue
      # if ZZ(-1 + k * modulus).is_prime(proof=False):
      #   prs += nonzero_pr_facs_of_rat_pr_idl(modulus, nf, -1 + k * modulus)
      if ZZ(1 + k * modulus).is_prime(proof=False):
        prs += nonzero_pr_facs_of_rat_pr_idl(modulus, nf, 1 + k * modulus)
      if len(prs) > num_characters:
        break
  elif modulus % 4 == 2:
    for k in range(start_value, max_mult):
      if ZZ(1 + k * modulus).is_prime(proof=False):
        prs += nonzero_pr_facs_of_rat_pr_idl(modulus, nf, 1 + k * modulus)
      if len(prs) > num_characters:
        break
  dgs = [_.residue_class_degree() for _ in prs]
  if len(prs) >= num_characters:
    return prs[:num_characters]
  else:
    raise ValueError("did not find enough primes for modulus pr={} with max_mult={}".format(modulus, max_mult))

def idl_facs(nf, elem):
  '''
  :param elem: can be an integer
  :return:
  '''
  return [_[0] for _ in nf.ideal(elem).factor()]

def inerts_old(nf, num_characters, Delta):
  max_trials_inerts = 10000*num_characters ** 2
  pr = 0
  poly = nf.polynomial()
  inerts_found = 0
  modulus = ZZ(2)
  inerts_start_time = time.process_time()
  inert_prs=[None]*num_characters
  for z in range(0, max_trials_inerts):
    pr = next_prime(pr,proof=False)
    if pr % modulus != 1:
      continue
    gf_polys = PolynomialRing(GF(pr), names='Z')
    # irreducible poly of degree 4 => poly.discriminant() is not a sqr (mod pr) => delta is not a square (mod pr)
    if jacobi_symbol(Delta, pr) != -1:
      continue
    if poly.change_ring(gf_polys).is_irreducible():
      inert_prs[inerts_found] = pr
      inerts_found += 1
    if inerts_found == num_characters:
      break
  if inerts_found < num_characters:
    raise ValueError("Found {} inerts but need {}".format(inerts_found, num_characters))
  inerts_total_time = time.process_time() - inerts_start_time
  print("Precomputation for character matrix: found logD inerts in {:.2f}. Time per tested prime {:.4f}".format(
    inerts_total_time, inerts_total_time / (z+1)))
  return inert_prs

def hnf_add_rows_with_transformation(hnf_data, new_relns):
  # if hnf_data[0].rank() != hnf_data[0].ncols():
  #   raise ValueError("hnf must have full rank")

  # expand the hnf matrix

  nn=1 # number of new rows to add at a time
  start_time=time.process_time()
  print("Starting to add {} rows to hnf".format(new_relns.nrows()))
  total_time = 0
  for i in range(new_relns.nrows()):
    if i % 5 == 0:
      delta=time.process_time() - start_time
      total_time+=delta
      print("Added 5 rows in {} total rows {}, {}. Transformation height {}".format(delta, i, total_time, hnf_data[1].height()))
      start_time=time.process_time()

    new_mat = hnf_data[0].stack(new_relns[i:i+nn,:])
    U0 = hnf_data[1]
    # expanded original transformation matrix
    U_old = block_matrix([ [U0, zero_matrix(U0.nrows(), nn)], [zero_matrix(nn, U0.ncols()), identity_matrix(nn)]]).change_ring(ZZ)
    # get new transformation matrix U_new by calling hermite_form
    new_hnf, U_new = new_mat.hermite_form(algorithm='flint', transformation=True, include_zero_rows=True)
    # multiply old transformation by new transformation
    hnf_data = [new_hnf, U_new * U_old]
  return hnf_data

def ntl_modular_hnf(mat, include_zero_rows, det_multiple):
  '''
  Based on the Cython code sage/matrix/matrix_integer_dense.pyx. That code contians a bug that makes it not work if the number of rows is greater than the number of cols
  NTL docs say
  // The input matrix A is an n x m matrix of rank m (so n >= m), and D
  // is a multiple of the determinant of the lattice L spanned by the
  // rows of A.  W is computed as the Hermite Normal Form of A; that is,
  // W is the unique m x m matrix whose rows span L, such that
  :param mat:
  :param include_zero_rows:
  :param det_multiple: a multiple of the determinant of the input mat.
  :return:
  '''
  nr = mat.nrows()
  nc = mat.ncols()
  import sage.libs.ntl.ntl_mat_ZZ
  v =  sage.libs.ntl.ntl_mat_ZZ.ntl_mat_ZZ(nr, nc)
  for i in range(0, nr):
    for j in range(0, nc):
      v[i,j] = mat[nr-i-1,nc-j-1]

  try:
    start_hnf_add_rows_time = time.process_time()
    w1 = v.HNF(D=det_multiple)
    timers.linalg_hnf_add_rows += time.process_time() - start_hnf_add_rows_time
  except RuntimeError as e: # HNF may fail if a nxm matrix has rank < m
    raise RuntimeError("NTL bad input for det {}, nrows {}, rank {}".format(det_multiple, v.nrows(), mat.rank()))
    # raise ValueError("ntl only computes HNF for square matrices of full rank.")

  if include_zero_rows:
    H_m = mat.new_matrix()
  else:
    H_m = mat.new_matrix(nrows=w1.nrows())

  nr = w1.nrows()
  nc = w1.ncols()

  for i in range(0, w1.nrows()):
    for j in range(0, w1.ncols()):
      H_m[i,j] = w1[nr-i-1,nc-j-1]
  return H_m

def col_inds_of_dg_two_primes_above_a_rational_with_two_facs(fb, rat_prs, alg_prime_str_dict):
  '''
  Returns the smaller column index of the pair of col inds corresponding to primes with residue_class_degree()==2 above a rational prime with two distinct prime ideal factors.
  Since the fb is sorted the larger index can easily be found if its needed.
  :param fb:
  :param rat_prs:
  :param alg_prime_str_dict:
  :return:
  '''
  dg_two_pair_col_inds =[]
  for pr in rat_prs:
    if pr == rat_prs[-1]:
      rng_max = len(fb)
    else:
      rng_max = alg_prime_str_dict[next_prime(pr,proof=False)]
    if (rng_max - alg_prime_str_dict[pr]) == 2 and fb[alg_prime_str_dict[pr]].residue_class_degree() == 2 and fb[alg_prime_str_dict[pr]+1].residue_class_degree() == 2:
        dg_two_pair_col_inds.append(alg_prime_str_dict[pr])
  return dg_two_pair_col_inds

def abs_std_embedding(elems, field):
  '''
  Source for this function is Sage's number_field code.
  That code has a bug that requires the number of elems to have the same length as the degree
  :param elems:
  :param prec:
  :return:
  '''
  nf = elems[0].parent()
  n = len(elems)
  r,s = nf.signature()
  places = nf.places(prec=field.precision())
  d = {}
  for col in range(n):
    for row in range(r):
      d[(row,col)] = (places[row](elems[col])).abs()
    for row in range(r,r+s):
      d[(row,col)] = (places[row](elems[col])).abs()**2
  return sage.matrix.all.matrix(d)

def class_num_formula_approx(fb, alg_prime_str_dict, prec):
  nf = fb[0].number_field()
  reals = RealField(prec)
  norm_bound = reals(nf.discriminant().abs())**0.5
  r1,r2 = nf.signature()
  euler_approx = euler_prod_approximation(fb, alg_prime_str_dict, prec)
  num_roots = len(nf.roots_of_unity())
  return reals(2)**(-r1)*reals(2*pi)**(-r2)*norm_bound*euler_approx*num_roots

def euler_prod_approximation(fb, alg_prime_str_dict, prec):
  Rprec = RealField(prec)
  prod = Rprec(1)
  rat_primes = list(alg_prime_str_dict.keys())

  i=0
  for i in range(len(rat_primes)):
    pr = rat_primes[i]
    prod *= (1-1/pr)
    if pr == rat_primes[-1]:
      rng_max = len(fb)
    else:
      # get the index of the next prime so we can determine the range of prime indexs to look at
      next_pr = next_prime(pr,proof=False)
      rng_max = alg_prime_str_dict[next_pr]
    for j in range(alg_prime_str_dict[pr], rng_max):
      alg_pr = fb[j]
      prod *= Rprec(1-1/alg_pr.norm())**(-1)

  return prod

def nf_elem_to_vec(nf, elem, ems=None):
  if ems == None:
    return vector([e(elem) for e in nf.embeddings(CC)])
  return vector([e(elem) for e in ems])

def log_T2(elem, prec):
  return RR(T2_elem_real(elem, prec)).log(2)

def T2_elem_real(elem, prec):
  nf = elem.parent()
  try:
    lat_vec = nf.minkowski_embedding([elem] + [0]*3, prec).transpose()[0]
  except Exception:
    raise 'Converting elem to lattice vector failed'
  return lat_vec.norm(2)**2

def T2_sqrt_elems_complex(elems, prec):
  '''
  :param elems: elems of the maximal order or of the field itself
  :param ii: can be a real number greater than 1, infinity ("oo" or "Infinity"), or a symbolic expression.
  :return: a tuple of lengths of the elems under t2
  '''
  if prec < 10:
    raise ValueError("precision too small")
  if len(elems) == 0:
    raise ValueError('must include at least 1 elem')
  lengths = []
  nf = elems[0].parent()
  complex = ComplexField(prec, names='I')
  ems = nf.embeddings(complex)
  for elem in elems:
    try:
      lat_vec = nf_elem_to_vec(nf, elem, ems)
    except Exception:
      raise 'Converting elem to lattice vector failed'
    lengths.append(lat_vec.norm(2))
  return lengths

def rgcd(a,b,regulator_lower_bound):
  '''
  Cohen 5.9.3
  :param a:
  :param b:
  :param regulator_lower_bound: pass 0.2 if no bound is known
  :return:
  '''
  a = abs(a)
  b = abs(b)
  #Cohen doesn't normalize this way
  if a<b:
    return rgcd(b,a, regulator_lower_bound)
  regulator_lower_bound=max(regulator_lower_bound, 0.2)
  while b > regulator_lower_bound:
    tmp_b = b
    b = a - b*floor(a/b)
    a = tmp_b
  return a

def robust_rgcd(a,b, regulator_lower_bound):
  # See May 30 p1-p4 for why this needed. Basically the problem is that rgcd can overshoot the value to return by a single step
  gcds = []
  for x in range(1,4):
    result = rgcd(x*a,x*b, regulator_lower_bound)
    result /= x
    gcds.append(result)
  return max(gcds)

def debugger_is_active() -> bool:
  """Return if the debugger is currently active"""
  return hasattr(sys, 'gettrace') and sys.gettrace() is not None

def asserts_on():
  return __debug__

def is_prod():
  return not asserts_on()

def vec_dense_rep_from_vec_sparse_rep(length, sparse_rep):
  dense_exp_vec=[0]*length
  for i in range(len(sparse_rep)):
    dense_exp_vec[sparse_rep[i][0]] = sparse_rep[i][1]
  return vector(dense_exp_vec)

def str_pairs_list_to_int_list(str):
  if str=='':
    return []
  # a python function that takes a string of comma separated ints as input and returns a list
  strs = str.split("), (")
  strs[0] = strs[0].split("(")[1]
  strs[-1] = strs[-1].split(")")[0]
  return [str_to_int_list(_) for _ in strs]

def str_to_int_list(str):
  if str=='':
    return []
  # a python function that takes a string of comma separated ints as input and returns a list
  strs = str.split(",")
  return [int(_.lstrip()) for _ in strs]

def str_to_rat_list(str):
  # a python function that takes a string of comma separated ints as input and returns a list
  strs = str.split(",")
  return [QQ(_.lstrip()) for _ in strs]

def str_to_ZZ_list(str):
  # a python function that takes a string of comma separated ints as input and returns a list
  strs = str.split(",")
  return [ZZ(_.lstrip()) for _ in strs]

def remove_first_element(lst):
  return lst[1:]

def is_nf_elem_int(elem):
  '''
  Returns true if elem is an int
  :param elem:
  :return:
  '''
  dg = elem.parent().degree()
  return remove_first_element(list(elem)) == [0]*(dg-1)

def conjugates(nf, prec, elem=None):
  if elem==None:
    elem = nf.gen()
  return [_.abs() for _ in elem.complex_embeddings(prec)]
  # return list(nf.minkowski_embedding([nf.gen(),0,0,0], prec).columns()[0])

def log_conjugates(nf, prec, print_prec, elem=None):
  return [RR(max(_,2)).log(2).n(print_prec) for _ in conjugates(nf, prec,elem)]

def max_conjugate(nf, prec, elem):
  return max([_.abs() for _ in conjugates(nf, prec, elem)])


def height(pol):
  return max([ZZ(_).abs() for _ in pol.coefficients(sparse=False)])

def log_height(pol, prec=17):
  return log(height(pol),2).n(prec)

def Delta(pol):
  return NumberField(pol,names='t').discriminant().abs()

def log_Delta(pol, prec=15):
  return log(Delta(pol),2).n(prec)

def log_Delta_nf(nf, prec=15):
  return log(nf.discriminant().abs(),2).n(prec)


def height_as_pow_of_discrim(pol, prec):
  disc = NumberField(pol,names='t').discriminant().abs()
  largest_coeff = max([ZZ(_).abs() for _ in pol.coefficients(sparse=False)])
  log_disc = log(disc, 2)
  return [log(largest_coeff, 2).n(prec) / log_disc.n(prec), log_disc.n(prec)]

def all_zero(bool_list):
  for x in bool_list:
    if x != 0:
      return False
  return True

def prod_all_elems(elems):
  if len(elems) == 0:
    return None

  prod = elems[0]
  for i in range(1, len(elems)):
    prod = prod * elems[i]
  return prod

def pr_idl_not_inert(pr):
  '''
  Returns True if the prime idl is not inert
  :param pr:
  :return:
  '''
  dg=pr.number_field().degree()
  return not pr.residue_class_degree() == dg

def rat_pr_not_inert(rat_pr, def_poly_coeffs):
  '''
  Returns True if the prime idl generated by rat_pr is not inert in Q[X]/(f) where f is the monic poly with coeffs def_poly_coeffs
  :param pr:
  :return:
  '''
  if def_poly_coeffs[-1] != 1:
    raise ValueError("poly must be monic")
  gf_poly = PolynomialRing(GF(rat_pr), names='Z')(def_poly_coeffs)
  return not gf_poly.is_irreducible()

def prime_facs(num):
  return [pair[0] for pair in list(factor(num))]

def prime_idl_facs(numfield, prime):
  '''
  This should be used instead of Sage's nf.primes_above(prime) function which isn't deterministic
  :param numfield:
  :param prime:
  :return:
  '''
  return [numfield.ideal([_[0], _[1]]) for _ in numfield.pari_nf().idealprimedec(prime)]

def index_of_elem_multiplied(elem, elems, mul):
  '''
  In typical usage elems is an index calculus factor base
  :param elem:
  :param elems:
  :param mul:
  :return:
  '''
  if len(elems) == 0:
    raise ValueError('elems should be a list')
  return mul*elems.index(elem)

def count_generator(reader):
  b = reader(1024 * 1024)
  while b:
    yield b
    b = reader(1024 * 1024)

def intersection(lst1, lst2):
  """
  If there are repetitions in one of the inputs then the order of the output depends on the order of the input
  :param lst1:
  :param lst2:
  :return:
  """
  return [element for element in lst1 if element in lst2]

def is_bound_smooth(elem, smoothnessBound, ret_facs=False):
  """
  DEF bound-smooth is an integer that has all prime factors less than smoothness bound
  Result is undefined if smoothnessBound > 2^32

  :param elem:
  :param smoothnessBound:
  :return:
  """
  # for my application it does not make sense to include these as smooths
  if elem == 0:
    raise ValueError("0 does not have a factorization")
  if elem == -1 or elem == 0 or elem == 1:
    if ret_facs == False:
      return True
    return [True, [1]]
  # assumption last element of factor_trial_division is the "unfactored" part
  facs = [pair[0] for pair in list(factor_trial_division(elem, smoothnessBound))]
  # facs = [pair[0] for pair in list(factor(elem))]
  facs = [_ for _ in facs if _ != 1]
  # facs = [pair[0] for pair in list(ZZ(elem).factor())]
  if ret_facs == False:
    return facs[-1] <= smoothnessBound
  return [facs[-1] <= smoothnessBound, facs]

def is_bound_smooth_rat(elem, smoothnessBound, ret_facs=False):
  '''
  :param elem: assumed to be in lowest terms
  :param smoothnessBound:
  :param ret_facs:
  :return:
  '''
  if ret_facs == False:
    return is_bound_smooth(elem.numerator(), smoothnessBound) and is_bound_smooth(elem.denominator(), smoothnessBound)
  [is_num_smooth, num_facs] = is_bound_smooth(elem.numerator(), smoothnessBound, True)
  [is_dnm_smooth, dnm_facs] = is_bound_smooth(elem.denominator(), smoothnessBound, True)
  facs = [_ for _ in (num_facs+dnm_facs) if _ != 1]
  return [is_num_smooth and is_dnm_smooth, facs]

def mat_weight(mat):
  weight=0
  for ri in range(mat.nrows()):
    weight+=mat[ri].hamming_weight()
  return weight

def print_T2_elem_size_and_poly_size(crt_elem, li, lll_elem, nf, sieve_pol):
  print("CRT elem T2={}\nLLL elem T2={}".format(T2_elem_real(crt_elem, 15), T2_elem_real(lll_elem, 15)))
  norm_frm = 1/li * norm_form([nf(li), crt_elem])
  # norm_frm_minus = norm_form([nf(li), crt_elem-nf(li)])
  print(crt_elem)
  print(lll_elem)
  print("max conj of CRT norm form {}".format(conjugates(nf, 15, crt_elem)))
  # print("max conj of CRT norm form {}".format([_[0].abs().n(15) for _ in norm_frm_minus.change_ring(CC).roots()]))
  print("max conj of LLL norm form {}".format(conjugates(nf, 15, lll_elem)))
  # print("max conj of LLL norm form {}".format(max_conjugate(nf_sieve_pol, 15)))
  # print("height of CRT norm form {}".format(log_height(norm_frm, 15)))
  # print("height of LLL norm form {}".format(log_height(sieve_pol, 15)))
  print("size of CRT norm form {}".format(shi_bai_size(CadoQuartic(norm_frm)).n(15)))
  print("size of LLL norm form {}".format(shi_bai_size(CadoQuartic(sieve_pol)).n(15)))

