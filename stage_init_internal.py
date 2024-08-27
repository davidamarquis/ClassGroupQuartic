from sage.all import *
from src.utils_sage import sorted_primes_above
from src.utils import prod_all_elems
from src.determ_rng import *
from src.utils import *

def get_sieve_primes(rat_prs, nf, norm_bound, disc, restrict_index_prs_and_restrict_deg_to_one):
  index=ZZ(nf.polynomial().discriminant().abs()/disc)
  sieve_primes = []
  sieve_rat_prs = []
  other_prs=[]
  for pr in rat_prs:
    if gcd(pr, index)>1 and restrict_index_prs_and_restrict_deg_to_one:
      continue
    prs_above_all_ram_ind = prime_idl_facs(nf, pr)
    if norm_bound==None:
      prs_above = [Pr for Pr in prs_above_all_ram_ind]
    else:
      prs_above = [Pr for Pr in prs_above_all_ram_ind if Pr.norm() <= norm_bound]

    if len(prs_above) >= 2 and prs_above[0].norm()==prs_above[1].norm():
      if restrict_index_prs_and_restrict_deg_to_one and (prs_above[0].residue_class_degree()>1 or prs_above[1].residue_class_degree()>1):
          continue
      first_pr = prs_above[0]
      sieve_primes.append(first_pr)
      sieve_primes.append(prs_above[1])
      sieve_rat_prs.append(pr)
      other_prs.append(prs_above[2:])
  return sieve_primes, sieve_rat_prs, other_prs

def get_sieve_primes_less_than_bound(rat_prs, range_upper_bnd, nf, pr_deg, smoothness_bound, disc, restrict_index_prs_and_restrict_deg_to_one):
  sieve_data = get_sieve_primes(rat_prs, nf, smoothness_bound, disc, restrict_index_prs_and_restrict_deg_to_one)
  return sieve_data

def mat_col_to_elem(col_ind, idl_mat, ww, force_integer=True):
  '''
  :param col_ind:
  :param idl_mat: a square ZZ_matrix representing an idl
  :param ww:
  :param force_integer:
  :return: sage.rings.number_field.number_field_element.OrderElement_absolute (you may want to cast to an nf element)
  '''
  if idl_mat.ncols() == 0:
    raise ValueError('no cols')
  if idl_mat.ncols() < col_ind:
    raise ValueError('not enough columns')
  nf = ww[0].parent()
  elem = nf(0)
  for ri in range(idl_mat.nrows()):
    if force_integer:
      try:
        elem += ZZ(idl_mat[ri][col_ind]) * ww[ri]
      except IndexError as e:
        raise IndexError("tried to access {},{} in a matrix with dims {},{}".format(ri, col_ind, idl_mat.nrows(), idl_mat.ncols()))
    else:
      elem += idl_mat[ri][col_ind] * ww[ri]
  return elem

def alg_primes_w_norm_in_range(rmin, rmax, rat_condition, alg_condition, nf, pari_idl=True, get_one=False, num_primes_required=None):
  """
  Returns a list of pairs [p, p_i] where p is a rational prime and p_i is a pari prime above (p) st N(p_i) in [rmin, rmax)
  :param rmin:
  :param rmax:
  :param rat_condition: a function for adding extra conditions to the rational prime p
  :param alg_condition: a function for adding extra conditions to the algebraic prime p_i
  Note in the alg_condition if pari_idl=True then you can use
  pr_get_f() # the degree
  pr_get_e() # the ramification index
  :return:
  """
  # get a set of rational primes that includes all the possible rational primes
  cand_rat_primes = [ZZ(num) for num in range(1, rmax) if is_prime(num) and rat_condition(num)]
  primes = []
  if num_primes_required == None:
    pr_num = nf.degree()
  else:
    pr_num = num_primes_required
  for prime in cand_rat_primes:
    if pari_idl:
      new_primes = [[prime, k[0].pari_prime()] for k in sorted_primes_above(nf, prime)[:pr_num] if k[0].norm() in range(rmin, rmax) and alg_condition(k[0].pari_prime())]
      if num_primes_required is not None:
        if num_primes_required == len(new_primes):
          primes += new_primes
      else:
        primes += new_primes
    else:
      # primes_above is not deterministic
      # factor is not deterministic
      new_primes = [[prime, k] for k in sorted_primes_above(nf, prime)[:pr_num] if k.norm() in range(rmin, rmax) and alg_condition(k)]
      # new_primes = [[prime, k[0]] for k in list(nf.ideal(prime).factor())[:pr_num] if k[0].norm() in range(rmin, rmax) and alg_condition(k[0])]
      if num_primes_required is not None:
        if num_primes_required == len(new_primes):
          primes += new_primes
      else:
        primes += new_primes
  filtered_primes = []
  if not get_one:
    filtered_primes = primes
  if get_one:
    for pr in primes:
      [a,b] = pr
      if alg_condition(b):
        filtered_primes = [pr]

  return filtered_primes

def primes_in_range(rmin, rmax, condition):
  """
  :param rmin:
  :param rmax:
  :param condition: a function for adding extra conditions to the prime
  :return:
  """
  return [ZZ(num) for num in range(rmin, rmax) if is_prime(num) and condition(num)]

def sparse_vec_of_factor_base_inds(idl, fb):
  if idl.number_field() != fb[0].number_field():
    raise ValueError("nfs should be the same")
  facs = idl.factor()
  indexs=[fb.index(_) for _ in [_[0] for _ in facs]]
  pows=[_[1] for _ in facs]
  sparse_rep=list(zip(indexs, pows))
  return sparse_rep

def fb_of_idls_of_bounded_norm_bounded(fb, smoothness_bound):
  log_smooth_bound = RR(smoothness_bound).log(2)
  new_fb = []
  for Pr in fb:
    pr = Pr.smallest_integer()
    pr_dg = Pr.residue_class_degree()
    pr_rd = Pr.ramification_index()
    # if norm is larger than smoothness bound then continue
    log_nrm = (pr_rd*pr_dg)*RR(pr).log(2)
    if log_nrm > log_smooth_bound or pr_dg==4:
      continue
    new_fb.append(Pr)
  return new_fb

def fac_base_data(numfield, bound=None):
  '''
  Returns a factor basis of prime ideals so that the rational primes below them satisfy pr <= bound
  A latex string is also generated that has information about all the factors
  :param coeffs:
  :return:
  '''
  if bound == None:
    primes_bound = floor(numfield.minkowski_bound())
  else:
    primes_bound = bound

  sanity_check_primes = []
  fb = []
  rat_primes = primes_in_range(0, primes_bound+1, lambda x: True)
  alg_prime_dict = dict.fromkeys(rat_primes)
  prime_numfacs_dict = dict.fromkeys(rat_primes)
  temp_dict = dict.fromkeys(rat_primes)
  temp_dict2 = dict.fromkeys(rat_primes)

  row_str = ''
  alg_pr_start_ind = 0
  for prime in rat_primes:
    prime_facs = prime_idl_facs(numfield, prime)
    sanity_check_primes.append(prime_facs[0]) # use the first prime as a way to sanity check

    for nf_prime in prime_facs:
      fb.append(nf_prime)
    temp_dict[prime] = alg_pr_start_ind
    temp_dict2[prime] = len(prime_facs)

    alg_pr_start_ind += len(prime_facs)

    row_str = row_str + "&(" + str(prime) + ") = " + latex(prime_facs) + str('\\\\')

  alg_prime_dict |= temp_dict
  prime_numfacs_dict |= temp_dict2

  return [row_str, sanity_check_primes, fb, rat_primes, alg_prime_dict, prime_numfacs_dict, numfield]

def reduced_basis_in_random_dir(nf, idl_basis, prec, pow_bound):
  '''
  :param nf:
  :param direction: direction vector
  :param bound: Note: Cohen suggests a pow_bound of 20
  :param num_elems: if 0 we return all the elems that are found
  :param ideal: [N, M, V] where N is num found, M is max T_2, V is a matrix with columns equal to the found vecs
  :param prec:
  :return:
  '''
  if pow_bound < 0:
    raise ValueError("only pow_bound of 0 is allowed")
  r1,r2=nf.signature()
  if pow_bound > 0 and r1+r2-1 == 1:
    raise ValueError("Its pointless to have a pow_bound>0 b/c the unit rank is 1")

  # direction = get_rand_lst(-pow_bound, pow_bound + 1, r1+r2)
  # Cohen suggests only using non-negative powers leq 20 (p.355 Ch 6.5)
  direction = get_rand_lst(0, pow_bound + 1, r1+r2)

  direction += direction[r1:]
  return reduced_basis_in_direction(nf, idl_basis, prec, direction)

def gmat(i,j, ww, ss, vec):
  '''
  Gram matrix defined in Cohen's book Ch6
  :param i:
  :param j:
  :return:
  '''
  n = len(ss)
  total = 0
  for k in range(0, n):
    s_part = (ss[k](ww[i])).conjugate() * ss[k](ww[j])
    total += exp(vec[k]) * s_part
  return total.real_part()

def sorted_embeddings(nf):
  complex_nums = ComplexField(100)
  roots = sorted(nf.polynomial().roots(complex_nums, multiplicities=False),
                 key=lambda x: 1000*x.imag().abs()+x.real().abs())
  v = [nf.hom([e], check=False) for e in roots]
  return v

def fill(nf, vec, basis, real_prec=300):
  # comps = ComplexField(300);
  ss = sorted_embeddings(nf)
  # ss = nf.embeddings(comps)

  deg = nf.polynomial().degree()
  mat = Matrix(RealField(real_prec), deg, deg)
  for i in range(0, mat.nrows()):
    for j in range(0, mat.ncols()):
      mat[i,j] = gmat(i, j, basis, ss, vec)
  return mat

def reduced_basis_in_direction(nf, idl_basis, real_prec, direction):
  '''
  :param nf:
  :param direction: direction vector
  :param bound:
  :param num_elems: if 0 we return all the elems that are found
  :param ideal: [N, M, V] where N is num found, M is max T_2, V is a matrix with columns equal to the found vecs
  :param real_prec:
  :return:
  '''
  #TODO we should check that this direction is legitimate. This means checking that the complex components are equal

  dg=nf.degree()
  if len(idl_basis) < dg:
    raise ValueError('Somehow Zk is not an integral basis of the ideal')
  rr_mat = fill(nf, direction, idl_basis, real_prec)
  #TODO is this used
  pari.set_real_precision(real_prec)
  T = pari(rr_mat).qflllgram()
  if T.ncols() < dg:
    raise ValueError('Somehow minkowski embedding matrix is not full rank. Rank={}'.format(T.ncols()))
  return [sum([ZZ(T[i][j]) * idl_basis[j] for j in range(dg)]) for i in range(dg)]

# Can also compare with FpLLL_reduced_basis
def closest_num(target, nums):
  if len(nums) == 0:
    raise ValueError("nums cant be empty")
  closest_num = None
  min_distance = target + max(nums)
  for num in nums:
    if abs(target - num) < min_distance:
      closest_num = num
      min_distance = abs(target - num)
  return closest_num, min_distance

def consecutive_pairs(iterable):
  '''

  :param iterable: [s0,s1,s2,s3,...]
  :return: [(s0, s1), (s2, s3), (s4, s5)]
  '''
  a = iter(iterable)
  return list(zip(a, a))

def divisor_combinations(primes, num_primes_in_the_divisor):
  '''
  :param primes:
  :param num_primes_in_the_divisor:
  :return: a list of divisors of prod_{p in primes}{p}
  '''
  if num_primes_in_the_divisor > len(primes):
    raise ValueError('num_primes_in_the_divisor cant be larger than len(primes)')
  combs = list(Combinations(primes, num_primes_in_the_divisor))
  return inner_all_prod(combs)

def inner_all_prod(list_of_lists):
  # a function that takes a list of lists as input and returns a list of the product of every element in the inner lists
  return [prod_all_elems(_) for _ in list_of_lists]

def closest_divisor(target, primes, num_primes_in_the_divisor):
  '''
  :param target: the target you want the divisor to be close to
  :param primes:
  :param num_primes_in_the_divisor:
  :return: closest num to divisor, min distance to divisor
  '''
  divs = divisor_combinations(primes, num_primes_in_the_divisor)
  return closest_num(target, divs)