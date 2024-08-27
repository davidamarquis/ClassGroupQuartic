import os
import subprocess
import time

from src.poly_utils import shi_bai_size
from src.cado_poly import *
from src.io_utils import get_poly_file
from src.paths import cado_makefb_binary
from src.utils import norm_form, log_height
from src.utils_arithmetic import construct_idl_with_close_norm
from src.constants import *
from src.stage_init_internal import *
from src.utils_arithmetic import construct_idl_with_prescribed_divisor_and_close_norm
from src.utils import lll_reduced_sieve_poly
import src.timers as timers
from statistics import mean

def roots_diffs(nf, sieve_primes):
  '''
  Assumes that there are at least 2 roots for each of the primes
  :param nf:
  :param sieve_primes:
  :return:
  '''
  poly = nf.polynomial()
  root_diffs = []
  r0s = []
  r1s = []
  for pr in sieve_primes:
    gf_polys = PolynomialRing(GF(pr), names='Y')
    gf_poly = gf_polys(poly)
    roots = gf_poly.roots()
    roots.sort()
    r0 = roots[0][0]
    r1 = roots[1][0]
    if len(roots) < 2:
      raise ValueError("not enough roots")
    r0s.append(r0)
    r1s.append(r1)
    root_diffs.append(r0-r1)
  return [[ZZ(_) for _ in root_diffs], [ZZ(_) for _ in r0s]]

def ideal_generator_differences_and_initial_elems(sieve_primes, def_pol):
  '''
  Assumes that all primes in sieve_primes factor
  :param sieve_primes: a list of pairs of prime ideals that are above the same prime p
  :return:
  '''
  elem_diffs = []
  initial_elems = []
  for pr_pair in sieve_primes:
    if gcd(pr_pair[0].norm(), pr_pair[1].norm())==1:
      raise ValueError("primes must be the same. They are {} and {}".format(pr_pair[0].norm(), pr_pair[1].norm()))
    elem0 = map_pr_idls_to_dedekind_kummer_elements([[pr_pair[0].smallest_integer(), pr_pair[0]]], def_pol)[0]
    elem1 = map_pr_idls_to_dedekind_kummer_elements([[pr_pair[1].smallest_integer(), pr_pair[1]]], def_pol)[0]
    # [_, elem0] = pr_pair[0].gens_two()
    # [_, elem1] = pr_pair[1].gens_two()
    elem_diffs.append(elem0 - elem1)
    initial_elems.append(elem0)
  init_elem_norms = [_.norm().factor() for _ in initial_elems]
  return [elem_diffs, initial_elems]

def lll_reduced_sieve_poly_for_self_init(idl, smallest_int, lll_prec, buchmann_pow_bound):
  nf=idl.number_field()
  nrm=idl.norm()
  lll_time_prev = time.process_time()
  ww = reduced_basis(nf, idl, prec=lll_prec)
  # ww = reduced_basis_in_random_dir(nf, idl.integral_basis(), lll_prec, buchmann_pow_bound)
  timers.init_lll += time.process_time() - lll_time_prev
  count_number_of_times_first_lll_elem_was_not_li = 0
  if smallest_int != ww[0]:
    count_number_of_times_first_lll_elem_was_not_li += 1
  if not is_nf_elem_int(ww[1]):
    lll_elem = ww[1]
  else:
    lll_elem = ww[0]

  nrm=lll_elem.norm()
  if nrm.is_square():
    # an elem of a degree 4 field thats not in Z but has square norm is very likely to be in a quadratic subfield
    lll_elem=ww[2]

  module_gens = [nf(smallest_int), lll_elem]
  # initial_sieve_pol = 1/li * sieving_poly(module_gens, prod)
  resultant_prev_time = time.process_time()
  form = norm_form(module_gens)
  try:
    coeff_gcd = GCD_list(form.coefficients(sparse=False))
  except TypeError as e:
    print(form)
    raise e
  assert nrm % smallest_int == 0
  assert coeff_gcd % nrm == 0
  assert coeff_gcd % nf.ideal(nf(smallest_int), lll_elem).norm() == 0 # they don't need to be equal though
  sieve_pol = 1 / smallest_int * form
  timers.init_resultant += time.process_time() - resultant_prev_time
  if asserts_on():
    smallest_coeffs_pol = 1 / coeff_gcd * form
    height_data = log_height(sieve_pol, 15)
  polys = PolynomialRing(QQ, names='X')
  assert polys(sieve_pol.coefficients(sparse=False)).is_irreducible(),"sieve poly for {} is not irreducible in Q.\nFactors {}".format(module_gens, sieve_pol.factor())
  return [lll_elem, sieve_pol, count_number_of_times_first_lll_elem_was_not_li]

def validate_idl(hnf, idl_gens_two, pari_zk):
  if len(idl_gens_two) != 2:
    raise ValueError("idl_gens_two should have 2 elements")
  nf = idl_gens_two[1].parent()
  idl = nf.ideal(idl_gens_two)
  for x in pari_zk * hnf:
    elem = nf(x)
    assert elem in idl, "expected {} to be in the ideal with factorization {}\n  gens={}\n  HNF={}".format(elem, idl.factor(), idl_gens_two, hnf)

def gray_vec_to_sparse_fb_vec(inds_vec, gray_vec):
  if len(inds_vec)!=2*len(gray_vec):
    raise ValueError("wrong length")
  sparse_vec=[]
  for i in range(len(gray_vec)):
    if gray_vec[i]==0:
      sparse_vec.append((inds_vec[2*i+1],1))
    elif gray_vec[i]==1:
      sparse_vec.append((inds_vec[2*i],1))
    else:
      raise ValueError("gray_vec should only have values 0 or 1")
  return sparse_vec

from src.gray_codes import *
from src.poly_utils import *
def self_init_data(algo, fb, pr_idl_pairs, bounded_rat_prs, disc, is_SI):
  '''
  if algo=sieve_algo_pairs_of_deg1 the Z-module's we sieve are ZZ[li,Y-elem] with elem an integer smaller than li. They are represented as the pair [int, int] for doing IO>
  if algo=sieve_algo_all_prs then the Z-module's we sieve are ZZ[li, elem] with elem in Ok and they are represented as [int, list of rationals] for doing IO.
  Either of these IO reps can be turned into string.
  Python can initialize ints from strings. Sage can initialize rationals from strings.
  Technically, theres little difference but the rep with a pair of integers is analogous to the quadratic case so it needs to stay.
  :param algo:
  :param fb: a list of prime ideals
  :return: [polys, idl_dat]
  '''
  nf = fb[0].number_field()
  polys = []
  idls = []
  idl_gens = []
  fb_exponent_vecs = []
  poly_construct_start_time = time.process_time()

  #TODO remove
  if len(pr_idl_pairs) > 0:
    SI_poly_gen_start=time.process_time()
    rat_sieve_primes = [_[0].smallest_integer() for _ in pr_idl_pairs]
    if len(rat_sieve_primes) != len(pr_idl_pairs):
      raise ValueError("rat_sieve_primes and pr_idl_pairs should have the same length")
    lll_prec = 200
    #TODO add a check that none of the primes divide the discrim
    if len(rat_sieve_primes) == 0 and len(bounded_rat_prs)==0:
      raise ValueError("must be at least one sieve prime")
    # li is short for least_integer
    li = prod_all_elems(rat_sieve_primes)
    li_nf = nf(li)
    initial_roots_list = None
    roots_diff_list = None
    if algo==sieve_algo_pairs_of_deg_1:
      #TODO change root_diffs so prime ideals are passed in
      roots_diff_list, initial_roots_list = roots_diffs(nf, rat_sieve_primes)
    elif algo==sieve_algo_all_prs:
      roots_diff_list, initial_roots_list = ideal_generator_differences_and_initial_elems(pr_idl_pairs, nf.polynomial())
      if len(roots_diff_list) != len(initial_roots_list):
        raise ValueError("lists returned by ideal_generator_differences_and_initial_elems should have the same length")
    print("Num sieve prs {}".format(len(rat_sieve_primes)))
    print("Sieve prs {}\nSecond sieving elems {}".format(rat_sieve_primes, initial_roots_list))

    gray_vec = [1]*len(rat_sieve_primes)
    gray_vecs=[]
    if asserts_on():
      gray_vecs.append(copy(gray_vec))

    Y = nf.gen()
    rat_pr_idls = [nf.ideal(_) for _ in rat_sieve_primes]
    # crt_elem is used to construct the non-integer generator of the Z-module we are sieving.
    crt_elem = None
    coeff_gcds = []
    indexs_of_prs = [fb.index(pr) for pr in [item for sublist in pr_idl_pairs for item in sublist]]

    use_crt=True
    crt_elems=[]

    if algo==sieve_algo_pairs_of_deg_1:
      # elem in ZZ
      crt_elem = CRT_list(initial_roots_list, rat_sieve_primes)
      if crt_elem >= li:
        raise ValueError("elem={} but should be smaller than norm={}".format(crt_elem, li))
      sieve_pol = 1/li*multivar_norm_form([li_nf, Y-crt_elem])
      idl_gens.append([li, crt_elem])
    elif algo==sieve_algo_all_prs:
      crt_elem = nf.idealchinese(rat_pr_idls, initial_roots_list) #an alternative would be to use crt_elems. I like idealchinese more b/c it normalizes the coefficents so they're in the interval [-modulus/2,modulus/2]
      module_gens, sieve_pol, sparse_rep = get_module_gens_and_sieve_poly(use_crt, crt_elem, li, lll_prec, nf, fb,
                                                                          gray_vec, indexs_of_prs, is_SI)
      idl_gens.append([module_gens[0], list(module_gens[1])])

      fb_exponent_vecs.append(sparse_rep)
    crt_elems.append(crt_elem)
    polys.append(sieve_pol)
    try:
      coeff_gcds.append(GCD_list(sieve_pol.coefficients(sparse=False)))
    except TypeError as e:
      raise TypeError("with lll_prec={} the sieve pol was not in ZZ[X]. idl factorization is {}".format(lll_prec, idl.factor()))

    tt = len(rat_sieve_primes)
    crt_coeffs = [CRT_list([0,1],[ZZ(li/_), ZZ(_)]) for _ in rat_sieve_primes]

    for i in range(len(crt_coeffs)):
      assert crt_coeffs[i] ** 2 % rat_sieve_primes[i] == 1
      assert crt_coeffs[i] ** 2 % (li / rat_sieve_primes[i]) == 0

    for expected_gray_vec, changed_ind, sign in GrayIterator(tt):
      if changed_ind >= tt:
        # this means that the iteration is finished
        break
      gray_vec[changed_ind] += sign
      if asserts_on():
        gray_vecs.append(copy(gray_vec))
      if algo==sieve_algo_pairs_of_deg_1:
        crt_elem += sign * crt_coeffs[changed_ind] * roots_diff_list[changed_ind]
        crt_elem = crt_elem % li # this is not done on SPE p343 (TODO compare)
        if crt_elem in crt_elems:
          continue
        norm_form=multivar_norm_form([li_nf, Y-crt_elem])
        sieve_pol = 1/li*norm_form

        idl_gens.append([li, crt_elem])
      elif algo==sieve_algo_all_prs:
        crt_elem += sign * crt_coeffs[changed_ind] * roots_diff_list[changed_ind]
        crt_elem = reduce_coordinates_of_nf_elem(crt_elem, li)
        if crt_elem in crt_elems:
          continue
        module_gens, sieve_pol, sparse_rep = get_module_gens_and_sieve_poly(use_crt, crt_elem, li, lll_prec, nf, fb,
                                                                            gray_vec, indexs_of_prs, is_SI)
        idl_gens.append([module_gens[0], list(module_gens[1])])
        fb_exponent_vecs.append(sparse_rep)
      crt_elems.append(crt_elem)
      polys.append(sieve_pol)
      coeff_gcds.append(GCD_list(sieve_pol.coefficients(sparse=False)))

  if not is_SI:
    roots_for_sieve_data(list(set(bounded_rat_prs)), polys)
    print("Time for SI is {}".format(time.process_time()-SI_poly_gen_start))

  timers.init_poly_construct = time.process_time()-poly_construct_start_time
  return polys, idl_gens, fb_exponent_vecs

def target_sieve_data(bounded_fb, disc, fb, fb_exponent_vecs, idl_gens, nf, polys, pr_idls_to_target, is_SI, bounded_rat_prs):
  #TODO remove code for fb_exponent_vecs. It makes no difference

  prec_for_knapsack = 6
  prec = 100
  target = lattice_sieving_bound(nf.polynomial())
  other_prs_sizes = []
  sieve_extra_poly_gen_start = time.process_time()
  indx = index(nf.polynomial())

  new_fb = [_ for _ in fb if _.residue_class_degree() == 1 and gcd(_.norm(), indx) == 1]
  fb_nrms = [_.norm() for _ in new_fb]
  fb_rat_prs = [ZZ(_.smallest_integer()) for _ in new_fb]
  dede_kummer_time_start = time.process_time()
  second_gen_residues = map_pr_idls_to_dedekind_kummer_elements(
    list(zip(fb_rat_prs, new_fb)), nf.polynomial())
  print("FB subset we may sieve on has length {} out of {} prime ideals above primes less than B".format(len(new_fb), len(fb)) + os.linesep)
  print("Time for Dede kummer is {:.2f}".format(time.process_time() - dede_kummer_time_start))
  LLL_time = 0
  shi_time = 0
  pr_idls_to_target = [_ for _ in pr_idls_to_target if _.residue_class_degree() == 1 and gcd(_.norm(), indx) == 1]
  for pr in pr_idls_to_target:

    binary = False
    max_num_primes = None

    MP_knapsack_start = time.process_time()
    try:
      [prescribed_idl, idl_norm, prs_used, exp, rat_prs_used] = construct_idl_with_prescribed_divisor_and_close_norm(
      new_fb, fb_nrms, fb_rat_prs, target, prec_for_knapsack, pr, binary, max_num_primes)
    except sage.numerical.mip.MIPSolverException:
      raise ValueError("Failed to find targeted solution for prime pr={}".format(pr))
    timers.MP_knapsack_total += time.process_time()-MP_knapsack_start

    assert len(new_fb) == len(exp)
    list_of_pairs = [(i, exp[i]) for i in range(len(exp)) if exp[i] != 0]

    LLL_sieve_poly_start = time.process_time()
    second_gen_residues_used = [second_gen_residues[_] for _ in range(len(second_gen_residues)) if exp[_]!=0]
    try:
      sieve_pol, common_factor_of_poly_coeffs, mod_gens = crt_sieve_poly(prs_used, second_gen_residues_used)
    except ZeroDivisionError:
      print("crt of Dedekind elems failed. Probably b/c knapsack selected the same prime")
      continue
    LLL_time += time.process_time()-LLL_sieve_poly_start

    sieve_pol_with_common_fac = common_factor_of_poly_coeffs*sieve_pol  # write_cado handles the common factor
    polys.append(sieve_pol_with_common_fac)
    idl_gens.append([list(_) for _ in mod_gens])
    fb_exponent_vecs.append(list_of_pairs)

  if not is_SI:
    roots_for_sieve_data(bounded_rat_prs, polys)

  print("Time for sieve extra is {}\nLLL time {}\nShi time {}".format(time.process_time()-sieve_extra_poly_gen_start,
                                                                      LLL_time, shi_time))

def roots_for_sieve_data(bounded_rat_prs, polys):
  root_time_start = time.process_time()
  for prime in bounded_rat_prs:
    int_mod_polys = PolynomialRing(GF(prime), names='Z')
    for poly in polys:
      mod_poly = int_mod_polys(poly)
      root_poly = mod_poly.roots()
  root_time = time.process_time()-root_time_start
  sys.stderr.write("MP root time {:.2f} ({:.2f} per poly)".format(root_time, root_time/(len(polys)))+os.linesep)


def get_module_gens_and_sieve_poly(use_crt, crt_elem, li, lll_prec, nf, fb, gray_vec, indexs_of_prs, is_self_init):
  '''
  I use 'sieve_poly' to mean the norm form of I divided by N(I)
  :param use_crt:
  :param crt_elem:
  :param li:
  :param lll_prec:
  :param nf:
  :param fb:
  :param gray_vec:
  :param indexs_of_prs:
  :return:
  '''
  if use_crt:
    module_gens = [nf(li), crt_elem]
    allow_deg_2=True
    norm_form_start=time.process_time()

    if is_self_init:
      # crt_elem has the form X-r so we can calculate it using the formula from my thesis
      coeffs=nf.polynomial().coefficients(sparse=False)
      form = deg4_norm_form_linear(coeffs, ZZ(li), -ZZ(list(crt_elem)[0]))
    else:
      form = norm_form(module_gens, allow_deg_2)
    sieve_pol = 1/li * form
    #note: it appears resultants for the same defining poly are cached so there's minimal time savings for using the formula
    idl = nf.ideal(li, crt_elem)
  else:
    [lll_elem, lll_sieve_pol, _] = lll_reduced_sieve_poly_for_self_init(nf.ideal(li, crt_elem), li, lll_prec, None)
    module_gens = [nf(li), lll_elem]
    sieve_pol = lll_sieve_pol
    idl = nf.ideal(li, lll_elem)
  # print_T2_elem_size_and_poly_size(crt_elem, li, lll_elem, nf, lll_sieve_pol)
  sparse_rep = sparse_rep_of_new_idl_on_factor_base(fb, gray_vec, idl, indexs_of_prs, li, crt_elem, nf, 0)
  # sparse_rep = gray_vec_to_sparse_fb_vec(indexs_of_prs, gray_vec)
  # sparse_rep = sparse_vec_of_factor_base_inds(idl, fb)
  return module_gens, sieve_pol, sparse_rep

def sparse_rep_of_new_idl_on_factor_base(fb, gray_vec, original_idl, indexs_of_prs, li, lll_elem, nf, index):
  new_idl = nf.ideal(li, lll_elem)
  sparse_rep_of_idl = gray_vec_to_sparse_fb_vec(indexs_of_prs, gray_vec)
  sparse_rep = sparse_vec_of_factor_base_inds(new_idl/original_idl, fb)+sparse_rep_of_idl
  if set(sparse_rep) != set(sparse_vec_of_factor_base_inds(new_idl, fb)):
    print("Factorization of this ideal is different from the original ideal at FB index {}. Factorization of their quotient is {}".format(index, new_idl/original_idl))
  return sparse_rep_of_idl

def roots_matrix_lambda(pr, mat, rat_fb, col_ind, nf_deg):
  '''
  :param pr: ZZ
  :param mat: mat_ZZ assumes the matrix uses the root-rows layout
  :param rat_fb: [ZZ]
  :param col_ind: int, specifies which of the ideals in the basis will be used
  :param nf_deg
  :return:
  '''
  # Roots matrix is given in root-rows layout
#
#             Col1  Col2
#            (Idl1) (Idl2)
#            _____________....
# Pr1 Root1  |     |     |
# Pr1 Root2  |     |     |
#            |     |     |
#            |_____|_____|
# Pr2 Root1  |     |     |
#            |     |     |
#            |     |     |
#            |_____|_____|
#            .
#            .
  row_start_ind = index_of_elem_multiplied(pr, rat_fb, nf_deg)
  ind = None
  for i in range(nf_deg):
    if mat[row_start_ind + i, col_ind] == -1:
      ind = row_start_ind + i
  roots_list = mat[row_start_ind: ind, col_ind].transpose().list()
  return lambda line_param : line_param % pr in roots_list

def get_fb_args(path_to_params_dir, params_spec_str, verbose, smoothness_bound):
  '''
  Parse a parameters file to get a parameters:
  qrange, I lim0 lpb0 mfb0
  This is used to get a full list parameters for running las
  :return:
  '''
  params_file = get_params_file(path_to_params_dir, params_spec_str, verbose)
  params_full_path = path_to_params_dir + params_file

  # Step 1: look in the params file for name = NAME if it doesn't exist then return
  # if there's only one file with '_params' then just use that
  with open(params_full_path, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      newline = line.replace(" ", '')
      name_str = "name="
      if newline.startswith(name_str):
        job_name = newline.split(name_str,1)[1]
        job_name = job_name.rstrip()
        if verbose:
          print("Job name {}".format(job_name))

  poly_path = get_poly_file(path_to_params_dir, False, params_spec_str, verbose)

  if params_spec_str == None:
    out_name = path_to_params_dir + job_name + ".roots0.gz".format()
  else:
    out_name = path_to_params_dir + params_spec_str + ".roots0.gz".format()

  args = " -poly {}".format(poly_path)
  args += " -lim {}".format(smoothness_bound)
  args += " -out {}".format(out_name)
  args += " -t 2"
  args += " -side 0"
  return args

def fb_cmd(path_to_workdir, param_spec_str, verbose, smoothness_bound):
  args = get_fb_args(path_to_workdir, param_spec_str, verbose, smoothness_bound)
  cmd = cado_makefb_binary + args
  if verbose:
    print(cmd)
  return cmd

def get_bound(path_to_params_dir, params_spec_str, sieve_rad, verbose):
  '''
  :param path_to_params_dir:
  :param params_spec_str:
  :param verbose:
  :return:
  '''
  poly_file = get_poly_file(path_to_params_dir, False, params_spec_str, verbose)
  with open(poly_file, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      newline = line.replace(" ", "")

      max_conj_str = "max_conj:"
      if newline.startswith(max_conj_str):
        max_conj = float(newline.split(max_conj_str)[1].rstrip())
        if verbose:
          print("max_conj {}".format(max_conj))
      elif "conj" in newline:
        raise ValueError("found substring {} in a line but formatting on the line was incorrect".format("conj"))

  return floor(abs(max_conj)/(4*sieve_rad))

def lattice_sieving_bound(poly):
  '''
  :param Delta:
  :param dg: the degree of the number field
  :return:
  '''
  return cohen_size(poly)

def check_coefficients_satisfy_Cohen_bound(poly, bound_multiplier, raise_err):
  '''
  Tests that the poly satisfies $|c_{n-k}| <= bound_multiplier * bin(n,k)(T_2(g)/n)^(k/2)
  :param nf:
  :param bound_multiplier: a fudge factor larger than 1
  :return:
  '''
  if bound_multiplier < 1:
    raise ValueError("fudge factor should be at least 1")
  nf = NumberField(poly, names='Y')
  Y = nf.gen()
  dg = nf.degree()
  dsc = nf.discriminant().abs()
  # t2_poly = ti_lengths([Y] + [0]*(dg-1), ii=2)[0]
  t2_poly_if_poly_is_close_to_polred_abs = (bound_multiplier * dsc/4)**(1/(2*(dg-1))).n()
  SS = t2_poly_if_poly_is_close_to_polred_abs**2

  bounds = [bound_multiplier * binomial(dg, _) * (SS / dg) ** ((dg - _)/2.0) for _ in range(dg + 1)]
  coeffs = nf.polynomial().coefficients(sparse = False)
  bound_minus_coeff = [(bounds[_] - coeffs[_]).n(15) for _ in range(dg+1)]
  for i in range(dg+1):
    if coeffs[i] > bounds[i]:
      if raise_err:
        raise ValueError("Coeff i={} did not satisfy c_i < bound_i".format(i))
      else:
        print("Warning: Coeff i={} did not satisfy c_i < bound_i".format(i))
        return False
  return True

def run_batch_fb(path_to_workdir, path_to_params_file, verbose, pol, max_num_primes, algo, fb, smoothness_bound,
                 sieve_radius, job_name, disc_div, disc, is_SI):
  '''
  Does the following
  1. Generates sieving polys
  2. Constructs "_params" and ".poly" files for all the polys
  3. Runs las on all the files

  :param is_SI:
  :param path_to_workdir: the path to the directory that the defining poly's params file is in
  :param path_to_params_file:
  :param verbose:
  :param pol: the poly is assumed to be polredabs
  :return:[procs, params_paths] these values are useful when you intend to run the sieving stage right after all the procs in the fb stage are done.
  For larger computations the plan is to run these stages separately
  '''
  init_time_start = time.process_time()
  # throws error if poly doesn't have the coeffs that are the right shape
  if max_num_primes == None:
    raise ValueError("max_num_primes must be set")
  if algo != sieve_algo_pairs_of_deg_1 and algo != sieve_algo_all_prs:
    raise ValueError('algo parameter must be one of {} or {}'.format(sieve_algo_pairs_of_deg_1, sieve_algo_all_prs))

  check_coefficients_satisfy_Cohen_bound(pol, 2, False)
  disc=Delta(pol)
  nf = fb[0].number_field()
  skew_control_constant=0.0

  print("Starting get_sieve_primes")
  prec_for_knapsack=32
  SI_knapsack_start = time.process_time()

  restrict_index_prs_and_restrict_deg_to_one=True

  bounded_fb = fb_of_idls_of_bounded_norm_bounded(fb, smoothness_bound)
  bounded_nrms_start=time.process_time()
  bounded_nrms = [_.norm() for _ in bounded_fb]
  bounded_rat_prs = [_.norm()**(_.residue_class_degree()**(-1)) for _ in bounded_fb]
  bounded_nrms_end=time.process_time()-bounded_nrms_start
  try:
    alg_sieve_prs, nrm, unused_pr_idls_above_same_rat_prs_as_sieve_prs = get_alg_prs_by_knapsack(max_num_primes, nf,
                                                                                                 smoothness_bound,
                                                                                                 prec_for_knapsack,
                                                                                                 disc,
                                                                                                 restrict_index_prs_and_restrict_deg_to_one,
                                                                                                 list(set(
                                                                                                   bounded_rat_prs)))
  except sage.numerical.mip.MIPSolverException:
    raise ValueError("Failed to find solution for max_num_primes={}".format(max_num_primes))
  timers.SI_knapsack = time.process_time()-SI_knapsack_start

  print("Time to compute bounded norms {}".format(bounded_nrms_end))

  max_num_extra_polys = ZZ(2)**8
  unused_pr_idls_above_same_rat_prs_as_sieve_prs = unused_pr_idls_above_same_rat_prs_as_sieve_prs[:min(max_num_extra_polys, len(unused_pr_idls_above_same_rat_prs_as_sieve_prs))]
  sys.stderr.write("Number of extra polys that _could_ be sieved is {} out of a max of {}".format(len(unused_pr_idls_above_same_rat_prs_as_sieve_prs), max_num_extra_polys) + os.linesep)
  polys=[]
  idls=[]
  fb_exponent_vecs=[]

  [polys, idls, fb_exponent_vecs] = self_init_data(algo, fb, consecutive_pairs(alg_sieve_prs),
                                                   bounded_rat_prs, disc, is_SI)
  sys.stderr.write("Def poly {}".format(nf.polynomial()) + os.linesep)
  sys.stderr.write("First sieving poly is {}".format(polys[0]) + os.linesep)
  sys.stderr.write("|initial idls|={}".format(len(idls)) + os.linesep)
  sieve_extra=len(idls) < ZZ(2)**8 and not is_pure_quartic(nf.polynomial())
  if sieve_extra:
    sys.stderr.write("Sieve extra: enabled" + os.linesep)
    print("Sieve extra: enabled {}".format(unused_pr_idls_above_same_rat_prs_as_sieve_prs))
    target_sieve_data(fb, disc, fb, fb_exponent_vecs, idls, nf, polys, unused_pr_idls_above_same_rat_prs_as_sieve_prs,
                      is_SI, bounded_rat_prs)

  if len(fb_exponent_vecs) != len(polys):
    raise ValueError("fb_exp_vecs should have same length as idls")
  if fb_exponent_vecs[0] == None:
    raise ValueError("fb_exp_vecs[0] should not be none")

  num_deg_two=0
  for i in range(len(polys)):
    poly = polys[i]
    if poly.degree()==2:
      num_deg_two+=1
    common_fac = GCD_list(poly.coefficients(sparse=False))
  print("num degree two polys={} out of {}".format(num_deg_two, len(polys)))

  print("Starting to write poly files")
  poly_paths = [0]*len(polys)
  params_paths = [0]*len(polys)
  params_files = [0]*len(polys)

  # make a temp params file
  path_to_temp = path_to_workdir + "temp_params"
  print("copying from {} to {}".format(path_to_params_file, path_to_temp))
  subprocess.run(["cp", path_to_params_file, path_to_temp])

  # remove tasks.polyselect.import value from temp file
  polyselect_str = "tasks.polyselect.import"
  found_polyselect = False
  with open(path_to_temp, "r") as f:
    lines = f.readlines()
  # this overwrites the existing file and writes all the lines except the polyselect string
  str_written = ''
  with open(path_to_temp, "w") as f:
    for line in lines:
      newline = line.replace(" ","")
      if not newline.startswith(polyselect_str):
        str_written += newline
        f.write(newline)
      else:
        found_polyselect = True
  if found_polyselect == False:
    raise ValueError("did not find path to polynomial in the file")

  for ind in range(len(polys)):
    cado_pol = CadoQuartic(polys[ind])
    path_to_poly = path_to_workdir + job_name + "_Q{}i{}.poly".format(nrm, ind)
    poly_paths[ind] = path_to_poly
    write_cado(cado_pol, poly_paths[ind], True, idls[ind], fb_exponent_vecs[ind], log_output_bitlength=None, disc_div=disc_div)

    param_filename = job_name + "_Q{}i{}_params".format(nrm,ind)
    path_to_new_params = path_to_workdir + param_filename
    params_paths[ind] = path_to_new_params
    params_files[ind] = param_filename

    subprocess.run(["cp", path_to_temp, path_to_new_params])

    if cado_pol.poly.degree() == 2:
      # sys.stderr.write("Found degree 2 sieve poly. Writing params file with sieve radius I={}.".format(sieve_radius) + os.linesep)
      output=[]
      with open(path_to_new_params, "r") as fp:
        lines = fp.readlines()
        for line in lines:
          if line.startswith("tasks.I"):
            output.append("tasks.I = {}\n".format(sieve_radius))
            continue
          output.append(line)
      with open(path_to_new_params, "w") as fp:
        fp.writelines(output)

    params_poly_str = "\n"
    params_poly_str += polyselect_str
    params_poly_str += " = " + path_to_poly
    with open(path_to_new_params, "a") as fp:
      fp.write(params_poly_str)

  start_start_time = time.process_time()

  heights=[log_height(_,15) for _ in polys]
  num_prs=int(len(alg_sieve_prs)/2.0)
  print_stats(fb, heights, nf, pol, prec_for_knapsack, smoothness_bound, num_prs)

  timers.stats_time = time.process_time() - start_start_time

  cmds = []
  makefb_time_start = time.process_time()
  (params_dir, _) = os.path.split(path_to_params_file)

  for path_to_params_file in params_files:
    cmds.append(fb_cmd(path_to_workdir, path_to_params_file, True, smoothness_bound))
  pool_size=1
  print("Running makefb with pool size={}".format(pool_size))
  limited_concurrent_commands(cmds, pool_size)
  timers.init_makefb = time.process_time() - makefb_time_start
  timers.init_total = time.process_time() - init_time_start
  timers.init_total -= timers.stats_time
  return [cmds, params_files]

def print_stats(fb, heights, nf, pol, prec_for_knapsack, smoothness_bound, num_sieving_prs):
  num_self_init_pols = ZZ(2)**num_sieving_prs
  cado_orig_poly = CadoQuartic(pol)
  size_value_nbits = 10
  heights=heights[:num_self_init_pols]

  min_height = min(heights)
  mean_height = mean(heights)
  height_def = log_height(nf.polynomial(), 15)
  print("Knapsack prec={}".format(prec_for_knapsack))
  norm_bnd_fb = fb_of_idls_of_bounded_norm_bounded(fb, smoothness_bound)
  print("Truncated heights list to length {}=2^{}. SIQS smallest height={}\nDefining polynomial height {}\nSIQS mean height {}".format(len(heights), num_sieving_prs, min_height, height_def, mean_height))
  print("heights {}".format(heights[:]))
  print("Factor base contains {} prime ideals".format(len(norm_bnd_fb)))

def get_alg_prs_by_knapsack(max_num_primes, nf, smoothness_bound, prec_for_knapsack, disc,
                            restrict_index_prs_and_restrict_deg_to_one, rat_prs):
  target_product=floor(lattice_sieving_bound(nf.polynomial()))
  NN=target_product
  pr_deg=None
  S_Z = get_sieve_primes(rat_prs, nf, smoothness_bound, disc, restrict_index_prs_and_restrict_deg_to_one)[1]

  candidate_sieve_prs_data = get_sieve_primes_less_than_bound(S_Z, floor(NN**(1/5)), nf, pr_deg, smoothness_bound, disc, restrict_index_prs_and_restrict_deg_to_one)
  tt=min(max_num_primes, len(S_Z))
  S_N=candidate_sieve_prs_data[1]
  for i in range(100):
    if len(S_N) >= 5:
      break
    NN*=2
    candidate_sieve_prs_data = get_sieve_primes_less_than_bound(S_Z, floor(NN**(1/5)), nf, pr_deg, smoothness_bound, disc, restrict_index_prs_and_restrict_deg_to_one)
    S_N=candidate_sieve_prs_data[1]

  candidate_sieve_prs=candidate_sieve_prs_data[0]
  if len(candidate_sieve_prs)==0:
    raise ValueError("WARNING: No candidate sieve prs found.")

  other_prs=candidate_sieve_prs_data[2]
  binary=True

  initial_pr_selection_start=time.process_time()
  [_, nrm, sieve_prs, exps, _] = construct_idl_with_close_norm(candidate_sieve_prs_data[0], target_product, prec_for_knapsack, binary, tt)

  solve_time = time.process_time() - initial_pr_selection_start
  sys.stderr.write("sieving stats\ntarget={}\n"
                   "t={}\n"
                   "|S_Z|={}\n"
                   "|S_N|={}\n"
                   "Num idls above primes in S_N={}\n"
                   "Sieve prs {}\n"
                   "solve time {}".format(NN, tt, len(S_Z), len(S_N), len(candidate_sieve_prs_data[0]), sieve_prs, solve_time) + os.linesep)

  full_sieve_prs=[]
  for i in range(len(exps)):
    if exps[i]!=0:
      first_pr=candidate_sieve_prs[i]
      full_sieve_prs.append(first_pr)
      second_pr=candidate_sieve_prs[i+1]
      full_sieve_prs.append(second_pr)
  return full_sieve_prs, nrm, flatten(other_prs)

def get_alg_prs_by_searching_divisors(algo, bnd, max_num_primes, nf, smoothness_bound, raise_if_nearest_divisor_is_far_from_target):
  if algo == sieve_algo_pairs_of_deg_1:
    prime_size_parameter = 10
    prime_bound = min(floor(prime_size_parameter*bnd**(max_num_primes**(-1))), smoothness_bound)
    [_, candidate_sieve_prs] = get_sieve_primes_less_than_bound(prime_bound, nf, pr_deg=1)
  elif algo == sieve_algo_all_prs:
    prime_size_parameter = 10
    prime_bound = min(floor(prime_size_parameter*bnd**(max_num_primes**(-1))), smoothness_bound)
    [_, candidate_sieve_prs] = get_sieve_primes_less_than_bound(prime_bound, nf,
                                                                pr_deg=None)  # probably a smaller bound would work
    # filter the list of rational primes
    _, candidate_sieve_prs = get_sieve_primes(candidate_sieve_prs, nf, )
    for i in range(3):
      if len(candidate_sieve_prs) > 0:
        break
      print("No sieving primes found. Doubling the size of the sieving bound")
      prime_size_parameter *= 2
      prime_bound = min(floor(prime_size_parameter*bnd**(max_num_primes**(-1))), smoothness_bound)
      [_, candidate_sieve_prs] = get_sieve_primes_less_than_bound(prime_bound, nf,
                                                                  pr_deg=None)  # probably a smaller bound would work
      _, candidate_sieve_prs = get_sieve_primes(candidate_sieve_prs, nf, )
  print("the number of possible sieve primes less than bound={} is {}".format(bnd, len(candidate_sieve_prs)))
  if len(candidate_sieve_prs) < max_num_primes:
    print("Warning: truncating the number of polynomials. If more polynomials are needed, set the prime_size_parameter")
  num_primes = min(len(candidate_sieve_prs), max_num_primes)
  num_divs = binomial(len(candidate_sieve_prs), num_primes)
  if len(candidate_sieve_prs) > 10:
    print("Starting to compute bin({} {})={} divisors".format(len(candidate_sieve_prs), num_primes, num_divs))
  closest_divisor_start = time.process_time()
  (nrm, distance) = closest_divisor(bnd, candidate_sieve_prs, num_primes)
  print("time for closest divisor {}".format(time.process_time()-closest_divisor_start))
  if raise_if_nearest_divisor_is_far_from_target and distance > 5000000:
    raise ValueError(
      "WARNING: distance of closest divisor is large {}. Number of primes may be too large. Possible sieving primes are {}".format(
        distance, candidate_sieve_prs))
  print("Num availble primes for sieving (given how the set of possible primes has been restricted) {}".format(
    candidate_sieve_prs))
  print(
    "Rational sieving primes {} after rejecting those primes with fewer than two degree 1 prime factors, and of ramification index == 1.".format(
      candidate_sieve_prs))
  actual_sieve_prs = prime_facs(nrm)
  alg_sieve_prs, _, other_prs = get_sieve_primes(actual_sieve_prs, nf, )
  return alg_sieve_prs, nrm, other_prs
