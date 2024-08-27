import os
import shlex
import argparse

from src.stage_init import *
from src.utils import *
from src.utils_arithmetic import *
import src.timers as timers

def get_las_args(path_to_params_dir, out_dir, params_spec_str, verbose, smoothness_bound, expected_cg_generator_bound):
  '''
  path_to_params_dir: should be terminated with '/'
  Open the parameters file to get ```qrange, I lim0 lpb0 mfb0``` and construct two input strings for las for the SI and pad steps respectively.
  :return:
  '''
  params_file = get_params_file(path_to_params_dir, params_spec_str, verbose)
  params_full_path = path_to_params_dir + params_file

  # Step 1: look in the params file for name = NAME if it doesn't exist then return
  # if there's only one file with '_params' then just use that
  job_name = get_job_name(path_to_params_dir, params_spec_str, verbose)

  # Step 2: look for NAME.roots0.gz
  roots0_name = params_spec_str + ".roots0.gz"
  roots0_full_path = path_to_params_dir + roots0_name
  try:
    open(roots0_full_path, "r")
  except:
    print("Las Warning: polynomial file with no roots found but I am continuing")
    pass

  # Step 3: cp NAME.roots0.gz > NAME.roots1.gz
  roots1_name = params_spec_str + ".roots1.gz"
  roots1_full_path = path_to_params_dir + roots1_name

  ret_code = subprocess.run(["cp", roots0_full_path, roots1_full_path], stdout=subprocess.DEVNULL)
  lim0 = get_smooth_bound(path_to_params_dir, params_spec_str, verbose)

  # Step 4: get the other params
  with open(params_full_path, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      newline = line.replace(' ','')

      poly_str = "tasks.polyselect.import="
      I_str = "tasks.I="
      las_num_threads = "tasks.las.threads="
      qrange_str = "tasks.sieve.SI.qrange="
      lpb0_str = "tasks.sieve.SI.lpb0="
      mfb0_str = "tasks.sieve.SI.mfb0="
      if newline.startswith(qrange_str):
        qrange = int(newline.split(qrange_str)[1].rstrip())
        if verbose:
          print("qrange {}".format(qrange))
      if newline.startswith(I_str):
        II = int(newline.split(I_str)[1].rstrip())
        if verbose:
          print("I {}".format(II))
      if newline.startswith(lpb0_str):
        lpb0 = int(newline.split(lpb0_str)[1].rstrip())
        lpb1 = lpb0
        if verbose:
          print("lpb0 {}".format(lpb0))
      if newline.startswith(mfb0_str):
        mfb0 = int(newline.split(mfb0_str)[1].rstrip())
        mfb1 = mfb0
        if verbose:
          print("mfb0 {}".format(mfb0))
      if newline.startswith(poly_str):
        poly_path = newline.split(poly_str)[1].rstrip()
        if verbose:
          print("poly path {}".format(poly_path))
      if newline.startswith(las_num_threads):
        las_num_threads_from_params= newline.split(las_num_threads)[1].rstrip()
  lim1=lim0

  ap = argparse.ArgumentParser('parse las arguments')
  ap.add_argument('--las_t', nargs='?')
  las_parse = ap.parse_args()

  las_args = ""
  # las_args += " -adjust-strategy 2".format(poly_path)
  las_args += " -poly {}".format(poly_path)
  las_args += " -lim0 {}".format(lim0)
  las_args += " -lim1 {}".format(lim1)
  las_args += " -fb0 {}".format(roots0_full_path)
  las_args += " -fb1 {}".format(roots1_full_path)

  ### set the parts of the arguments that are different
  q0 = 10
  q1 = qrange + q0
  out_filename = out_dir + params_spec_str + ".{}-{}.gz".format(q0, q1)
  las_args += " -I {}".format(II)
  las_args += " -q0 {}".format(q0)
  las_args += " -q1 {}".format(q1)
  # las_args += " -tdthresh 1024 "
  # las_args += " -lambda0 0.1"
  # las_args += " -lambda1 0.1"
  # las_args += " -batch "
  # las_args += " -batchlpb0 {}".format(lpb0)
  # las_args += " -batchlpb1 {}".format(lpb1)
  # las_args += " -batchmfb0 {}".format(mfb0)
  # las_args += " -batchmfb1 {}".format(mfb1)
  # las_args += " -dup -dup-qmin q0 -dup-qmax q1" # '2024 Jan experiment showed it slowed down sieving on an nf of about 100 bits
  las_args += " -lpb0 {}".format(lpb0)
  las_args += " -lpb1 {}".format(lpb1)
  las_args += " -mfb0 {}".format(mfb0)
  las_args += " -mfb1 {}".format(mfb1)
  las_args += " -out {}".format(out_filename)
  val = las_parse.las_t
  if val != None and int(val) != las_num_threads_from_params:
    print("WARNING: command line las_t was given and is overriding the value t={} in the params file with {}".format(las_num_threads_from_params, val))
    las_args += " -t {}".format(val)
  else:
    las_args += " -t {}".format(las_num_threads_from_params)

  if not asserts_on():
    las_args += " -production"

  return las_args

def target_poly(verbose, smoothness_bound, path_to_target_dir, fb_filename, poly_filename, prs_used, fb, rat_fb,
                alg_prime_str_dict, max_num_relns, log_sieve_radius, sieving_poly, num_threads, disc, mod_gens,
                exit_early, special_q_filepath, is_special_q, disc_div):
  '''
  Writes the poly to a file and then calls las on it
  :param verbose:
  :param smoothness_bound:
  :param path_to_target_dir:
  :param fb_filename:
  :param rat_pr:
  :param fb:
  :param rat_fb:
  :param alg_prime_str_dict:
  :param max_num_relns:
  :param log_sieve_radius:
  :return:
  '''
  nf = fb[0].number_field()
  cado_quartic = CadoQuartic(sieving_poly)
  # sys.stderr.write("{}".format(prs_used) + os.linesep)
  # sys.stderr.write("log height {}\nshi bai size {}\nskew {}\nindex facs {}".format(log_height(sieving_poly), shi_bai_size(CadoQuartic(sieving_poly))) + os.linesep)
  if prs_used != [] and max(prs_used) > smoothness_bound:
    raise ValueError("either the ideal or the common factor of the sieving poly's coeffs are too large")
  poly_filepath = path_to_target_dir + poly_filename
  write_cado(cado_quartic, poly_filepath, False, None, None, None, disc_div)

  args = " -poly {}".format(poly_filepath)
  args += " -lim {}".format(smoothness_bound)
  fb_file_path = path_to_target_dir+fb_filename
  args += " -out {}".format(fb_file_path)
  args += " -t {}".format(num_threads)
  args += " -side 0"
  fb_cmd = cado_makefb_binary + args
  fb_args = shlex.split(fb_cmd)
  if verbose:
    print("Running makefb for target pr idl {}".format(fb_filename))
    print(fb_cmd)

  my_env = set_ld_lib_path()
  if my_env!=None:
    subprocess.run(fb_args, capture_output=True, env=my_env)
  else:
    subprocess.run(fb_args, capture_output=True)

  las_args = " -poly {}".format(poly_filepath)
  #according to the docs this won't do anything
  las_args += " -seed 0"
  las_args += " -lim0 {}".format(smoothness_bound)
  las_args += " -lim1 {}".format(smoothness_bound)
  las_args += " -fb0 {}".format(fb_file_path)
  las_args += " -fb1 {}".format(fb_file_path)

  las_args += " -I {}".format(log_sieve_radius)
  #q0 and q1 are overridden by special_q_filepath if that is not None
  las_args += " -q0 {}".format(10) #if this is too small for Cado then we may get ```code BUG() : condition (((p)->_mp_size != 0) & (static_cast<int> ((p)->_mp_d[0])))```
  las_args += " -q1 {}".format(60)
  if special_q_filepath != None:
    las_args += " -todo {}".format(special_q_filepath)
  lpb_bound = ceil(RR(smoothness_bound).log(2))
  sys.stderr.write("target I {}".format(log_sieve_radius) + os.linesep)
  las_args += " -lpb0 {}".format(lpb_bound+1)
  las_args += " -lpb1 {}".format(lpb_bound+1)
  las_args += " -mfb0 {}".format(lpb_bound)
  las_args += " -mfb1 {}".format(lpb_bound)
  las_args += " -t {}".format(num_threads)

  # from las docs "-exit-early: once a relation has been found, go to next special-q (value==1), or exit (value==2)"
  if exit_early != None:
    las_args += " -exit-early {}".format(exit_early)
  if is_prod():
    las_args += " -production"

  target_cmd = cado_las_binary + las_args
  if verbose:
    print("Running las for target pr idl {}".format(fb_filename))
    print(target_cmd)
  sys.stderr.write(target_cmd + os.linesep)

  output_start = time.process_time()
  stdout = subprocess.check_output(target_cmd, shell=True)
  lines = stdout.decode('utf-8').split("\n")
  #TODO add timer for special q
  timers.filter_target_new_polys_process_las_output += time.process_time()-output_start
  if lines[-5].startswith('# Total cpu time'):
    las_cpu_time = find_float(lines[-5])
    if is_special_q:
      timers.filter_target_special_q_las += las_cpu_time
    else:
      timers.filter_target_new_polys_las_time += las_cpu_time
  if 'reports' in lines[-2]:
    sys.stderr.write(lines[-2] + os.linesep)

  relns_count=0
  for line in lines:
    if line.startswith("#"):
      continue
    relns_count+=1
    if relns_count > max_num_relns:
      break

  #TODO min not needed ?
  num_relns = min(max_num_relns, relns_count)
  relns_mat = Matrix(ZZ, nrows=num_relns, ncols=len(fb), sparse=True)
  elems_mat = Matrix(nf, nrows=num_relns, ncols=1)

  reln_ind=0
  for i in range(len(lines)-1):
    line=lines[i]
    if line.startswith("#"):
      continue
    pairs_str, reln_primes0_str, reln_primes1_str = line.split(":", 2)
    cc = tuple([int(_) for _ in pairs_str.split(",", 2)])
    #TODO could add a duplicate check
    if cc[0]==0 and cc[1]==0:
      continue
    norm_prime_facs = [int(_, 16) for _ in reln_primes1_str.split(",")]
    max_pr=max(norm_prime_facs)
    if max_pr > smoothness_bound:
      print("found large prime. Logs are {}>{}. Skipping".format(RR(max_pr).log(2).n(30), RR(smoothness_bound).log(2).n(30)))
      continue
    norm_prime_facs += prs_used
    # sort the pr facs
    # norm_prime_facs = list(dict.fromkeys(norm_prime_facs))
    norm_prime_facs = sorted(set(norm_prime_facs))
    elem = cc[0]*mod_gens[0] + cc[1]*mod_gens[1]
    if elem==nf(0):
      continue
    elems_mat[reln_ind,0] = elem
    sieve_idl=nf.ideal(elem)
    if norm_prime_facs == []:
      print("Warning: found unit when reading relations at line index {}".format(i))
      pass
    else:
      for i in range(len(norm_prime_facs)):
        pr = norm_prime_facs[i]
        if pr == rat_fb[-1]:
          rng_max = len(fb)
        else:
          # get the index of the next prime so we can determine the range of prime indexs to look at
          next_pr = next_prime(pr, proof=False)
          rng_max = alg_prime_str_dict[next_pr]
        vals_at_primes_above_p = []
        for j in range(alg_prime_str_dict[pr], rng_max):
          alg_pr = fb[j]
          start_valuation_time_sp=time.process_time()
          vv = sieve_idl.valuation(alg_pr)
          #TODO make own timer
          # timers.filter_valuation += time.process_time() - start_valuation_time_sp
          vals_at_primes_above_p += [vv]
          if vv != 0:
            relns_mat[reln_ind, j] = vv
    reln_ind +=1
    if reln_ind == num_relns:
      break
  relns_mat = copy(relns_mat[:reln_ind])
  elems_mat = copy(elems_mat[:reln_ind])
  return relns_mat, elems_mat

def compute_poly_of_idl_and_target_idl(verbose, smoothness_bound, prescribed_div_idl, path_to_target_dir,
                                       filename_without_extension, prescribed_rat_pr, fb, rat_fb, alg_prime_str_dict,
                                       max_num_relns, log_sieve_radius, disc, bounded_fb, second_gen_residues,
                                       nrms_bounded_fb, prs_bounded_fb, disc_div, pol_ind):
  from src.stage_init import lattice_sieving_bound
  if len(bounded_fb)!=len(second_gen_residues):
    raise ValueError("second_gen_residues has the wrong number of elems")

  prec_for_knapsack=1
  binary=False
  max_num_primes=None
  # print("Target knapsack settings: prec {}\tnum norms to choose from {}".format(prec_for_knapsack, len(nrms_bounded_fb)))
  nf=fb[0].number_field()
  bnd= lattice_sieving_bound(nf.polynomial())

  #Oct 17, testing crt_elem
  if prescribed_div_idl.residue_class_degree() > 1 or gcd(prescribed_div_idl.norm(), pol_ind)>1:
    filter_target_new_polys_knapsack_start=time.process_time()
    [prescribed_idl, idl_norm, prs_used, exp, rat_prs_used]= construct_idl_with_prescribed_divisor_and_close_norm(
      bounded_fb, nrms_bounded_fb, prs_bounded_fb, bnd, prec_for_knapsack, prescribed_div_idl, binary, max_num_primes)
    timers.filter_target_new_polys_knapsack+=time.process_time()-filter_target_new_polys_knapsack_start

    prec=200
    sieving_poly, common_factor_of_poly_coeffs, mod_gens, idl_norm = lll_reduced_sieve_poly(prescribed_idl, prec)
    lll_elem=mod_gens[1]
    poly=nf.polynomial()
    dg=poly.degree()
    expected_T2_size_of_lll_reduced_elem=1.0/(dg-1) * (d_log2(prescribed_idl.norm())+0.5*d_log2(disc)) # b/c what you get from taking log of (N(I) * D^0.5)^(1/(d-2)) is this. The power (d-2) comes from the worst case which is that the field has a quadratic subfield. The 2nd elem in the LLL reduced basis is usually in the subfield. So we need to take the 3rd elem
    sys.stderr.write("For Pr idl of degree > 1. size of lll_elem is {} vs what 'might be expected' from the LLL bound: the log of (N(I)D^0.5)^(1/(d-1)={}".format(0.5*log_T2(lll_elem, 15).n(15), expected_T2_size_of_lll_reduced_elem) + os.linesep)
    sys.stderr.write("For Pr idl of degree > 1. height of poly is {}".format(log_height(sieving_poly)) + os.linesep)
  else:
    filter_target_new_polys_knapsack_start=time.process_time()
    [prescribed_idl, idl_norm, prs_used, exp, rat_prs_used]= construct_idl_with_prescribed_divisor_and_close_norm(
      bounded_fb, nrms_bounded_fb, prs_bounded_fb, bnd, prec_for_knapsack, prescribed_div_idl, binary, max_num_primes)
    timers.filter_target_new_polys_knapsack+=time.process_time()-filter_target_new_polys_knapsack_start
    second_gen_residues_used = [second_gen_residues[_] for _ in range(len(second_gen_residues)) if exp[_]!=0]

    try:
      sieving_poly, common_factor_of_poly_coeffs, mod_gens = crt_sieve_poly(prs_used, second_gen_residues_used)
    except ZeroDivisionError:
      print("Warning: crt of Dedekind elems failed. Probably b/c knapsack selected the same prime. Falling back to using reduced sieve poly")
      sieving_poly, common_factor_of_poly_coeffs, mod_gens, _ = lll_reduced_sieve_poly(prescribed_idl, prec)

  if common_factor_of_poly_coeffs != 1:
    rat_prs_used += prime_factors(common_factor_of_poly_coeffs)

  poly_filename = filename_without_extension + ".poly"
  fb_filename = filename_without_extension + ".roots0.gz"

  num_threads=5
  exit_early=2
  # print("Target las settings: I={}\tnum_threads {}\texit_early {}".format(log_sieve_radius, num_threads, exit_early))
  special_q_filepath=None
  is_special_q=False
  relns_mat, elems_mat = target_poly(verbose, smoothness_bound, path_to_target_dir, fb_filename, poly_filename, rat_prs_used, fb, rat_fb, alg_prime_str_dict, max_num_relns, log_sieve_radius, sieving_poly, num_threads, disc, mod_gens, exit_early, special_q_filepath, is_special_q, disc_div)
  return relns_mat, elems_mat, sieving_poly

def las_cmd(path_to_workdir, path_to_relns_dir, params_spec_str, fb_proc, verbose, smoothness_bound, expected_cg_generator_bound):
  if fb_proc != None:
    print("Waiting for fb")
    # fb_proc.wait()

  print("fb proc finished. Running siever")
  # get the relations dir path
  las_sieve_args= get_las_args(path_to_workdir, path_to_relns_dir, params_spec_str, verbose, smoothness_bound, expected_cg_generator_bound)
  sieve_cmd = cado_las_binary + las_sieve_args
  return sieve_cmd

def run_batch_las(path_to_workdir, path_to_relns_dir, smoothness_bound, job_name, fb_procs=None, filenames=None, expected_cg_generator_bound=None):
  '''
  Two use cases:
  1) This function is run immediately after the factor base is computed.
  It needs the process ids so it can wait for them to finish. The filenames are also provided
  2) This function is run separately from factor base construction.
  Can assume all the fb processes have finished.
  This function finds the filenames itself
  :param path_to_workdir:
  :param path_to_relns_dir:
  :param fb_procs: a list of processes running the fb script
  :param filenames:
  :return:
  '''
  if filenames == None:
    filenames = []
    for filename in os.listdir(path_to_workdir):
      poly_str = '.poly'
      if not filename.endswith(poly_str):
        continue
      spec_str = filename.rstrip(poly_str) + "_params"
      filenames.append(spec_str)
    fb_procs = [None]*len(filenames)
  if filenames == None or len(filenames) == 0:
    raise ValueError("Error getting filenames")
  else:
    print("Starting sieving on {} polynomial files".format(filenames))

  #Case 1: if the dir does not exist then create it
  #Case 2: if the dir exists then empty it
  procs = []
  mk_proc = subprocess.run(["mkdir", path_to_relns_dir])
  procs.append(mk_proc)
  for file in os.listdir(path_to_relns_dir):
    print("deleting file {}".format(file))
    del_path = path_to_relns_dir + file
    proc = subprocess.run(["rm", del_path])
    procs.append(proc)

  cmds=[]
  if len(filenames)==0:
    raise ValueError("sieving filenames tuple is empty")
  las_sieve_args= get_las_args(path_to_workdir, path_to_relns_dir, filenames[0], False, smoothness_bound, expected_cg_generator_bound)
  sys.stderr.write("Las running with args {}".format(las_sieve_args) + os.linesep)
  for i in range(len(filenames)):
    filename = filenames[i]
    if filename == path_to_relns_dir:
      continue
    verbose=True
    cmds.append(las_cmd(path_to_workdir, path_to_relns_dir, filenames[i], None, verbose, smoothness_bound, expected_cg_generator_bound))

  pool_size=10
  print("Running las with pool size={}".format(pool_size))
  user_sieve_start=time.perf_counter()
  limited_concurrent_commands(cmds, pool_size)
  timers.user_sieve=time.perf_counter() - user_sieve_start
  return []

