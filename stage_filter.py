import os
import shutil
import sys
import time

from tqdm import *
from src.stage_filter_internal import *
from src.stage_linalg_internal import submat_full_rank_at_columns
from src.stage_sieve import *
from src.magma_utils import *
from src.constants import *

class FilterState():
  global task_filter_concat
  global task_filter_load_mat
  global task_filter_dup1
  global task_filter_merge
  global task_filter_dup2
  global task_filter_triv
  global task_filter_singleton
  global task_filter_zero
  global task_filter_write
  global task_filter_targeted_relations
  global task_filter_save_mat

  def __init__(self, log_level):
    self.states = [task_filter_concat, task_filter_load_mat, task_filter_dup1, task_filter_merge, task_filter_dup2, task_filter_triv, task_filter_singleton, task_filter_zero, task_filter_write, task_filter_targeted_relations, task_filter_save_mat]
    self.log_level = log_level
    # set the initial state
    self.cur_state = task_filter_concat

  def update_filter_state(self, new_state):
    if not new_state in self.states:
      raise ValueError('not a valid state')
    self.cur_state = new_state

  def print_all_states(self):
    for state in self.states:
      self.update_filter_state(state)
      print(self.log_state())

  def log_state(self):
    return "Task " + self.cur_state

def unzip_and_concat(path_to_params_dir, path_to_reln_dir, params_spec_str, verbose):
  '''
  :param path_to_params_dir:
  :param params_spec_str:
  :param verbose:
  :return:
  '''
  # Make this function idempotent

  if verbose:
    print("Unzipping relations files in {}".format(path_to_reln_dir))

  # find each file in dir with a .gz suffix
  procs = []
  #TODO I need to sort these files
  files = os.listdir(path_to_reln_dir)
  #the ordering doesnt really matter but be aware that this will sort the strings ['i0','i1','i2',...,'i10',...] as ['i0', 'i10', 'i1', 'i2',..]
  sorted_files = sorted(files)
  non_gz_filenames=[]
  for filename in sorted_files:
    if not filename.endswith('.gz'):
      continue
    if 'tmp' in filename:
      print("WARNING (Unzip): filename={} contains 'tmp' which probably means that finding the roots failed".format(filename))
      continue
    root, ext = os.path.splitext(filename)
    non_gz_filenames.append(root)
    file_path = path_to_reln_dir + filename
    if verbose:
      print("Unzipping {}".format(file_path))
    proc = subprocess.call(["gunzip", file_path]) # use call b/c it waits for the subprocess to finish
    procs.append(proc)

  idl_datas = dict()
  for filename in os.listdir(path_to_params_dir):
    if not filename.endswith('.poly'):
      continue
    poly_name = filename.split(".poly",1)[0]
    with open(path_to_params_dir + filename, 'r') as poly_fp:
      lines = poly_fp.readlines()

      #TODO I shouldnt be parsing here. I should just put the whole line into idl_data
      idl_data = [0]*4
      for line in lines:
        if line.startswith(idl_least_int_str):
          idl_data[0] = line
        elif line.startswith(idl_poly_gen_str):
          idl_data[1] = line
        elif line.startswith(idl_coeff_gcd_str):
          idl_data[2] = line
        elif line.startswith(idl_fb_exp_vec_str):
          idl_data[3] = line

      idl_datas[poly_name] = idl_data

  # concat all the files
  input_files = []
  # for this file split the name at '_params' and append the [least_int, root] that was found in the poly file
  sorted_files_no_gz = []
  print("files in this directory are {}".format(sorted_files))
  contains_gz = False
  for _ in sorted_files:
    if not _.endswith('.gz'):
      raise ValueError("upload directory not clean. Path {} contains {}. Delete it and restart".format(path_to_reln_dir, _))
    else:
      contains_gz = True
    sorted_files_no_gz.append(_[:-3])
  if not contains_gz:
    raise ValueError("upload dir does not contain a gz file")

  total_uniques=0
  sorted_files_no_gz_dedup=[]
  for filename in sorted_files_no_gz:
    file_path = path_to_reln_dir + filename
    new_name=filename+"_dedup"
    new_path=path_to_reln_dir+new_name
    sorted_files_no_gz_dedup.append(new_name)
    cmd_str="cat {} | sort | uniq > {}".format(file_path, new_path)

    cmd_str_count="grep -v '^#' {} | sort | uniq | wc -l".format(file_path) # to count how many unique there are
    num_uniques_result=subprocess.check_output([cmd_str_count, file_path], shell=True)
    num_uniques=int(num_uniques_result)
    total_uniques+=num_uniques
    if not is_prod():
      cmd_str_count_total="grep -v '^#' {} | wc -l".format(file_path) # to count how many unique there are
      num_total_result=subprocess.check_output([cmd_str_count_total, file_path], shell=True)
      num_total=int(num_total_result)
      print("Deduplicating reln file {}\nFound this many unique relations: {} out of {}".format(filename, num_uniques, num_total))
    subprocess.call([cmd_str, file_path], shell=True)

  print("Concatting these files (in this order) {}".format(sorted_files_no_gz_dedup))
  for filename in sorted_files_no_gz_dedup:
    if filename.endswith('.gz'):
      raise ValueError("Gunzip failed on file {}".format(filename))
    if 'tmp' in filename:
      print("WARNING (Concat): filename={} contains 'tmp' which probably means that finding the roots failed".format(filename))
      continue
    params_str = "_params"
    if not params_str in filename:
      continue
    file_path = path_to_reln_dir + filename
    with open(file_path, "rb") as reln_fp_check:
      # throws Unicode error if invalid
      utf8_decode(reln_fp_check.read(), filename)

      reln_fp_check.seek(0)
      last_line = reln_fp_check.readlines()[-1].decode("utf-8")
      if idl_poly_gen_str in last_line:
        raise ValueError("Found {} on last line of file {}.\nAborting because its likely that some files in this directory are invalid".format(idl_poly_gen_str, filename))

    poly_name = filename.split(params_str, 1)[0]
    str = ""
    with open(file_path, "ab") as reln_fp:
      str += "{}".format( idl_datas[poly_name][0])
      str += "{}".format( idl_datas[poly_name][1])
      str += "{}".format( idl_datas[poly_name][2])
      str += "{}".format( idl_datas[poly_name][3])
      reln_fp.write(str.encode("utf-8"))

    input_files.append(file_path)

  dest_filename = "all_reln"
  dest_path = path_to_reln_dir + dest_filename
  with open(dest_path,'wb') as wfd:
    for f in input_files:
      with open(f,'rb') as fd:
        shutil.copyfileobj(fd, wfd)
  return total_uniques

def filter(log_level, poly, smooth_bound, path_to_reln_file, write_rel_path, sieving_stats_rel_path,
           filtering_ratio_of_relns_to_load, sieve_radius, timing_test_for_full_rank, disc, target_path, disc_div,
           verbose, sanity_check, time_per_relation_timing, total_reln_available):
  '''
  :param log_level:
  :param poly:
  :param smooth_bound:
  :param path_to_reln_file: the absolute path to the relations file
  :param write_rel_path:
  :return:
  '''
  global task_filter_concat
  global task_filter_load_mat
  global task_filter_dup1
  global task_filter_merge
  global task_filter_dup2
  global task_filter_triv
  global task_filter_singleton
  global task_filter_zero
  global task_filter_write
  global task_filter_targeted_relations
  global task_filter_save_mat
  print("\nStage: Filter")
  if total_reln_available == 0:
    raise ValueError("0 reports in all_reln file {}. Has the file been written yet ?".format(path_to_relns_file))

  ### TASK: load
  filt_state = FilterState(log_level)
  filt_state.update_filter_state(task_filter_load_mat)
  print(filt_state.log_state())

  start_filter_time = time.process_time()
  start_setup_time = time.process_time()
  start_load_time = time.process_time()

  path_to_relns_file = path_to_reln_file
  relns_mat_path = write_rel_path + "_relns"
  sp_rat_reln_mat_path = write_rel_path + "_sp_rat_relns"
  flint_relns_mat_path_initial = relns_mat_path + "_initial_flint"
  flint_relns_mat_path_merged = relns_mat_path + "_merged_flint"
  flint_relns_mat_path_final = relns_mat_path + "_final_flint"
  elems_mat_path = write_rel_path + "_elems"
  num_rels = 0

  dg = poly.degree()
  nf = NumberField(poly, names='YY')
  [_, _, fb, distinct_rat_prs_of_fb, alg_prime_str_dict, _, _] = fac_base_data(nf, smooth_bound)
  rat_pr_fb = [_.smallest_integer() for _ in fb]
  norm_bnd_fb = fb_of_idls_of_bounded_norm_bounded(fb, smooth_bound)
  seen = set()

  print("Starting to load self-initialization relations")
  timers.las_time = 0
  if not os.path.isfile(path_to_relns_file):
    raise ValueError("all_relns file does not exist at {}".format(path_to_relns_file))
  with open(path_to_relns_file, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      if line.startswith('# Total cpu time'):
        las_cpu_time = find_float(line)
        timers.las_time += las_cpu_time
      # count relns using the lines of the form '# Total <X> reports..'
      words = line.split(" ")
  sys.stderr.write("all_reln took {:.3f}".format(timers.las_time) + os.linesep)
  sys.stderr.write("{} {}".format(sieve_log_str_num_relns, total_reln_available) + os.linesep)
  print("all_reln took {:.3f}".format(timers.las_time))
  print("{} {}".format(sieve_log_str_num_relns, total_reln_available))

  if time_per_relation_timing:
    min_num_relns=0.9*len(norm_bnd_fb)
    c=1
    if disc < ZZ(2)**16:
      c=8
    if disc < ZZ(2)**32:
      c=4
    max_num_relns=c*len(fb)

    if total_reln_available < min_num_relns:
      raise ValueError("too few relations collected for time_per_relation_timing. Found {} and should be more than 0.9*len(norm_bnd_fb)={}".format(total_reln_available, min_num_relns))
    elif max_num_relns < total_reln_available:
      raise ValueError("too many relations collected for time_per_relation_timing. Found {} and should be less than {}*len(fb)={}".format(total_reln_available, c, max_num_relns))

  if time_per_relation_timing:
    return

  sys.stderr.write("Size of FB is {} vs number of prime ideals above rat prs {}".format(len(norm_bnd_fb), len(fb)) + os.linesep)
  num_relns_requested = floor(filtering_ratio_of_relns_to_load*len(norm_bnd_fb))
  if num_relns_requested > total_reln_available:
    raise ValueError("Requested more relations {}x|FB|={} than were found {}. "
                     "Aborting b/c we should always get as many as are requested".format(filtering_ratio_of_relns_to_load, num_relns_requested, total_reln_available))
  else:
    # get a random 0-1 tuple that's guaranteed to be at least as long as the number of relations
    initial_reln_inds = get_random_zero_one_list_and_increment_seed(total_reln_available, num_relns_requested)
    sys.stderr.write("\nWARNING Reading a maximum of {} relns. This is {}% of the available relations.\n".format(num_relns_requested, (float(num_relns_requested)/total_reln_available)*100) + os.linesep)
    print("\nWARNING Reading a maximum of {} relns. This is {}% of the available relations.\n".format(num_relns_requested, (float(num_relns_requested)/total_reln_available)*100))

  YY = nf.gen()
  sp_reln_mat = Matrix(ZZ, nrows=num_relns_requested, ncols=len(fb), sparse=True)
  lp_reln_mat = Matrix(ZZ, nrows=sp_reln_mat.nrows(), ncols=dg+1, sparse=True)
  elems_mat = Matrix(nf, nrows=sp_reln_mat.nrows(), ncols=1, sparse=True)

  print("Num of rows={} and cols={}".format(sp_reln_mat.nrows(), sp_reln_mat.ncols()))
  #TODO if nrows < ncols shouldnt we just abort here

  loaded_reln_ind = 0
  reln_ind = 0
  count_relns_with_too_many_large_primes = 0
  count_dups = 0

  with open(path_to_relns_file, 'r') as cado_file:
    first_elem = None
    idl_least_ints = []
    idl_poly_gens = []
    idl_coeff_gcds = []
    idl_fb_exp_vecs = []
    for line in cado_file.readlines():
      if line.startswith(idl_least_int_str):
        if "[" and "]" in line:
          #TODO the leftmost split does not seem to be required
          no_brackets_str = line.split(idl_least_int_str, 1)[1].split("[", 1)[1].split("]", 1)[0]
          idl_least_ints.append(str_to_rat_list(no_brackets_str))
        else:
          idl_least_ints.append(int(line.split(idl_least_int_str + ": ",1)[1]))
      elif line.startswith(idl_poly_gen_str):
        if "[" and "]" in line:
          #TODO the leftmost split does not seem to be required
          no_brackets_str = line.split(idl_poly_gen_str, 1)[1].split("[", 1)[1].split("]", 1)[0]
          idl_poly_gens.append(str_to_rat_list(no_brackets_str))
        else:
          idl_poly_gens.append(int(line.split(idl_poly_gen_str + ": ",1)[1]))
      elif line.startswith(idl_coeff_gcd_str):
        idl_coeff_gcds.append(int(line.split(idl_coeff_gcd_str + ": ",1)[1]))
      elif line.startswith(idl_fb_exp_vec_str):
        no_brackets_str = line.split(idl_fb_exp_vec_str, 1)[1].split("[", 1)[1].split("]", 1)[0]
        idl_fb_exp_vecs.append(str_pairs_list_to_int_list(no_brackets_str))
      else:
        continue
    if len(idl_fb_exp_vecs) != len(idl_coeff_gcds):
      raise ValueError("idl_fb_exp_vecs has length {}. idl_coeff_gcds has length {}".format(len(idl_fb_exp_vecs), len(idl_coeff_gcds)))

    cado_file.seek(0)

    # Set data for the first sieving ideal
    idl_idx = 0
    initial_idl, facs_of_com, initial_idl_rat_prs = get_initial_ideal_and_common_factor_of_sieve_poly_coeffs(fb, idl_coeff_gcds, idl_fb_exp_vecs, idl_idx, nf, rat_pr_fb)
    initial_idl_index_pairs = idl_fb_exp_vecs[idl_idx]
    idl_one = nf.ideal(1)

    for line in cado_file.readlines():
      if line.startswith('#') or line.startswith(idl_poly_gen_str) or line.startswith(idl_coeff_gcd_str) or line.startswith(idl_fb_exp_vec_str):
        # Case: skip. This line does not trigger a change of ideal and is not a relation in an ideal
        continue
      elif line.startswith(idl_least_int_str):
        # Case: switch to a new ideal
        idl_idx += 1
        if idl_idx == len(idl_least_ints):
          # Take this to mean that we're at the end of the file
          break
        initial_idl, facs_of_com, initial_idl_rat_prs = get_initial_ideal_and_common_factor_of_sieve_poly_coeffs(fb, idl_coeff_gcds, idl_fb_exp_vecs, idl_idx, nf, rat_pr_fb)
        initial_idl_index_pairs = idl_fb_exp_vecs[idl_idx]
        continue

      # initial_reln_inds is an indicator vector. 1 means this relation should be included
      if initial_reln_inds[loaded_reln_ind] == 0:
        loaded_reln_ind += 1
        continue
      loaded_reln_ind += 1

      if len(line.split(":", 2)) == 3:
        pairs_str, reln_primes0_str, reln_primes1_str = line.split(":", 2)
      else:
        raise ValueError("wrong number of entries in line {}".format(line))

      cc = tuple([int(_) for _ in pairs_str.split(",", 2)])

      # this algo to updates only those entries of the matrix that match the norms found by line sieving
      # up to this point I had just been doing trial division
      norm_prime_facs = [int(_, 16) for _ in reln_primes1_str.split(",")]

      # if idl_least_ints is an integer then the exact norm is elem_norm = com * idl_least_ints[idl_idx] * ZZ(prod_all_elems(norm_prime_facs))
      elem_norm_divisor=ZZ(prod_all_elems(norm_prime_facs))

      try:
        elem = cc[0]*nf(idl_least_ints[idl_idx]) + cc[1]*nf(idl_poly_gens[idl_idx])
      except TypeError:
        elem = cc[0]*idl_least_ints[idl_idx] + cc[1]*(YY - idl_poly_gens[idl_idx])

      nrm=elem.norm()
      missed_norm_factor = nrm/gcd(nrm, elem_norm_divisor)
      missed_prs = prime_facs(missed_norm_factor)

      # las sometimes seems to miss small primes dividing the leading coefficient of the sieving polynomial so we include initial_idl_rat_prs here even though it shouldn't be required
      norm_prime_facs = sorted(list(set(norm_prime_facs + facs_of_com + initial_idl_rat_prs + missed_prs)))

      coords_on_pow_basis=tuple([int(_) for _ in list(elem)])
      if coords_on_pow_basis in seen:
        count_dups += 1
        continue
      else:
        seen.add(coords_on_pow_basis)

      if first_elem == None:
        first_elem = copy(elem)
      elems_mat[reln_ind,0] = elem

      # step 1: fill in lp mat (valuations)
      small_primes = norm_prime_facs

      last_pr = small_primes[-1]
      if last_pr > smooth_bound:
        small_primes.remove(last_pr)
        norm_prime_facs = [last_pr]
      else:
        norm_prime_facs = []
      new_last_pr = small_primes[-1] # last pr was removed
      if new_last_pr > smooth_bound:
        print("There's more than 1 large prime. Right now I dont want to handle that. Skip")
        count_relns_with_too_many_large_primes += 1
        continue

      #constructing the principal ideal generated by elem is needed to compute the valuation efficiently. It is not used for anything else so we add the ideal construction time to valuation
      idl_convert_start_time = time.process_time()
      idl = nf.ideal(elem)
      timers.filter_valuation += time.process_time()-idl_convert_start_time

      # only a single large prime is allowed
      if len(norm_prime_facs) == 1:
        lp = ZZ(norm_prime_facs[0]) # the large prime
        lp_reln_mat[reln_ind, 0] = lp
        large_prime_facs=prime_idl_facs(nf, norm_prime_facs[0])
        valn_sum = 0
        if not is_prod():
          nrm_valn = elem_norm_divisor.valuation(lp)
        for i in range(0,len(large_prime_facs)):
          alg_pr = large_prime_facs[i]
          start_valuation_time_lp=time.process_time()
          vv = idl.valuation(alg_pr)
          timers.filter_valuation += time.process_time() - start_valuation_time_lp
          if vv != 0:
            lp_reln_mat[reln_ind, i+1] += vv

      # step 2: fill in small primes mat and elems mat
      # step 2a: fill in valuations for rational primes up to and including the 2nd last rational prime

      for i in range(len(small_primes)):
        pr = small_primes[i]
        if pr == distinct_rat_prs_of_fb[-1]:
          rng_max = len(fb)
        else:
          # get the index of the next prime so we can determine the range of prime indexs to look at
          next_pr = next_prime(pr,proof=False)
          rng_max = alg_prime_str_dict[next_pr]

        vals_at_primes_above_p = []
        #TODO move the belabas valuation code. I need it to work but I dont want to use it in production
        if not is_prod():
          nrm_valn = elem_norm_divisor.valuation(pr)
          rat_pr_ind = distinct_rat_prs_of_fb.index(pr)
          # sp_rat_reln_mat[reln_ind, rat_pr_ind] = nrm_valn

        for j in range(alg_prime_str_dict[pr], rng_max):
          alg_pr = fb[j]
          start_valuation_time_sp=time.process_time()
          vv = idl.valuation(alg_pr) #
          timers.filter_valuation += time.process_time() - start_valuation_time_sp
          vals_at_primes_above_p += [vv]
          if vv != 0:
            sp_reln_mat[reln_ind, j] += vv
        if log_level >= 3:
          print("val of ideal at primes above p={}".format(vals_at_primes_above_p))
      assert assert_row_multiplies_to_correct_idl(reln_ind, sp_reln_mat, elems_mat, nf, fb, lp_reln_mat)
      reln_ind +=1

  # Validation
  if (reln_ind + count_relns_with_too_many_large_primes + count_dups) != num_relns_requested:
    raise ValueError("Wrong number of relns")

  print("Relations loaded")
  if count_relns_with_too_many_large_primes > 0:
    sys.stderr.write("{} relns had too many large primes. Ratio of rows with too many large primes was {:.2f}".format(count_relns_with_too_many_large_primes, float(count_relns_with_too_many_large_primes)/reln_ind) + os.linesep)
  print("Elem set function excluded {} dups".format(count_dups))
  print("Shrinking matrix to fit actual number of relations {:.2f}".format(reln_ind/float(num_relns_requested)))
  copy_start=time.process_time()
  sp_reln_mat = copy(sp_reln_mat[:reln_ind])
  lp_reln_mat = copy(lp_reln_mat[:reln_ind])
  # sp_rat_reln_mat = copy(sp_rat_reln_mat[:reln_ind])
  elems_mat = copy(elems_mat[:reln_ind])
  print("time to copy {}".format(time.process_time() - copy_start))

  if asserts_on():
    print("asserts after load")
    assert assert_rows_multiply_to_correct_ideal(sp_reln_mat, elems_mat, nf, fb, lp_reln_mat)

  timers.filter_load_mat = time.process_time() - start_load_time

  ###TASK: Trivial Relations
  filt_state.update_filter_state(task_filter_triv)
  print(filt_state.log_state())
  trivial_relns_start_time = time.process_time()

  triv_mat, triv_elems_mat, rat_prs_without_triv_reln, triv_rat_pr_mat = triv_rels(fb, alg_prime_str_dict, smooth_bound)
  sp_reln_mat = triv_mat.stack(sp_reln_mat)
  elems_mat = triv_elems_mat.stack(elems_mat)
  lp_zm = zero_matrix(triv_mat.nrows(), lp_reln_mat.ncols())
  lp_reln_mat = lp_zm.stack(lp_reln_mat)
  assert lp_reln_mat.nrows() == sp_reln_mat.nrows()
  print("Num triv relations {}".format(triv_mat.nrows()))
  timers.filter_trivial_relns = time.process_time() - trivial_relns_start_time

  if log_level >= 1:
    print("time spent on loading/parsing operations was {}\n'Merging' means any non-IO and non-parsing operation".format(time.time() - start_load_time))

  print("starting timer for merging")
  start_merge_time = time.process_time()

  if asserts_on():
    for row in sp_reln_mat:
      if list(row) == [0]*sp_reln_mat.ncols():
        raise ValueError("there are 0 rows")

  filt_state.update_filter_state(task_filter_merge) #Remove large primes
  print(filt_state.cur_state)
  # let A be a tuple of the row indexs in the mat to delete after 'merging' (1st unique occurence of a nonzero exponent) and B be a tuple of tuple of indexes of rows to merge them into
  visited, inds_of_rows_to_modify, inds_of_rows_to_del = inds_for_merge_one(lp_reln_mat)

  sp_reln_mat, elems_mat, sp_rat_reln_mat = merge_large_primes_and_slice_mat(visited, inds_of_rows_to_modify, inds_of_rows_to_del, sp_reln_mat, elems_mat, None)
  nrms_bounded_fb = [_.norm() for _ in norm_bnd_fb]
  prs_bounded_fb = [_.norm()**(_.residue_class_degree()**(-1)) for _ in norm_bnd_fb]

  print("Timers: merge time {}".format(timers.filter_merge))
  print(task_filter_merge)
  if log_level >= 1:
    print("Mat dimensions post-merge are {}x{}".format(sp_reln_mat.nrows(), sp_reln_mat.ncols()))
    print("Size of bounded FB is {}".format(len(nrms_bounded_fb)))
  if asserts_on():
    print("Asserts post-merge:")
  if not is_running_on_remote():
    print("Stopping, can't go further without Magma")
  # assert assert_rat_pr_mat_rows_multiply_to_correct_ideal(sp_rat_reln_mat, elems_mat, distinct_rat_prs_of_fb)
  assert assert_elems_smooth_and_nonzero(elems_mat, nf, smooth_bound)
  assert assert_rows_multiply_to_correct_ideal(sp_reln_mat, elems_mat, nf, fb)

  timers.filter_merge = time.process_time() - start_merge_time
  timers.filter_setup = time.process_time() - start_setup_time

  ###TASK target Spec-Q - Off
  # elems_mat, sp_reln_mat, spec_q_target_relns_mat, target_dir, target_log_sieve_radius = target_with_spec_q_on_def_poly(alg_prime_str_dict, disc, distinct_rat_prs_of_fb, elems_mat, fb, nf, smooth_bound, sp_reln_mat, target_path)

  ###TASK target New Polys
  start_target_time = time.process_time()
  target_dir = basepath_artifacts+target_path

  mag_relns = magma_mat_from_sparse_sage_mat(sp_reln_mat)
  print("\nStarting Target New Polys")

  cand_class_num_before = compute_class_num(mag_relns)
  print("cand class num {}".format(cand_class_num_before))

  target_verbose=False
  # a bigger sieve radius is needed for targetting
  log_sieve_radius_target=2 + floor(RR(sieve_radius).log(2))
  sys.stderr.write("radius is {}".format(log_sieve_radius_target) + os.linesep)
  rank_needed, target_newpolys_elems_mat, target_newpolys_relns_mat = target_nonpivot_cols(alg_prime_str_dict, disc, fb,
                                                                                           nf, fb, target_dir,
                                                                                           distinct_rat_prs_of_fb,
                                                                                           smooth_bound, sp_reln_mat,
                                                                                           target_verbose,
                                                                                           log_sieve_radius_target,
                                                                                           disc_div)

  stack_start = time.process_time()
  sp_reln_mat = sp_reln_mat.stack(target_newpolys_relns_mat)
  elems_mat = elems_mat.stack(target_newpolys_elems_mat)
  timers.filter_target_new_polys_stack += time.process_time()-stack_start

  if not is_prod():
    print("asserts on targeted relations")
    assert assert_rows_multiply_to_correct_ideal(target_newpolys_relns_mat, target_newpolys_elems_mat, nf, fb)
    assert assert_elems_smooth_and_nonzero(elems_mat, nf, smooth_bound)
    assert assert_rows_multiply_to_correct_ideal(sp_reln_mat, elems_mat, nf, fb)

  mag_relns = magma(sp_reln_mat)
  print("Mat dimensions post-target {}x{}".format(sp_reln_mat.nrows(), sp_reln_mat.ncols()))
  cand_class_num = compute_class_num(mag_relns)
  print("cand class num {}".format(cand_class_num))

  _, remaining_pr_inds = prs_to_target(fb, mag_relns, smooth_bound)

  print("Num prs to target after Target phase {}/{}\n".format(len(remaining_pr_inds), rank_needed))
  miss_ratio=len(remaining_pr_inds)/rank_needed
  if miss_ratio > 0.01:
    raise ValueError("Column miss ratio {} is too high".format(miss_ratio))

  timers.filter_target_new_polys_total = time.process_time()-start_target_time
  if timing_test_for_full_rank:
    timers.filter_total = time.process_time() - start_filter_time
    return

  ###TASK: Remove Singleton Cols
  remove_singletons=False
  if remove_singletons:
    filt_state.update_filter_state(task_filter_singleton)
    print(filt_state.log_state())
    # Removing singleton cols is done in 2 parts. First: remove the row corresponding to a singleton. Second: remove the singleton column
    start_singleton_cols_time = time.process_time()

    singleton_row_inds, non_singleton_row_inds = get_rows_of_singleton_columns(sp_reln_mat)

    # TODO verify this swapping process
    num_singleton = len(singleton_row_inds)
    for i in singleton_row_inds:
      if i < num_singleton:
        continue
      # now swap with the first non-singleton row and delete the index of that row from the non_singleton_row_inds
      ind_non_singleton = non_singleton_row_inds[0]
      sp_reln_mat.swap_rows(i, ind_non_singleton)

      # manual swap b/c calling swap_rows on elems_mat throws TypeError:object cannot be interpreted as an integer
      elems_mat[i], elems_mat[ind_non_singleton] = elems_mat[ind_non_singleton], elems_mat[i]

      non_singleton_row_inds.remove(non_singleton_row_inds[0])

    if asserts_on():
      # check the singleton rows of the new matrix are sorted
      singleton_inds, _ = get_rows_of_singleton_columns(sp_reln_mat)
      expected_singleton_inds = [_ for _ in range(len(singleton_row_inds))]
      if singleton_inds != expected_singleton_inds:
        raise ValueError("expected singleton indexs to be {}\n\ngot\n{}".format(expected_singleton_inds, singleton_inds))

    sp_reln_mat = copy(sp_reln_mat[num_singleton:])
    elems_mat = copy(elems_mat[num_singleton:])
    if log_level >= 1:
      if asserts_on():
        print("Doing assertions for task {}".format(task_filter_singleton))
        assert assert_rows_multiply_to_correct_ideal(sp_reln_mat, elems_mat, nf, fb)
      print("Found {} singleton columns".format(num_singleton))
      print("Mat dimensions post-singleton {}x{}".format(sp_reln_mat.nrows(), sp_reln_mat.ncols()))
    timers.filter_singleton_columns = time.process_time() - start_singleton_cols_time
    print("Timers: singleton time {}".format(timers.filter_singleton_columns))

  use_purge=False
  if use_purge:
    print("Purge step 1: remove largest prime ideal above each rational prime p")
    #If I've added a trivial relation to the matrix for this rational prime then don't use it to eliminate the largest prime. The HNF algorithm can decide the best way to use the trivial relation.
    for pr in rat_prs_without_triv_reln:
      if pr == distinct_rat_prs_of_fb[-1]:
        rng_max = len(fb)
      else:
        # get the index of the next prime so we can determine the range of prime indexs to look at
        rng_max = alg_prime_str_dict[next_prime(pr,proof=False)]
      # fb primes are sorted so rng_max-1 is the column index of a prime of largest residue class degree above the rational prime pr
      rat_ind = distinct_rat_prs_of_fb.index(pr)
      for rr in range(sp_reln_mat.nrows()):
        vv = sp_reln_mat[rr,rng_max-1]
        if vv != 0:
          # next loop subtracts vv*(trivial relation) from this row
          for cc in range(alg_prime_str_dict[pr], rng_max):
            sp_reln_mat[rr, cc] -= vv
          elems_mat[rr,0] /= pr**vv

  # print("Purge step 2: deleted {} relations".format(num_rows_to_del))
  # print("Weight after {}".format(mat_weight(sp_reln_mat)))
  if log_level >= 1:
    print("Asserts for purge and delete".format("column elimination"))
    # assert assert_rat_pr_mat_rows_multiply_to_correct_ideal(sp_rat_reln_mat, elems_mat, distinct_rat_prs_of_fb)
    assert assert_rows_multiply_to_correct_ideal(sp_reln_mat, elems_mat, nf, fb)

  filt_state.update_filter_state(task_filter_targeted_relations)
  print(filt_state.cur_state)

  remove_nonpivots_start_time = time.process_time()
  remove_nonpiv_rows=True
  if remove_nonpiv_rows:
    probable_piv_cols,_ = probable_nonpivot_columns_fastest_available(sp_reln_mat)
    row_inds_to_keep = list(range(sp_reln_mat.nrows()))
    for r in range(sp_reln_mat.nrows()):
      for c in probable_piv_cols:
        if sp_reln_mat[r,c] != 0:
          try:
            row_inds_to_keep.remove(r)
          except ValueError:
            continue
    num_rows_deleted = sp_reln_mat.nrows() - len(row_inds_to_keep)
    print("Removing {} nonpivot rows. Fraction of all rows is {:.2f}".format(num_rows_deleted, float(num_rows_deleted)/sp_reln_mat.nrows()))
    sp_reln_mat = sp_reln_mat.matrix_from_rows(row_inds_to_keep)
    elems_mat = elems_mat.matrix_from_rows(row_inds_to_keep)
  timers.filter_remove_nonpivot_rows = time.process_time() - remove_nonpivots_start_time

  eliminate_cols_start_time = time.process_time()

  #Removing singleton columns. This is cheap and valuable
  remove_cols= True
  if remove_cols:
    filt_state.update_filter_state(task_filter_zero)
    print(filt_state.log_state())
    indexs_of_zero_cols, nonzero_indexs = col_indexs_with_weight_below(sp_reln_mat, 1)
    sp_reln_mat = sp_reln_mat.matrix_from_columns(nonzero_indexs)
    fb = [fb[_] for _ in range(len(fb)) if _ in nonzero_indexs]

    if asserts_on():
      indexs_of_zero_cols, nonzero_indexs = col_indexs_with_weight_below(sp_reln_mat, 1)
      if len(indexs_of_zero_cols) > 0:
        raise ValueError("failed at removing 0 cols")
  cols_after_prune=sp_reln_mat.ncols()
  filter_dim_str="Mat dimensions after pruning FB are {}x{}".format(sp_reln_mat.nrows(), sp_reln_mat.ncols())
  sys.stderr.write(filter_dim_str + os.linesep)
  print(filter_dim_str)

  time_per_time=timers.las_time + timers.filter_target_new_polys_las_time
  print("{} {:.2f}".format(sieve_log_str_time_per_reln, 1000*time_per_time/sp_reln_mat.nrows()))

  print("{} {:.2f}".format(sieve_log_str_time_to_generate_initial_relns, timers.las_time))

  if log_level >= 1:
    if asserts_on():
      print("Doing assertions for task {}".format("column elimination"))
    # assert assert_rat_pr_mat_rows_multiply_to_correct_ideal(sp_rat_reln_mat, elems_mat, distinct_rat_prs_of_fb)
    assert assert_rows_multiply_to_correct_ideal(sp_reln_mat, elems_mat, nf, fb)
  timers.filter_eliminate_cols = time.process_time() - eliminate_cols_start_time

  print("Task: bring submat to top")
  num_rows=floor(1.2*sp_reln_mat.ncols())
  # get_full_rank_submat()

  timers.filter_total = time.process_time() - start_filter_time
  sp_reln_mat = sp_reln_mat[:num_rows]
  elems_mat = elems_mat[:num_rows]

  filter_saved_mat_dims_str="Starting Linalg on mat of dimensions {}x{}".format(sp_reln_mat.nrows(), sp_reln_mat.ncols())
  sys.stderr.write(filter_saved_mat_dims_str + os.linesep)
  return sp_reln_mat, elems_mat

def remove_dups(elems_mat, fb, filt_state, log_level, nf, sp_reln_mat, task_filter_dup2):
  filt_state.update_filter_state(task_filter_dup2)
  print(filt_state.log_state())
  dups_post_merge_start_time = time.process_time()
  inds_of_dups, inds_of_non_dups, units = inds_of_dup_rows_and_units(sp_reln_mat, elems_mat)
  sp_rat_reln_mat = None
  sp_reln_mat, lp_reln_mat, elems_mat, sp_rat_reln_mat = swap_rows_to_top(sp_reln_mat, None, elems_mat, sp_rat_reln_mat,
                                                                          inds_of_dups, inds_of_non_dups)
  print(
    "Found {} duplicate rows in the relations matrix and {} nontrivial units (units that are not 4th roots of unity)".format(
      len(inds_of_dups), len(units)))
  sp_reln_mat = copy(sp_reln_mat[len(inds_of_dups):])
  elems_mat = copy(elems_mat[len(inds_of_dups):])
  timers.filter_dups_post_merge = time.process_time()-dups_post_merge_start_time
  print("Timers: de-dup time {}".format(timers.filter_dups_post_merge))
  assert len(inds_of_dup_rows_and_units(sp_reln_mat, elems_mat)[0]) == 0
  if log_level >= 1 and asserts_on():
    print("Doing assertions for task {}".format(task_filter_dup2))
    assert assert_rows_multiply_to_correct_ideal(sp_reln_mat, elems_mat, nf, fb)
  return elems_mat, sp_reln_mat

def write_matrices(elems_mat, elems_mat_path, filt_state, flint_relns_mat_path_final, indexs_of_nonzero_cols, log_level,
                   relns_mat_path, sp_reln_mat, task_filter_write, write_rel_path):
  filt_state.update_filter_state(task_filter_write)
  print(filt_state.log_state())
  write_start_time = time.process_time()
  ############## Filtering Steps Over
  print("Dimensions at end of filter stage {}x{}".format(sp_reln_mat.nrows(), sp_reln_mat.ncols()))
  if asserts_on():
    print("writing final relns mat in flint format")
    write_int_mat(flint_relns_mat_path_final, sp_reln_mat, log_level, formatter=mat_to_flint_mat_str)
  # assert sp_rat_reln_mat.nrows() == sp_reln_mat.nrows()
  # write_int_mat(sp_rat_reln_mat_path, sp_rat_reln_mat, log_level)
  write_int_mat(relns_mat_path, sp_reln_mat, log_level)
  write_nf_elem_mat(elems_mat_path, elems_mat, log_level)
  # TODO give these a path
  write_fb_inds = basepath_artifacts+write_rel_path+"_fb_inds"
  with open(write_fb_inds, "w") as fp:
    fp.write(str(indexs_of_nonzero_cols))
  timers.filter_write_mats = time.process_time()-write_start_time

def set_up_mat_for_Vollmer_approach(elems_mat, full_rank_submat_start, row_inds_used_submat0, sp_reln_mat):
  other_inds = sorted(list(set(range(sp_reln_mat.nrows()))-set(row_inds_used_submat0)))
  inds = list(row_inds_used_submat0)+other_inds  # inds is a permutation of the row indices of sp_reln_mat
  sp_reln_mat = sp_reln_mat.matrix_from_rows(inds)
  elems_mat = elems_mat.matrix_from_rows(inds)
  return elems_mat, sp_reln_mat

def target_with_spec_q_on_def_poly(alg_prime_str_dict, disc, distinct_rat_prs_of_fb, elems_mat, fb, nf, smooth_bound,
                                   sp_reln_mat, target_path):
  print("Starting Special-q Target")
  target_dir = basepath_artifacts+target_path
  mk_proc = subprocess.run(["mkdir", "-p", target_dir])
  target_special_q_total_start = time.process_time()
  spec_q_log_sieve_radius = 10
  las_num_threads = 16
  spec_q_num_threads = las_num_threads
  spec_q_target_relns_mat, spec_q_target_elems_mat = target_special_q(alg_prime_str_dict, disc, fb, nf, target_dir,
                                                                      distinct_rat_prs_of_fb, smooth_bound, sp_reln_mat,
                                                                      spec_q_log_sieve_radius, spec_q_num_threads)
  print("asserts on special-q relations")
  assert assert_rows_multiply_to_correct_ideal(spec_q_target_relns_mat, spec_q_target_elems_mat, nf, fb)
  print("Done special q. Found {}".format(spec_q_target_elems_mat.nrows()))
  sp_reln_mat = sp_reln_mat.stack(spec_q_target_relns_mat)
  elems_mat = elems_mat.stack(spec_q_target_elems_mat)
  timers.filter_target_special_q_total = time.process_time()-target_special_q_total_start
  return elems_mat, sp_reln_mat, spec_q_target_relns_mat

def get_initial_ideal_and_common_factor_of_sieve_poly_coeffs(fb, idl_coeff_gcds, sparse_factor_base_vecs, idl_idx, nf, rat_pr_fb):
  com = idl_coeff_gcds[idl_idx]
  exp = sparse_factor_base_vecs[idl_idx]
  facs_of_com = [int(_) for _ in prime_factors(com)]
  initial_idl = nf.ideal(1)
  pr_facs=[]
  for i in range(len(exp)):
    ind_of_pr_idl=exp[i][0]
    initial_idl *= fb[ind_of_pr_idl]**exp[i][1]
    pr_facs.append(rat_pr_fb[ind_of_pr_idl])
  return initial_idl, facs_of_com, pr_facs

def target_special_q(alg_prime_str_dict, disc, fb, nf, path_to_target_dir, rat_fb, smooth_bound, sp_reln_mat,
                     log_sieve_radius, las_num_threads):
  if is_running_on_remote():
    mag_relns = magma_mat_from_sparse_sage_mat(sp_reln_mat)
    pivs = compute_piv_cols_fastest_available_and_update_time(mag_relns)

  log_smooth_bound = RR(smooth_bound).log(2)
  inds_that_must_be_made_pivots = []
  #TODO: this loop could be restricted to column index of the bounded fb rather than the full fb
  rat_prs_to_target=[]
  for _ in range(sp_reln_mat.ncols()):
    if _ in pivs:
      continue
    Pr = fb[_]
    pr_dg = Pr.residue_class_degree()
    if pr_dg > 1:
      continue
    pr = Pr.smallest_integer()
    pr_rd = Pr.ramification_index()
    # if norm is larger than smoothness bound then continue
    log_nrm = (pr_rd*pr_dg)*RR(pr).log(2)
    if log_nrm > log_smooth_bound or pr_dg == 4:
      continue
    rat_prs_to_target.append(pr)
    inds_that_must_be_made_pivots.append(_)
  rat_prs_to_target=sorted(list(set(rat_prs_to_target)))

  sieving_poly=nf.polynomial()
  mod_gens=1,-nf.gen()
  prs_used=[]

  special_q_filepath=path_to_target_dir+"/special_q"
  #example of todo format '1 29 14', '1 37 27'
  roots_start=time.process_time()
  roots_str=''
  print("Finding roots of {} primes".format(len(rat_prs_to_target)))
  for pr in rat_prs_to_target:
    gf_poly = PolynomialRing(GF(pr), names='W')(list(sieving_poly))
    roots=[_[0] for _ in list(gf_poly.roots())]
    for root in roots:
      roots_str += "1 {} {}\n".format(str(pr), str(root))
  timers.filter_target_special_q_roots=time.process_time()-roots_start

  with open(special_q_filepath, 'w+') as fp:
    fp.write(roots_str)

  verbose=False
  max_num_relns=10*len(fb)
  filename_without_extension='special_q'
  poly_filename = filename_without_extension + ".poly"
  fb_filename = filename_without_extension + ".roots0.gz"

  exit_early=None
  is_special_q=True
  disc_div=None
  target_relns_mat, target_elems_mat = target_poly(verbose, smooth_bound, path_to_target_dir, fb_filename,
                                                   poly_filename, prs_used, fb, rat_fb, alg_prime_str_dict,
                                                   max_num_relns, log_sieve_radius, sieving_poly,las_num_threads, disc,
                                                   mod_gens, exit_early, special_q_filepath, is_special_q, disc_div)
  return target_relns_mat, target_elems_mat

def compute_piv_cols_fastest_available_and_update_time(mat):
  piv_start_time = time.process_time()
  if is_running_on_remote():
    mag_piv_start_time = magma.cputime()

  pivs = pivot_cols_fastest_available(mat)

  if is_running_on_remote():
    timers.filter_pivots += time.process_time() - piv_start_time + magma.cputime(mag_piv_start_time)
  else:
    timers.filter_pivots += time.process_time()
  return pivs

def target_nonpivot_cols(alg_prime_str_dict, disc, fb, nf, norm_bnd_fb, poly_dir, rat_fb, smooth_bound, sp_reln_mat,
                         target_verbose, log_sieve_radius, disc_div):
  rank_needed = len(norm_bnd_fb)
  mag_relns = magma_mat_from_sparse_sage_mat(sp_reln_mat)
  target_prs,inds_that_must_be_made_pivots = prs_to_target(fb, mag_relns, smooth_bound)

  target_relns_mat = Matrix(ZZ, nrows=0, ncols=len(fb), sparse=True)
  target_elems_mat = Matrix(nf, nrows=0, ncols=1)

  #TODO this should be computed in a single place
  pol_ind=index(nf.polynomial())
  nrms_bounded_fb = [_.norm() for _ in norm_bnd_fb]
  rat_prs_bounded_fb = [_.norm()**(_.residue_class_degree()**(-1)) for _ in norm_bnd_fb] #TODO make a comparison of this with using the smallest_integer function
  second_gen_residues = map_pr_idls_to_dedekind_kummer_elements(list(zip(rat_prs_bounded_fb, norm_bnd_fb)), nf.polynomial())

  if target_verbose:
    rng = trange(len(target_prs))
  else:
    rng = trange(len(target_prs))

  log_heights_and_yields=[]
  for i in rng:
    Pr_pair = target_prs[i]
    max_num_relns = 10
    idl=Pr_pair[0]

    target_relns, target_elems, sieving_poly = compute_poly_of_idl_and_target_idl(target_verbose, smooth_bound, idl,
                                                                                  poly_dir, Pr_pair[1], Pr_pair[2], fb,
                                                                                  rat_fb, alg_prime_str_dict,
                                                                                  max_num_relns, log_sieve_radius, disc,
                                                                                  norm_bnd_fb, second_gen_residues,
                                                                                  nrms_bounded_fb, rat_prs_bounded_fb,
                                                                                  disc_div, pol_ind)
    log_heights_and_yields.append([log_height(sieving_poly, 15), target_relns.nrows()])
    #TODO temp moved the timer so that I could see the elem stack time
    target_relns_mat = target_relns_mat.stack(target_relns)

    stack_start = time.process_time()
    target_elems_mat = target_elems_mat.stack(target_elems)
    timers.filter_target_new_polys_stack += time.process_time()-stack_start
  print("Log heights vs relation yields {}".format(log_heights_and_yields))
  return rank_needed, target_elems_mat, target_relns_mat

def prs_to_target(fb, mag_relns, smooth_bound):
  ncols=mag_relns.NumberOfColumns().sage()

  pivs=compute_piv_cols_fastest_available_and_update_time(mag_relns)

  print("smooth bound {}".format(smooth_bound))
  log_smooth_bound = RR(smooth_bound).log(2)
  inds_that_must_be_made_pivots = []
  # TODO: this loop could be restricted to column index of the bounded fb rather than the full fb
  prs_to_target = []
  for _ in range(ncols):
    if _ in pivs:
      continue
    Pr = fb[_]
    pr = Pr.smallest_integer()
    pr_dg = Pr.residue_class_degree()
    pr_rd = Pr.ramification_index()
    # if norm is larger than smoothness bound then continue
    log_nrm = (pr_rd*pr_dg)*RR(pr).log(2)
    if log_nrm > log_smooth_bound or pr_dg == 4:
      continue
    idl_name = "p"+str(pr)+"_dg"+str(pr_dg)+"_rm"+str(pr_rd)
    prs_to_target.append([Pr, idl_name, pr])
    inds_that_must_be_made_pivots.append(_)
  return prs_to_target, inds_that_must_be_made_pivots

def factor_least_int_and_common_fac(YY, idl_coeff_gcds, idl_idx, idl_least_ints, idl_poly_gens, nf):
  com = idl_coeff_gcds[idl_idx]
  facs_of_com = [int(_) for _ in prime_factors(com)]
  return com, facs_of_com

def merge_large_primes_and_slice_mat(visited, inds_of_rows_to_modify, inds_of_rows_to_del, sp_reln_mat, elems_mat, sp_rat_reln_mat):
  inds_of_post_merge_dups = []
  for i in range(len(visited)):
    ind_to_del = inds_of_rows_to_del[i]
    sp_reln_mat.swap_rows(ind_to_del, i)
    if sp_rat_reln_mat != None:
      sp_rat_reln_mat.swap_rows(ind_to_del, i)
    if elems_mat != None:
      # manual swap b/c calling swap_rows on elems_mat throws TypeError:object cannot be interpreted as an integer
      elems_mat[i], elems_mat[ind_to_del] = elems_mat[ind_to_del], elems_mat[i]
    for j in range(len(inds_of_rows_to_modify[i])):
      # this is equivalent to making the row equal mat[inds_of_rows_to_modify[i][j]] - mat[inds_of_rows_to_del[i]]
      ind = inds_of_rows_to_modify[i][j]
      sp_reln_mat.add_multiple_of_row(ind, i, -1)  # the add_multiple_of_row function is optimized
      if sp_rat_reln_mat != None:
        sp_rat_reln_mat.add_multiple_of_row(ind, i, -1)
      if elems_mat != None:
        elems_mat[ind, 0] = elems_mat[ind, 0] / elems_mat[i, 0]  # this line is just for testing. In production the vector of logarithms will be used

  nm_del = len(inds_of_rows_to_del)
  sp_reln_mat = copy(sp_reln_mat[nm_del:])
  if sp_rat_reln_mat != None:
    sp_rat_reln_mat = copy(sp_rat_reln_mat[nm_del:])
  if elems_mat != None:
    elems_mat = copy(elems_mat[nm_del:])
  return sp_reln_mat, elems_mat, sp_rat_reln_mat

def assert_correct_norm(sieve_idl, expected_norm):
  '''
  :param sieve_idl: an ideal with norm that is not equal to 1
  :param expected_norm:
  :return:
  '''
  #TODO rewrite this function to use assertions instead of raising errors
  if not asserts_on():
    return
  nf = sieve_idl.number_field()
  if sieve_idl == nf.ideal(1):
    raise ValueError("Ideal must not have norm equal to 1")
  fac_data = sieve_idl.factor()
  facs = [_[0] for _ in fac_data]
  pows = [_[1] for _ in fac_data]
  for _ in pows:
    if _ != 1:
      raise ValueError("Sieve idl is divisible by a power larger than 1. Facs={}".format(fac_data))
  sieve_norm = prod_all_elems(facs).norm()
  if sieve_norm != expected_norm:
    raise ValueError("Sieve idl should have norm equal to the least integer that was stored")
  return True

def assert_elems_smooth_and_nonzero(elems_mat, nf, smooth_bound):
  if not asserts_on():
    return
  count_non_smooth_numerator = 0
  count_elems_that_are_0 = 0
  n_elems = 0
  for i in range(0, elems_mat.nrows()):
    n_elems += 1
    if elems_mat[i,0] == 0:
      count_elems_that_are_0 += 1
      continue
    nrm = nf.ideal(elems_mat[i,0]).norm()
    [is_smooth, fac] = is_bound_smooth(nrm.numerator(), smooth_bound, True)
    if not is_smooth:
      print("Found an elem that's not B-smooth, largest fac is {}\nFound at row ind={}".format(fac,i))
      count_non_smooth_numerator += 1
  #TODO (Sep 29, 2022) after I do a pass to eliminate duplicates before the first merge replace this print with ValueError
  if count_elems_that_are_0 > 0:
    raise ValueError("Validation: found {} elems were 0 out of {}".format(count_elems_that_are_0, n_elems))
  if count_non_smooth_numerator > 0:
    raise ValueError("Validation: after filtering number of non-smooths elements is {}".format(count_non_smooth_numerator))
  return True

def assert_rat_pr_mat_rows_multiply_to_correct_ideal(relns_mat, elems_mat, rat_prs):
  '''
  :param relns_mat:
  :param elems_mat:
  :param rat_prs:
  :return:
  '''
  if not asserts_on():
    return
  nrows= elems_mat.nrows()
  for j in range(nrows):
    if j % 500 == 0:
      print("Progress on asserts {:2f}%".format(100*j/nrows))

    row = relns_mat[j:j + 1, :].list()
    nrm = elems_mat[j,0].norm()

    actual_nrm = ZZ(1)
    for i in range(0, len(row)):
      ee = row[i]
      if ee != 0:
        new_term = rat_prs[i]**ee
        actual_nrm *= new_term
    if actual_nrm.abs() != nrm.abs():
      raise ValueError("The norm of row j={} is wrong\nActual facs={}\nExpected facs={}\nActual/expected {}".format(j, actual_nrm.factor(), nrm.factor(), actual_nrm.abs()/nrm.abs()))
  return True

def assert_rows_multiply_to_correct_ideal(relns_mat, elems_mat, nf, fb, large_primes_mat=None):
  '''
  :param relns_mat:
  :param elems_mat:
  :param nf:
  :param fb:
  :param large_primes_mat: None or a matrix of 5 columns. First column holds the norm
  :return:
  '''
  if not asserts_on():
    return
  nrows= elems_mat.nrows()
  print("starting to assert rows multiply to right ideal")
  for j in trange(nrows):
    assert_row_multiplies_to_correct_idl(j, relns_mat, elems_mat, nf, fb, large_primes_mat)
  return True


def assert_row_multiplies_to_correct_idl(j, relns_mat, elems_mat, nf, fb, large_primes_mat=None):
  row = relns_mat[j:j+1, :].list()
  idl = nf.ideal(elems_mat[j, 0])
  if idl == nf.ideal(0):
    raise ValueError("elems mat contains the element 0 at index {}".format(j));
  prod = nf.ideal(1)
  for i in range(0, len(row)):
    ee = row[i]
    if ee != 0:
      new_term = fb[i]**ee
      prod *= new_term
  err_str = "At index j={} (fb product corresponding to jth row of relns_mat) != <elems_mat[j]>\nrow is {}\nfactorization of quo {}\nelem{}\nelem norm={}\nidl={}\nprod={}".format(
    j, row, (idl/prod).factor(), elems_mat[j, 0], elems_mat[j, 0].norm().factor(), idl.factor(), prod.factor())
  if prod == nf.ideal(0):
    raise ValueError("prod is the 0 idl")
  if large_primes_mat == None:
    if not prod/idl == 1:
      raise ValueError(err_str)
  else:
    # if prod does not divide idl then (idl/prod).norm() is (very likely) in QQ not ZZ
    try:
      nrm = ZZ((idl/prod).norm())
    except TypeError:
      raise ValueError(err_str)
    if large_primes_mat[j, 0] == 0:
      if not prod/idl == 1:
        raise ValueError(err_str)
    elif not nrm.divides(large_primes_mat[j, 0]):
      raise ValueError(err_str)
  return True
