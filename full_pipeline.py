import random
import sys
import os
import glob
from datetime import datetime
from contextlib import redirect_stdout
from pathlib import Path

from sage.all import *

from src.stage_init import *
from src.stage_sieve import run_batch_las
from src.stage_filter import *
from src.stage_linalg import *
from src.polredabs_polys import *
from src.paths import basepath_artifacts
from src.paths import basepath_save
from src.cado_poly import *
from src.utils import class_num_formula_approx
import src.timers as timers
from src.json_io import *
from src.data_workflow.sane_parameter_defaults import *

def default_params_path(p_str):
  return "{0}/workdir/{0}_params".format(p_str)

def run(p_str, start_at_str, full_rank_timing, time_per_relation_timing, is_SI):
  '''
  Uses the parameters file to run the full algorithm
  :param is_SI:
  :param dirname: expected usage is that the dirname is just the polynomial's name
  :return:
  '''
  params_path=default_params_path(p_str)
  test_data = basepath_save + "{}/test_data.txt".format(p_str)
  set_seed(0)

  #TODO consider a tuple of stages instead of all these flags?
  if start_at_str=='filter':
    start_at_filter=True
    start_at_linalg=False
  elif start_at_str=='linalg':
    start_at_filter=False
    start_at_linalg=True
  else:
    start_at_filter=False
    start_at_linalg=False

  (rel_dirname, params_filename) = os.path.split(params_path)
  rel_dirname += '/'
  path_to_workdir = basepath_artifacts + rel_dirname
  if is_running_on_remote():
    path_to_params_dir = basepath_save_remote + rel_dirname
  else:
    path_to_params_dir = basepath_save + rel_dirname
  pol = read_cado(path_to_params_dir, params_filename, verbose=False, dg=4)
  smoothness_bound, max_num_primes, sieve_radius, name= get_stage_init_params(path_to_params_dir, params_filename, verbose=False)

  start_dt = datetime.now()
  timestamp_str = "Starting {}\t{}".format(params_filename, start_dt)
  file_timestamp_str = "{}_{}_{}_{}".format(datetime.today().month, datetime.today().day, start_dt.hour, start_dt.minute)
  log_path = "{}full_pipeline_log_{}.txt".format(path_to_params_dir, file_timestamp_str)
  #readable_log_path is a path that other programs can read from
  if is_SI:
    readable_log_path = "{}/sieve.log".format(path_to_params_dir)
  else:
    # use a different log for multiple polys
    readable_log_path = "{}/sieve_mp.log".format(path_to_params_dir)

  def run_tasks():
    print("Running")
    timers.start_time_user = time.perf_counter()
    timers.start_time_cpu = time.process_time()
    start_time = time.process_time()
    print(timestamp_str)
    print("SEED {}".format(get_seed()))
    verbose = True

    # write setup data
    nf = NumberField(pol, names='YY')
    [_, _, fb, rat_fb, alg_prime_str_dict, _, _] = fac_base_data(nf, smoothness_bound)
    if not is_prod():
      pari_zk = nf.pari_zk()
      for idl in fb:
        validate_idl(idl.pari_hnf(), idl.gens_two(), pari_zk)

    print("Write basis - make sure this doesn't take too long b/c this can just be a precomputation")
    ic_path = basepath_artifacts + rel_dirname + 'ind_calc.json'
    nf_path = basepath_artifacts + rel_dirname + 'nf.json'
    fb_path = basepath_artifacts + rel_dirname + 'fb.json'
    write_ic_and_fb_json(fb, rat_fb, smoothness_bound, ic_path, nf_path, fb_path)

    verbose=True

    # clean pad directory
    target_dirpath = basepath_save + rel_dirname + 'target'
    if os.path.exists(target_dirpath) and len(os.listdir(target_dirpath)) > 0:
      destroy_files_path = basepath_save + rel_dirname + 'target'
      files = glob.glob(destroy_files_path + "/*")
      for f in files:
        os.remove(f)

    params_filename = Path(path_to_params_dir).parts[-2] + '_params'
    # params_filename = get_params_file(path_to_params_dir, params_spec_str=, verbose=False)
    job_name = get_job_name(path_to_params_dir, params_filename, False)
    params_path = path_to_params_dir + params_filename
    with open(params_path, "r") as fp:
      params=fp.read()
      # this makes every log file start with the list of params
      print(params)

    path_to_relns_dir = get_relations_dir(path_to_workdir, path_to_params_dir, params_filename)
    upload_rel_filepath = path_to_relns_dir + '/all_reln'
    path_to_params_file = path_to_params_dir + get_params_file(path_to_params_dir, params_filename, verbose)

    Delta=nf.discriminant().abs()
    log_Del = floor(Delta.log(2)) # must be larger than unit rank. Larger = fewer primes to saturate AND greater liklihood saturation works.
    disc_div = disc_divisor(pol)
    if disc_div==None:
      disc_div=Delta

    if start_at_filter or start_at_linalg:
      print("WARNING - start_at_filter=True")
      if os.path.exists(path_to_relns_dir) and len(os.listdir(path_to_relns_dir)) == 0:
        raise ValueError('start_at_filter or start_at_linalg=True but upload dir is empty')
    else:
      if os.path.exists(path_to_relns_dir) and len(os.listdir(path_to_relns_dir)) > 0:
        print("WARNING - I am deleting the files in your upload directory (start_at_filter=False)")
        # this is the same as path_to_relns_dir but I want to make this name as explicit as possible so I don't delete files in the wrong place
        destroy_files_path = basepath_artifacts + rel_dirname + job_name + '.upload'
        files = glob.glob(destroy_files_path + "/*")
        for f in files:
          os.remove(f)
      smoothness_bound_las_and_fb=smoothness_bound # quick experiment suggests it's not worth changing this parameter
      [_, params_files] = run_batch_fb(path_to_workdir, path_to_params_file, verbose, pol, max_num_primes,
                                       sieve_algo_all_prs, fb, smoothness_bound_las_and_fb, sieve_radius, name,
                                       disc_div, disc, is_SI)
      time_per_poly = timers.init_poly_construct / 2 ** max_num_primes
      run_batch_las(path_to_workdir, path_to_relns_dir, smoothness_bound_las_and_fb, job_name, None, params_files, None)
      num_relns=unzip_and_concat(path_to_workdir, path_to_relns_dir, params_filename, verbose=False)

    log_level=1
    sieve_filtered_mat_path = rel_dirname + 'filtered_mat'
    time_per_poly=0
    target_path=rel_dirname + "target/"

    if not start_at_linalg:
      # stage_filter stuff
      sanity_check=True

      filter_relns_to_load_ratio = get_stage_filter_params(path_to_params_dir, params_filename, False)
      sieve_path = rel_dirname + 'sieve_stats.txt'

      try:
        relns_mat, elems_mat = filter(log_level, pol, smoothness_bound, upload_rel_filepath, sieve_filtered_mat_path,
                                      sieve_path, filter_relns_to_load_ratio, sieve_radius, full_rank_timing, Delta,
                                      target_path, disc_div, verbose, sanity_check, time_per_relation_timing, num_relns)

        if time_per_relation_timing:
          return
        if full_rank_timing:
          print("---------\nTime to full rank timing profile:")
          print_profile(time_per_poly, max_num_primes, full_rank_timing, is_SI)
          return
      except ValueError as e:
        print_profile(time_per_poly, max_num_primes, full_rank_timing, is_SI)
        print("FAIL in stage {}: ValueError {}".format("Filter", e))
        raise e
      except AssertionError as e:
        print_profile(time_per_poly, max_num_primes, full_rank_timing, is_SI)
        print("FAIL in stage {}: Assertion failed {}".format("Filter", e))
        raise e

    stage_linalg(elems_mat, nf, relns_mat, False)

    print_profile(time_per_poly, max_num_primes, full_rank_timing, is_SI)
    print_stdout_and_err("PASS")

  if not is_prod():
    run_tasks()
  else:
    with open(log_path, "w") as ff:
      with redirect_stdout(ff):
        run_tasks()
    sleep(1)
    shutil.copyfile(log_path, readable_log_path)

def print_stdout_and_err(str):
  sys.stderr.write(str)
  print(str)

def print_profile(time_per_poly, max_num_primes, full_rank_timing_test, is_SI):
  timers.filter_target_new_polys_total += timers.filter_target_new_polys_las_time
  timers.filter_total += timers.filter_target_new_polys_las_time
  external_time = timers.las_time + timers.filter_target_new_polys_las_time + timers.linalg_magma_hnf_time
  excluded_time = timers.stats_time
  overhead_time = 0.1*timers.linalg_total_time
  timers.linalg_total_time += overhead_time
  timers.total_time_cpu = time.process_time() - timers.start_time_cpu
  timers.total_time_cpu -= excluded_time
  timers.total_time_cpu += external_time
  timers.total_time_cpu += overhead_time
  print("Excluded time from CPU time: {:.1f}".format(excluded_time))
  print_stdout_and_err("CPU time:")
  if asserts_on():
    print("WARNING: Asserts were ON for these timings")
  level0_spacer = ' '*1
  level1_spacer = ' '*2
  level2_spacer = ' '*3
  level3_spacer = ' '*4
  stats_str = '- '
  user_str = ''
  if is_SI:
    total_str="{}{:.2f} total time for whole SI algorithm, this should roughly be the same as (INIT + SIEVE + FILTER + LINALG)={:.1f}\n---------\n".format("", timers.total_time_cpu, timers.init_total + timers.las_time + timers.filter_total + timers.linalg_total_time)
  else:
    total_str="{}{:.2f} total time for whole MP algorithm, this should roughly be the same as (INIT + SIEVE + FILTER + LINALG)={:.1f}\n---------\n".format("", timers.total_time_cpu, timers.init_total + timers.las_time + timers.filter_total + timers.linalg_total_time)
  cpu_str_brief=total_str
  cpu_str=total_str
  init_str="{}{:.1f} INIT\n".format(level0_spacer, timers.init_total)
  cpu_str+=init_str
  cpu_str_brief+=init_str
  cpu_str+="{}{:.1f} SI knapsack\n".format(level1_spacer, timers.SI_knapsack)
  cpu_str+="{}{:.2f} polys construct, for 2^{} polys.\n".format(level1_spacer, timers.init_poly_construct, max_num_primes)
  cpu_str+="{}{:.1f} total time applying lll to sieving ideals\n".format(level2_spacer, timers.init_lll)
  cpu_str+="{}{:.1f} total time computing resultant of sieving ideals\n".format(level2_spacer, timers.init_resultant)
  cpu_str+="{}{:.1f} MP total time solving knapsack\n".format(level2_spacer, timers.MP_knapsack_total)
  cpu_str+="{}{}Time per poly {:.1f}\n".format(level1_spacer, stats_str, time_per_poly)
  cpu_str+="{}{:.1f} makefb calls\n".format(level1_spacer, timers.init_makefb)
  sieve_str="{}{:.1f} SIEVE\n".format(level0_spacer, timers.las_time)
  cpu_str+=sieve_str
  cpu_str_brief+=sieve_str
  user_str+="{}{:.2f} SIEVE\n".format(level0_spacer, timers.user_sieve)
  filter_str="{}{:.1f} FILTER\tsetup+newpolys={:.1f}\n".format(level0_spacer, timers.filter_total, timers.filter_setup + timers.filter_target_new_polys_total + timers.filter_convert_time)
  cpu_str+=filter_str
  cpu_str_brief+=filter_str
  cpu_str+="{}{:.1f} SETUP\n".format(level1_spacer, timers.filter_setup)
  cpu_str+="{}{:.1f} load mat. Loaded {} relations\n".format(level2_spacer, timers.filter_load_mat, timers.filter_load_mat_num_relns)
  cpu_str+="{}{:.1f} valuation\n".format(level3_spacer, timers.filter_valuation)
  cpu_str+="{}{:.1f} merge\n".format(level2_spacer, timers.filter_merge)
  cpu_str+="{}{:.1f} trivial reln\n".format(level2_spacer, timers.filter_trivial_relns)
  cpu_str+="{}{:.1f} FILTER - TARGET NEW POLYS\n".format(level1_spacer, timers.filter_target_new_polys_total)
  cpu_str+="{}{:.1f} knapsack\n".format(level2_spacer, timers.filter_target_new_polys_knapsack)
  cpu_str+="{}{:.1f} lll total\n".format(level2_spacer, timers.filter_target_new_polys_lll)
  cpu_str+="{}{:.1f} form total\n".format(level2_spacer, timers.filter_target_new_polys_form)
  cpu_str+="{}{:.1f} las total\n".format(level2_spacer, timers.filter_target_new_polys_las_time)
  cpu_str+="{}{:.1f} process las output\n".format(level2_spacer, timers.filter_target_new_polys_process_las_output)
  cpu_str+="{}{:.1f} constructing matrices\n".format(level2_spacer, timers.filter_target_new_polys_stack)
  cpu_str+="{}{:.1f} submat setup \n".format(level3_spacer, timers.filter_submat_setup)
  cpu_str+="{}{:.1f} remove singleton cols\n".format(level1_spacer, timers.filter_singleton_columns)
  cpu_str+="{}{:.1f} remove nonpivot rows (uses Magma)\n".format(level1_spacer, timers.filter_remove_nonpivot_rows)
  cpu_str+="{}{:.1f} eliminate all-zero cols\n".format(level1_spacer, timers.filter_eliminate_cols)
  cpu_str+="{}{:.1f} full rank submat (uses Magma) \n".format(level1_spacer, timers.filter_full_rank_submat)
  cpu_str+="{}{:.1f} FILTER Multiple Call Timers\n".format(level0_spacer, timers.filter_convert_time + timers.filter_pivots + timers.filter_class_num)
  cpu_str+="{}{:.1f} Convert sparse mat into Magma sparse\n".format(level1_spacer, timers.filter_convert_time)
  cpu_str+="{}{:.1f} Pivot columns\n".format(level1_spacer, timers.filter_pivots)
  cpu_str+="{}{:.1f} Class num\n".format(level1_spacer, timers.filter_class_num)
  if full_rank_timing_test:
    user_str="----\nUser time: {:.2f}\n{}".format(time.perf_counter() - timers.start_time_user, user_str)
    print(cpu_str)
    print(user_str)
    return
  linalg_str="{}{:.1f} LINALG\n".format(level0_spacer, timers.linalg_total_time)
  cpu_str+=linalg_str
  cpu_str_brief+=linalg_str
  cpu_str+="{}{:.2f} magma convert \n".format(level1_spacer, timers.linalg_magma_convert)
  cpu_str+="{}{:.2f} magma HNF \n".format(level1_spacer, timers.linalg_magma_hnf_time)
  print(cpu_str)
  sys.stderr.write(cpu_str_brief)
  user_time_str="----\nUser time:\n{:.2f}\n".format(time.perf_counter()-timers.start_time_user)
  print_stdout_and_err(user_time_str)

def timers_for_Vollmer_HNF(cpu_str, level0_spacer, level1_spacer, level2_spacer, stats_str):
  cpu_str += "{}{:.1f} singleton columns \n".format(level1_spacer, timers.filter_singleton_columns)
  cpu_str += "{}{:.1f} rank check \n".format(level1_spacer, timers.filter_rank_check)
  cpu_str += "{}{:.1f} eliminate zero columns \n".format(level1_spacer, timers.filter_eliminate_cols)
  cpu_str += "{}{:.1f} write mats \n".format(level1_spacer, timers.filter_write_mats)
  cpu_str += "{}{}Should be that filter's subtimers sum close to filter's time (within 1s)\n".format(level1_spacer,
                                                                                                     stats_str)
  cpu_str += "{}{:.1f} LINALG\n".format(level0_spacer, timers.linalg_total_time)
  cpu_str += "{}{:.1f} setup time\n".format(level1_spacer, timers.linalg_setup_time)
  cpu_str += "{}{:.1f} Step0 initial hnf\n".format(level1_spacer, timers.linalg_step0_cg_initial_relns)
  cpu_str += "{}{:.1f} Step1 additional relns\n".format(level1_spacer, timers.linalg_step1_cg_add_relns)
  cpu_str += "{}{:.1f} Step2 initial relns\n".format(level1_spacer, timers.linalg_step2_ug_initial_relns)
  cpu_str += "{}{:.1f} project submodule \n".format(level2_spacer, timers.linalg_project_submodule)
  cpu_str += "{}{:.1f} ug add rows (vollmer) \n".format(level2_spacer, timers.linalg_ug_system_solving)
  cpu_str += "{}{:.1f} Step3 saturate\n".format(level1_spacer, timers.linalg_step3_ug_saturate)
  cpu_str += "{}{:.1f} 3a) constructing inerts and character mats \n".format(level2_spacer,
                                                                             timers.linalg_saturate_setup_prs_and_character_mat)
  cpu_str += "{}{:.1f} 3b) finding a solution to saturation \n".format(level2_spacer,
                                                                       timers.linalg_saturate_individual_mats)
  cpu_str += "{}{:.1f} 3c) constructing the combined matrix \n".format(level2_spacer,
                                                                       timers.linalg_saturate_intersection_reg)
  cpu_str += "{}{:.1f} Step4 complete\n".format(level1_spacer, timers.linalg_step4_cg_complete)
  cpu_str += "{} LINALG COMPONENTS \n".format(level1_spacer)
  cpu_str += "{}{:.1f} time spent in add_rows (HNF)\n".format(level1_spacer, timers.linalg_add_rows)
  cpu_str += "{}{:.1f} time add_rows spends just doing HNF calculations\n".format(level2_spacer,
                                                                                  timers.linalg_hnf_add_rows)
  cpu_str += "{}{:.1f} time spent in left_kernel\n".format(level1_spacer, timers.linalg_kern)
  cpu_str += "{}{:.1f} regulator \n".format(level1_spacer, timers.linalg_regulator)
  return cpu_str

def profile(func_to_profile):
  pr = cProfile.Profile()
  pr.enable()
  func_to_profile()
  pr.disable()
  stats = Stats(pr)
  stats.sort_stats('tottime').print_stats(10)
