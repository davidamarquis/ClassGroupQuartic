import time
import json
from sage.all import *
from src.utils import count_generator
from src.paths import *
from src.utils import str_to_int_list
from src.utils import asserts_on
from src.obj_string_converters import *
import re

def find_float(str):
  m = re.search(r'(\d+\.\d+)', str)
  if m:
    return float(m.group(1))
  raise ValueError("string should have a float")

def get_invariants_from_number_field_test_data(filepath):
  #TODO maybe using RDF precision for the regulator is too low
  with open(filepath, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      newline = line.replace(" ", "")
      regulator_str = "regulator="
      class_group_str = "group="
      class_number_str = "number="
      if newline.startswith(regulator_str):
        regulator = RDF(newline.split(regulator_str,1)[1].rstrip())
      elif newline.startswith(class_group_str):
        class_group = str_to_int_list(newline.split(class_group_str,1)[1].rstrip())
      elif newline.startswith(class_number_str):
        class_number = ZZ(newline.split(class_number_str,1)[1].rstrip())
  return [regulator, class_group, class_number]

def write_stats(small_primes_mat, fb, rel_path, sieve_idl_norm, num_times_targetted):
  '''
  :param small_primes_mat:
  :param fb:
  :param rel_path:
  :param sieve_idl_norm: set to None if using random reduction
  :param num_times_targetted: set to None if using sieving
  :return:
  '''
  if sieve_idl_norm == None:
    sieve_idl_nrm = 1
  if num_times_targetted == None:
    num_times_targetted = 0

  frequency_pr_idl_divides_reln = [0]*len(fb)
  for _ in range(len(fb)):
    for rr in range(small_primes_mat.nrows()):
      if small_primes_mat[rr,_] != 0:
        frequency_pr_idl_divides_reln[_] += 1

  write_path = basepath_artifacts + rel_path
  number_reln_actual = small_primes_mat.nrows()
  with open(write_path, 'w+') as fp:
    fp.write("fb size = {}\n".format(len(fb)))
    fp.write("number of relns = {}\n".format(number_reln_actual))
    fp.write("{:<16}|{:<18}|{:<8}|{:<8}\n".format("actual freq", "expected freq", "dg", "sieve prime"))
    for _ in range(len(fb)):
      nrm = fb[_].norm()
      is_sieve_idl = gcd(nrm, sieve_idl_norm) > 1
      adjusted_freq = frequency_pr_idl_divides_reln[_] - num_times_targetted
      stats_str ="{:<16.5f}|{:<18.5f}|{:<8}|{:<8}\n".format(adjusted_freq/number_reln_actual, 1/float(nrm), fb[_].residue_class_degree(), is_sieve_idl)
      fp.write(stats_str)

def write_int_mat(write_rel_path, mat, log_level, formatter=mat_to_str):
  '''
  Writes a flint matrix to the file given by combining a basepath stored in constants.py with the relative path given as input
  :param write_rel_path:
  :param mat:
  :param log_level:
  :param formatter: either mat_to_str or mat_to_flint_str
  :return:
  '''
  write_path = basepath_artifacts + write_rel_path
  start_write_time = time.process_time()
  if formatter == mat_to_str:
    with open(write_path, 'wb+') as fp:
      #I tested fp.write(str(small_primes_mat).encode("utf-8")) and found it was no faster than writing each line with this loop
      str = formatter(mat)
      fp.write(str.encode("utf-8"))
  elif formatter == mat_to_flint_mat_str:
    with open(write_path, 'w+') as fp:
      fp.write(formatter(mat))
  else:
    raise ValueError("Invalid choice of formatter")
  if log_level >= 1:
    diff = time.process_time() - start_write_time
    print("time spent on writing the matrix (TODO lines) was {}".format(diff))
    time_per_line = diff/mat.nrows()
    print("time spent per line is {}".format(time_per_line))
  return

def write_nf_elem_mat(write_rel_path, mat, log_level):
  '''
  Writes a matrix to the file given by combining a basepath stored in constants.py with the relative path given as input
  :param write_rel_path:
  :param mat:
  :param log_level:
  :return:
  '''
  if log_level >= 1:
    start_write_time = time.process_time()
  write_path = basepath_artifacts + write_rel_path
  with open(write_path, 'wb+') as fp:
    #I tested fp.write(str(small_primes_mat).encode("utf-8")) and found it was no faster than writing each line with this loop
    n_rows = mat.nrows()
    for i in range(n_rows):
      # the next line is the only one that differs from the function for writing integer matrices
      line_str = nf_elem_to_str(mat[i,0])
      if i < n_rows-1:
        line_str = line_str + '\n'

      fp.write(line_str.encode("utf-8"))
  if log_level >= 1:
    diff = time.process_time() - start_write_time
    print("time spent on writing the matrix (TODO lines) was {}".format(diff))
    time_per_line = diff/mat.nrows()
    print("time spent per line is {}".format(time_per_line))
  return

def read_int_mat(log_level, path, basepath_override, rows_to_cols_ratio_for_load):
  '''
  loads a matrix stored in a binary file encoded in utf-8
  :param log_level:
  :param path:
  :param basepath_override:
  :param rows_to_cols_ratio_for_load:
  :return:
  '''
  if basepath_override == None:
    path = basepath_artifacts + path
  else:
    path = basepath_override + path

  num_rows = 0
  with open(path, 'rb') as fp:
    c_generator = count_generator(fp.raw.read)
    num_rows = sum(buffer.count(b'\n') for buffer in c_generator) # count each \n
  num_rows += 1 # last line does not end with '\n'

  if log_level > 1:
    print("Loaded a matrix with nrows={}".format(num_rows))

  num_cols = 0
  # read one line and get the number of cols from that
  with open(path, 'rb') as fp:
    line = fp.readline();
    str_line = line.decode('utf-8')
    str_line.rstrip()
    str_list = str_line.split(',')
    num_cols = len(str_list)

  num_rows_to_load = min(num_rows, floor(rows_to_cols_ratio_for_load * num_cols))

  mat = Matrix(ZZ, num_rows_to_load, num_cols, sparse=True)

  if log_level > 1:
    print("Created a matrix with dims {}x{}".format(num_rows_to_load, num_cols))

  line_num = 0
  with open(path, 'rb') as fp:
    for _ in range(num_rows_to_load):
      line = fp.readline()
      str_line = line.decode('utf-8')
      str_line.rstrip()
      str_list = str_line.split(',')
      try:
        # looping over columns is much faster than setting the entire row by assigning mat[i]
        for ci in range(len(str_list)):
          mat[_, ci] = int(str_list[ci])
      except IndexError as e:
        print("IndexError at line num {}".format(line_num))
        raise IndexError(e)
      line_num+=1
  return mat

def read_nf_elem_mat(log_level, path, basepath_override, num_rows_to_load, nf):
  '''
  loads a matrix stored in a binary file encoded in utf-8
  :param log_level:
  :param path:
  :param basepath_override:
  :param num_rows_to_load:
  :return:
  '''
  if basepath_override == None:
    path = basepath_artifacts + path
  else:
    path = basepath_override + path

  num_rows = 0
  with open(path, 'rb') as fp:
    c_generator = count_generator(fp.raw.read)
    num_rows = sum(buffer.count(b'\n') for buffer in c_generator) # count each \n

  if log_level > 1:
    print("Loaded a matrix with nrows={}".format(num_rows))

  mat = Matrix(nf, num_rows_to_load, 1)

  if log_level > 1:
    print("Created a matrix with dims {}x{}".format(num_rows_to_load, 1))

  line_num = 0
  with open(path, 'rb') as fp:
    for _ in range(num_rows_to_load):
      line = fp.readline()
      str_line = line.decode('utf-8')
      str_line.rstrip()
      str_list = str_line.split(',')
      mat[_, 0] = nf([QQ(_) for _ in str_list])
      line_num+=1
  return mat

def read_fb_inds(path, basepath_override, fb):
  '''
  :param path:
  :param basepath_override:
  :param fb:
  :return: a factorbase that's reordered to match the columns
  '''
  if not asserts_on():
    raise ValueError("this function should never be called when asserts are off")
  if basepath_override == None:
    path = basepath_artifacts + path
  else:
    path = basepath_override + path

  with open(path, "r") as fp:
    lines = fp.readlines()
    fb_inds = [0]*(len(lines)-1)
    for i in range(len(lines)-1):
      fb_inds[i] = [int(_) for _ in lines[i].split(",")]
    num_to_delete = int(lines[len(lines)-1].rstrip())

  fb_in_same_order_as_cols = copy(fb) # assignment doesn't copy data
  for i in range(len(fb_inds)):
    i0,i1 = fb_inds[i]
    fb_in_same_order_as_cols[i0], fb_in_same_order_as_cols[i1] = fb_in_same_order_as_cols[i1], fb_in_same_order_as_cols[i0]
  fb_in_same_order_as_cols = fb_in_same_order_as_cols[num_to_delete:]
  return fb_in_same_order_as_cols

def utf8_decode(data, path):
  "Returns a Unicode object on success, or None on failure"
  # https://stackoverflow.com/a/3269387
  try:
    return data.decode('utf-8')
  except UnicodeDecodeError as e:
    e.reason = "data in file at {} is not valid unicode".format(path)
    raise e

def get_relations_dir(path_to_workdir, path_to_params_dir, params_spec_str):
  job_name = get_job_name(path_to_params_dir, params_spec_str)
  return path_to_workdir + job_name + ".upload/"

def get_params_file(path_to_params_dir, params_spec_str=None, verbose=False):
  count_params = 0
  last_params_file_found = None
  try:
    files=os.listdir(path_to_params_dir)
  except FileNotFoundError as e:
    print("Error found while running on remote: {}".format(is_running_on_remote()))
    raise(e)

  for file in files:
    filename = os.fsdecode(file)
    if filename.endswith("_params") and params_spec_str == None:
      if verbose:
        print("Found params filename {}".format(filename))
      count_params += 1
      last_params_file_found = copy(filename)
    elif params_spec_str == filename:
      if verbose:
        print("Found params filename {}".format(filename))
      count_params += 1
      return filename

  if params_spec_str == None:
    if count_params == 1 and params_spec_str == None:
      params_filename = last_params_file_found
    elif count_params == 0:
      raise ValueError("A params file is required in the directory you've specified")
    else:
      raise ValueError("this case is not supported")
  elif last_params_file_found == None:
    raise ValueError("no params file found in {}".format(path_to_params_dir))

  return params_filename

def get_job_name(path_to_params_dir, params_spec_str=None, verbose=False):
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
  job_name = None
  with open(params_full_path, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      newline = line.replace(" ", "")
      name_str = "name="
      if newline.startswith(name_str):
        job_name = newline.split(name_str,1)[1]
        job_name = job_name.rstrip()
        if verbose:
          print("Job name {}".format(job_name))
  return job_name

def get_smooth_bound(path_to_params_dir, params_spec_str=None, verbose=False):
  '''
  Parse a parameters file to get a parameters:
  qrange, I lim0 lpb0 mfb0
  This is used to get a full list parameters for running las
  :return:
  '''
  params_file = get_params_file(path_to_params_dir, params_spec_str, verbose)
  params_full_path = path_to_params_dir + params_file

  with open(params_full_path, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      newline = line.replace(" ", "")
      lim0_str = "tasks.lim0="
      lim0_str2 = "tasks.sieve.lim0="
      if newline.startswith(lim0_str):
        lim0 = newline.split(lim0_str,1)[1]
        lim0 = int(lim0.rstrip())
      if newline.startswith(lim0_str2):
        lim0 = newline.split(lim0_str2,1)[1]
        lim0 = int(lim0.rstrip())
  return lim0

def get_stage_filter_params(path_to_params_dir, params_spec_str=None, verbose=False):
  '''
  :return:
  '''
  params_file = get_params_file(path_to_params_dir, params_spec_str, verbose)
  params_full_path = path_to_params_dir + params_file
  with open(params_full_path, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      newline = line.replace(" ", "")
      relns_str="tasks.filter.relns_to_load_ratio"
      if line.startswith(relns_str):
        relns_to_load = float(newline.split("{}=".format(relns_str),1)[1].rstrip())
  return relns_to_load

def get_linalg_relns_to_load_ratio(path_to_params_dir, params_spec_str=None, verbose=False):
  '''
  :return:
  '''
  params_file = get_params_file(path_to_params_dir, params_spec_str, verbose)
  params_full_path = path_to_params_dir + params_file

  with open(params_full_path, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      if line.startswith("tasks.linalg.relns_to_load_ratio"):
        newline = line.replace(" ", "")
        ret_str = newline.split("tasks.linalg.relns_to_load_ratio=",1)[1]
        ret_str = float(ret_str.rstrip())
  if ret_str > 1:
    print("###\nWARNING: filter.relns_to_load > 1\n###")
  return ret_str

def get_large_prime_bound(path_to_params_dir, params_spec_str=None, verbose=False):
  '''
  Parse a parameters file to get a parameters:
  qrange, I lim0 lpb0 mfb0
  This is used to get a full list parameters for running las
  :return:
  '''
  params_file = get_params_file(path_to_params_dir, params_spec_str, verbose)
  params_full_path = path_to_params_dir + params_file

  with open(params_full_path, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      newline = line.replace(" ", "")
      lpb0_str = "tasks.sieve.SI.lpb0="
      if newline.startswith(lpb0_str):
        lpb0 = newline.split(lpb0_str,1)[1]
        lpb0 = int(lpb0.rstrip())
  return ZZ(2)**lpb0

def get_poly_file(path_to_params_dir, can_read_local_poly, params_spec_str=None, verbose=False):
  '''
  Returns the path to the polynomial file
  :param path_to_params_dir:
  :param params_spec_str:
  :param verbose:
  :return:
  '''
  params_file = get_params_file(path_to_params_dir, params_spec_str, verbose)
  params_full_path = path_to_params_dir + params_file
  with open(params_full_path, 'r') as fp:
    lines = fp.readlines()
    if is_running_on_remote() or not can_read_local_poly:
      poly_str = "tasks.polyselect.import="
    elif can_read_local_poly:
      poly_str = "tasks.polyselect.local_import="
    else:
      raise ValueError("could not find path to polynomial")

    for line in lines:
      newline = line.replace(" ", "")
      if newline.startswith(poly_str):
        poly_path = newline.split(poly_str)[1].rstrip()
        if verbose:
          print("poly path {}".format(poly_path))
        break
  return poly_path

def get_stage_init_params(path_to_params_dir, params_spec_str=None, verbose=False):
  '''
  Parse parameters file to get all parameters for the init stage
  :return:
  '''
  params_file = get_params_file(path_to_params_dir, params_spec_str, verbose)
  params_full_path = path_to_params_dir + params_file

  with open(params_full_path, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      newline = line.replace(" ", "")
      name_str = "name="
      lim0_str = "tasks.lim0="
      max_num_primes_str = "tasks.init.max_num_primes="
      sieve_rad_str = "tasks.I="
      if newline.startswith(name_str):
        name = newline.split(name_str,1)[1].rstrip()
      elif newline.startswith(lim0_str):
        lim0 = int(newline.split(lim0_str,1)[1].rstrip())
      elif newline.startswith(max_num_primes_str):
        max_num_primes = int(newline.split(max_num_primes_str,1)[1].rstrip())
      elif newline.startswith(sieve_rad_str):
        sieve_rad = int(newline.split(sieve_rad_str,1)[1].rstrip())
  return lim0, max_num_primes, ZZ(2)**sieve_rad, name

def str_to_bool(mystr):
  return mystr == "True"