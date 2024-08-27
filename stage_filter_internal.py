import subprocess
import argparse
import shlex
from sage.all import *
from src.io_utils import *
from src.paths import cado_las_binary
from src.utils import *
import src.timers as timers

def inds_for_merge_one(mat):
  '''
  The expected usage is that the rows of mat are associated with another matrix other_mat
  A=inds_of_rows_to_del be a tuple of the row indexs in other_mat to delete after. Ive chosen these to be the 1st occurrence of a nonzero row. Zero rows are ignored
  B=inds_of_rows_to_modify be a tuple of tuple of indexes of rows
  In normal usage this modification is a 'merge' of the rows of other_mat so that large primes are eliminated.
  Specifically we do this
  other_mat[B[i][j]] = other_mat[B[i][j]] - other_mat[A[i]]
  and then delete the rows of other_mat with indexs in A
  :param mat:
  :return:
  '''
  visited = []
  inds_of_rows_to_del = []
  inds_of_rows_to_modify = []
  for i in range(mat.nrows()):
    large_prime_row = list(mat[i])
    if large_prime_row == [0]*mat.ncols():
      continue
    if large_prime_row in visited:
      ind = visited.index(large_prime_row)
      inds_of_rows_to_modify[ind] += [i]
    else:
      visited += [large_prime_row]
      # grow by one the inds_of_rows_to_modify so that visited[i] corresponds to inds_of_rows_to_modify[i] for any i
      inds_of_rows_to_modify += [[]]
      inds_of_rows_to_del += [i]
  return visited, inds_of_rows_to_modify, inds_of_rows_to_del

def triv_rels(fb, alg_prime_str_dict, smooth_bound):
  '''
  :param fb:
  :param alg_prime_str_dict:
  :return:
  '''
  rat_primes = list(alg_prime_str_dict.keys())
  triv_rels_mat = Matrix(ZZ, len(rat_primes), len(fb), sparse=True)
  rat_prs_mat = Matrix(ZZ, len(rat_primes), len(rat_primes), sparse=True)
  elems_mat = Matrix(ZZ, len(rat_primes), 1)
  rat_prs_without_triv_reln = []
  nf = fb[0].number_field()
  disc = nf.discriminant()
  reln_ind = 0
  for i in range(len(rat_primes)):
    pr = rat_primes[i]
    if pr == rat_primes[-1]:
      rng_max = len(fb)
    else:
      # get the index of the next prime so we can determine the range of prime indexs to look at
      rng_max = alg_prime_str_dict[next_prime(pr,proof=False)]
    if (rng_max - alg_prime_str_dict[pr] == 1):
      # skip inert primes
      continue
    for _ in range(alg_prime_str_dict[pr], rng_max):
      if fb[_].norm() > smooth_bound:
        rat_prs_without_triv_reln.append(pr)
        # the norm of this algebraic prime is too large. Skip this rational prime
        break
    else:
      rat_prs_mat[reln_ind,i]=4
      for j in range(alg_prime_str_dict[pr], rng_max):
        if disc % pr != 0:
          triv_rels_mat[reln_ind,j] = 1
        else:
          # get the actual fb element and compute its valuation
          triv_rels_mat[reln_ind,j] = (fb[j]).ramification_index()
      elems_mat[reln_ind,0] = pr
      reln_ind += 1
  rat_prs_without_triv_reln = sorted(list(set(rat_prs_without_triv_reln)))
  triv_rels_mat = copy(triv_rels_mat[:reln_ind])
  rat_prs_mat = copy(rat_prs_mat[:reln_ind])
  elems_mat = copy(elems_mat[:reln_ind])
  return triv_rels_mat, elems_mat, rat_prs_without_triv_reln, rat_prs_mat

def get_rows_of_singleton_columns(mat):
  '''
  A singleton column c has c[i] is nonzero and c[k] = 0 for all i!=k.
  This function returns rows[i] for each such column.
  :param mat:
  :return:
  '''

  non_singleton_row_inds = [_ for _ in range(mat.nrows())]
  singleton_row_inds = []
  # loop over cols
  for j in range(mat.ncols()):
    nonzero_indexs = []
    # loop over rows
    for i in range(mat.nrows()):
      if mat[i,j] != 0:
        nonzero_indexs.append(i)
    if len(nonzero_indexs) == 1:
      ind = nonzero_indexs[0]
      if ind in non_singleton_row_inds:
        non_singleton_row_inds.remove(ind)
        singleton_row_inds.append(ind)
  return sorted(singleton_row_inds), non_singleton_row_inds

def inds_of_dup_rows_and_units(mat, elems_mat):
  '''
  Computes a partition [inds_of_dup_rows, inds_of_non_dup_rows] of [0,...,nrows-1] so that for every i in
   inds_of_dup_rows there is a j in inds_of_non_dup_rows so that mat[i]=mat[j]
  :param mat:
  :return:
  '''
  visited_rows = []
  inds_of_dup_rows = []
  nrows = mat.nrows()
  inds_of_non_dup_rows = [_ for _ in range(nrows)]
  units = []
  original_indexs = [] # indexs of the visited rows in the original matrix
  for i in range(nrows):
    row = mat[i]
    if mat[i] in visited_rows:
      index = visited_rows.index(mat[i])
      if elems_mat != None:
        unit = elems_mat[i]/elems_mat[original_indexs[index]]
        nrm_abs=unit.norm().abs()
        if nrm_abs!= 1:
          raise ValueError("supposed unit has norm {} at index {}. N(elems_mat[i])={} and N(elems_mat[original_indexs[index]])={}".format(nrm_abs, i, elems_mat[i].norm(), elems_mat[original_indexs[index]].norm()))
        #hardcode to degree 4
        if unit**2 != 1:
          units.append(unit)
      inds_of_dup_rows.append(i)
      inds_of_non_dup_rows.remove(i)
    else:
      visited_rows.append(row)
      original_indexs.append(i)
  return inds_of_dup_rows, inds_of_non_dup_rows, units

def find_row_inds_of_deg_two_pairs(alg_prime_str_dict, fb, idl_least_ints, rat_fb, sp_reln_mat):
  inds_of_rows_with_nonzero_at_prime_to_eliminate_from_fb = []
  inds_of_rows_to_keep = list(range(sp_reln_mat.nrows()))
  dg_two_pair_col_inds = col_inds_of_dg_two_primes_above_a_rational_with_two_facs(fb, rat_fb, alg_prime_str_dict)
  print("Purge step 2: remove relations that are nonzero at certain primes of degree larger than 1")
  for cc in dg_two_pair_col_inds:
    if gcd(fb[cc].norm(), idl_least_ints[0]) > 1:  # TODO assumes all values in idl_least_ints are the same
      print(
        "Sieving prime with norm {} was not purged because we probably have enough relations such that this prime divides the principal ideal generated by the relation's element".format(
          fb[cc].norm()))
      continue
    for rr in range(sp_reln_mat.nrows()):
      vv = sp_reln_mat[rr, cc]
      if vv != 0:
        inds_of_rows_with_nonzero_at_prime_to_eliminate_from_fb.append(rr)
        if rr in inds_of_rows_to_keep:
          inds_of_rows_to_keep.remove(rr)
  inds_of_rows_with_nonzero_at_prime_to_eliminate_from_fb = sorted(
    inds_of_rows_with_nonzero_at_prime_to_eliminate_from_fb)
  return inds_of_rows_to_keep, inds_of_rows_with_nonzero_at_prime_to_eliminate_from_fb

def test_swaps():
  mat = Matrix([[1,2],[3,4],[5,6]])
  for _ in range(mat.nrows()):
    mat.swap_rows(1,_)
    print("R1->R{}".format(_))
    print(mat)

def swap_rows_to_top(sp_reln_mat, lp_reln_mat, elems_mat, rat_pr_mat, inds_to_move_to_top, other_inds):
  '''
  Warning: destructively changes other_inds
  :param sp_reln_mat:
  :param lp_reln_mat: can be a matrix or None
  :param elems_mat:
  :param rat_pr_mat:
  :param inds_to_move_to_top:
  :param other_inds:
  :return:
  '''
  mat = copy(sp_reln_mat)
  for _ in inds_to_move_to_top:
    if _ < len(inds_to_move_to_top):
      continue
    ind_non_dup = other_inds[0]
    mat.swap_rows(ind_non_dup, _)
    if rat_pr_mat != None:
      rat_pr_mat.swap_rows(ind_non_dup, _)
    if lp_reln_mat != None:
      lp_reln_mat.swap_rows(ind_non_dup, _)
    if elems_mat != None:
      # manual swap b/c calling swap_rows on elems_mat throws TypeError:object cannot be interpreted as an integer
      elems_mat[_], elems_mat[ind_non_dup] = elems_mat[ind_non_dup], elems_mat[_]

    other_inds.remove(other_inds[0])
  return [mat, lp_reln_mat, elems_mat, rat_pr_mat]

def col_indexs_with_weight_below(mat, ham_weight):
  return filtered_col_indexs_with_weight_below(mat, ham_weight, lambda _:True)

def nonzero_col_indexs_with_weight_below(mat, ham_weight):
  if ham_weight > mat.nrows() or ham_weight < 0:
    raise ValueError("ham_weight is too large or small")
  return filtered_col_indexs_with_weight_below(mat, ham_weight, lambda _:not _.is_zero())

def filtered_col_indexs_with_weight_below(mat, ham_weight, filter):
  ncols = mat.ncols()
  inds_with_weight_geq_bound = []
  for j in range(ncols):
    col = mat.column(j)
    if filter(col):
      inds_with_weight_geq_bound.append(j)

  filtered_inds = copy(inds_with_weight_geq_bound)
  inds_with_weight_below_bound = []
  for j in filtered_inds:
    col = mat.column(j)
    if col.hamming_weight() < ham_weight:
      inds_with_weight_below_bound.append(j)
      inds_with_weight_geq_bound.remove(j)
  return inds_with_weight_below_bound, inds_with_weight_geq_bound

def pad_prepare_params_files(path_to_pad_parameters, rat_primes_to_sieve):
  '''
  Add a function that takes a filepath as an argument and a list of primes q. It opens every file at that path. It splits every line by spaces and if the second argument is in the list of primes it saves the file to a file with the same but with “for_pad” appended.
  If no rows remain after this filtering then delete the params file.
  Then calls a function pad_run_las_on_primes which calls pad_get_las_arg.
  '''
  lines_to_write = []
  for filename in os.listdir(path_to_pad_parameters):
    with open(filename, 'w') as file:
      for line in file.readlines():
        line_pr = int(line.split(" ")[1])
        if line_pr in rat_primes_to_sieve:
          lines_to_write.append(line)
    with open(filename + 'for_pad', 'w+') as write_file:
      for line in lines_to_write:
        write_file.write(line)

def pad_step_las_with_early_exit(path_to_workdir, smoothness_bound):
  '''
  Two use cases:
  1) This function is run immediately after the factor base is computed.
  It needs the process ids so it can wait for them to finish. The filenames are also provided
  2) This function is run separately from factor base construction.
  Can assume all the fb processes have finished.
  This function finds the filenames itself
  :param path_to_workdir:
  :return:
  '''
  filenames = []
  for filename in os.listdir(path_to_workdir):
    poly_str = '.poly'
    if not filename.endswith(poly_str):
      continue
    spec_str = filename.rstrip(poly_str) + "_params"
    filenames.append(spec_str)

  #Case 1: if the dir does not exist then create it
  #Case 2: if the dir exists then empty it
  out_dir = path_to_workdir + '/pad'
  procs = []
  mk_proc = subprocess.run(["mkdir", out_dir])
  procs.append(mk_proc)

  for i in range(len(filenames)):
    verbose=True
    proc = pad_run_las(path_to_workdir, filenames[i], verbose)

def pad_run_las(path_to_workdir, params_spec_str, verbose):
  params_full_path = path_to_workdir + params_spec_str
  [las_num_threads_from_params, poly_path, pad_I, pad_mfb, smoothness_bound] = read_sieve_params(params_full_path)
  las_sieve_args = pad_get_las_args(path_to_workdir, params_spec_str, path_to_workdir + '/pad', smoothness_bound,
                                    las_num_threads_from_params, poly_path, pad_I, pad_mfb)
  sieve_args = shlex.split(cado_las_binary + las_sieve_args)
  if verbose:
    print("Running las to pad:")
    print(sieve_args)
  return subprocess.run(sieve_args, capture_output=True)

def read_sieve_params(params_full_path):
  with open(params_full_path, 'r') as fp:
    lines = fp.readlines()
    for line in lines:
      newline = line.replace(' ','')

      smoothness_bound_str="tasks.lim0="
      las_num_threads = "tasks.las.threads="
      poly_str = "tasks.polyselect.import="
      pad_I_str = "tasks.filter.pad.I="
      pad_mfb_str="tasks.filter.pad.mfb0="
      if newline.startswith(las_num_threads):
        las_num_threads_from_params= newline.split(las_num_threads)[1].rstrip()
      if newline.startswith(poly_str):
        poly_path = newline.split(poly_str)[1].rstrip()
      if newline.startswith(pad_I_str):
        pad_I=int(newline.split(pad_I_str)[1].rstrip())
      if newline.startswith(pad_mfb_str):
        pad_mfb=int(newline.split(pad_mfb_str)[1].rstrip())
      if newline.startswith(smoothness_bound_str):
        smoothness_bound=int(newline.split(smoothness_bound_str)[1].rstrip())
  return las_num_threads_from_params, poly_path, pad_I, pad_mfb, smoothness_bound

def pad_get_las_args(path_to_params_dir, params_spec_str, pad_out_dir, smoothness_bound, las_num_threads_from_params,
                     poly_path, pad_I, pad_mfb, path_to_todo):
  '''
  path_to_params_dir: should be terminated with '/'
  Open the parameters file to get ```qrange, I lim0 lpb0 mfb0``` and construct two input strings for las for the SI and pad steps respectively.
  :return:
  '''
  roots0_name = params_spec_str + ".roots0.gz"
  roots0_full_path = path_to_params_dir + roots0_name
  try:
    open(roots0_full_path, "r")
  except:
    raise ValueError("pad_get_las_args: polynomial file with no roots found but I am continuing")
  roots1_name = params_spec_str + ".roots1.gz"
  roots1_full_path = path_to_params_dir + roots1_name
  ap = argparse.ArgumentParser('parse las arguments')
  ap.add_argument('--las_t', nargs='?')
  if pad_out_dir != None:
    pad_out_filename = path_to_params_dir + pad_out_dir + params_spec_str + ".gz".format()

  las_parse = ap.parse_args()

  las_args = ""
  las_args += " -poly {}".format(poly_path)
  las_args += " -lim0 {}".format(smoothness_bound)
  las_args += " -lim1 {}".format(smoothness_bound)
  las_args += " -fb0 {}".format(roots0_full_path)
  las_args += " -fb1 {}".format(roots1_full_path)

  val = las_parse.las_t
  if is_prod():
    las_args += " -t {}".format(las_num_threads_from_params)
  else:
    if val != None and int(val) != las_num_threads_from_params:
      print("WARNING: command line las_t was given and is overriding {}".format(las_num_threads_from_params))
    las_args += " -t {}".format(val)
  if not asserts_on():
    las_args += " -production"

  las_args += " -I {}".format(pad_I)
  las_args += " -lpb0 {}".format(pad_mfb+1) # add 1 to ensure large primes don't happen
  las_args += " -lpb1 {}".format(pad_mfb+1)
  las_args += " -mfb0 {}".format(pad_mfb)
  las_args += " -mfb1 {}".format(pad_mfb)
  las_args += " -todo {}".format(path_to_todo)
  las_args += " -exit-early"
  if pad_out_dir != None:
    las_args += " -out {}".format(pad_out_filename)
  return las_args
