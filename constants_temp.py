# setup constants by copying this file to `constants.py` and fill in

# constants to configure for each machine
is_magma_available =
is_gp_available =

# fixed constants
sieve_algo_pairs_of_deg_1 = 'sieve_algo_pairs_of_deg_1'
sieve_algo_all_prs = 'sieve_algo_all_primes'
hnf_algo_sage = "sage"
hnf_algo_magma = "magma"
hnf_algo_gp = "gp"

idl_coeff_gcd_str = "idl_coeff_gcd"
idl_least_int_str = "idl_least_int"
idl_poly_gen_str = "idl_poly_gen"
idl_fb_exp_vec_str = "fb_exp_vec"

task_setup_gen_polys = "gen polys"
task_setup_factor_base = "gen polys"
task_sieve_sieve = "sieve"
task_filter_concat = "TASK: gzip and concat files"
task_filter_load_mat = "TASK: load matrix"
task_filter_dup1 = "TASK: remove dups / divide out divisors"
task_filter_merge = "TASK: merging large primes"
task_filter_dup2 = "TASK: remove duplicates"
task_filter_triv = "TASK: add triv rels"
task_filter_singleton = "TASK: remove singletons"
task_filter_zero = "TASK: remove zero cols"
task_filter_write = "TASK: check the matrix is full rank"
task_filter_targeted_relations = "TASK: generate more relations"
task_filter_save_mat = "TASK: save matrix"
task_linalg_class_num = "class num"
task_linalg_modhnf = "mod hnf"
task_linalg_unit_lattice_basis = "unit lattice basis"
task_linalg_check_hR = "check if we have generates of extended relations lattice"
task_linalg_targeted_relations = "generate more relations"
task_finalize = "smith form"
task_finalize = "prove alg primes in (smooth_bound, GRH_bound] are in the subgroup"
