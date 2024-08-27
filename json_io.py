from src.obj_string_converters import *
from src.paths import basepath_artifacts

mat_idl_dict_key = "mat"
gen0_idl_dict_key = "gen0"
gen1_idl_dict_key = "gen1"
residue_class_degree_idl_dict_key = "residue_class_degree"
ramification_index_idl_dict_key = "ramification_index"

def idl_to_idl_str_dict(idl):
  '''
  Creates a dict with the following keys. This dict is intended to be serialized to json format
  gen0: ZZ, this is the first elem in the 2-elem rep
  gen1: [ZZ], this is the second elem in the 2-elem rep
  mat: ["[ZZ]"] these are the rows of the ideal's basis matrix
  :param idl:
  :return:
  '''
  gens = idl.gens_two()
  return {mat_idl_dict_key : mat_to_str(idl.pari_hnf().sage().transpose()),
          gen0_idl_dict_key : str(gens[0]),
          gen1_idl_dict_key : nf_elem_to_str(gens[1]),
          residue_class_degree_idl_dict_key : int(idl.residue_class_degree()),
          ramification_index_idl_dict_key : int(idl.ramification_index())}

def fb_data_to_str_dict(idls, rat_primes, num_facs_above_prime, smoothness_bound):
  '''
  Intended use: writing json
  :param idls:
  :param rat_primes:
  :param num_facs_above_prime:
  :param smoothness_bound:
  :return:
  '''
  idls_dicts = [idl_to_idl_str_dict(_) for _ in idls]
  return {"prime_idls" : idls_dicts,
          "rat_primes" : [str(_) for _ in rat_primes],
          "num_facs_above_prime" : [str(_) for _ in num_facs_above_prime],
          "smoothness_bound" : smoothness_bound}

def nf_to_str_dict(nf):
  '''
  Intended use: writing json
  :param idls:
  :return:
  '''
  pari_basis = nf.pari_zk()
  precomputed_int_bas_inv = Matrix(QQ, [list(nf(_)) for _ in pari_basis])**(-1)
  precomputed_int_bas_inv_num, precomputed_int_bas_inv_dnm = mat_to_num_denom(precomputed_int_bas_inv)
  return {"def_poly" : poly_to_str(nf.polynomial()),
          "Delta" : str(nf.discriminant().abs()),
          "int_basis_mat_inv_num" : mat_to_str(precomputed_int_bas_inv_num),
          "int_basis_mat_inv_dnm" : str(precomputed_int_bas_inv_dnm),
          "pari_int_bas" : {"basis" : [nf_elem_to_str(elem) for elem in pari_basis]}}

def write_ic_and_fb_json(idls, sp_rat_prs, smoothness_bound, ind_calc_path, nf_path, fb_path):
  '''
  :param idls: a list of idls
  :param ind_calc_path TODO implement if needed
  :param nf_path
  :param fb_path:
  :return:
  '''
  if len(idls) == 0:
    raise ValueError("list of ideals must contain at least one ideal")
  nf = idls[0].number_field()
  with open(nf_path, 'w') as ic_fp:
    nf_dict = nf_to_str_dict(nf)
    json.dump(nf_dict, ic_fp)
  with open(fb_path, 'w') as fb_fp:
    num_primes_above = [len(nf.ideal(_).factor()) for _ in sp_rat_prs]
    fb_dict = fb_data_to_str_dict(idls, sp_rat_prs, num_primes_above, smoothness_bound)
    json.dump(fb_dict, fb_fp)

def read_ic_and_fb_json(read_ic_rel_path, read_fb_rel_path, nf):
  '''
  This function returns a list of pairs ```gens```, a list of mats ```basis_mats```, and a list of nf_elems ```pari_int_bas```
  :param read_ic_rel_path:
  :param read_fb_rel_path:
  :param nf:
  :return:
  '''
  is_int_bas_loaded = False
  dg = nf.degree()
  ic_path = basepath_artifacts + read_ic_rel_path
  while not os.path.exists(ic_path):
    time.sleep(1)
  with open(ic_path, 'r') as fp:
    data = fp.read()
    strs_dict = json.loads(data)
    pari_int_bas = str_dict_to_ind_calc_data(strs_dict['pari_int_bas'], nf)
  with open(basepath_artifacts + read_fb_rel_path, 'r') as fp_fb:
    fb_strs_dict = json.loads(fp_fb.read())
    gens_and_mats_list = [str_dict_to_idl_data(_, nf) for _ in fb_strs_dict]
  gens_list = [_[0] for _ in gens_and_mats_list]
  bas_mats_list = [_[1] for _ in gens_and_mats_list]

  return gens_list, bas_mats_list, pari_int_bas
