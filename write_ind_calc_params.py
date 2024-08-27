from src.polredabs_polys import *
from src.obj_string_converters import *
from src.stage_init_internal import *
from src.json_io import *
'''
The purpose of this file is to write files that can be used by the cpp version of the code
'''

# a list of all predefined polys
polys = [b97, b102, b103, b104, b105, b109, b110, b111, b112, b113, b114]

def b103_B10000():
  smoothness_bound=10000
  nf = NumberField(b103, names='Y')
  [_, _, fb, rat_prs, _, _, _] = fac_base_data(nf, smoothness_bound)
  basepath = "/Users/davidmarquis/git/thesis/ANTL/tests/DataForTests/IndCalcParams/"
  fb_path = "{}{}_B{}_fb.json".format(basepath, "b103", smoothness_bound)
  nf_path = "{}{}_B{}_nf.json".format(basepath, "b103", smoothness_bound)
  ind_calc_path = None
  write_ic_and_fb_json(fb, rat_prs, smoothness_bound, ind_calc_path, nf_path, fb_path)

b103_B10000()