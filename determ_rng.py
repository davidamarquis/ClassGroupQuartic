from sage.all import * #Dont move this, Sage clobbers the seed function in Python with its own
from random import seed
from random import sample
from random import randrange
import src.globals

def increment_seed():
  src.globals.curr_seed += 1
  seed(src.globals.curr_seed)

def set_seed(new_seed):
  src.globals.curr_seed = new_seed
  seed(src.globals.curr_seed)

def get_seed():
  return src.globals.curr_seed

def get_randrange_and_increment_seed(start, stop, step=1):
  increment_seed()
  # Return a randomly selected element from range(start, stop, step) and increment the global seed
  r_int = randrange(start, stop, step)
  return r_int

def get_sample_and_increment_seed(population, k):
  increment_seed()
  samp = sample(population, k)
  return samp

def get_rand_lst(start, stop, len_lst):
  increment_seed()
  lst = []
  for i in range(len_lst):
    lst.append(get_randrange_and_increment_seed(start, stop))
  return lst

def get_random_and_increment_seed():
  increment_seed()
  return randrange(0, 2**32-1)

def get_random_zero_one_list_and_increment_seed(vec_len, num_ones):
  increment_seed()
  samp = sample(range(0,vec_len), num_ones)
  lst = [0]*vec_len
  for _ in samp:
    lst[_]=1
  return lst