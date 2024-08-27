from sage.all import *
def gray_code_shift(tt):
  '''
  For testing purposes this is the Gray code that starts from 0.
  :param tt:
  :return:
  '''
  gray_vec=[0]*tt
  gray_vecs=[[0]*tt] # in python 0.bit_length()=0 so we handle it separately
  inds_to_change=[]
  sign=[]
  for ii in range(1, 2**tt):
    next_int=ii ^ (ii >> 1)
    gray_vec=[int(d) for d in reversed(format(next_int, 'b'))] + [0]*(tt-next_int.bit_length())
    changed_ind=ZZ(ii+1).valuation(2)
    inds_to_change.append(changed_ind)
    gray_vecs.append(gray_vec)
    pow=(ii+1)/2**(changed_ind+1) # could also be written v_2(2*ii)
    sign.append((-1)**ceil(pow+1))
  return gray_vecs, inds_to_change, sign

class GrayIterator:
  '''An iterator of binary Gray codes as tuples that starts with the [1] vector and also returns the index and the sign to obtain the next tuple.
  The caller should break out of the loop without using the very last index and sign values because this index value is larger than the list.
  Proof of correctness in Contini's thesis p68.
  Note that sign and changed_ind are constructed independently of gray_vec'''
  def __init__(self, tt):
    '''
    :param tt: length of the Gray code
    '''
    self.ii = -1
    self.tt = tt
    self.gray_vec =None
    self.changed_ind=None
    self.sign=None

  def __iter__(self):
    return self

  def __next__(self):
    self.ii += 1
    if self.ii == 0:
      self.gray_vec=[1]*self.tt
    elif self.ii < 2**self.tt:
      next_int=self.ii ^ (self.ii >> 1)
      self.gray_vec=[int(d) for d in reversed(format(next_int, 'b'))] + [0]*(self.tt-next_int.bit_length())
      self.gray_vec=[(_+1)%2 for _ in self.gray_vec]
    else:
      raise StopIteration
    self.changed_ind=ZZ(self.ii+1).valuation(2) # could also be written v_2(2*(ii+1))
    pow=(self.ii+1)/2**(self.changed_ind+1)
    self.sign=(-1)**ceil(pow) # Gray code starting from the 0 vector would use pow+1 instead of pow
    return self.gray_vec, self.changed_ind, self.sign

