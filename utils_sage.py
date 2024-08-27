from functools import cmp_to_key
def sort_idl(idl, other_idl):
  '''
  :param mat: cypari mat
  :param other_mat: cypari mat
  :return:
  '''
  mat = idl.pari_hnf()
  other_mat = other_idl.pari_hnf()
  ncols=mat.length()
  nrows=mat[0].length()
  ncols_other=mat.length()
  nrows_other=mat[0].length()
  if ncols != ncols_other:
    raise ValueError("Number of cols should be equal")
  if nrows != nrows_other:
    raise ValueError("Number of rows should be equal")
  for rr in range(nrows):
    for cc in range(ncols):
      if mat[rr,cc] < other_mat[rr,cc]:
        return True
      elif mat[rr,cc] > other_mat[rr,cc]:
        return False
  return True

def sorted_primes_above(nf, pr):
  return sorted(nf.primes_above(pr), key=cmp_to_key(sort_idl))
