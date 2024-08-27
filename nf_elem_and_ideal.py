from sage.all import *
from collections import defaultdict
from src.utils import *

def perm_elem(perm, elem, primes):
  '''
  For a1 feb 2022
  Permutes an elem of the number field
  TODO validate that the length of the permutation is the same as length of the action
  :param perm: sage.groups.perm_gps.permgroup_element.PermutationGroupElement
  :param elem:
  :param primes:
  :return:
  '''
  if len(perm.domain()) != len(primes):
    raise ValueError('{}\nNumber of gens of permutations {} must equal number of primes {}'.format(perm.domain(), perm.parent().ngens(), len(primes)))
  bas = idl_crt_basis(primes)
  perm_bas = vector(Permutation(perm).action(bas))
  # is reduce_at_primes here analogous to what SolvePell does by reducing at B ? I dont think so
  # return perm_bas.dot_product(reduce_at_primes(elem, primes))
  return perm_bas.dot_product(vector([elem]*len(primes)))

def cyc_pow(grp_size, pow):
  '''
  For a1 feb 2022
  :param grp_size:
  :param pow:
  :return: CyclicPermutationGroup(grp_size) direct producted with itself pow times
  '''
  G = CyclicPermutationGroup(grp_size)
  for _ in range(pow-1):
    G = CyclicPermutationGroup(grp_size).direct_product(G, False)
  return G

def enumerate_image_of_action_of_the_power_of_C2_group(elem, primes):
  '''
  For a1 feb 2022
  enumerates the image of the action of the power of C2 group
  :param elem:
  :param primes: it is assumed that the list has even length
  :return: [sigma_j(elem) (mod Z)]_{j=0}^{2^tt-1} where Z is the product of all primes in primes and sigma_i is in Aut(dirprod_{i=0}^{t-1} Ok/qfrak_i / Z/q_i)
  '''
  if not len(primes) % 2 == 0:
    raise ValueError('primes must have even length and the primes at 2i and 2i+1 should be conjugate')
  tt = len(primes)//2
  return [perm_elem(perm, elem, primes) for perm in cyc_pow(2, tt).list()]

def enumerate_image_of_action_of_the_power_of_C2_group_idl(idl_bas, primes):
  '''
  For a1 feb 2022
  enumerates the image of the action of the power of C2 group
  :param idl_bas:
  :param primes: it is assumed that the list has even length
  :return: [sigma_j(elem) (mod Z)]_{j=0}^{2^tt-1} where Z is the product of all primes in primes and sigma_i is in Aut(dirprod_{i=0}^{t-1} Ok/qfrak_i / Z/q_i)
  '''
  return [enumerate_image_of_action_of_the_power_of_C2_group(_, primes) for _ in idl_bas]

def prod_even(elems):
  '''
  For a1 feb 2022
  Returns the prod of every 2nd entry
  :param elems:
  :return:
  '''
  if len(elems) == 0:
    return None

  prod = 1
  for i in range(0, len(elems), 2):
    prod = prod * elems[i]
  return prod

def bas_elems_list_to_mat(nf, idl_bas):
  '''
  For a1 feb 2022
  :param nf:
  :param idl_bas:
  :return:
  '''
  RR = nf.maximal_order()
  return Matrix(ZZ, [list(RR.coordinates(y)) for y in idl_bas]).transpose()

def reduce_at_primes(elem, primes):
  '''
  For a1 feb 2022
  :param elem:
  :param primes:
  :return:
  '''
  return vector([pr.reduce(elem) for pr in primes])

def idl_crt_basis(primes):
  '''
  For a1 feb 2022
  :param primes:
  :return:
  '''
  nn = len(primes)
  if nn == 0:
    raise ValueError('len must be longer than 0')
  nf = primes[0].number_field()
  std_bas = VectorSpace(QQ, nn).basis()
  return vector([nf.idealchinese(primes, vv) for vv in std_bas])