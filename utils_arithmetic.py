import time
from sage.all import *
import src.timers as timers

def nonnegative_knapsack(seq, sum_max, binary, max_num_primes, solver=None):
  '''
  sage.numerical.knapsack.knapsack gives wrong values when value_only is False so I had to write this.
  Ex for
  seq=[(1,1),(2,2),(3,3),(4,4)]
  binary=False
  value_only=False
  max=5
  verbose=True
  sage.numerical.knapsack.knapsack gives
  [5.0, [(2, 2), (4, 4)]]

  Maximizes sum{b_i*u_i} subject to sum{b_i*w_i}<=max
  :param seq: list of pairs (w_i, u_i)
  :param sum_max: positive integer
  :param solver: if None then the default solver is used
  :return:
  '''
  from sage.numerical.mip import MixedIntegerLinearProgram
  tolerance=1e-3 # If the variable value differs from the nearest integer by more than ``tolerance``,raise a ``RuntimeError``.
  p = MixedIntegerLinearProgram(solver=solver, maximization=True)
  if binary:
    present = p.new_variable(binary = True, nonnegative=True)
  else:
    present = p.new_variable(integer = True, nonnegative=True)
  p.set_objective(p.sum([present[i] * seq[i][1] for i in range(len(seq))]))
  p.add_constraint(p.sum([present[i] * seq[i][0] for i in range(len(seq))]), max=sum_max)
  if max_num_primes != None:
    p.add_constraint(p.sum([present[i] for i in range(len(seq))]), max=max_num_primes, min=max_num_primes)
  objective = p.solve()
  present = p.get_values(present, convert=ZZ, tolerance=tolerance)
  return [objective,list(present.values())]

def closest_integer_power_prod(pr_pows, target, prec, binary, max_num_primes):
  '''
  :param pr_pows:
  :param target:
  :param prec:
  :param binary:
  :param max_num_primes:
  :return:
  '''
  short_pr_pows_list = sorted(list(set(pr_pows)))
  start_seq=[]
  reals=RealField(prec)
  rr_pr_pows=[reals(_) for _ in short_pr_pows_list]
  for pr in rr_pr_pows:
    start_seq.append(pr.log(2))
  start_seq_round=[up_round(_,prec) for _ in start_seq]
  knapsack_pow=0.2
  pr_pow_utilities = [_**knapsack_pow for _ in start_seq_round]
  start_seq_round_pairs=list(zip(start_seq_round, pr_pow_utilities))

  bnd = reals(max(target,2))
  log_bnd=up_round(bnd.log(2), prec)
  [_, exps] = nonnegative_knapsack(start_seq_round_pairs, log_bnd, binary, max_num_primes)
  exps_in_full_list = [0]*len(pr_pows)
  for i in range(len(exps)):
    if exps[i] != 0:
      pr_to_find=short_pr_pows_list[i]
      ind=pr_pows.index(pr_to_find)
      exps_in_full_list[ind]=exps[i]

  prod=1
  for i in range(len(pr_pows)):
    prod*= pr_pows[i]**exps_in_full_list[i]
  prs_used=[pr_pows[_] for _ in range(len(pr_pows)) if exps_in_full_list[_]!=0]
  return [prod, exps_in_full_list, prs_used]

def up_round(x, prec):
  pow= ZZ(2)**prec
  return ZZ(pow * x)

def down_round(x, prec):
  return x/ZZ(2)**prec

def construct_idl_with_close_norm(bounded_fb, target, prec, binary, max_num_primes):
  '''
  Convenience function for construct_idl_with_prescribed_divisor_and_close_norm()
  :param bounded_fb:
  :param target:
  :param prec:
  :return:
  '''
  nf=bounded_fb[0].number_field()
  prescribed_div=nf.ideal(1)
  bounded_nrms = [_.norm() for _ in bounded_fb]
  bounded_rat_prs = [_.norm()**(_.residue_class_degree()**(-1)) for _ in bounded_fb]
  return construct_idl_with_prescribed_divisor_and_close_norm(bounded_fb, bounded_nrms, bounded_rat_prs, target, prec,
                                                              prescribed_div, binary, max_num_primes)

def construct_idl_with_prescribed_divisor_and_close_norm(pr_idls, pr_idls_nrms, pr_idls_rat_prs, target, prec, prescribed_div, binary, max_num_primes):
  '''
  Finds an ideal that's a product of the given list of ideals and has norm close the target.
  This ideal is found by solving the knapsack problem.
  :param pr_idls:
  :param pr_idls_nrms:
  :param pr_idls_rat_prs:
  :param target:
  :param prec: the precision of the knapsack problem
  :param prescribed_div: the ideal that must divide the product
  :return:
  '''
  nf=pr_idls[0].number_field()
  # the list of nrms will have repetitions
  idl_one=nf.ideal(1)
  if prescribed_div != idl_one:
    prescribed_div_nrm=prescribed_div.norm()
    # prescribed_div_rat_pr=prescribed_div.smallest_integer()
  else:
    prescribed_div_nrm=ZZ(1)
    prescribed_div_rat_pr=ZZ(1)

  [idl_nrm, exps, _] = closest_integer_power_prod(pr_idls_nrms, target/prescribed_div_nrm, prec, binary, max_num_primes)

  idl = nf.ideal(1)
  for i in range(len(pr_idls)):
    if exps[i]!=0:
      idl *= pr_idls[i]**exps[i]
  idl_with_prescribed_divisor=idl*prescribed_div
  idl_nrm_with_prescribed=idl_nrm*prescribed_div_nrm

  if prescribed_div != idl_one:
    try:
      ind=pr_idls.index(prescribed_div)
    except ValueError:
      print("prescribed_div was not found in the list of prime ideals. Exponent vector should not be used.")
    exps[ind]=1
  rat_prs_used=set(pr_idls_rat_prs[_] for _ in range(len(pr_idls_rat_prs)) if exps[_] != 0)
  rat_prs_used=sorted(list(rat_prs_used))
  # prs_used=set(pr_idls[_] for _ in range(len(pr_idls)) if exps[_] != 0)
  prs_used=[pr_idls[_] for _ in range(len(pr_idls)) if exps[_] != 0]
  return [idl_with_prescribed_divisor, idl_nrm_with_prescribed, prs_used, exps, rat_prs_used]
