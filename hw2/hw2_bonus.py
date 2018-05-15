import numpy as np
import sys
import os
from scipy.stats import norm

sys.setrecursionlimit(6000)

from hw2_util import *


def combinatorial_pricing(S0, K, r, q, sigma, T, n):
    """ Implementation of bonus 2."""
    if DEBUG:
        print("Number of periods: {}".format(n))
        print("S0 = {}, K = {}, r = {}, q = {}, sigma = {}, T = {}".format(S0, K, r, q, sigma, T))

    t = 0.0 # always starts at t = 0

    dT = T/n
    u = np.exp(sigma * np.sqrt(dT))
    d = 1 / u
    p = (np.exp((r-q) * dT) - d)/(u - d)
    #print("p: , u: , d: , dT: \n".format(p, u, d, dT)

    ret_call = 0.0
    ret_put = 0.0

   
    for i in range(n + 1): # i = 0, 1, ..., n
        lnSi = np.log(S0) + (n - i) * np.log(u) + i * np.log(d) # Si at t = i if S0 goes up i times.
        Si = np.exp(lnSi)
        call_payoff_i = calc_call_payoff(Si, K)  # payoff of call at t = i 
        put_payoff_i = calc_put_payoff(Si, K)    # payoff of put at t = i 
        #print("i is: %d, Si is: , call_payoff_i is: , put_payoff_i is: \n".format(i, Si, call_payoff_i, put_payoff_i)

        if call_payoff_i != 0:
                #print("ln_combination is: \n".format(ln_combination(n, i))
                # S0 * u ** (n - i * d ** i
                ret_call += np.exp( ln_combination(n, i) + (n - i) * np.log(p) + i * np.log(1 - p) + np.log(call_payoff_i) )
                #print("ret_call: \n".format(ret_call)
        if put_payoff_i != 0:
                ret_put  += np.exp( ln_combination(n, i) + (n - i) * np.log(p) + i * np.log(1 - p) + np.log(put_payoff_i)  ) 
        #print("ret_call now is: , ret_put now is: \n".format(ret_call, ret_put)
    

    C0 = ret_call * np.exp(-r * T)
    P0 = ret_put * np.exp(-r * T)

    return (C0, P0)


if __name__ == "__main__":

    # read input from stdin

    _input = input("Getting input: S0, K, r, q, sigma, T, num of simulations, num of repetitions, n \n")

    _input = [float(i) for i in _input.split(" ")]

    S0, K, r, q, sigma, T, simulations, repetitions, n = _input
    simulations = int(simulations)
    repetitions = int(repetitions)
    n = int(n)

    print("S0          = {}".format(S0))
    print("K           = {}".format(K))
    print("r           = {}".format(r))
    print("q           = {}".format(q))
    print("sigma       = {}".format(sigma))
    print("T           = {}".format(T))
    print("simulations = {}".format(simulations))
    print("repetitions = {}".format(repetitions))
    print("n           = {}".format(n))

    # Bonus 2
    print_line()
    print("Bonus 2: combinatorial method to price European options\n")
    bonus_call, bonus_put = combinatorial_pricing(S0, K, r, q, sigma, T, n) 
    print("Call price = {0:.10f}\nPut price  = {1:.10f}".format(bonus_call, bonus_put))

