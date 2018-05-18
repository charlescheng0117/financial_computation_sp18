import numpy as np
import sys
import os
from scipy.stats import norm

DEBUG = False
#DEBUG = True

def N(x):
    return norm.cdf(x)

def calc_d1(S0, K, r, q, sigma, T): 
    d1 = (np.log(S0/K) + (r - q + np.power(sigma, 2.0) / 2.0) * T) / (sigma * np.sqrt(T))
    return d1

def calc_d2(S0, K, r, q, sigma, T): 
    d2 = (np.log(S0/K) + (r - q - np.power(sigma, 2.0) / 2.0) * T) / (sigma * np.sqrt(T))
    return d2

def BS_call(S0, K, r, q, sigma, T):
    d1 = calc_d1(S0, K, r, q, sigma, T)
    d2 = calc_d2(S0, K, r, q, sigma, T)
    return S0 * np.exp(-q * T) * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)

def BS_put(S0, K, r, q, sigma, T):
    d1 = calc_d1(S0, K, r, q, sigma, T);
    d2 = calc_d2(S0, K, r, q, sigma, T);
    return K * np.exp(-r * T) * norm.cdf(-d2) - S0 * np.exp(-q * T) * norm.cdf(-d1)

def calc_payoff(ST, K1, K2, K3, K4):
    if ST <= K1 or ST >= K4:
        return 0
    elif K2 <= ST and ST <= K3:
        return K2 - K1
    elif K1 < ST and ST < K2:
        return ST - K1
    else: # K3 < ST < K4
        return K4 - ST

def get_option_value(S0, r, q, sigma, T, K1, K2, K3, K4):

    # calculate d1, d2 for K1, K2, K3, K4
    d1_K1, d2_K1 = calc_d1(S0, K1, r, q, sigma, T), calc_d2(S0, K1, r, q, sigma, T)
    d1_K2, d2_K2 = calc_d1(S0, K2, r, q, sigma, T), calc_d2(S0, K2, r, q, sigma, T)
    d1_K3, d2_K3 = calc_d1(S0, K3, r, q, sigma, T), calc_d2(S0, K3, r, q, sigma, T)
    d1_K4, d2_K4 = calc_d1(S0, K4, r, q, sigma, T), calc_d2(S0, K4, r, q, sigma, T)
    
    value = S0 * np.exp(-q * T) * ( norm.cdf(d1_K1) - norm.cdf(d1_K2) ) \
            - K1 * np.exp(-r * T) * ( norm.cdf(d2_K1) - norm.cdf(d2_K2) ) \
            + (K2 - K1) * np.exp(-r * T) *  ( norm.cdf(d2_K2) - norm.cdf(d2_K3) ) \
            + ( (K1 - K2)/(K4 - K3) ) * ( S0 * np.exp(-q * T) * (norm.cdf(d1_K3) - norm.cdf(d1_K4)) \
                                          - K4 * np.exp(-r * T) * (norm.cdf(d2_K3) - norm.cdf(d2_K4) )
                                        )

    return value

def calc_payoff(ST, K1, K2, K3, K4):
    if ST <= K1 or ST >= K4:
        return 0
    elif K2 <= ST and ST <= K3:
        return K2 - K1
    elif K1 < ST and ST < K2:
        return ST - K1
    else: # K3 < ST < K4
        return ( (K2 - K1)/(K4 - K3) ) * (K4 - ST)

def monte_carlo(S0, K1, K2, K3, K4, r, q, sigma, T, simulations, repetitions):
    print("--------------------------------------------------------------------------------")
    print("Monte Carlo Simution:")
    if DEBUG:
        print("Number of simulations: {}, Number of repetitions: {}".format(simulations, repetitions))
        print("S0 = {}, K1 = {}, K2 = {}, K3 = {}, K4 = {}, r = {}, q = {}, sigma = {}, T = {}".format(S0, K1, K2, K3, K4, r, q, sigma, T))

    # mean, sigma of lnSt
    mean_lnST = np.log(S0) + (r - q - np.power(sigma, 2.0) / 2.0) * T 
    sigma_lnST = sigma * np.sqrt(T) #// sigma of lnST
    
    if DEBUG:
        print("mean_lnST = {}, sigma_lnST = {}".format(mean_lnST, sigma_lnST))

    # draw samples from N(lnST, sigmaST)


    option_results = []
    
    for i in range(repetitions):
        # draw samples from N(lnST, sigmaST)
        samples = np.random.normal(mean_lnST, sigma_lnST, simulations)
        
        #print("samples: {}".format(samples[:5]))
        
        option_vals = []
    
        #for j in range(simulations):
        for s in samples:
            option_vals.append(calc_payoff(np.exp(s), K1, K2, K3, K4) * np.exp(-r * T))

        option_expected_val = np.mean(option_vals)
        option_results.append(option_expected_val)

    option_mean = np.mean(option_results)
    option_std = np.std(option_results)

    print("Option's mean = {}, std = {}".format(option_mean, option_std))
    print("0.95 C.I. for the MonteCarlo Simulation: [{}, {}]".format(option_mean - 2 * option_std, option_mean + 2 * option_std))

def testBS():
    C1, P1 = BS_call(S0, K1, r, q, sigma, T), BS_put(S0, K1, r, q, sigma, T)
    C2, P2 = BS_call(S0, K2, r, q, sigma, T), BS_put(S0, K2, r, q, sigma, T)
    C3, P3 = BS_call(S0, K3, r, q, sigma, T), BS_put(S0, K3, r, q, sigma, T)
    C4, P4 = BS_call(S0, K4, r, q, sigma, T), BS_put(S0, K4, r, q, sigma, T)
    
    if DEBUG:
        print("S0 = {}, K1 = {}, r = {}, q = {}, sigma = {}, T = {}".format(S0, K1, r, q, sigma, T))
        print("call: {}, put: {}".format(C1, P1))

        print("S0 = {}, K2 = {}, r = {}, q = {}, sigma = {}, T = {}".format(S0, K2, r, q, sigma, T))
        print("call: {}, put: {}".format(C2, P2))
        
        print("S0 = {}, K3 = {}, r = {}, q = {}, sigma = {}, T = {}".format(S0, K3, r, q, sigma, T))
        print("call: {}, put: {}".format(C3, P3))
    
        print("S0 = {}, K4 = {}, r = {}, q = {}, sigma = {}, T = {}".format(S0, K4, r, q, sigma, T))
        print("call: {}, put: {}".format(C4, P4))

    option = C1 - C2 - C3 + C4

used_sys_arg = False

if __name__ == "__main__":
    #S0, K1, K2, K3, K4, r, q, sigma, T = 50, 40, 50, 60, 70, 0.05, 0.08, 0.2, 1
    #S0, r, q, sigma, T, K1, K2, K3, K4= 50, 0.05, 0.08, 0.2, 1, 40, 50, 60, 70
    
    if used_sys_arg:
        _input = sys.argv[1:]
        _input = [float(args) for args in _input]
        S0, r, q, sigma, T, K1, K2, K3, K4 = _input
        if DEBUG:
            print(_input)

        if len(_input) != 9:
            print("usage: python hw1.py S0 r q sigma T K1 K2 K3 K4")
            sys.exit()
    
    _input = input("Getting input...\n") 
    
    _input = _input.split(" ")
    _input = [float(args) for args in _input]

    S0, r, q, sigma, T, K1, K2, K3, K4 = _input
    
    print("S0    = {}".format(S0))
    print("r     = {}".format(r))
    print("q     = {}".format(q))
    print("sigma = {}".format(sigma))
    print("T     = {}".format(T))
    print("K1    = {}".format(K1))
    print("K2    = {}".format(K2))
    print("K3    = {}".format(K3))
    print("K4    = {}".format(K4))

    if DEBUG:
        C1, P1 = BS_call(S0, K1, r, q, sigma, T), BS_put(S0, K1, r, q, sigma, T)
        C2, P2 = BS_call(S0, K2, r, q, sigma, T), BS_put(S0, K2, r, q, sigma, T)
        C3, P3 = BS_call(S0, K3, r, q, sigma, T), BS_put(S0, K3, r, q, sigma, T)
        C4, P4 = BS_call(S0, K4, r, q, sigma, T), BS_put(S0, K4, r, q, sigma, T)
    
        print("S0 = {}, K1 = {}, r = {}, q = {}, sigma = {}, T = {}".format(S0, K1, r, q, sigma, T))
        print("call: {}, put: {}".format(C1, P1))

        print("S0 = {}, K2 = {}, r = {}, q = {}, sigma = {}, T = {}".format(S0, K2, r, q, sigma, T))
        print("call: {}, put: {}".format(C2, P2))
        
        print("S0 = {}, K3 = {}, r = {}, q = {}, sigma = {}, T = {}".format(S0, K3, r, q, sigma, T))
        print("call: {}, put: {}".format(C3, P3))
    
        print("S0 = {}, K4 = {}, r = {}, q = {}, sigma = {}, T = {}".format(S0, K4, r, q, sigma, T))
        print("call: {}, put: {}".format(C4, P4))

        value_BS = C1 - C2 - C3 + C4
        print("BS value = {}".format(value_BS))
    
    
    value = get_option_value(S0, r, q, sigma, T, K1, K2, K3, K4)
    print("option value = {}".format(value))

    simulations = 10000
    repetitions = 20

    monte_carlo(S0, K1, K2, K3, K4, r, q, sigma, T, simulations, repetitions)




    
