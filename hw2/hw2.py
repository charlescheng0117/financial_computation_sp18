import numpy as np
import sys
import os
from scipy.stats import norm

DEBUG = True

def print_line():
    print("--------------------------------------------------------------------------------")


def N(x):
    return norm.cdf(x)

def calc_d1(S0, K, r, q, sigma, T): 
    d1 = (np.log(S0/K) + (r - q + np.power(sigma, 2.0) / 2.0) * T) / (sigma * np.sqrt(T))
    return d1

def calc_d2(S0, K, r, q, sigma, T): 
    d2 = (np.log(S0/K) + (r - q - np.power(sigma, 2.0) / 2.0) * T) / (sigma * np.sqrt(T))
    return d2

def Black_Scholes_call(S0, K, r, q, sigma, T):
    d1 = calc_d1(S0, K, r, q, sigma, T)
    d2 = calc_d2(S0, K, r, q, sigma, T)
    return S0 * np.exp(-q * T) * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)

def Black_Scholes_put(S0, K, r, q, sigma, T):
    d1 = calc_d1(S0, K, r, q, sigma, T)
    d2 = calc_d2(S0, K, r, q, sigma, T)
    return K * np.exp(-r * T) * norm.cdf(-d2) - S0 * np.exp(-q * T) * norm.cdf(-d1)

def calc_call_payoff_payoff(ST, K): 
    return max(ST - K, 0.0)

def calc_put_payoff(ST, K): 
    return max(K - ST, 0.0)

def monte_carlo(S0, K, r, q, sigma, T, simulations, repetitions):
    
    if DEBUG:
        print("Number of simulations: , Number of repetitions: ".format(simulations, repetitions))
        print("S0 = , K = , r = , q = , sigma = , T = ".format(S0, K, r, q, sigma, T))

    # mean, sigma of lnST
    mean_lnST = np.log(S0) + (r - q - np.power(sigma, 2.0) / 2.0) * T 
    sigma_lnST = sigma * np.sqrt(T)
    
    if DEBUG:
        print("mean_lnST = , sigma_lnST = ".format(mean_lnST, sigma_lnST))

    # draw samples from N(lnST, sigmaST)
    call_option_results = []
    put_option_results = []
    
    for i in range(repetitions):
        # draw samples from N(lnST, sigmaST)
        samples = np.random.normal(mean_lnST, sigma_lnST, simulations)
        
        print("samples: ".format(samples[:5]))
        
        call_option_vals = []
        put_option_vals = []

        #for j in range(simulations):
        for s in samples:
            ST = np.exp(s)
            call_option_vals.append(calc_call_payoff_payoff(ST, K) * np.exp(-r * T))
            put_option_vals.append(calc_put_payoff(ST, K) * np.exp(-r * T))

        call_option_np.expected_val = np.mean(call_option_vals)
        put_option_np.expected_val = np.mean(put_option_vals)
        
        call_option_results.append(call_option_np.expected_val)
        put_option_results.append(put_option_np.expected_val)

    call_option_mean = np.mean(call_option_results)
    call_option_std = np.std(call_option_results)
    
    put_option_mean = np.mean(put_option_results)
    put_option_std = np.std(put_option_results)

    print("Call Option's mean = , std = ".format(call_option_mean, call_option_std))
    print("0.95 C.I. for the Monte Carlo Simulation: [, ]\n".format(call_option_mean - 2 * call_option_std,
                                                                      call_option_mean + 2 * call_option_std))
    print("Put Option's mean = , std = ".format(put_option_mean, put_option_std))
    print("0.95 C.I. for the Monte Carlo Simulation: [, ]".format(put_option_mean - 2 * put_option_std,
                                                                      put_option_mean + 2 * put_option_std))


def factorial(n):
    """ return n! """
    if (n == 1 or n == 0):
        return 1
    return np.exp( np.log(n) + np.log( factorial(n - 1) ) )

def combination(n , i):
    """ return nCi """
    if (n < i):
        return -1
    return factorial(n) / (factorial(i) * factorial(n - i))

def ln_factorial(n):
    """ return ln(n!)"""
    if (n == 1 or n == 0):
        return 0
    return np.log(n) + ln_factorial(n - 1)


def ln_combination(n, i):
    """ return ln(nCi)"""
    if (n < i):
        return -1
    return ln_factorial(n) - (ln_factorial(i) + ln_factorial(n - i))

def CRR_binomial(S0, K, r, q, sigma, T, n):
    # n: # of periods
    print("\nCRR binomial tree model.\n")
    
    if DEBUG:
        print("Number of periods: ".format(n))
        print("S0 = , K = , r = , q = , sigma = , T = ".format(S0, K, r, q, sigma, T))

    t = 0.0 # always starts at t = 0

    dT = T/n
    u = np.exp(sigma * np.sqrt(dT))
    d = 1/u
    p = ( np.exp((r-q) * dT) - d )/(u - d)
    
    if DEBUG:
        print("p: , u: , d: , dT: ".format(p, u, d, dT))

    ret_call = 0.0
    ret_put = 0.0
    results = []

    temp_Si = np.zeros(n + 1) 
    
    # european Ci, Pi
    temp_Ci_euro = np.zeros(n + 1) 
    temp_Pi_euro = np.zeros(n + 1) 
    
    # american Ci, Pi
    temp_Ci_amer = np.zeros(n + 1) 
    temp_Pi_amer = np.zeros(n + 1) 

    for i in range(n + 1):
    #for (int i = 0 i <= n i += 1)  #leafs values
        Si = S0 * np.power(u, n - i) * np.power(d, i) #Si at t = i if S0 goes up n - i times.
        temp_Si[i] = Si
        
        Ci = calc_call_payoff(Si, K)
        Pi = calc_put(Si, K)

        temp_Ci_euro[i] = Ci
        temp_Pi_euro[i] = Pi

        temp_Ci_amer[i] = Ci
        temp_Pi_amer[i] = Pi
    
    
    # from i = n-1, n-2. ..., 0
    for i in range(n - 1, -1, -1):
        # from j = 0, 1, ..., i
        for j in range(i + 1):
    #for (int i = n - 1 i >= 0 i -= 1)  # from n - 1, n - 2, ..., 0
        #for (int j = 0 j <= i j += 1)  # c(i, j) = e^(-r dt) * (p * c(i + 1, j) + (1 - p) * c(i + 1, j + 1))
            # Si = p * temp[j] + (1 - p) * temp[j + 1]
            # temp_Si[j] = Si * np.exp(-r * T/n)
            Sij = S0 * pow(u, i - j) * pow(d, j)
            Cu_euro = temp_Ci_euro[j]
            Cd_euro = temp_Ci_euro[j + 1]
            
            Pu_euro = temp_Pi_euro[j]
            Pd_euro = temp_Pi_euro[j + 1]		

            temp_Ci_euro[j] = np.exp(-r * dT) * (p * Cu_euro + (1 - p) * Cd_euro) 
            temp_Pi_euro[j] = np.exp(-r * dT) * (p * Pu_euro + (1 - p) * Pd_euro)

            Cu_amer = temp_Ci_amer[j]
            Cd_amer = temp_Ci_amer[j + 1]
            
            Pu_amer = temp_Pi_amer[j]
            Pd_amer = temp_Pi_amer[j + 1]	

            rv_Ci_amer = np.exp(-r * dT) * (p * Cu_amer + (1 - p) * Cd_amer)  # retension value of Call
            rv_Pi_amer = np.exp(-r * dT) * (p * Pu_amer + (1 - p) * Pd_amer)  # retension value of Put
            
            if DEBUG:
                print("Exercise value of call: ".format(calc_call_payoff(Sij, K)))
                print("Current value of call:  ".format(temp_Ci_euro[j]))
            
            temp_Ci_amer[j] = max(rv_Ci_amer, calc_call_payoff(Sij, K))
            temp_Pi_amer[j] = max(rv_Pi_amer, calc_put_payoff(Sij, K))
            
            if DEBUG:
                print("Value of American call:  ".format(temp_Ci_amer[j]))
                print("Value of American call:  ".format(temp_Ci_amer[j]))
            
    
    C0_euro = temp_Ci_euro[0]
    P0_euro = temp_Pi_euro[0]
    C0_amer = temp_Ci_amer[0]
    P0_amer = temp_Pi_amer[0]

    return (C0_euro, P0_euro, C0_amer, P0_amer)


def bonusCRRBinomial(S0, K, r, q, sigma, T, n): 
    print("\nBonus CRRBinomial pricing.\n")
    print("Number of periods: {}\n".format(n))
    print("S0 = {}, K = {}, r = {}, q = {}, sigma = {}, T = {}\n".format(S0, K, r, q, sigma, T))

    t = 0.0 # always starts at t = 0

    dT = T/n
    u = np.exp(sigma * sqrt(dT))
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
                #print("lnComb is: \n".format(lnComb(n, i))
                ret_call += np.exp( lnComb(n, i) + (n - i) * np.log(p) + i * np.log(1 - p) + np.log(call_payoff_i) ) # S0 * u^(n - i) * d^(i)
                #print("ret_call: \n".format(ret_call)
        if put_payoff_i != 0:
                ret_put  += np.exp( lnComb(n, i) + (n - i) * np.log(p) + i * np.log(1 - p) + np.log(put_payoff_i)  ) 
        #print("ret_call now is: , ret_put now is: \n".format(ret_call, ret_put)
    

    C0 = ret_call * np.exp(-r * T)
    P0 = ret_put * np.exp(-r * T)

    return (C0, P0)





if __name__ == "__main__":

    # read input from stdin

    _input = input("Getting input: S0, K, r, q, sigma, T, num of simulations, num of repetitions, n \n")

    _input = [float(i) for i in _input.split(" ")]
    print(_input)

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

    print_line()
    # Basic requirement
    print("Basic requirement:\n")
    print("Black Scholes formulas: ")
    call_price = Black_Scholes_call(S0, K, r, q, sigma, T)
    put_price  = Black_Scholes_put(S0, K, r, q, sigma, T)
    print("Call price = {}\nPut price  = {}".format(call_price, put_price))
    
    print_line()
    print("Monte Carlo simulation:\n")
    monte_carlo(S0, K, r, q, sigma, T, simulations, repetitions) 

    print_line()
    print("CRR binomial tree model:\n")
    call_euro, put_euro, call_amer, put_amer = CRR_binomial(S0, K, r, q, sigma, T, n)

    print("European Call = {}, Put = {}".format(call_euro, put_euro))
    print("American Call = {}, Put = {}".format(put_amer, put_amer))




