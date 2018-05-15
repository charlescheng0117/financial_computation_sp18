import numpy as np
import sys
import os
from scipy.stats import norm

DEBUG = False

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

def calc_call_payoff(ST, K): 
    return max(ST - K, 0.0)

def calc_put_payoff(ST, K): 
    return max(K - ST, 0.0)

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
