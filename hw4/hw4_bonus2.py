import numpy as np
import os
import sys

if len(sys.argv) != 2:
    print("Usage: python3 hw4_bonus2.py <input>")
    sys.exit()

in_f = open(sys.argv[1], "r")
line = in_f.readline()
line = line.strip().split(" ")
St, r, q, sigma, t, T, S_max_t, n, num_simulations, num_repetitions = [float(ele) for ele in line]
n = int(n)
num_simulations = int(num_simulations)
num_repetitions = int(num_repetitions)

print("Bonus 2: Implement the method: Cheuk and Vorst (1997) ")
print("St          = " + str(St))
print("r           = " + str(r))
print("q           = " + str(q))
print("sigma       = " + str(sigma))
print("t           = " + str(t))
print("T           = " + str(T))
print("S_max_t     = " + str(S_max_t))
print("n           = " + str(n))
print("simulations = " + str(num_simulations))
print("repetitions = " + str(num_repetitions))

dT = (T - t) / n
u = np.exp(sigma * np.sqrt(dT))
d = 1/u
mu = np.exp(r * dT)

p_u = (mu * u - 1) / (mu * (u - d) )
#p_d = 1 - p_u

# last layer
list_payoff_eu = []
list_payoff_am = []

for j in range(n + 1): # calculate the payoff from top to bottom
    payoff_j = np.power(u, n - j)
    list_payoff_eu.append(np.max(payoff_j - 1, 0))
list_payoff_am = list_payoff_eu[:]
#print(list_payoff)

for i in range(n - 1, -1, -1): # period n - 1, n - 2, ..., 0

    for j in range(i): # depth j has i nodes, and again from top to bottom
         payoff_eu_j = ( (1 - p_u) * list_payoff_eu[j] + p_u * list_payoff_eu[j + 2] )
         list_payoff_eu[j] = payoff_eu_j

         payoff_am_j = ( (1 - p_u) * list_payoff_am[j] + p_u * list_payoff_am[j + 2] )
         #print(payoff_am_j)
         #print(np.max( np.power(u, n - j) - 1, 0))
         #print(max(  np.max( np.power(u, n - j) - 1, 0 ), payoff_am_j))
         list_payoff_am[j] = max( payoff_am_j, np.max( np.power(u, i - j) - 1, 0) )
    
    j = i
    payoff_eu_j = ( (1 - p_u) * list_payoff_eu[j] + p_u * list_payoff_eu[j + 1] )
    list_payoff_eu[j] = payoff_eu_j
    payoff_am_j = ( (1 - p_u) * list_payoff_am[j] + p_u * list_payoff_am[j + 1] )
    list_payoff_am[j] = max( payoff_am_j, np.max( np.power(u, i - j) - 1, 0) )

put_eu = list_payoff_eu[0] * St
put_am = list_payoff_am[0] * St

print("### Answer for Bonus 2 ###")
print("European: {}".format(put_eu))
print("American: {}".format(put_am))
