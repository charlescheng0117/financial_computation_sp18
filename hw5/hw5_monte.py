import numpy as np

import numpy as np
import random

def exp(x):
	return np.e ** x

if __name__ == '__main__':
	file_name = "hw5_input.txt"
	inputs = list(open(file_name))[0]
	inputs = [float(x) for x in inputs.split()]

	S_t, K, r, q, sigma, T_minus_t, M, n, S_ave_t, passing_time, simulations, repetitions = inputs
	M, n, simulations, repetitions = int(M), int(n), int(simulations), int(repetitions)

	print("Monte Carlo simulation to arithmetic average call.")
	print("S_t = {}, K = {}, r = {}, q = {}, sigma = {}, T_minus_t = {}, M = {}, n = {}".format(S_t, K, r, q, sigma, T_minus_t, M, n))
	print("S_ave_t = {}, passing_time = {}, simulations = {}, repetitions = {}".format(S_ave_t, passing_time, simulations, repetitions))

	days = np.ceil((T_minus_t) * 365)  # in days
	days = int(days)
	dT = 1.0 / 365.0

	expected_asian_vals = []
	for rep in range(repetitions):
		print("#{} iteration.".format(rep))
		sim_results = []

		for i in range(simulations):
			S_beg = S_t
			sim_S_vals = []

			S_prev, S_next = S_beg, 0
			remained_days = days
			while (remained_days > 0):
				mean_lnS = np.log(S_prev) + (r - q - ( ( sigma ** 2.0 ) / 2.0) ) * dT # mean of lnST
				sigma_lnS = sigma * np.sqrt(dT)								   	    # sigma of lnST
				
				#if remained_days == days:
				#	print("mean = {}, sigmd = {}".format(mean_lnS, sigma_lnS))

				S_next = exp( random.gauss(mean_lnS, sigma_lnS) )

				sim_S_vals.append(S_next)
				S_prev = S_next
				remained_days -= 1
			
			#print("len of S_u_list = {}".format(len(S_u_list)))
			total_days = np.ceil( (passing_time + T_minus_t) * 365 ) + 1
			sim_S_ave = ( (passing_time * 365 + 1) * S_ave_t + sum(sim_S_vals) ) / (total_days)   # days = 0, 1, 2, ..., T * 365
			#print(max_S_u)
			payoff = max(sim_S_ave - K, 0)
			payoff = payoff * exp(- r * (T_minus_t) )
			sim_results.append(payoff)

		expected_asian = np.mean(sim_results)
		expected_asian_vals.append(expected_asian)

	sim_mean = np.mean(expected_asian_vals)
	sim_std  = np.std(expected_asian_vals)

	print("\nMonte Carlo results.")
	print("mean: {}, std: {}".format(sim_mean, sim_std))
	print("0.95 C.I. of a asian option: [{}, {}]".format(sim_mean - 2 * sim_std, sim_mean + 2 * sim_std))

			

