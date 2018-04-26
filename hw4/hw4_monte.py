import numpy as np
import random

def exp(x):
	return np.e ** x

if __name__ == '__main__':
	file_name = "hw4_input.txt"
	inputs = list(open(file_name))[0]
	inputs = [float(x) for x in inputs.split()]

	St, r, q, sigma, t, T, S_max_t, n, simulations, repetitions = inputs
	n, simulations, repetitions = int(n), int(simulations), int(repetitions)

	print("Monte Carlo simulation to price European lookback puts.")
	print("St = {}, r = {}, q = {}, sigma = {}, t = {}, T = {}, S_max_t = {}".format(St, r, q, sigma, t, T, S_max_t))
	print("n  = {}, simulations = {}, repetitions = {}".format(n, simulations, repetitions))

	S_beg = St
	days = np.ceil((T - t) * 365)
	#dT = 1/days
	dT = 1 / 365

	print("S_beg = {}, days = {}, dT = {}".format(S_beg, days, dT))
	#print("S_beg")

	expected_put_vals = []
	for rep in range(repetitions):
		print("#{} iteration.".format(rep))
		sim_results = []

		for i in range(simulations):
			S_u = S_beg
			S_u_list = [S_max_t]

			S_prev, S_next = S_u, 0
			remained_days = days
			while (remained_days > 0):
				mean_lnS_u = np.log(S_prev) + (r - q - ( ( sigma ** 2.0 ) / 2.0) ) * dT # mean of lnST
				sigma_lnS_u = sigma * np.sqrt(dT)								   	    # sigma of lnST

				#mean_lnS_u = np.log(S_prev) + (r - q - ( ( sigma ** 2.0 ) / 2.0) ) * (1 / 365) # mean of lnST
				#sigma_lnS_u = sigma * np.sqrt((1 / 365))								   	    # sigma of lnST

				S_next = exp( random.gauss(mean_lnS_u, sigma_lnS_u) )
				#S_next = exp( np.random.normal(mean_lnS_u, sigma_lnS_u) )
				#print("S_next = {}".format(S_next))

				S_u_list.append(S_next)
				S_prev = S_next
				remained_days -= 1
			
			#print("len of S_u_list = {}".format(len(S_u_list)))
			max_S_u = max(S_u_list)
			#print(max_S_u)
			payoff = max(max_S_u - S_beg, 0)
			payoff = payoff * exp(- r * (T - t) )
			sim_results.append(payoff)

		expected_put = np.mean(sim_results)
		expected_put_vals.append(expected_put)

	sim_mean = np.mean(expected_put_vals)
	sim_std  = np.std(expected_put_vals)

	print("Monte Carlo results.")
	print("mean: {}, std: {}".format(sim_mean, sim_std))
	print("0.95 C.I. of a lookback put option: [{}, {}]".format(sim_mean - 2 * sim_std, sim_mean + 2 * sim_std))

			

