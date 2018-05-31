import numpy as np
import numpy.linalg as la
import math
import sys

def sample_mat_Z(simulations, n):
	mat_Z = []
	mu, sigma = 0, 1
	for i in range(simulations):
		row_i = np.random.normal(mu, sigma, n)
		mat_Z.append(row_i)
	return np.array(mat_Z)

def calc_mat_Z_tilde(mat_Z):
	mat_Z_T = np.transpose(mat_Z)
	mat_Z_tilde_T = []

	for i in range(len(mat_Z_T)):
		zi = mat_Z_T[i]				# mat_Z_T[i] is the row i of mat_Z_T, which is exactly column i of mat_Z
		ui = np.mean(zi)			# mean of column i of mat_Z
		z_tilde = zi - ui
		'''
		print("zi: {}".format(zi))
		print("ui: {}".format(ui))
		print("z^~: {}".format(z_tilde))
		'''
		mat_Z_tilde_T.append(z_tilde)
	mat_Z_tilde_T = np.array(mat_Z_tilde_T)
	#print("mat_Z_tilde's transpose:\n{}".format(mat_Z_tilde_T))
	return np.transpose(mat_Z_tilde_T)

def calc_mat_C(mat_Z):
	mat_Z_T = np.transpose(mat_Z)
	mat_Z_tilde_T = []

	for i in range(len(mat_Z_T)):
		zi = mat_Z_T[i]				# mat_Z_T[i] is the row i of mat_Z_T, which is exactly column i of mat_Z
		ui = np.mean(zi)			# mean of column i of mat_Z
		z_tilde = zi - ui
		'''
		print("zi: {}".format(zi))
		print("ui: {}".format(ui))
		print("z^~: {}".format(z_tilde))
		'''
		mat_Z_tilde_T.append(z_tilde)
	mat_Z_tilde_T = np.array(mat_Z_tilde_T)
	#print("mat_Z_tilde's transpose:\n{}".format(mat_Z_tilde_T))

	#mat_Z_tilde = np.transpose(mat_Z_tilde)
	#print("mat_Z_tilde:\n{}".format(mat_Z_tilde))
	mat_C = np.cov(mat_Z_tilde_T)

	return mat_C

def sumSquareOfColPrevEle(mat, i):
	# return the sum of square of ith columns' element previous to entry i
	ret = 0.0
	for k in range(i):
		ret += pow(mat[k][i], 2.0)
	return ret

def sumProdOfTwoColPrevEle(mat, i, j): 
	# return the sum of square of ith columns' element previous to entry i
	ret = 0.0
	for k in range(i):
		ret += (mat[k][i] * mat[k][j])
	return ret

def CholeskyDecomp(cov_mat):
	n = len(cov_mat)
	mat_A = [ [0 for i in range(n)] for j in range(n)]
	c11 = cov_mat[0][0]
	mat_A[0][0] = np.sqrt(c11); # a11 = sqrt(c11)

	a11 = mat_A[0][0];
	for j in range(1, n):    # a1j = c1j / a11, j = 2, ..., n
		mat_A[0][j] = cov_mat[0][j] / a11
	
	# loop of step 2 and step 3
	for i in range(1, n - 1): #for (int i = 1; i < n - 1; i += 1) { 

		# step 2: aii = sqrt(cii - \sum_{k = 1}^{i - 1} a_{k,i}^2)
		# printf("sum_{k = 1}^{i - 1} a_{k, %d}^2 = %f\n", i, sumSquareOfColPrevEle(mat_A, i) );
		mat_A[i][i] = np.sqrt(cov_mat[i][i] - sumSquareOfColPrevEle(mat_A, i) )

		# step 3: aij = 1/aii (cij - \sum_{k = 1}^{i - 1} a_{k, i} * a_{k, j}), j = i + 1, ..., n
		for j in range(i + 1, n):
			#printf("sum_{k = 1}^{i - 1} a_{k, %d} * a_{k, %d} = %f\n", i, j, sumProdOfTwoColPrevEle(mat_A, i, j));
			mat_A[i][j] = 1 / mat_A[i][i] * ( cov_mat[i][j] - sumProdOfTwoColPrevEle(mat_A, i, j) )
	
	# step 4
	mat_A[n - 1][n - 1] = np.sqrt(cov_mat[n - 1][n - 1] - sumSquareOfColPrevEle(mat_A, n - 1) )
	return np.array(mat_A)

def calcPayoffRainbow(Si_values, K):
	# payoff = max(max(Si_values) - K, 0)
	max_Si = -10e7
	max_index = -1
	n = len(Si_values)
	ret = {}

	#print(Si_values)

	for i in range(n):
		if Si_values[i] > max_Si:
			max_Si = Si_values[i]
			max_index = i
			#cout << "current Si is: " << Si_values[i] << "\n";
			#cout << "max_Si is: " << max_Si << "\n";
	#print("Si_values: {}\n max_Si: {}".format(Si_values, max_Si))	

	ret['payoff'] = max(max_Si - K, 0.0)
	ret[i] = max_index
	#print("payoff = {}".format(ret['payoff']))
	#cout<< "Max value: " << ret.payoff << endl;
	#cout<< "Corresponding i: " << ret.i << endl;
	return ret


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 hw3_bonus2.py <input_file>")
        sys.exit()
    file_name = sys.argv[1]
    inputs = list(open(file_name))

    line1 = [float(x) for x in inputs[0].split()]
    K, r, T, simulations, repetitions, n = line1
    simulations, repetitions, n = int(simulations), int(repetitions), int(n)

    line2 = [float(x) for x in inputs[1].split()]  # S0_vals
    S0_vals = np.array(line2)

    line3 = [float(x) for x in inputs[2].split()]  # qi_vals
    qi_vals = np.array(line3)

    line4 = [float(x) for x in inputs[3].split()]  # sigma_vals
    sigma_vals = np.array(line4)

    rho_mat = []
    for i in range(4, len(inputs)):				   # rho_matrix
            line = [float(x) for x in inputs[i].split()]
            rho_mat.append(line)
    rho_mat = np.array(rho_mat)

    cov_mat = [ [ 0.0 for i in range(n) ] for j in range(n) ]
    for i in range(n):
            for j in range(n):
                    cov_mat[i][j] = rho_mat[i][j] * sigma_vals[i] * sigma_vals[j]
    cov_mat = np.array(cov_mat)

    mat_A = CholeskyDecomp(cov_mat)
    
    DEBUG = False
    print("----------     Bonus 2     ----------")
    print("Monte Carlo Simulation - Bonus 2")
    print("Number of Simulations: {}\nNumber of Repetitions: {}\nn: {}".format(simulations, repetitions, n))
    print("K: {}\nr: {}\nT: {}".format(K, r, T))
    if DEBUG:
        print("S0_vals = {}".format(S0_vals))
        print("qi_vals = {}".format(qi_vals))
        print("sigma_vals = {}".format(sigma_vals))
        print("rho_mat = \n{}".format(rho_mat))
        print("cov_mat = \n{}".format(cov_mat))
        print("mat_A = \n{}".format(mat_A))

    rainbow_results = []

    for rep in range(repetitions):
            mat_Z = sample_mat_Z(simulations, n)     # step 1 of new method
            #print("mat_Z = \n{}".format(mat_Z))
            #print("mat_Z^T = \n{}".format(np.transpose(mat_Z)))

            mat_Z_tilde = calc_mat_Z_tilde(mat_Z)    # step 2 of new method
            #print("mat_Z_tilde = \n{}".format(mat_Z_tilde))
            mat_C = np.cov(mat_Z_tilde.T)
            #print("mat_C = \n{}".format(mat_C))

            mat_A_tilde = CholeskyDecomp(mat_C)		 # step 3 of new method
            #print("mat_A_tilde = \n{}".format(mat_A_tilde))
            inv_mat_A_tilde = la.inv(mat_A_tilde)

            mat_Z_prime = np.dot(mat_Z_tilde, inv_mat_A_tilde)
            #print("mat_Z_prime = \n{}".format(mat_Z_prime))

            rainbow_vals = []
            for i in range(simulations):
                    r_vals = np.dot(mat_Z_prime[i], mat_A)
                    '''
                    print("Z'_i = {}".format(mat_Z_prime[i]))
                    print("mat_A = \n{}".format(mat_A))
                    print("r_vals = {}".format(r_vals))
                    '''
                    mean_lnST_vals, sigma_lnST_vals = [], []

                    for j in range(n):

                            S0, sigma, q = S0_vals[j], sigma_vals[j], qi_vals[j]
                            mean_lnST = np.log(S0) + (r - q - ( ( sigma ** 2.0 )/ 2.0) ) * T # mean of lnST
                            sigma_lnST = sigma * np.sqrt(T)								     # sigma of lnST

                            mean_lnST_vals.append(mean_lnST)
                            sigma_lnST_vals.append(sigma_lnST)

                    # get S_1T, S_2T, ..., S_nT
                    SiT_vals = []

                    for j in range(n):
                            lnST = r_vals[j] * np.sqrt(T) + mean_lnST_vals[j] # // ri ~ N(0, sigma^2), but
                                                                                                                              # // lnST ~ N(log(S0) + (r - q - pow(sigma, 2.0) / 2.0) * T, sigma^2 T)
                            SiT_vals.append(np.e ** lnST)
                    #print("#sim: {}, SiT_vals: {}".format(i, SiT_vals))
                            
                    rainbow_payoff = calcPayoffRainbow(SiT_vals, K)
                    payoff = rainbow_payoff['payoff']
                    #qi = qi_vals[rainbow_result.i]
                    
                    #double pv_rainbow_val = payoff * exp(-(r - qi) * T);
                    pv_rainbow_val = payoff * (np.e ** (-r * T))
                    rainbow_vals.append(pv_rainbow_val)

            expected_rainbow_val = np.mean(rainbow_vals)
            #print("rainbow_vals = {}, expected_rainbow_val = {}".format(rainbow_vals, expected_rainbow_val))
            rainbow_results.append(expected_rainbow_val)

    #print("rainbow_results = {}".format(rainbow_results))
    rainbow_val_mean = np.mean(rainbow_results)
    rainbow_val_std = np.std(rainbow_results)
    
    print("------------------------------")
    print("### Answer for Bonus 2 ###")
    print("mean: {:.6f}, std: {:.6f}".format(rainbow_val_mean, rainbow_val_std))
    print("0.95 C.I.: [{:.6f}, {:.6f}]".format(rainbow_val_mean - 2 * rainbow_val_std, rainbow_val_mean + 2 * rainbow_val_std))
