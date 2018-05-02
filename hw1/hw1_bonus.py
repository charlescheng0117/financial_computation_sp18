import numpy as np

def calc_payoff(ST, K1, K2, K3, K4):
    if ST <= K1 or ST >= K4:
        return 0
    elif K2 <= ST and ST <= K3:
        return K2 - K1
    elif K1 < ST and ST < K2:
        return ST - K1
    else: # K3 < ST < K4
        return ( (K2 - K1)/(K4 - K3) ) * (K4 - ST)

def monteCarlo(S0, K1, K2, K3, K4, r, q, sigma, T, simulations, repetitions):
    print("A MonteCarlo Simution.")
    print("Number of simulations: {}, Number of repetitions: {}".format(simulations, repetitions))
    print("S0 = {}, K1 = {}, K2 = {}, K3 = {}, K4 = {}, r = {}, q = {}, sigma = {}, T = {}".format(S0, K1, K2, K3, K4, r, q, sigma, T))

    # mean, sigma of lnSt
    mean_lnST = np.log(S0) + (r - q - np.power(sigma, 2.0) / 2.0) * T 
    sigma_lnST = sigma * np.sqrt(T) #// sigma of lnST

    print("mean_lnST = {}, sigma_lnST = {}".format(mean_lnST, sigma_lnST))

    # draw samples from N(lnST, sigmaST)


    option_results = []
    # set random seed    
    #generator.seed(time(NULL)) // reset seed
    
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


    print("option: mean = {}, std = {}".format(option_mean, option_std))
    print("0.95 C.I. of Option value for MonteCarlo simulations: [{}, {}]".format(option_mean - 2 * option_std, option_mean + 2 * option_std))


if __name__ == "__main__":
    S0, K1, K2, K3, K4, r, q, sigma, T = 50, 40, 50, 60, 70, 0.05, 0.08, 0.2, 1
    """
    S0 = 100
    K1, K2, K3, K4 = 85, 95, 105, 115
    r = 0.5
    q = 0.1
    sigma = 0.3
    T = 1
    """
    simulations = 10000
    repetitions = 20

    monteCarlo(S0, K1, K2, K3, K4, r, q, sigma, T, simulations, repetitions)
