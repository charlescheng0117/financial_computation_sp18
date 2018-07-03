#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <algorithm>
#include <random>
#include <stdlib.h>
#include <time.h>

using namespace std;

void display(vector<double>& v) {
	int size = v.size();
	for (int i = 0; i < size; i += 1) {
		printf("%f,  ", v[i]);
	}
	printf("\n");
}

double mean_vector(vector<double>& v) {
	double ret = 0.0;
	int size = v.size();
	for (int i = 0; i < size; i += 1) {
		ret += v[i];
	}
	return ret / size;
}

double std_vector(vector<double>& v) {
	double ret = 0.0;
	int size = v.size();
	double mean = mean_vector(v);

	for (int i = 0; i < size; i += 1) {
		ret += pow(v[i] - mean, 2.0);
	}
	return sqrt(ret / (size));
}

double max_vec(vector<double>& v) {
	double max_val = -10e07;
	int size = v.size();
	for (int i = 0; i < size; i += 1) {
		double cur = v[i];
		if (cur > max_val)
			max_val = cur;
	}
	return max_val;
}

int main(int argc, char **argv) {
	double S_t, K, r, q, sigma, T_minus_t, S_ave_t, passing_time;
    int M, n, n_sim, n_rep;
	double u, d, p, dT;
    double passing_period;

    scanf("%lf %lf %lf %lf %lf %lf %d %d %lf %lf %d %d", &S_t, &K, &r, &q, &sigma, &T_minus_t, &M, &n, &S_ave_t, &passing_time, &n_sim, &n_rep);

	printf("Basic requirement: Monte Carlo simulation\n");
    printf("S_t          = %f\n", S_t);
    printf("K            = %f\n", K);
    printf("r            = %f\n", r);
    printf("q            = %f\n", q);
    printf("sigma        = %f\n", sigma);
    printf("T_minus_t    = %f\n", T_minus_t);
    printf("M            = %d\n", M);
    printf("n            = %d\n", n);
    printf("S_ave_t      = %f\n", S_ave_t);
    printf("passing_time = %f\n", passing_time);
    printf("n_sim        = %d\n", n_sim);
    printf("n_rep        = %d\n", n_rep);

	dT = (T_minus_t) / (double) n;
    passing_period = passing_time / (dT);
	u = exp(sigma * sqrt(dT));
	d = 1/u;
	p = (exp((r-q) * dT) - d)/(u - d);

	printf("u: %f, d: %f, p: %f, dT: %f\n", u, d, p, dT);

	default_random_engine generator;
	generator.seed(time(NULL)); // reset seed

	vector<double> expected_option_vals;
	for (int rep = 0; rep < n_rep; rep++) {
		//printf("#%d iteration.\n", rep);
		vector<double> sim_results;

		for (int i = 0; i < n_sim; i++) {
			vector<double> S_t_list;
            //S_t_list.push_back(S_t);
            
            double S_prev = S_t, S_next;
            
			int remained_days = n; // my way
			while (remained_days > 0) {
				//printf("dt = %f\n", dT);
				double mean_lnS_t = log(S_prev) + (r - q - ( pow(sigma, 2.0) / 2.0) ) * dT; // mean of lnST
				double sigma_lnS_t = sigma * sqrt(dT);								   	    // sigma of lnST

				normal_distribution<double> normal(mean_lnS_t, sigma_lnS_t);
				S_next = exp( normal(generator) );

				S_t_list.push_back(S_next);
				S_prev = S_next;
				--remained_days;
			}

            double sum = S_ave_t * passing_period + 
                         accumulate(S_t_list.begin(), S_t_list.end(), 0.0);

			//double S_ave = sum / (passing_period + n + 1); // S_t + (remained_days) of S_next => n + 1 values
			//double S_ave = sum / (passing_period + n + 1); // S_t + (remained_days) of S_next => n + 1 values
			double S_ave = sum / (passing_period + n ); // S_t + (remained_days) of S_next => n + 1 values

			double payoff = max(S_ave - K, 0.0);
			payoff = payoff * exp(-r * (T_minus_t) );
			sim_results.push_back(payoff);
		}

		//printf("len of sim_results = %d\n", sim_results.size());
		double expected_option = mean_vector(sim_results);
		expected_option_vals.push_back(expected_option);
	}

	double sim_mean = mean_vector(expected_option_vals);
	double sim_std  = std_vector(expected_option_vals);
   
    printf("-----------------------------------------\n");
    printf("### Answer for Monte Carlo Simulation ###\n");
	printf("mean: %f, std: %f\n", sim_mean, sim_std);
	printf("0.95 C.I. [%f, %f]\n", sim_mean - 2 * sim_std, sim_mean + 2 * sim_std);
    printf("-----------------------------------------\n");
	return 0;
}
