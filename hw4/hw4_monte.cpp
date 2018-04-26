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
	// load data
	if (argc != 2) {
		cout << "Command arguments not correct.\n";
		//return -1;
	}
	char* in_file = argv[1];
	ifstream in(in_file);

	double St, r, q, sigma, t, T, S_max_t;
	int n, simulations, repetitions;

	/*
	in >> St >> r >> q >> sigma >> t >> T >> S_max_t;
	in >> n >> simulations >> repetitions;
	*/
	
	St = 50.0, r = 0.1, q = 0.05, sigma = 0.4, t = 1.0, T = 4.0, S_max_t = 60.0;
	n = 20, simulations = 1000, repetitions = 20;

	double S_beg = St;
	int days = (int) ceil((T - t) * 365);
	double dT = 1.0 / 365.0; // in terms of days

	printf("Monte Carlo simulation to price European lookback puts.\n");
	printf("St = %f, r = %f, q = %f, sigma = %f, t = %f, T = %f, S_max_t = %f\n", St, r, q, sigma, t, T, S_max_t);
	printf("n  = %d, simulations = %d, repetitions = %d\n", n, simulations, repetitions);

	printf("S_beg = %f, days = %d, dT = %f\n", S_beg, days, dT);

	default_random_engine generator;
	generator.seed(time(NULL)); // reset seed

	vector<double> expected_put_vals;
	for (int rep = 0; rep < repetitions; rep += 1) {
		printf("#%d iteration.\n", rep);
		vector<double> sim_results;

		for (int i = 0; i < simulations; i += 1) {
			vector<double> S_u_list;
			S_u_list.push_back(S_max_t);

			double S_prev, S_next;
			int remained_days = 0;
			S_prev = S_beg, S_next = 0.0;

			remained_days = days;
			//printf("runs = %d\n", runs);
			while (remained_days > 0) {
				//printf("dt = %f\n", dT);
				double mean_lnS_u = log(S_prev) + (r - q - ( pow(sigma, 2.0) / 2.0) ) * dT; // mean of lnST
				double sigma_lnS_u = sigma * sqrt(dT);								   	    // sigma of lnST

				//printf("mean_lnS_u = %f, sigma_lnS_u = %f\n", mean_lnS_u, sigma_lnS_u);

				normal_distribution<double> normal(mean_lnS_u, sigma_lnS_u);
				S_next = exp( normal(generator) );

				S_u_list.push_back(S_next);
				S_prev = S_next;
				remained_days -= 1;
			}

			//printf("Su list =\n");
			
			
			//printf("len of S_u_list = %d\n", S_u_list.size());

			double max_S_u = max_vec(S_u_list);
			double S_last = S_next;  // or = *(S_u_list.end() - 1)

			double payoff = max(max_S_u - S_last, 0.0);
			//printf("payoff = %f\n", payoff);
			payoff = payoff * exp(- r * (T - t) );
			//printf("discounted payoff = %f\n", payoff);
			sim_results.push_back(payoff);
		}
		//printf("len of sim_results = %d\n", sim_results.size());
		double expected_put = mean_vector(sim_results);
		expected_put_vals.push_back(expected_put);
	}

	//printf("len of expected_put_vals = %d\n", expected_put_vals.size());
	//display(expected_put_vals);
	double sim_mean = mean_vector(expected_put_vals);
	double sim_std  = std_vector(expected_put_vals);

	printf("Monte Carlo results.\n");
	printf("mean: %f, std: %f\n", sim_mean, sim_std);
	printf("0.95 C.I. of a lookback put option: [%f, %f]\n", sim_mean - 2 * sim_std, sim_mean + 2 * sim_std);
	return 0;
}