#include <iostream>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <algorithm>
#include <random>
#include <stdlib.h>
#include <time.h>
#include <fstream>

#define e exp(1)
#define PI M_PI

using namespace std;

struct rainbowPair {
	double payoff; // max payoff
	int i;	   // corresponding index i
};

void display(vector<double>& v) {
	int size = v.size();
	for (int i = 0; i < size; i += 1) {
		printf("%f  ", v[i]);
	}
	printf("\n");
}

void display(vector< vector<double> >& mat) {
	int m = mat.size(), n = mat[0].size();
	printf("m is %d, n is %d\n", m, n);

	for (int i = 0; i < m; i += 1) {
		for (int j = 0; j < n; j += 1) {
			printf("%f\t", mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

double dot(vector<double>& v1, vector<double>& v2) { // return the dot product of a vector
	if (v1.size() != v2.size()) {
		printf("Incorrect dimension. v1: %lu, v2: %lu\n", v1.size(), v2.size());
		return 1;
	}
	double result = 0;
	for (int i = 0; i < v1.size(); i += 1) {
		result += (v1[i] * v2[i]);
	}
	return result;
}

vector<double> getColumn(vector< vector<double> > mat, int k) { // return the kth column vector of a matrix
	int m = mat.size();  // mat: m x n
	vector<double> ret;

	for (int i = 0; i < m; i += 1) {
		ret.push_back(mat[i][k]);
	}
	return ret;
}

vector<double> dot(vector<double>& vec, vector< vector<double> >& mat) { // multiplication of a vector vs a matrix
	int m = 1, n_m1 = vec.size(); // dimension of m1: 1 x n
	int n_m2 = mat.size(), k = mat[0].size(); // dimension of m2: n x k

	if (n_m1 != n_m2) { // m1: m x n,  m2: n x k.   n = m1.size(),  n = m2[0].size()
		printf("Incorrect dimension. m1: %d x %d, m2: %d x %d\n", m, n_m1, n_m2, k);
		return ( vector<double> ) 0;
	}
	vector<double> ret(k);     // Dimension of ret is: 1 x k.

	for (int i = 0; i < k; i += 1) {
		vector<double> ci = getColumn(mat, i);
		ret[i] = dot(vec, ci);
	}
	return ret;
}

vector< vector<double> > dot(vector< vector<double> >& m1, vector< vector<double> >& m2) { // multiplication of a square matrix
	int m = m1.size(), n_m1 = m1[0].size(); // dimension of m1: m x n
	int n_m2 = m2.size(), k = m2[0].size(); // dimension of m2: n x k

	if (n_m1 != n_m2) { // m1: m x n,  m2: n x k.   n = m1.size(),  n = m2[0].size()
		printf("Incorrect dimension. m1: %d x %d, m2: %d x %d\n", m, n_m1, n_m2, k);
		return (vector< vector<double> >) 0;
	}
	vector< vector<double> > ret(m);     // Dimension of ret is: m x k.
	for (int i = 0; i < k; i += 1) {
		vector<double> row(k);
		ret[i] = row;
	}

	for (int i = 0; i < m; i += 1) {
		for (int j = 0; j < k; j += 1) {
			vector<double> ri = m1[i];
			vector<double> cj = getColumn(m2, j);
			ret[i][j] = dot(ri, cj);
		}
	}
	return ret;
}

vector<double> initVec(double* values, int N) {
	//printf("sizeof values is: %lu, sizeof double is: %lu \n", sizeof(values), sizeof(double));
	vector<double> v(values, values + N);
	return v;
}

//template <int N>
vector< vector<double> > initMat(int M, int N, double** values) { // initialize a matrix: M x N
	vector< vector<double> > ret(M);

	for (int i = 0; i < M; i += 1) {
		ret.push_back(initVec(values[i], N));
	}
	return ret;
}


double sumSquareOfColPrevEle(vector< vector<double> >& mat, int i) { 
	// return the sum of square of ith columns' element previous to entry i
	double ret = 0.0;
	for (int k = 0; k < i; k += 1) {
		ret += pow(mat[k][i], 2.0);
	}
	return ret;
}

double sumProdOfTwoColPrevEle(vector< vector<double> >& mat, int i, int j) { 
	// return the sum of square of ith columns' element previous to entry i
	double ret = 0.0;
	for (int k = 0; k < i; k += 1) {
		ret += (mat[k][i] * mat[k][j]);
	}
	return ret;
}

void CholeskyDecomposition(vector< vector<double> >& cov_matrix, vector< vector<double> >& mat_A) {
	// step 1
	int n = cov_matrix.size();
	double c11 = cov_matrix[0][0];
	mat_A[0][0] = sqrt(c11); // a11 = sqrt(c11)

	double a11 = mat_A[0][0];
	for (int j = 1; j < n; j += 1) { // a1j = c1j / a11, j = 2, ..., n
		mat_A[0][j] = cov_matrix[0][j] / a11;
	}

	// loop of step 2 and step 3
	for (int i = 1; i < n - 1; i += 1) { 

		// step 2: aii = sqrt(cii - \sum_{k = 1}^{i - 1} a_{k,i}^2)
        #ifdef DEBUG
        printf("sum_{k = 1}^{i - 1} a_{k, %d}^2 = %f\n", i, sumSquareOfColPrevEle(mat_A, i) );
        #endif
		mat_A[i][i] = sqrt(cov_matrix[i][i] - sumSquareOfColPrevEle(mat_A, i) );

		// step 3: aij = 1/aii (cij - \sum_{k = 1}^{i - 1} a_{k, i} * a_{k, j}), j = i + 1, ..., n
		for (int j = i + 1; j < n; j += 1) {
            #ifdef DEBUG
			printf("sum_{k = 1}^{i - 1} a_{k, %d} * a_{k, %d} = %f\n", i, j, sumProdOfTwoColPrevEle(mat_A, i, j));
            #endif
			mat_A[i][j] = 1 / mat_A[i][i] * ( cov_matrix[i][j] - sumProdOfTwoColPrevEle(mat_A, i, j) );
		}
	}

	// step 4
	mat_A[n - 1][n - 1] = sqrt(cov_matrix[n - 1][n - 1] - sumSquareOfColPrevEle(mat_A, n - 1));
	return;
}

rainbowPair calcPayoffRainbow(vector<double>& Si_values, double K) {
	// payoff = max(max(Si_values) - K, 0)
	double max_Si = -10e7;
	int max_index = -1;
	int n = Si_values.size();
	rainbowPair ret;

	for (int i = 0; i < n; i += 1) {
		if (Si_values[i] > max_Si) {
			max_Si = Si_values[i];
			max_index = i;
			//cout << "current Si is: " << Si_values[i] << "\n";
			//cout << "max_Si is: " << max_Si << "\n";
		}
	}
	ret.payoff = max(max_Si - K, 0.0), ret.i = max_index;
	//cout<< "Max value: " << ret.payoff << endl;
	//cout<< "Corresponding i: " << ret.i << endl;
	return ret;
}

double mean_vector(vector<double>& v) {
	double ret = 0.0;
	int size = v.size();
	for (int i = 0; i < size; i += 1) {
		ret += v[i];
	}
	return ret / size;
}

double std_vector(vector<double>& v, double mean) {
	double ret = 0.0;
	int size = v.size();

	for (int i = 0; i < size; i += 1) {
		ret += pow(v[i] - mean, 2.0);
	}
	return sqrt(ret / (size));
}

void moment_match(vector<double>& z) {
	double mean = mean_vector(z);
	double std = std_vector(z, mean);
	int size = z.size();

	for (int i = 0; i < size; i += 1) {
		double tmp = (z[i] - mean) / std;
		z[i] = tmp;
	}
	return;
}

double monteCarloRainbow(double K, double r, double T, int simulations, int repetitions, int n,
						 vector<double>& S0_vals, vector<double>& qi_vals, vector<double>& sigma_vals,
						 vector< vector<double> >& mat_A) {
	printf("Monte Carlo Simulation - Basic Requirement\n");
	printf("Number of simulations: %d\nNumber of repetitions: %d\nn: %d\n", simulations, repetitions, n);
	printf("K: %f\nr: %f\nT: %f\n", K, r, T);

    #ifdef DEBUG
	printf("S0_vals are \n");
	display(S0_vals);

	printf("qi_vals are \n");
	display(qi_vals);

	printf("sigma_vals are \n");
	display(sigma_vals);
    #endif

	// sampling z1, z2, ..., zn from iid standard normal.
	default_random_engine generator;
	normal_distribution<double> z_normal(0, 1);

	vector<double> rainbowResults;

	generator.seed(time(NULL)); // reset seed
	for (int rep = 0; rep < repetitions; rep += 1) {
		vector<double> rainbowVals;

		// Following is ``1 senario'' of a rainbow option.
		for (int k = 0; k < simulations; k += 1) {

			vector<double> z_vals;
			for (int i = 0; i < n; i += 1) {
				double zi = z_normal(generator);
				z_vals.push_back(zi);
			}

			// get r1, r2, ..., rn by dot(z_vals, A)
			vector<double> r_vals;
			r_vals = dot(z_vals, mat_A);

			vector<double> mean_lnST_vals, sigma_lnST_vals;
			vector< normal_distribution<double> > normals;
			for (int i = 0; i < n; i += 1) {
				double S0 = S0_vals[i], sigma = sigma_vals[i], q = qi_vals[i];
				double mean_lnST = log(S0) + (r - q - (pow(sigma, 2.0) / 2.0) ) * T; // mean of lnST
				double sigma_lnST = sigma * sqrt(T); // sigma of lnST

				mean_lnST_vals.push_back(mean_lnST);
				sigma_lnST_vals.push_back(sigma_lnST);

			}

			// get S_1T, S_2T, ..., S_nT
			vector<double> SiT_vals;
			for (int i = 0; i < n; i += 1) {
				double lnST = r_vals[i] * sqrt(T) + mean_lnST_vals[i]; // ri ~ N(0, sigma^2), but
																	   // lnST ~ N(log(S0) + (r - q - pow(sigma, 2.0) / 2.0) * T, sigma^2 T)
				SiT_vals.push_back(exp(lnST));
			}

			rainbowPair rainbow_result = calcPayoffRainbow(SiT_vals, K);
			double payoff = rainbow_result.payoff, qi = qi_vals[rainbow_result.i];
			//double pv_rainbow_val = payoff * exp(-(r - qi) * T);
			double pv_rainbow_val = payoff * exp(-r * T);

			rainbowVals.push_back(pv_rainbow_val);
		}
		double expected_rainbow_val = mean_vector(rainbowVals);
		rainbowResults.push_back(expected_rainbow_val);
	}

	double rainbow_val_mean = mean_vector(rainbowResults);
	double rainbow_val_std = std_vector(rainbowResults, rainbow_val_mean);
    
    printf("------------------------------\n");
	printf("### Answer for Basic Requirement ###\n");
	printf("mean: %f, std: %f\n", rainbow_val_mean, rainbow_val_std);
	printf("0.95 CI: [%f, %f]\n", rainbow_val_mean - 2 * rainbow_val_std, rainbow_val_mean + 2 * rainbow_val_std);

	return 0;
}

double bonus_monteCarloRainbow(double K, double r, double T, int simulations, int repetitions, int n,
						 vector<double>& S0_vals, vector<double>& qi_vals, vector<double>& sigma_vals,
						 vector< vector<double> >& mat_A) {

    // Bonus 1
	printf("Monte Carlo Simulation - Bonus 1\n");
	printf("Number of simulations: %d\nNumber of repetitions: %d\nn: %d\n", simulations, repetitions, n);
	printf("K: %f\nr: %f\nT: %f\n", K, r, T);
    
    #ifdef DEBUG    
    printf("S0_vals are \n");
	display(S0_vals);

	printf("qi_vals are \n");
	display(qi_vals);

	printf("sigma_vals are \n");
	display(sigma_vals);
    #endif

	// sampling z1, z2, ..., zn from iid standard normal.
	default_random_engine generator;
	normal_distribution<double> z_normal(0, 1);

	vector<double> rainbowResults;

	generator.seed(time(NULL)); // reset seed
	for (int rep = 0; rep < repetitions; rep += 1) {
		vector<double> rainbowVals;
		/* Sampling <n * simulations> zi's, then use the method of moment match and antithetic_variate. */
		// moment matching
		int total_sample = n * simulations;
		vector<double> vec_z_vals;
		for (int i = 0; i < total_sample; i += 1) {
			double zi = z_normal(generator);
			vec_z_vals.push_back(zi);
			//printf("%d. zi is: %f\n", i, zi);
		}	
		//printf("size of vec_z_vals: %d\n", vec_z_vals.size());
		//printf("Before\n");
		//display(vec_z_vals);
		moment_match(vec_z_vals);
		//printf("After\n");
		//display(vec_z_vals);
		//printf("### Mean is %f, std is %f\n", mean_vector(vec_z_vals), std_vector(vec_z_vals, mean_vector(vec_z_vals)));

		// antithetic_variate
		for (int i = 0; i < total_sample; i += 1) {
			double tmp = -vec_z_vals[i];
			vec_z_vals.push_back(tmp);
		}
		//printf("Hey size of vec_z_vals: %d\n", vec_z_vals.size());
		//printf("std = %f\n", std_vector(vec_z_vals, mean_vector(vec_z_vals)));

		vector< vector<double> > all_z_vals;
		vector<double> tmp_z_vals;
		total_sample = total_sample * 2;
		for (int i = 0; i < total_sample; i += 1) {
			double tmp_z = vec_z_vals[i]; 
			tmp_z_vals.push_back(tmp_z);
			if (tmp_z_vals.size() == n) {
				all_z_vals.push_back(tmp_z_vals);
				vector<double> new_vec;
				tmp_z_vals = new_vec;
			}
		}

		//printf("size of all_z_vals: %d\n", all_z_vals.size());

		// Following is ``1 senario'' of a rainbow option.
		for (int k = 0; k < simulations * 2; k += 1) {
			vector<double> z_vals = all_z_vals[k];
			//display(z_vals);

			// get r1, r2, ..., rn by dot(z_vals, A)
			vector<double> r_vals;
			r_vals = dot(z_vals, mat_A);

			vector<double> mean_lnST_vals, sigma_lnST_vals;
			vector< normal_distribution<double> > normals;
			for (int i = 0; i < n; i += 1) {
				double S0 = S0_vals[i], sigma = sigma_vals[i], q = qi_vals[i];
				double mean_lnST = log(S0) + (r - q - (pow(sigma, 2.0) / 2.0) ) * T; // mean of lnST
				double sigma_lnST = sigma * sqrt(T); // sigma of lnST

				mean_lnST_vals.push_back(mean_lnST);
				sigma_lnST_vals.push_back(sigma_lnST);

			}

			// get S_1T, S_2T, ..., S_nT
			vector<double> SiT_vals;
			for (int i = 0; i < n; i += 1) {
				double lnST = r_vals[i] * sqrt(T) + mean_lnST_vals[i]; // ri ~ N(0, sigma^2), but
																	   // lnST ~ N(log(S0) + (r - q - pow(sigma, 2.0) / 2.0) * T, sigma^2 T)
				SiT_vals.push_back(exp(lnST));
			}

			rainbowPair rainbow_result = calcPayoffRainbow(SiT_vals, K);
			double payoff = rainbow_result.payoff, qi = qi_vals[rainbow_result.i];
			//double pv_rainbow_val = payoff * exp(-(r - qi) * T);
			double pv_rainbow_val = payoff * exp(-r * T);

			rainbowVals.push_back(pv_rainbow_val);
		}
		double expected_rainbow_val = mean_vector(rainbowVals);
		rainbowResults.push_back(expected_rainbow_val);
	}
	
    #ifdef DEBUG
	printf("rainbowResults: \n");
	display(rainbowResults);
    #endif

	double rainbow_val_mean = mean_vector(rainbowResults);
	double rainbow_val_std = std_vector(rainbowResults, rainbow_val_mean);

    printf("------------------------------\n");
	printf("### Answer for Bonus 1 ###\n");
	printf("mean: %f, std: %f\n", rainbow_val_mean, rainbow_val_std);
	printf("0.95 CI: [%f, %f]\n", rainbow_val_mean - 2 * rainbow_val_std, rainbow_val_mean + 2 * rainbow_val_std);

	return 0;
}

int main(int argc, char const *argv[]) {
	
	if (argc != 2) {
		cout<< "usage: " << argv[0] << " <filename>\n";
		return 0;
	}
	ifstream in(argv[1]); 
	ofstream out(argv[2]);

	double K, r, T;
	int simulations, repetitions, n;
	vector<double> Si_values;
	vector<double> qi_values;
	vector<double> sigma_values;

	in >> K >> r >> T >> simulations >> repetitions >> n;

	vector< vector<double> > rho_matrix(n);
	vector< vector<double> > cov_matrix(n);
	vector< vector<double> > mat_A(n);

	for (int i = 0; i < n; i += 1) {
		vector<double> v1(n), v2(n), v3(n);

		rho_matrix[i] = v1;
		cov_matrix[i] = v2;
		mat_A[i] = v3;
	}

	//printf("K: %2f\nr: %2f\nT: %2f\nNumber of simulations: %2d\nNumer of repetitions: %2d\nn: %2d\n", K, r, T, simulations, repetitions, n);
	double si;
	for (int i = 0; i < n; i += 1) {
		in >> si;
		Si_values.push_back(si);
	}
    #ifdef DEBUG
	display(Si_values);
    #endif
    
    #ifdef DEBUG
	rainbowPair pairTest = calcPayoffRainbow(Si_values, K);
	printf("Payoff of the rainbow option: %f\n", pairTest.payoff);
    #endif

	double qi;
	for (int i = 0; i < n; i += 1) {
		in >> qi;
		qi_values.push_back(qi);
	}

	double sigma_i;
	for (int i = 0; i < n; i += 1) {
		in >> sigma_i;
		//printf("%f ", sigma_i);
		sigma_values.push_back(sigma_i);
	}

    #ifdef DEBUG
	display(sigma_values);
    #endif

	//printf("sigma_i is: %f ", sigma_i);

	double rho_ij;
	for (int i = 0; i < n; i += 1) {
		for (int j = 0; j < n; j += 1) {
			in >> rho_ij;
			//printf(" %f ", rho_ij);
			rho_matrix[i][j] = rho_ij;
			cov_matrix[i][j] = rho_ij * sigma_values[i] * sigma_values[j];
		}
	}

	CholeskyDecomposition(cov_matrix, mat_A);
    
    printf("----------Basic Requirement----------\n");
	monteCarloRainbow(K, r, T, simulations, repetitions, n, Si_values, qi_values, sigma_values, mat_A);
    
    printf("----------     Bonus 1     ----------\n");
	bonus_monteCarloRainbow(K, r, T, simulations, repetitions, n, Si_values, qi_values, sigma_values, mat_A);
	
    #ifdef DEBUG	
	cout << "rho_matrix\n";
	display(rho_matrix);
	cout << "cov_matrix\n"; 
	display(cov_matrix);
	cout << "matrix_A\n";
	display(mat_A);
    #endif
	
	
	return 0;
}
