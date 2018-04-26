#include <iostream>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <algorithm>
#include <random>
#include <stdlib.h>
#include <time.h>

#define e exp(1)
#define PI M_PI

using namespace std;

void displayVector(vector<double> v) {
	for (int i = 0; i < v.size(); i += 1) {
		cout << v[i] << " ";
	}
	cout << "\n";
}

double calcD1(double S0, double K, double r, double q, double sigma, double T) {
	double d1 = (log(S0/K) + (r - q + pow(sigma, 2.0) / 2.0) * T) / (sigma * sqrt(T));
	return d1;
}

double calcD2(double S0, double K, double r, double q, double sigma, double T) {
	double d2 = (log(S0/K) + (r - q - pow(sigma, 2.0) / 2.0) * T) / (sigma * sqrt(T));
	return d2;
}

double pdfNormal(double x, double mu, double sigma) { // Normal dist. pdf
	return (1 / (sigma * sqrt(2 * PI))) * exp(-pow((x - mu), 2) / (2 * pow(sigma, 2)));
}

double N(double x) { // Normal dist. cdf
	double gamma = 0.2316419;
	double a1 = 0.319381530;
	double a2 = -0.356563782;
	double a3 = 1.781477937;
	double a4 = -1.821255978;
	double a5 = 1.330274429;
	double pdf = pdfNormal(x, 0, 1); // standard normal

	double k, tmp;
	if (x >= 0) {
		k = 1 / (1 + gamma * x);
	} else {
		k = 1 / (1 + gamma * -x);
	}

	tmp = a1 * k + a2 * pow(k, 2) + a3 * pow(k, 3) + a4 * pow(k, 4) + a5 * pow(k, 5);

	if (x >= 0) {
		return 1 - pdf * tmp;
	} else { // x < 0
		return 1 - (1 - pdf * tmp);
	}
} 

double callBlackScholes(double S0, double K, double r, double q, double sigma, double T) {
	double d1 = calcD1(S0, K, r, q, sigma, T);
	double d2 = calcD2(S0, K, r, q, sigma, T);
	return S0 * exp(-q * T) * N(d1) - K * exp(-r * T) * N(d2);
}

double putBlackScholes(double S0, double K, double r, double q, double sigma, double T) {
	double d1 = calcD1(S0, K, r, q, sigma, T);
	double d2 = calcD2(S0, K, r, q, sigma, T);
	return K * exp(-r * T) * N(-d2) - S0 * exp(-q * T) * N(-d1);
}

int test() {
	// pdf 
	double x = pdfNormal(2.5, 12, 5);
	cout << "pdf is: " << x << '\n';

	// cdf 
	x = -1.96;
	printf("x = %f, cdf is: %f \n", x, N(x));

	// d1
	double S0 = 102, K = 100, r = 0.05;
	double q = 0, sigma = 0.15, T = 465.0/365.0;
	double d1 = calcD1(S0, K, r, q, sigma, T);
	printf("S0: %f, K: %f, r: %f, q: %f, sigma: %f, T: %f, d1: %f \n", S0, K, r, q, sigma, T, d1);
	double d2 = calcD2(S0, K, r, q, sigma, T);
	printf("S0: %f, K: %f, r: %f, q: %f, sigma: %f, T: %f, d2: %f \n", S0, K, r, q, sigma, T, d2);
	printf("call price: %f, put price: %f\n", callBlackScholes(S0, K, r, q, sigma, T), putBlackScholes(S0, K, r, q, sigma, T));
	return 0;
}

double calcCall(double ST, double K) {
	return max(ST - K, 0.0);
}

double calcPut(double ST, double K) {
	return max(K - ST, 0.0);
}

double mean_vector(vector<double> v) {
	double ret = 0.0;
	int size = v.size();
	for (int i = 0; i < size; i += 1) {
		ret += v[i];
	}
	return ret / size;
}

double std_vector(vector<double> v, double mean) {
	double ret = 0.0;
	int size = v.size();

	for (int i = 0; i < size; i += 1) {
		ret += pow(v[i] - mean, 2);
	}
	return sqrt(ret / (size));
}

double monteCarlo(double S0, double K, double r, double q, 
				  double sigma, double T, int simulations, int repetitions) {
	printf("\nA MonteCarlo Simution.\n");
	printf("Number of simulations: %d, Number of repetitions: %d\n", simulations, repetitions);
	printf("S0 = %f, K = %f, r = %f, q = %f, sigma = %f, T = %f\n", S0, K, r, q, sigma, T);

	double mean_lnST = log(S0) + (r - q - pow(sigma, 2.0) / 2.0) * T; // mean of lnST
	double sigma_lnST = sigma * sqrt(T); // sigma of lnST

	printf("mean_lnST = %f, sigma_lnST = %f\n", mean_lnST, sigma_lnST);

	default_random_engine generator;
	normal_distribution<double> normal(mean_lnST, sigma_lnST);

	vector<double> callRes;
	vector<double> putRes;

	generator.seed(time(NULL)); // reset seed
	for (int i = 0; i < repetitions; i += 1) {
		vector<double> callVal;
		vector<double> putVal;

		//cout << "ST is: ";
		for (int j = 0; j < simulations; j += 1) {
			//cout << "simulations is:" << simulations << '\n';
			double sample = normal(generator);
			//cout << exp(sample) << " ";
			callVal.push_back(calcCall(exp(sample), K) * exp(-r * T));
			putVal.push_back(calcPut(exp(sample), K) * exp(-r * T));
		}
		//cout << "\n";

		double call_expected_val = mean_vector(callVal);
		double put_expected_val = mean_vector(putVal);

		callRes.push_back(call_expected_val);
		putRes.push_back(put_expected_val);
	}

	double call_val_mean = mean_vector(callRes);
	double call_val_std = std_vector(callRes, call_val_mean);

	double put_val_mean = mean_vector(putRes);
	double put_val_std = std_vector(putRes, put_val_mean);
	//double put_val_std = 

	printf("call: mean = %f, std = %f\n", call_val_mean, call_val_std);
	printf("put: mean = %f, std = %f\n", put_val_mean, put_val_std);
	printf("0.95 C.I. of Call value for MonteCarlo simulations: [%f, %f]\n", call_val_mean - 2 * call_val_std, call_val_mean + 2 * call_val_std);
	printf("0.95 C.I. of Put value for MonteCarlo simulations: [%f, %f]\n\n", put_val_mean - 2 * put_val_std, put_val_mean + 2 * put_val_std);

	return 0;
}

double factorial(int n) { // return n!
	if (n == 1 || n == 0) 
		return 1;
	return exp( log(n) + log( factorial(n - 1) ) );
}

double combination(int n , int i) { // nCi
	if (n < i)
		return -1;  // bug
	return factorial(n) / (factorial(i) * factorial(n - i));
}

double lnFactorial(int n) { // return n!
	if (n == 1 || n == 0)
		return 0;
	return log(n) + lnFactorial(n - 1) ;
}

double lnComb(int n , int i) { // nCi
	if (n < i)
		return -1;  // bug
	return lnFactorial(n) - (lnFactorial(i) + lnFactorial(n - i));
}

vector<double> CRRBinomial(double S0, double K, double r, double q, double sigma, double T, int n) { // n: # of periods
	printf("\nCRRBinomial pricing.\n");
	printf("Number of periods: %d\n", n);
	printf("S0 = %f, K = %f, r = %f, q = %f, sigma = %f, T = %f\n", S0, K, r, q, sigma, T);

	long double u, d, p, dT;
	double t = 0.0; // always starts at t = 0

	dT = T/n;
	u = exp(sigma * sqrt(dT));
	d = 1/u;
	p = (exp((r-q) * dT) - d)/(u - d);
	printf("p: %Lf, u: %Lf, d: %Lf, dT: %Lf\n", p, u, d, dT);

	double ret_Call = 0.0, retPut = 0.0;
	vector<double> results;

	double tempSi[n + 1]; // compute Si
	double tempCi_eu[n + 1]; // compute european Ci
	double tempPi_eu[n + 1]; // compute european Pi
	double tempCi_am[n + 1]; // compute american Ci
	double tempPi_am[n + 1]; // compute american Pi

	//cout << "Ci is: ";
	for (int i = 0; i <= n; i += 1) { // leaf's values
		double Si = S0 * pow(u, n - i) * pow(d, i); // Si at t = i if S0 goes up n - i times.
		tempSi[i] = Si;
		//cout << Si << " ";
		double Ci = calcCall(Si, K), Pi = calcPut(Si, K);

		tempCi_eu[i] = Ci;
		tempPi_eu[i] = Pi;

		tempCi_am[i] = Ci;
		tempPi_am[i] = Pi;
	}
	//cout << '\n';

	for (int i = n - 1; i >= 0; i -= 1) { // from n - 1, n - 2, ..., 0
		for (int j = 0; j <= i; j += 1) { // c(i, j) = e^(-r dt) * (p * c(i + 1, j) + (1 - p) * c(i + 1, j + 1))
			//double Si = p * temp[j] + (1 - p) * temp[j + 1];
			//tempSi[j] = Si * exp(-r * T/n);
			double Sij = S0 * pow(u, i - j) * pow(d, j);
			double Cu_eu = tempCi_eu[j], Cd_eu = tempCi_eu[j + 1];
			double Pu_eu = tempPi_eu[j], Pd_eu = tempPi_eu[j + 1];		

			tempCi_eu[j] = exp(-r * dT) * (p * Cu_eu + (1 - p) * Cd_eu); 
			tempPi_eu[j] = exp(-r * dT) * (p * Pu_eu + (1 - p) * Pd_eu);

			double Cu_am = tempCi_am[j], Cd_am = tempCi_am[j + 1];
			double Pu_am = tempPi_am[j], Pd_am = tempPi_am[j + 1];	

			double rv_Ci_am = exp(-r * dT) * (p * Cu_am + (1 - p) * Cd_am);  // retension value of Call
			double rv_Pi_am = exp(-r * dT) * (p * Pu_am + (1 - p) * Pd_am);  // retension value of Put

			//printf("Exercise value of call: %f\n", calcCall(Sij, K));
			//printf("Current value of call:  %f\n", tempCi_eu[j]);
			tempCi_am[j] = max(rv_Ci_am, calcCall(Sij, K));
			//printf("Value of American call:  %f\n", tempCi_am[j]);
			tempPi_am[j] = max(rv_Pi_am, calcPut(Sij, K));
		}
	}
	double C0_eu = tempCi_eu[0];
	double P0_eu = tempPi_eu[0];
	double C0_am = tempCi_am[0];
	double P0_am = tempPi_am[0];

	results.push_back(C0_eu), results.push_back(P0_eu);
	results.push_back(C0_am), results.push_back(P0_am);

	return results;
}

vector<double> bonusCRRBinomial(double S0, double K, double r, double q, double sigma, double T, int n) {
	printf("\nBonus CRRBinomial pricing.\n");
	printf("Number of periods: %d\n", n);
	printf("S0 = %f, K = %f, r = %f, q = %f, sigma = %f, T = %f\n", S0, K, r, q, sigma, T);

	long double u, d, p;
	double dT;
	double t = 0.0; // always starts at t = 0

	dT = T/n;
	u = exp(sigma * sqrt(dT));
	d = 1 / u;
	p = (exp((r-q) * dT) - d)/(u - d);
	//printf("p: %Lf, u: %Lf, d: %Lf, dT: %f\n", p, u, d, dT);

	double retCall = 0.0, retPut = 0.0;
	vector<double> results;

	for (int i = 0; i <= n; i += 1) {
		double lnSi = log(S0) + (n - i) * log(u) + i * log(d); // Si at t = i if S0 goes up i times.
		double Si = exp(lnSi);
		double callXi = calcCall(Si, K);  // payoff of call at t = i 
		double putXi = calcPut(Si, K);    // payoff of put at t = i 
		//printf("i is: %d, Si is: %f, callXi is: %f, putXi is: %f\n", i, Si, callXi, putXi);

		if (callXi != 0)
			//printf("lnComb is: %f\n", lnComb(n, i));
			retCall += exp( lnComb(n, i) + (n - i) * log(p) + i * log(1 - p) + log(callXi) ); // S0 * u^(n - i) * d^(i)
			//printf("retCall: %f\n", retCall);
		if (putXi != 0)
			retPut  += exp( lnComb(n, i) + (n - i) * log(p) + i * log(1 - p) + log(putXi)  ); 
		//printf("retCall now is: %f, retPut now is: %f\n", retCall, retPut);
	}

	double C0 = retCall * exp(-r * T);
	double P0 = retPut * exp(-r * T);

	results.push_back(C0), results.push_back(P0);

	return results;

}

double recursiveCRRBinomial(double S0, double K, double r, double q, double p, double u, double d, double i) {
	if (i == 0) // leaf
		return S0;
	return 0.0;
}

void testComb() {
	printf("Testing bonus comb & fact.\n");
	printf("%f\n", lnComb(20, 10));
	printf("%f %f %f\n", lnFactorial(200), lnFactorial(0), lnComb(200, 0));
}

int main(int argc, char const *argv[])
{
	//test();


	double S0 = 30, K = 29, r = 0.05, q = 0.025, sigma = 0.3, T = 1;
	int simulations = 10000, repetitions = 20, n = 10;
	
	double callPrice = callBlackScholes(S0, K, r, q, sigma, T);
	double putPrice = putBlackScholes(S0, K, r, q, sigma, T);
	printf("Black-Scholes:\nCall = %f, Put = %f\n", callPrice, putPrice);

	vector<double> optionPrices = CRRBinomial(S0, K, r, q, sigma, T, n);
	printf("European Call = %f, Put = %f\n", optionPrices[0], optionPrices[1]);
	printf("American Call = %f, Put = %f\n", optionPrices[2], optionPrices[3]);

	monteCarlo(S0, K, r, q, sigma, T, simulations, repetitions);

	vector<double> bonusOptions = bonusCRRBinomial(S0, K, r, q, sigma, T, n);
	printf("Call = %f, Put = %f\n", bonusOptions[0], bonusOptions[1]);
	
	//printf("\n!!!Should add American Option Price!!!\n");

	//testComb();

	return 0;
}
