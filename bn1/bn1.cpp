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

#define BS 0
#define BINOM 1
#define CALL 0
#define PUT 1
#define EU 0
#define AM 1

using namespace std;
vector<double> ln_fac;

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

//double callBlackScholes(double S0, double K, double r, double q, double sigma, double T) 
//double putBlackScholes(double S0, double K, double r, double q, double sigma, double T) 
double callBlackScholes(double sigma);
double putBlackScholes(double sigma);

double calcCall(double ST, double K) {
	return max(ST - K, (double) 0.0);
}

double calcPut(double ST, double K) {
	return max(K - ST, (double) 0.0);
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
	//if (n < i)
	//	return -1;  // bug
	//return lnFactorial(n) - (lnFactorial(i) + lnFactorial(n - i));
    return ln_fac[n] - ln_fac[i] - ln_fac[n - i];
}

vector<double> CRRBinomial(double S0, double K, double r, double q, double sigma, double T, int n) { // n: # of periods
    #ifdef DEBUG
    printf("Number of periods: %d\n", n);
	printf("S0 = %lf, K = %lf, r = %lf, q = %lf, sigma = %lf, T = %lf\n", S0, K, r, q, sigma, T);
    #endif
	
    double u, d, p, dT;
	double t = 0.0; // always starts at t = 0

	dT = T/n;
	u = exp(sigma * sqrt(dT));
	d = 1/u;
	p = (exp((r-q) * dT) - d)/(u - d);
	printf("p: %lf, u: %lf, d: %lf, dT: %lf\n", p, u, d, dT);

	double ret_Call = 0.0, retPut = 0.0;
	vector<double> results;

	double tempSi[n + 1]; // compute Si
	double tempCi_eu[n + 1][n + 1]; // compute european Ci
	double tempPi_eu[n + 1][n + 1]; // compute european Pi
	double tempCi_am[n + 1][n + 1]; // compute american Ci
	double tempPi_am[n + 1][n + 1]; // compute american Pi

	//cout << "Ci is: ";
	for (int i = 0; i <= n; i += 1) { // leaf's values
		double Si = S0 * pow(u, n - i) * pow(d, i); // Si at t = i if S0 goes up n - i times.
		tempSi[i] = Si;
		//cout << Si << " ";
		double Ci = calcCall(Si, K), Pi = calcPut(Si, K);

		tempCi_eu[n][i] = Ci;
		tempPi_eu[n][i] = Pi;

		tempCi_am[n][i] = Ci;
		tempPi_am[n][i] = Pi;
	}
	//cout << '\n';

	for (int i = n - 1; i >= 0; i -= 1) { // from n - 1, n - 2, ..., 0
		for (int j = 0; j <= i; j += 1) { // c(i, j) = e^(-r dt) * (p * c(i + 1, j) + (1 - p) * c(i + 1, j + 1))
			//double Si = p * temp[j] + (1 - p) * temp[j + 1];
			//tempSi[j] = Si * exp(-r * T/n);
			double Sij = S0 * pow(u, i - j) * pow(d, j);
			double Cu_eu = tempCi_eu[i + 1][j], Cd_eu = tempCi_eu[i + 1][j + 1];
			double Pu_eu = tempPi_eu[i + 1][j], Pd_eu = tempPi_eu[i + 1][j + 1];		

			tempCi_eu[i][j] = exp(-r * dT) * (p * Cu_eu + (1 - p) * Cd_eu); 
			tempPi_eu[i][j] = exp(-r * dT) * (p * Pu_eu + (1 - p) * Pd_eu);

			double Cu_am = tempCi_am[i + 1][j], Cd_am = tempCi_am[i + 1][j + 1];
			double Pu_am = tempPi_am[i + 1][j], Pd_am = tempPi_am[i + 1][j + 1];	

			double rv_Ci_am = exp(-r * dT) * (p * Cu_am + (1 - p) * Cd_am);  // retension value of Call
			double rv_Pi_am = exp(-r * dT) * (p * Pu_am + (1 - p) * Pd_am);  // retension value of Put

			//printf("Exercise value of call: %f\n", calcCall(Sij, K));
			//printf("Current value of call:  %f\n", tempCi_eu[j]);
			tempCi_am[i][j] = max(rv_Ci_am, calcCall(Sij, K));
			//printf("Value of American call:  %f\n", tempCi_am[j]);
			tempPi_am[i][j] = max(rv_Pi_am, calcPut(Sij, K));
		}
	}
	double C0_eu = tempCi_eu[0][0];
	double P0_eu = tempPi_eu[0][0];
	double C0_am = tempCi_am[0][0];
	double P0_am = tempPi_am[0][0];

	results.push_back(C0_eu), results.push_back(P0_eu);
	results.push_back(C0_am), results.push_back(P0_am);

	return results;
}

double combinatorial_binom_call(double sigma);
double combinatorial_binom_put(double sigma);
double bisection(int model, int option_t, int eu_na);
double newton(int model, int option_t, int eu_na);

void print_line() {
    printf("--------------------------------------------------------------------------------\n");
}

double S0, K, r, q, T, option_price;
int n;
char criterion[100];

int main(int argc, char const *argv[])
{
	//test();

    /*
	double S0 = 30, K = 29, r = 0.05, q = 0.025, sigma = 0.3, T = 1;
	int simulations = 10000, repetitions = 20, n = 10;
	*/

    fscanf(stdin, "%lf %lf %lf %lf %lf %lf %d %s", &S0, &K, &r, &q, &T, &option_price, &n, criterion);
    
    printf("Getting input... \n");
    printf("S0           = %lf\n", S0);
    printf("K            = %lf\n", K);
    printf("r            = %lf\n", r);
    printf("q            = %lf\n", q);
    printf("T            = %lf\n", T);
    printf("option_price = %lf\n", option_price);
    printf("n            = %d\n", n);
    printf("convergence  = %s\n", criterion);


    bisection(BS, CALL, EU);

    // preprocess ln_fac
    ln_fac = vector<double>(n + 1);

    ln_fac[0] = ln_fac[1] = 0.0;

    for (int i = 2; i <= n; ++i) {
        ln_fac[i] = ln_fac[i - 1] + log(i);
        //printf("%f ", ln_fac[i]);
    }

    /*
    print_line();
    // Basic requirement
    printf("Basic requirement:\n");
    printf("Black Scholes formulas: \n");
    double call_price = callBlackScholes(S0, K, r, q, sigma, T);
    double put_price  = putBlackScholes(S0, K, r, q, sigma, T);
    printf("Call price = %f\nPut price  = %f\n", call_price, put_price);

    print_line();
    printf("Monte Carlo simulation:\n");
    monteCarlo(S0, K, r, q, sigma, T, simulations, repetitions); 

    vector<double> optionPrices;
    if (n <= 500) {
        print_line();
        printf("CRR binomial tree model:\n");
        optionPrices = CRRBinomial(S0, K, r, q, sigma, T, n);
        //call_euro, put_euro, call_amer, put_amer = CRR_binomial(S0, K, r, q, sigma, T, n)

        printf("Price for European call and put:\nCall = %f, Put = %f\n", optionPrices[0], optionPrices[1]);
        printf("Price for American call and put:\nCall = %f, Put = %f\n", optionPrices[2], optionPrices[3]);
    }

    if (n <= 500) {
        print_line();
        printf("Bonus 1: CRR binomial tree with one column vector.\n");
        optionPrices = oneColumnCRRBinomial(S0, K, r, q, sigma, T, n);

        printf("Price for European call and put:\nCall = %f, Put = %f\n", optionPrices[0], optionPrices[1]);
        printf("Price for American call and put:\nCall = %f, Put = %f\n", optionPrices[2], optionPrices[3]);
    }

    print_line();
    printf("Bonus 2: combinatorial method to price European options\n");
    vector<double> bonusOptions = combinatorialPricing(S0, K, r, q, sigma, T, n); 
	printf("Call = %f\nPut = %f\n", bonusOptions[0], bonusOptions[1]);
	

	//testComb();
    */

	return 0;
}

double callBlackScholes(double sigma) {
	double d1 = calcD1(S0, K, r, q, sigma, T);
	double d2 = calcD2(S0, K, r, q, sigma, T);
	return S0 * exp(-q * T) * N(d1) - K * exp(-r * T) * N(d2);
}

double putBlackScholes(double sigma) {
	double d1 = calcD1(S0, K, r, q, sigma, T);
	double d2 = calcD2(S0, K, r, q, sigma, T);
	return K * exp(-r * T) * N(-d2) - S0 * exp(-q * T) * N(-d1);
}

double combinatorial_binom(double sigma, int option_t) { // for european
	double u, d, p;
	double dT;
	double t = 0.0; // always starts at t = 0

	dT = T/n;
	u = exp(sigma * sqrt(dT));
	d = 1 / u;
	p = (exp((r-q) * dT) - d)/(u - d);

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

    if (option_t == CALL) {
        return C0;
    } else {
        return P0;
    }

}

double combinatorial_binom_call(double sigma) { // for european
	double u, d, p;
	double dT;
	double t = 0.0; // always starts at t = 0

	dT = T/n;
	u = exp(sigma * sqrt(dT));
	d = 1 / u;
	p = (exp((r-q) * dT) - d)/(u - d);

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

    return C0;
}

double combinatorial_binom_put(double sigma) { // for european
	double u, d, p;
	double dT;
	double t = 0.0; // always starts at t = 0

	dT = T/n;
	u = exp(sigma * sqrt(dT));
	d = 1 / u;
	p = (exp((r-q) * dT) - d)/(u - d);

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

    return P0;
}

double bisection_helper( double (*f) (double), double option_price) {
    double ret;

    double a0 = -0.1126, b0 = 1.126;
    double an, bn, xn;
    double fa, fb, fxn;

    an = a0;
    bn = b0;
    xn = a0 + (b0 - a0) / 2.0;
    
    int cnt = 0;
    while ( ! (abs( (*f)(xn) - option_price ) < 0.001 ) ) {
        printf("option_price: %lf\n", option_price);
        printf("%d, xn: %lf, fxn: %lf \n", cnt, xn, (*f)(xn));
        printf("%d, an: %lf, fan: %lf \n", cnt, an, (*f)(an));
        printf("%d, bn: %lf, fbn: %lf \n", cnt, bn, (*f)(bn));
        

        fa = (*f)(an) - option_price;
        fxn = (*f)(xn) - option_price;

        if ( (fa * fxn) < 0 ) {
            bn = xn;
        } else {
            an = xn;
        }
        ++cnt;    
        xn = an + (bn - an) / 2.0;
    
    }
    return xn;
}

double bisection(int model, int option_t, int eu_na) {
    // model: 0, BS, 1, Binom
    // option_t: 0, call, 1, put
    // eu_na: 0, eu, 1, na

    auto f_ptr;

    if (model == BS) {
        if (option_t == CALL) {
            model_ptr = callBlackScholes;

        } else { //PUT
            model_ptr = putBlackScholes;
        }
    } else { // Binomial tree
        if (option_t == CALL) {
            if (eu_na == EU) {
                mo
            }
        }
    
    }

   

    // BS, eu, call
    double ans = bisection_helper(callBlackScholes, option_price);
    cout << callBlackScholes(0.3) << "\n";
    cout << ans << "\n";
    cout << option_price << "\n";


    return 0.0;
}
