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

typedef vector< vector<double>* > Matrix;
typedef vector< vector< vector<double>* >* > Cube;

double roundDouble(double val, int d) {
	return round(val * pow(10.0, d) ) / pow(10.0, d);
}

double S_i_j(double S_0, double i, double j, double u, double d) {
	double res = S_0 * pow(u, i - j) * pow(d, j);
	return roundDouble(res, 4);
}

double A_max(double S_0, double i, double j, double u, double d) {
	double res = (S_0 * (1 - pow(u, i - j + 1)) / (1 - u) + S_0 * pow(u, i - j) * d * (1 - pow(d, j)) / (1 - d) ) / (i + 1.0);
	return roundDouble(res, 4);
}

double A_min(double S_0, int i, int j, double u, double d) {
	double res = (S_0 * (1 - pow(d, j + 1)) / (1 - d) + S_0 * pow(d, j) * u * (1 - pow(u, i - j)) / (1 - u) ) / (i + 1.0);
	return roundDouble(res, 4);
}	

double A_i_j_k(double S_0, double i, double j, double k, double M, double u, double d) {
	double res = (M - k) / M * A_max(S_0, i, j, u, d) + k / M * A_min(S_0, i, j, u, d);
	return roundDouble(res, 4);
}

double A_u(Cube* ptr_Aijk, double S_0, int i, int j, int k, double u, double d) {
	double Aijk = (* (* (* ptr_Aijk)[i])[j])[k];
	double res = ( (i + 1) * Aijk + S_0 * pow(u, i + 1 - j) * pow(d, j) ) / (i + 2);
	return roundDouble(res, 4);
}

double A_d(Cube* ptr_Aijk, double S_0, int i, int j, int k, double u, double d) {
	double Aijk = (* (* (* ptr_Aijk)[i])[j])[k];
	double res = ( (i + 1) * Aijk + S_0 * pow(u, i + 1 - (j + 1)) * pow(d, j + 1) ) / (i + 2);
	return roundDouble(res, 4);
}

double w_u(Cube* ptr_Aijk, double Au, int i, int j, int k) {
	double nominator = (* (* (* ptr_Aijk)[i + 1])[j])[k - 1] - Au;
	double denominator = (* (* (* ptr_Aijk)[i + 1])[j])[k - 1] - (* (* (* ptr_Aijk)[i + 1])[j])[k];
	double res = nominator / denominator;
	if (isnan(res)) {
		res = 0.0;
	}
	return roundDouble(res, 4);
}

double w_d(Cube* ptr_Aijk, double Ad, int i, int j, int k) {
	double nominator = (* (* (* ptr_Aijk)[i + 1])[j + 1])[k - 1] - Ad;
	double denominator = (* (* (* ptr_Aijk)[i + 1])[j + 1])[k - 1] - (* (* (* ptr_Aijk)[i + 1])[j + 1])[k];
	double res = nominator / denominator;
	if (isnan(res)) {
		res = 0.0;
	}
	return roundDouble(res, 4);
}

double Cnjk(Cube* ptr_Aijk, int n, int j, int k, double K) {  // Call value for terminal node
	double Aijk = (* (* (* ptr_Aijk)[n])[j])[k];
	return max(Aijk - K, 0.0);
}

Matrix* create_Sij_mat(double S_0, double u, double d, int n) {
	Matrix* ret = new Matrix(n + 1);

	for (int i = 0; i <= n; i += 1) {
		vector<double>* vec = new vector<double>(n + 1);
		(*ret)[i] = vec;
	}

	for (int i = 0; i <= n; i += 1) {
		for (int j = 0; j <= i; j += 1) {
			(*(*ret)[i])[j] = S_i_j(S_0, (double) i, (double) j, u, d);
		}
	}
	return ret;
}

void delete_Sij_mat(Matrix* Sij_mat) {
	int n = (*Sij_mat)[0]->size() - 1;
	for (int i = 0; i <= n; i += 1) {
		delete (*Sij_mat)[i];
	}
	return;
}

Cube* create_Aijk_cube(double S_0, int M, double u, double d, int n) {
	Cube* ret = new Cube(n + 1);

	for (int i = 0; i <= n; i += 1) {
		vector< vector<double>* >* tmp = new (vector< vector<double>* >)(n + 1);
		(*ret)[i] = tmp;

		for (int j = 0; j <= i; j += 1) {
			vector<double>* vec_k = new vector<double>(M + 2); // k = -1, 0, 1, ..., M.   M + 2 values. 
			(* (*ret)[i] ) [j] = vec_k;
			//printf("A_max(%d, %d): %f ", (int) i, (int) j, A_max(S_0, (double) i, (double) j, u, d));
			//printf("A_min(%d, %d): %f\n", (int) i, (int) j, A_min(S_0, (double) i, (double) j, u, d));
			for (int k = -1; k <= M; k += 1) {
				(*(*(*ret)[i])[j])[k] = A_i_j_k(S_0, i, j, k, M, u, d);
			}
		}
		//cout << "\n";
	}
	return ret;
}

void display(Matrix* Sij_mat) {
	int n = (*Sij_mat)[0]->size() - 1;
	for (int i = 0; i <= n; i += 1) {
		for (int j = 0; j <= i; j += 1) {
			printf("S_%d_%d: %f ", i, j, (*(*Sij_mat)[i])[j]);
		}
		cout << '\n';
	}
}

void display(Cube* Aijk_cube) {
	printf("Displaying a cube.\n");

	int n = (*Aijk_cube)[0]->size() - 1;
	int M = (*(*Aijk_cube)[0])[0]->size() - 2;
	for (int i = 0; i <= n; i += 1) {
		for (int j = 0; j <= i; j += 1) {
			printf("\n");
			for (int k = -1; k <= M; k += 1) {
				printf("A_%d_%d_%d: %f\n", i, j, k, (*(*(*Aijk_cube)[i])[j])[k]);
			}
		}
		//cout << "\n";
	}
	return;
}

int test() {
	double S_0 = 150, K = 140, r = 0.06, q = 0, sigma = 0.6, n = 3, T = 1;
	double u, d, p, dT;

	int M = 5;

	dT = (T - 0.0) / n;
	u = exp(sigma * sqrt(dT));
	d = 1/u;
	p = (exp((r-q) * dT) - d)/(u - d);

	dT = roundDouble(dT, 4);
	u = roundDouble(u, 4);
	d = roundDouble(d, 4);
	p = roundDouble(p, 4);

	printf("u: %f, d: %f, p: %f\n", u, d, p);

	double S_1_0 = roundDouble(S_0 * u, 4); 
	double S_1_1 = roundDouble(S_0 * d, 4); 
	printf("S_1_0 %f, S_1_1 %f\n", S_1_0, S_1_1);

	double S_2_0 = roundDouble(S_1_0 * u, 4); 
	double S_2_1 = roundDouble(S_1_0 * d, 4); 
	double S_2_2 = roundDouble(S_1_1 * d, 4); 
	printf("S_2_0 %f, S_2_1 %f, S_2_2 %f\n", S_2_0, S_2_1, S_2_2);

	double S_3_0 = roundDouble(S_2_0 * u, 4); 
	double S_3_1 = roundDouble(S_2_0 * d, 4); 
	double S_3_2 = roundDouble(S_2_1 * d, 4); 
	double S_3_3 = roundDouble(S_2_2 * d, 4); 
	printf("S_3_0 %f, S_3_1 %f, S_3_2 %f, S_3_3 %f\n", S_3_0, S_3_1, S_3_2, S_3_3);

	double S_4_0 = roundDouble(S_3_0 * u, 4); 
	double S_4_1 = roundDouble(S_3_0 * d, 4); 
	double S_4_2 = roundDouble(S_3_1 * d, 4); 
	double S_4_3 = roundDouble(S_3_2 * d, 4); 
	double S_4_4 = roundDouble(S_3_3 * d, 4); 
	printf("S_4_0 %f, S_4_1 %f, S_4_2 %f, S_4_3 %f, S_4_4 %f\n", S_4_0, S_4_1, S_4_2, S_4_3, S_4_4);

	for (int i = 0; i <= n; i += 1 ) {
		for (int j = 0; j <= i; j += 1 ) {
			printf("S_%d_%d: %f ", i, j, S_i_j(S_0, (double) i, (double) j, u, d));
		}
		cout << '\n';
	}

	for (int i = 0; i <= n; i += 1) {
		for (int j = 0; j <= i; j += 1) {
			printf("A_max(%d, %d): %f ", (int) i, (int) j, A_max(S_0, (double) i, (double) j, u, d));
			printf("A_min(%d, %d): %f\n", (int) i, (int) j, A_min(S_0, (double) i, (double) j, u, d));
			for (int k = 0; k <= M; k += 1) {
				printf("A_i_j_k(%d, %d, %d): %f\n", i, j, k, A_i_j_k(S_0, i, j, k, M, u, d));
			}
		}
		cout << "\n";
	}
	return 0;
}

int main(int argc, char const *argv[]) {
	
	double S_0 = 150, K = 140, r = 0.06, q = 0, sigma = 0.6, n = 5, T = 1;
	double u, d, p, dT;

	int M = 15;

	dT = (T - 0.0) / n;
	u = exp(sigma * sqrt(dT));
	d = 1/u;
	p = (exp((r-q) * dT) - d)/(u - d);

	dT = roundDouble(dT, 4);
	u = roundDouble(u, 4);    // 1.414
	d = roundDouble(d, 4);    // 0.7072
	p = roundDouble(p, 4);    // 0.4428

	printf("u: %f, d: %f, p: %f\n", u, d, p);

	Matrix* Sij_mat = create_Sij_mat(S_0, u, d, n);
	display(Sij_mat);

	Cube* Aijk_cube = create_Aijk_cube(S_0, M, u, d, n);
	display(Aijk_cube);

	// create cube Cijk
	Cube* Cijk_cube = new Cube(n + 1);
	for (int i = 0; i <= n; i += 1) {
		vector< vector<double>* >* tmp = new (vector< vector<double>* >)(n + 1);
		(*Cijk_cube)[i] = tmp;
		for (int j = 0; j <= i; j += 1) {
			vector<double>* vec_k = new vector<double>(M + 2); // k = -1, 0, 1, ..., M.   M + 2 values. 
			(* (* Cijk_cube)[i] ) [j] = vec_k;
		}
		//cout << "\n";
	}

	// calculate C(n, j, k) for terminal node.
	for (int j = 0; j <= n; j += 1) {
		for (int k = -1; k <= M; k += 1) {
			(* (* (* Cijk_cube)[n])[j])[k] = Cnjk(Aijk_cube, n, j, k, K);
		}
	}
	//display(Cijk_cube);
	for (int i = n - 1; i >= 0; i -= 1) {
		for (int j = 0; j <= i; j += 1) {
			for (int k = 0; k <= M; k += 1) {
				double Au = A_u(Aijk_cube, S_0, i, j, k, u, d);
				double Ad = A_d(Aijk_cube, S_0, i, j, k, u, d);
				double wu = w_u(Aijk_cube, Au, i, j, k);
				double wd = w_d(Aijk_cube, Ad, i, j, k);

				double Cu_up = (* (* (* Cijk_cube)[i + 1])[j])[k];
				double Cu_down = (* (* (* Cijk_cube)[i + 1])[j])[k - 1];
				double Cd_up = (* (* (* Cijk_cube)[i + 1])[j + 1])[k];
				double Cd_down = (* (* (* Cijk_cube)[i + 1])[j + 1])[k - 1];

				double Cu = wu * Cu_up + (1 - wu) * Cu_down;
				double Cd = wd * Cd_up + (1 - wd) * Cd_down;
				
				printf("i = %d, j = %d, k = %d\n", i, j, k);
				printf("Au = %f, Ad = %f, wu = %f, wd = %f\n", Au, Ad, wu, wd);
				printf("Cu_up = %f, Cu_down = %f, Cd_up = %f, Cd_down = %f\n", Cu_up, Cu_down, Cd_up, Cd_down);
				printf("Cu = %f, Cd = %f\n", Cu, Cd);

				(* (* (* Cijk_cube)[i])[j])[k] = exp(-r * dT) * (p * Cu + (1 - p) * Cd);
				printf("Value of Cijk = %f\n", exp(-r * dT) * (p * Cu + (1 - p) * Cd));
			}
		}
	}

	printf("\n\n\nResults start here. \n\n\n");

	display(Cijk_cube);
	delete_Sij_mat(Sij_mat);
	//delete_Aijk_cube(Aijk_cube);

	return 0;
}