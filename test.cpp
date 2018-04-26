#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

typedef vector< vector<double> > Matrix;
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

bool isBinary(int x) {
	return (x & (x + 1));
}

vector<double> initVec(double* values, int N) {
	vector<double> ret(N);
	for (int i = 0; i < N; i += 1) {
		//ret.push_back(values[i]);
		ret[i] = values[i];
	}
	return ret;
}

int fib(int n) {
	if (n == 1 || n == 0) {
		cout << " " << n << " \n"; 
		return n;
	} else {
		int fib_n = fib(n - 1) + fib(n - 2);
		cout << " " << fib_n << " ";
		return fib_n;
	}
}

int testDouble(double i) {
	static double val[] = {3, 4, 5};
	return val[(int) i];
}

int main(int argc, char const *argv[]) {
	double val[] = {1, 2, 3};
	vector<double> vec = initVec(val, 3);

	display(vec);

	int i = 1;
	cout << testDouble(i) << "\n";

	double u = 1.414;

	cout << pow(u, 2.0) << "\n";
	cout << (u == 'nan') << "\n";
	double m = 2 ** 3;
	cout << abs(5 - 6) << "\n";
	cout << m << "\n";

	return 0;
}