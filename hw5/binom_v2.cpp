#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <algorithm>
#include <random>
#include <stdlib.h>
#include <time.h>
#include <functional>
#include <chrono>

#define EPSILON 0.00000001

using namespace std;

bool greater_double(double x1, double x2) {
    return x1 > x2;
}

void print_vec(vector<double>& vec) {
    int n = vec.size();
    for (int i = 0; i < n; ++i) {
        printf("%lf ", vec[i]);
    }
    cout << "\n";
}

pair<int, int> sequential_search(vector<double>& vec, double target) {
    // find range
    // vec: from largest to smallest
    int left = 0, right = 0;
    int n = vec.size();
    for (int i = 0; i < n - 1; ++i) {
        if (abs(target - vec[i]) < EPSILON ) {
            return pair<int, int>(i, i);
        }
        if ( vec[i] > target && target > vec[i + 1] ) {
            return pair<int, int>(i, i+1);
        }

    }
    /*
    if (abs(target - vec[n-1]) < 10e-4 ) {
        return pair<int, int>(n -1, n - 1);
    }*/
    return pair<int, int>(n -1, n - 1);
}

pair<int, int> binary_search(vector<double>& arr,int left, int right, double& key) {
    int mid = left + (right - left)/2 ;
    if ( abs(arr[mid] - key) < 0.0000001 ) {
        return pair<int, int>(mid, mid);
    } else if ( arr[mid] > key ) {
       return binary_search ( arr, mid + 1, right, key);
    } else if ( arr[mid] < key ) {
        if (abs ( arr[mid-1] - key ) < 0.0000001) {
            return pair<int, int>(mid-1, mid-1) ;
        } else if (arr[mid-1] > key) {
            return pair<int, int>(mid-1, mid) ;
        } else {
            return binary_search(arr, left, mid - 1, key );
        }
    }
}

pair<int, int> interpolation_search(vector<double>& arr, int left, int right, double key) {
    int mid ;
    if ( left == right ) {
        mid = left ;
    } else if ( ( arr[left] - arr[right] ) < 0.0000001 ) {
        return pair<int, int>(left, left) ;
    } else {
        mid = left + ( right - left ) * ( key - arr[left] ) / ( arr[right] - arr[left] ) ;
    }

    if(abs( arr[mid] - key ) < 0.0000001 ) {
        return pair<int, int>(mid, mid);
    } else if ( arr[mid] > key ) {
        return interpolation_search( arr, mid + 1, mid + 1, key);
    } else if ( arr[mid] < key ) {
        if ( abs ( arr[mid-1] - key ) < 0.0000001 ) {
            return pair<int, int>(mid-1, mid-1);
        } else if ( arr[mid-1] > key ) {
            return pair<int, int>(mid-1, mid);
        } else {
            return interpolation_search(arr, mid - 1, mid-1, key);
        }
    }
}

struct Node {
    int M;
    vector<double> A_vec;
    vector<double> A_log_vec;
    vector<double> C_vec;
    vector<double> C_log_vec;
    vector<double> C_am_vec;
    vector<double> C_am_log_vec;

    Node() { }

    Node(int new_M) {
        M = new_M;
        A_vec = vector<double>(M + 1, 0.0);
        A_log_vec = vector<double>(M + 1, 0.0);
        C_vec = vector<double>(M + 1, 0.0);
        C_log_vec = vector<double>(M + 1, 0.0);
        C_am_vec = vector<double>(M + 1, 0.0);
        C_am_log_vec = vector<double>(M + 1, 0.0);
    }
    
    void print() {
        printf("A    : ");
        for (int i = 0; i <= M; ++i) {
            printf("%lf ", A_vec[i]);
        }
        printf("\n");
        printf("C    : ");
        for (int i = 0; i <= M; ++i) {
            printf("%lf ", C_vec[i]);
        }
        printf("\n");
        printf("A_log: ");
        for (int i = 0; i <= M; ++i) {
            printf("%lf ", A_log_vec[i]);
        }
        printf("\n");
        printf("C_log: ");
        for (int i = 0; i <= M; ++i) {
            printf("%lf ", C_log_vec[i]);
        }
        printf("\n");
    }
};

double compute_A_max_ij(double S_0, double S_ave_t, double passing_period, double i, double j, double u, double d) {
	//double res = (S_0 * (1 - pow(u, i - j + 1) ) / (1 - u) + S_0 * pow(u, i - j) * d * (1 - pow(d, j)) / (1 - d) + S_ave_t * passing_period ) / (i + passing_period + 1.0);
	//double res = (S_0 * (1 - pow(u, i - j + 1) ) / (1 - u) + S_0 * pow(u, i - j) * d * (1 - pow(d, j)) / (1 - d) + S_ave_t * passing_period - S_0 ) / (i + passing_period);
	return (S_0 * (1 - pow(u, i - j + 1) ) / (1 - u) + S_0 * pow(u, i - j) * d * (1 - pow(d, j)) / (1 - d) + S_ave_t * passing_period - S_0 ) / (i + passing_period);
    //return res;
	//return roundDouble(res, 4);
}

double compute_A_min_ij(double S_0, double S_ave_t, double passing_period, double i, double j, double u, double d) {
	//double res = (S_0 * (1 - pow(d, j + 1)) / (1 - d) + S_0 * pow(d, j) * u * (1 - pow(u, i - j)) / (1 - u) + S_ave_t * passing_period ) / (i + passing_period + 1.0);
	//double res = (S_0 * (1 - pow(d, j +1 )) / (1 - d) + S_0 * pow(d, j) * u * (1 - pow(u, i - j)) / (1 - u) + S_ave_t * passing_period - S_0 ) / (i + passing_period );
	//return res;
	return (S_0 * (1 - pow(d, j +1 )) / (1 - d) + S_0 * pow(d, j) * u * (1 - pow(u, i - j)) / (1 - u) + S_ave_t * passing_period - S_0 ) / (i + passing_period );
    //return roundDouble(res, 4);
}	


void print_mat(vector< vector<double> >& mat) {
    // print upper triangle
    int n = mat[0].size() - 1;
    for (int j = 0; j <= n; ++j) {

        for (int k = 0; k < j; ++k) {
            printf("%.7f ", 0.0);
        }

        for (int i = j; i <= n; ++i) {

            printf("%lf ", mat[i][j]);
        }
        printf("\n");
    }
}

void print_tree(vector< vector<Node> >& tree) {
    int n = tree[0].size() - 1;
    for ( int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            printf("tree[%d][%d]:\n", i, j);
            tree[i][j].print(); 
        }
    }

}

void binomial(vector< vector<Node> >& tree, double S_t, double K, double r, double q, double sigma, double T_minus_t, double S_ave_t, 
              int M, int n, double u, double d, double p, double dT, double passing_period, int method);

int main(int argc, char const *argv[]) {

    vector< vector<Node> > tree;
    vector< vector<double> > mat_A_max;
    vector< vector<double> > mat_A_min;

    double S_t, K, r, q, sigma, T_minus_t, S_ave_t, passing_time;
    int M, n, n_sim, n_rep;
    double u, d, p, dT;
    double passing_period;

    scanf("%lf %lf %lf %lf %lf %lf %d %d %lf %lf %d %d", &S_t, &K, &r, &q, &sigma, &T_minus_t, &M, &n, &S_ave_t, &passing_time, &n_sim, &n_rep);

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

	dT = (T_minus_t - 0.0) / n;
    passing_period = passing_time / (T_minus_t / n) + 1;
	u = exp(sigma * sqrt(dT));
	d = 1/u;
	p = (exp((r-q) * dT) - d)/(u - d);

	printf("u: %f, d: %f, p: %f, dT: %f, passing_period: %f\n", u, d, p, dT, passing_period);
   
//    vector< vector<Node> >tree(n + 1, vector<Node>(n + 1, Node(M)) );


    tree = vector< vector<Node> >(n + 1, vector<Node>(n + 1, Node(M)) );

    /*
    // (1) calculate A_max, A_min first
    mat_A_max = vector< vector<double> >(n + 1, vector<double>(n + 1, 0.0) );
    mat_A_min = vector< vector<double> >(n + 1, vector<double>(n + 1, 0.0) );
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= n; ++j) {
            mat_A_max[i][j] =  compute_A_max_ij(S_t, S_ave_t, passing_period, i, j, u, d);
            mat_A_min[i][j] =  compute_A_min_ij(S_t, S_ave_t, passing_period, i, j, u, d);
        }
    }*/

    //Node node_ij;
    // (2) calculate A(i, j, k)
    double A_max_ij, A_min_ij;
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= i; ++j) {
            A_max_ij =  compute_A_max_ij(S_t, S_ave_t, passing_period, i, j, u, d);
            A_min_ij =  compute_A_min_ij(S_t, S_ave_t, passing_period, i, j, u, d);
            for (int k = 0; k <= M; ++k) {
                // A_ijk
                //tree[i][j].A_vec[k] = ( (M - k) / (double) M ) * mat_A_max[i][j] + (k / (double) M) * mat_A_min[i][j];
                // A_ijk log version
                //tree[i][j].A_log_vec[k] = exp( ( (M - k) / (double) M ) * log(mat_A_max[i][j]) + (k / (double) M) * log(mat_A_min[i][j]) );
                tree[i][j].A_vec[k] = ( (M - k) / (double) M ) * A_max_ij + (k / (double) M) * A_min_ij;
                // A_ijk log version
                tree[i][j].A_log_vec[k] = exp( ( (M - k) / (double) M ) * log(A_max_ij) + (k / (double) M) * log(A_min_ij) );
            }
        }
    }

    // (3) for terminal node(n, j), decide the payoff
    //Node node_nj;
    for (int j = 0; j <= n; ++j) {
        //node_nj = tree[n][j];
        for (int k = 0; k <= M; ++k) {
            /*
            double A_njk = tree[n][j].A_vec[k];
            double A_log_njk = tree[n][j].A_log_vec[k];
            tree[n][j].C_vec[k] = max(A_njk - K, 0.0);
            tree[n][j].C_log_vec[k] = max(A_log_njk - K, 0.0);
            tree[n][j].C_am_vec[k] = max(A_njk - K, 0.0);
            tree[n][j].C_am_log_vec[k] = max(A_log_njk - K, 0.0);
            */
            tree[n][j].C_vec[k] = max(tree[n][j].A_vec[k]- K, 0.0);
            tree[n][j].C_log_vec[k] = max(tree[n][j].A_log_vec[k] - K, 0.0);
            tree[n][j].C_am_vec[k] = max(tree[n][j].A_vec[k]- K, 0.0);
            tree[n][j].C_am_log_vec[k] = max(tree[n][j].A_log_vec[k] - K, 0.0);
        }
    }

    // (4) Backward
    // Record start time
    binomial(tree, S_t, K, r, q, sigma, T_minus_t, S_ave_t, 
              M, n, u, d, p, dT, passing_period, 0);
    binomial(tree, S_t, K, r, q, sigma, T_minus_t, S_ave_t, 
              M, n, u, d, p, dT, passing_period, 1);
    binomial(tree, S_t, K, r, q, sigma, T_minus_t, S_ave_t, 
              M, n, u, d, p, dT, passing_period, 2);

	return 0;
}

void binomial(vector< vector<Node> >& tree, double S_t, double K, double r, double q, double sigma, double T_minus_t, double S_ave_t, 
              int M, int n, double u, double d, double p, double dT, double passing_period, int method) {
//void binomial(int method) {

    double w_u, w_d;
    double C_u, C_d;
    int k_u, k_d;
    //double A_ijk, A_log_ijk;
    //double A_u, A_u_log;
    //double A_d, A_d_log;
    double A_u;
    double A_d;

    pair<int, int> range;

    auto start = chrono::high_resolution_clock::now();
    for (int i = n - 1; i >= 0; --i) {
        for (int j = 0; j <= i; ++j) {
            for (int k = 0; k <= M; ++k) {
                //A_ijk = tree[i][j].A_vec[k];
                //A_log_ijk = tree[i][j].A_log_vec[k];
                A_u = ( (i + passing_period ) * tree[i][j].A_vec[k] + S_t * pow(u, i + 1 - j) * pow(d, j)  ) / (i + passing_period + 1);
                A_d = ( (i + passing_period ) * tree[i][j].A_vec[k] + S_t * pow(u, i - j) * pow(d, j + 1) ) / (i + passing_period + 1);
                
                // linear
                if (method == 1) {
                    range = sequential_search(tree[i+1][j].A_vec, A_u);
                } else if (method == 0) {
                    range = binary_search(tree[i+1][j].A_vec, 0, M, A_u);
                } else {
                    range = interpolation_search(tree[i+1][j].A_vec, 0, M, A_u);
                }

                if ( range.first == range.second ) {
                    k_u = range.first;
                    w_u = 1;
                } else {
                    k_u = range.second;

                    if (abs(tree[i+1][j].A_vec[k_u - 1] - tree[i+1][j].A_vec[k_u]) < EPSILON) {
                        w_u = 1;
                    } else {
                        w_u = ( tree[i+1][j].A_vec[k_u -1] - A_u ) / ( tree[i+1][j].A_vec[k_u - 1] - tree[i+1][j].A_vec[k_u]);
                    }
                }
                C_u = w_u * tree[i+1][j].C_vec[k_u] + (1 - w_u) * tree[i+1][j].C_vec[k_u - 1];
               
                // linear
                if (method == 1) {
                    range = sequential_search(tree[i+1][j+1].A_vec, A_d);
                } else if (method == 0) {
                    range = binary_search(tree[i+1][j+1].A_vec, 0, M, A_d);
                } else {
                    range = interpolation_search(tree[i+1][j+1].A_vec, 0, M, A_d);
                }

                if ( range.first == range.second ) {
                    k_d = range.first;
                    w_d = 1;
                } else {
                    k_d = range.second;
                    
                    if (abs(tree[i+1][j+1].A_vec[k_d - 1] - tree[i+1][j+1].A_vec[k_d]) < EPSILON) {
                        w_d = 1;
                    } else {
                        w_d = ( tree[i+1][j+1].A_vec[k_d -1] - A_d ) / ( tree[i+1][j+1].A_vec[k_d - 1] - tree[i+1][j+1].A_vec[k_d] );
                    }
                }
                C_d = w_d * tree[i+1][j+1].C_vec[k_d] + (1 - w_d) * tree[i+1][j+1].C_vec[k_d - 1];
                tree[i][j].C_vec[k] = (p * C_u + (1 - p) * C_d) * exp( -r * dT );

                /* ### american ### */
                C_u = w_u * tree[i+1][j].C_am_vec[k_u] + (1 - w_u) * tree[i+1][j].C_am_vec[k_u - 1];
                C_d = w_d * tree[i+1][j+1].C_am_vec[k_d] + (1 - w_d) * tree[i+1][j+1].C_am_vec[k_d - 1];
                tree[i][j].C_am_vec[k] = max( tree[i][j].A_vec[k] - K, (p * C_u + (1-p) * C_d ) * exp( -r * dT) );

                /* ------------------------------------------------------------------------------------------------------------------------------------------*/
                /* log */

                A_u = ( (i + passing_period ) * tree[i][j].A_log_vec[k] + S_t * pow(u, i + 1 - j) * pow(d, j)  ) / (i + passing_period + 1);
                A_d = ( (i + passing_period ) * tree[i][j].A_log_vec[k] + S_t * pow(u, i - j) * pow(d, j + 1) ) / (i + passing_period + 1);
                //A_u_log = ( (i + passing_period ) * A_log_ijk + S_t * pow(u, i + 1 - j) * pow(d, j)  ) / (double) (i + passing_period + 1);
                //A_d_log = ( (i + passing_period ) * A_log_ijk + S_t * pow(u, i - j) * pow(d, j + 1) ) / (i + passing_period + 1);
                if (method == 1) {
                    range = sequential_search(tree[i+1][j].A_log_vec, A_u);
                } else if (method == 0) {
                    range = binary_search(tree[i+1][j].A_log_vec, 0, M, A_u);
                } else {
                    range = interpolation_search(tree[i+1][j].A_log_vec, 0, M, A_u);
                }

                if ( range.first == range.second ) {
                    k_u = range.first;
                    //C_u_log = tree[i+1][j].C_log_vec[k_u_log];
                    w_u = 1;
                } else {
                    k_u = range.second;
                    if (abs(tree[i+1][j].A_log_vec[k_u - 1] - tree[i+1][j].A_log_vec[k_u]) < EPSILON) {
                        w_u = 1;
                    } else {
                        w_u = ( tree[i+1][j].A_log_vec[k_u -1] - A_u) / ( tree[i+1][j].A_log_vec[k_u - 1] - tree[i+1][j].A_log_vec[k_u]);
                    }
                }
                C_u = w_u * tree[i+1][j].C_log_vec[k_u] + (1 - w_u) * tree[i+1][j].C_log_vec[k_u - 1];

                // log
                if (method == 1) {
                    range = sequential_search(tree[i+1][j+1].A_log_vec, A_d);
                } else if (method == 0) {
                    range = binary_search(tree[i+1][j+1].A_log_vec, 0, M, A_d);
                } else {
                    range = interpolation_search(tree[i+1][j+1].A_log_vec, 0, M, A_d);
                    //;
                }

                if ( range.first == range.second ) {
                    k_d = range.first;
                    //C_d_log = tree[i+1][j+1].C_log_vec[k_d_log];
                    w_d = 1;
                } else {
                    k_d = range.second;
                    if (abs(tree[i+1][j+1].A_log_vec[k_d - 1] - tree[i+1][j+1].A_log_vec[k_d]) < EPSILON) {
                        w_d = 1;
                    } else {
                        w_d = ( tree[i+1][j+1].A_log_vec[k_d -1] - A_d ) / ( tree[i+1][j+1].A_log_vec[k_d- 1] - tree[i+1][j+1].A_log_vec[k_d] );
                    }
                }
                C_d = w_d * tree[i+1][j+1].C_log_vec[k_d] + (1 - w_d) * tree[i+1][j+1].C_log_vec[k_d - 1];

                // update C(i, j, k)
                tree[i][j].C_log_vec[k] = (p * C_u + (1 - p) * C_d) * exp( -r * dT );
                // American: max( A(i, j, k) - K, (P * C_u + (1 - P) * C_d) * e^-r * dT
                C_u = w_u * tree[i+1][j].C_am_log_vec[k_u] + (1 - w_u) * tree[i+1][j].C_am_log_vec[k_u - 1];
                C_d = w_d * tree[i+1][j+1].C_am_log_vec[k_d] + (1 - w_d) * tree[i+1][j+1].C_am_log_vec[k_d - 1];

                tree[i][j].C_am_log_vec[k] = max( tree[i][j].A_log_vec[k] - K, (p * C_u + (1-p) * C_d) * exp(-r * dT) );
            }
        }
    }
    // Record end time
    //
    //print_tree(tree);
    
    string search_method;
    if (method == 1) {
        search_method = "Sequential search";
    } else if (method == 0) {
        search_method = "Binary search";
    } else {
        search_method = "Interpolation search";
    }
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    double elapsed_t = elapsed.count();

    printf("-----------------------------------------------------------------------------\n");
	printf("###Binomial Tree Model for European and American arithmetic average calls.###\n");
    printf("Computational time(%s):        %.4lf s\n", search_method.c_str(), elapsed_t);
    printf("European(Linear):    %lf\n", tree[0][0].C_vec[0]);
    printf("European(Logarithm): %lf\n", tree[0][0].C_log_vec[0]);
    printf("American(Linear):    %lf\n", tree[0][0].C_am_vec[0]);
    printf("American(Logarithm): %lf\n", tree[0][0].C_am_log_vec[0]);
    printf("-----------------------------------------------------------------------------\n");

} 
