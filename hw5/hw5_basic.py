import numpy as np

EPSILON = 0.00000001

def sequential_search(vec, target):
    left, right = 0, 0 
    n = len(vec)
    for i in range(n-1):  
        if (abs(target - vec[i]) < EPSILON ):
            return i, i 
        if ( vec[i] > target and target > vec[i + 1] ):
            return (i, i+1) 




    if (abs(target - vec[n-1]) < 10e-4 ) : 
        return (n -1, n-1)

    return (n -1, n - 1)

def binary_search(arr, left,  right, key):  
    mid = int(left + ( right - left ) / 2)  
    if ( abs ( arr[mid] - key ) < 0.0000001 ):
        return (mid, mid) 
    elif ( arr[mid] > key ):  
        return binary_search ( arr, mid + 1, right, key) 
    elif ( arr[mid] < key ):
        if (abs ( arr[mid-1] - key ) < 0.0000001):
            return (mid-1, mid-1)  
        elif (arr[mid-1] > key):
            return (mid-1, mid)  
        else: 
            return binary_search(arr, left, mid - 1, key ) 

def interpolation_search(arr,  left,  right, key): 
    mid = 0 
    if ( left == right ):
        mid = left  
    elif ( ( arr[left] - arr[right] ) < 0.0000001 ):
        return (left, left)  
    else:
        mid = left + ( right - left ) * ( key - arr[left] ) / ( arr[right] - arr[left] )  


    if(abs( arr[mid] - key ) < 0.0000001 ):
        return (mid, mid) 
    elif ( arr[mid] > key ):
        return interpolation_search( arr, mid + 1, mid + 1, key) 
    elif ( arr[mid] < key ): 
        if ( abs ( arr[mid-1] - key ) < 0.0000001 ):
            return (mid-1, mid-1) 
        elif ( arr[mid-1] > key ):
            return (mid-1, mid) 
        else:
            return interpolation_search(arr, mid - 1, mid-1, key) 

class Node(object):  
    M = 0 
    A_vec = []
    A_log_vec = []
    C_vec =[]
    C_log_vec = []
    C_am_vec  = []
    C_am_log_vec =[]

 
    def __init__(self, new_M):  
        self.M = new_M
        self.A_vec = [0.0 for i in range(new_M+1)]
        self.A_log_vec = [0.0 for i in range(new_M+1)]
        self.C_vec = [0.0 for i in range(new_M+1)] 
        self.C_log_vec = [0.0 for i in range(new_M+1)] 
        self.C_am_vec = [0.0 for i in range(new_M+1)]
        self.C_am_log_vec = [0.0 for i in range(new_M+1)]


def compute_A_max_ij( S_0,  S_ave_t,  passing_period,  i,  j,  u,  d):
# res = (S_0 * (1 - np.power(u, i - j + 1) ) / (1 - u) + S_0 * np.power(u, i - j) * d * (1 - np.power(d, j)) / (1 - d) + S_ave_t * passing_period ) / (i + passing_period + 1.0) 
	res = (S_0 * (1 - np.power(u, i - j + 1) ) / (1 - u) + S_0 * np.power(u, i - j) * d * (1 - np.power(d, j)) / (1 - d) + S_ave_t * passing_period - S_0 ) / (i + passing_period) 
	return res 
#return roundDouble(res, 4) 


def compute_A_min_ij( S_0,  S_ave_t,  passing_period,  i,  j,  u,  d):  
# res = (S_0 * (1 - np.power(d, j + 1)) / (1 - d) + S_0 * np.power(d, j) * u * (1 - np.power(u, i - j)) / (1 - u) + S_ave_t * passing_period ) / (i + passing_period + 1.0) 
	res = (S_0 * (1 - np.power(d, j +1 )) / (1 - d) + S_0 * np.power(d, j) * u * (1 - np.power(u, i - j)) / (1 - u) + S_ave_t * passing_period - S_0 ) / (i + passing_period ) 
	return res 
#return roundDouble(res, 4) 



#vector< vector<Node> > tree 
#vector<  > mat_A_max = 
#vector<  > mat_A_min =

in_f = open("1.in", "r")

for line in in_f:
	line = line.strip().split(" ")
	line = [float(ele) for ele in line]
	S_t, K, r, q, sigma, T_minus_t, M, n, S_ave_t, passing_time, n_sim, n_rep = line

M = int(M)
n = int(n)
n_sim = int(n_sim)
n_rep = int(n_rep)
	#S_t, K, &r, &q, &sigma, &T_minus_t, &M, &n, &S_ave_t, &passing_time, &n_sim, &n_rep = 
 #M, n, n_sim, n_rep 
 #    u, d, p, dT 
 #passing_period 

#scanf("%lf %lf %lf %lf %lf %lf %d %d %lf %lf %d %d", &S_t, &K, &r, &q, &sigma, &T_minus_t, &M, &n, &S_ave_t, &passing_time, &n_sim, &n_rep) 

print("S_t          = %f" %S_t) 
print("K            = %f" %K) 
print("r            = %f" %r) 
print("q            = %f" %q) 
print("sigma        = %f" %sigma) 
print("T_minus_t    = %f" % T_minus_t) 
print("M            = %d" % M) 
print("n            = %d" % n) 
print("S_ave_t      = %f" %S_ave_t) 
print("passing_time = %f" %passing_time) 
print("n_sim        = %d" %n_sim) 
print("n_rep        = %d" %n_rep) 

dT = (T_minus_t - 0.0) / n 
passing_period = passing_time / (T_minus_t / n) + 1 
u = np.exp(sigma * np.sqrt(dT)) 
d = 1/u 
p = (np.exp((r-q) * dT) - d)/(u - d) 

#vector< vector<Node> >tree(n + 1, vector<Node>(n + 1, Node(M)) ) 
tree = [ [Node(M) for i in range(n+1)] for j in range(n+1) ]

#tree = vector< vector<Node> >(n + 1, vector<Node>(n + 1, Node(M)) ) 
mat_A_max = [ [ 0.0 for i in range(n+1)] for j in range(n+1) ]
#mat_A_max = vector<  >(n + 1, (n + 1, 0.0) ) 
mat_A_min = [ [ 0.0 for i in range(n+1)] for j in range(n+1) ]

# (1) calculate A_max, A_min first
for i in range(n+1):
    for j in range(n+1):
        mat_A_max[i][j] =  compute_A_max_ij(S_t, S_ave_t, passing_period, i, j, u, d) 
        mat_A_min[i][j] =  compute_A_min_ij(S_t, S_ave_t, passing_period, i, j, u, d) 



#Node node_ij 
# (2) calculate A(i, j, k)
for i in range(n+1):
    for j in range(i+1):
        for k in range(M+1):
            # A_ijk
            tree[i][j].A_vec[k] = ( (M - k) /  M ) * mat_A_max[i][j] + (k /  M) * mat_A_min[i][j] 
            # A_ijk log version
            tree[i][j].A_log_vec[k] = np.exp( ( (M - k) /  M ) * np.log(mat_A_max[i][j]) + (k /  M) * np.log(mat_A_min[i][j]) ) 




# (3) for terminal node(n, j), decide the payoff
#Node node_nj 
#for ( j = 0  j <= n  ++j)  
for j in range(n+1):
    for k in range(M+1):
    #node_nj = tree[n][j] 
    #for ( k = 0  k <= M  ++k)  
        A_njk = tree[n][j].A_vec[k] 
        A_log_njk = tree[n][j].A_log_vec[k] 
        tree[n][j].C_vec[k] = max(A_njk - K, 0.0) 
        tree[n][j].C_log_vec[k] = max(A_log_njk - K, 0.0) 
        tree[n][j].C_am_vec[k] = max(A_njk - K, 0.0) 
        tree[n][j].C_am_log_vec[k] = max(A_log_njk - K, 0.0) 



# (4) Backward
# Record start time

range_, log_range_ = 0, 0

# computation_time 
#for ( t = 0  t < 3  ++t)   # three method
for i in range(n-1, -1, -1):
#for ( i = n - 1  i >= 0  --i)  
    for j in range(i+1):
    #for ( j = 0  j <= i  ++j)  
        for k in range(M+1):
        #for ( k = 0  k <= M  ++k)  
            A_ijk = tree[i][j].A_vec[k] 
            A_log_ijk = tree[i][j].A_log_vec[k] 
            # A_u = ( (i + passing_period + 1) * A_ijk + S_t * np.power(u, i + 1 - j) * np.power(d, j)  ) / (i + passing_period + 2) 
            # A_u_log = ( (i + passing_period + 1) * A_log_ijk + S_t * np.power(u, i + 1 - j) * np.power(d, j)  ) / (i + passing_period + 2) 
            A_u = ( (i + passing_period ) * A_ijk + S_t * np.power(u, i + 1 - j) * np.power(d, j)  ) /  (i + passing_period + 1) 
            A_u_log = ( (i + passing_period ) * A_log_ijk + S_t * np.power(u, i + 1 - j) * np.power(d, j)  ) /  (i + passing_period + 1) 

            
            range_ = binary_search(tree[i+1][j].A_vec, 0, M, A_u) 
            


            if ( range_[0] == range_[1] ): 
                k_u = range_[0] 
                #C_u = tree[i+1][j].C_vec[k_u] 
                w_u = 1 
            else:
                k_u = range_[1] 

                if (abs(tree[i+1][j].A_vec[k_u - 1] - tree[i+1][j].A_vec[k_u]) < EPSILON):
                    w_u = 1 
                else:
                    w_u = ( tree[i+1][j].A_vec[k_u -1] - A_u ) / ( tree[i+1][j].A_vec[k_u - 1] - tree[i+1][j].A_vec[k_u]) 



            C_u = w_u * tree[i+1][j].C_vec[k_u] + (1 - w_u) * tree[i+1][j].C_vec[k_u - 1] 

            
            log_range_ = binary_search(tree[i+1][j].A_log_vec, 0, M, A_u_log) 
            


            if ( log_range_[0] == log_range_[1] ):
                k_u_log = log_range_[0] 
                #C_u_log = tree[i+1][j].C_log_vec[k_u_log] 
                w_u_log = 1 
            else:
                k_u_log = log_range_[1] 
                if (abs(tree[i+1][j].A_log_vec[k_u_log - 1] - tree[i+1][j].A_log_vec[k_u_log]) < EPSILON): 
                    w_u_log = 1 
                else:
                    w_u_log = ( tree[i+1][j].A_log_vec[k_u_log -1] - A_u_log ) / ( tree[i+1][j].A_log_vec[k_u_log - 1] - tree[i+1][j].A_log_vec[k_u_log]) 

                #C_u_log = w_u_log * tree[i+1][j].C_log_vec[k_u_log] + (1 - w_u_log) * tree[i+1][j].C_log_vec[k_u_log - 1] 

            C_u_log = w_u_log * tree[i+1][j].C_log_vec[k_u_log] + (1 - w_u_log) * tree[i+1][j].C_log_vec[k_u_log - 1] 
            A_d = ( (i + passing_period ) * A_ijk + S_t * np.power(u, i - j) * np.power(d, j + 1) ) / (i + passing_period + 1) 
            A_d_log = ( (i + passing_period ) * A_log_ijk + S_t * np.power(u, i - j) * np.power(d, j + 1) ) / (i + passing_period + 1) 

        
            range_ = binary_search(tree[i+1][j+1].A_vec, 0, M, A_d) 
            


            if ( range_[0] == range_[1] ):
                k_d = range_[0] 
                #C_d = tree[i+1][j+1].C_vec[k_d] 
                w_d = 1 
            else:
                k_d = range_[1] 

                if (abs(tree[i+1][j+1].A_vec[k_d - 1] - tree[i+1][j+1].A_vec[k_d]) < EPSILON):
                    w_d = 1 
                else:
                    w_d = ( tree[i+1][j+1].A_vec[k_d -1] - A_d ) / ( tree[i+1][j+1].A_vec[k_d - 1] - tree[i+1][j+1].A_vec[k_d] ) 

                #C_d = w_d * tree[i+1][j+1].C_vec[k_d] + (1 - w_d) * tree[i+1][j+1].C_vec[k_d - 1] 

            C_d = w_d * tree[i+1][j+1].C_vec[k_d] + (1 - w_d) * tree[i+1][j+1].C_vec[k_d - 1] 

            
            log_range_ = binary_search(tree[i+1][j+1].A_log_vec, 0, M, A_d_log) 
            


            if ( log_range_[0] == log_range_[1] ):
                k_d_log = log_range_[0] 
                #C_d_log = tree[i+1][j+1].C_log_vec[k_d_log] 
                w_d_log = 1 
            else:
                k_d_log = log_range_[1] 
                if (abs(tree[i+1][j+1].A_log_vec[k_d_log - 1] - tree[i+1][j+1].A_log_vec[k_d_log]) < EPSILON):
                    w_d_log = 1 
                else:  
                    w_d_log = ( tree[i+1][j+1].A_log_vec[k_d_log -1] - A_d_log ) / ( tree[i+1][j+1].A_log_vec[k_d_log - 1] - tree[i+1][j+1].A_log_vec[k_d_log] ) 

                #C_d_log = w_d_log * tree[i+1][j+1].C_log_vec[k_d_log] + (1 - w_d_log) * tree[i+1][j+1].C_log_vec[k_d_log - 1] 

            C_d_log = w_d_log * tree[i+1][j+1].C_log_vec[k_d_log] + (1 - w_d_log) * tree[i+1][j+1].C_log_vec[k_d_log - 1] 




            # update C(i, j, k)
            tree[i][j].C_vec[k] = (p * C_u + (1 - p) * C_d) * np.exp( -r * dT ) 
            tree[i][j].C_log_vec[k] = (p * C_u_log + (1 - p) * C_d_log) * np.exp( -r * dT ) 
            # American: max( A(i, j, k) - K, (P * C_u + (1 - P) * C_d) * e^-r * dT

            C_u_am = w_u * tree[i+1][j].C_am_vec[k_u] + (1 - w_u) * tree[i+1][j].C_am_vec[k_u - 1] 
            C_d_am = w_d * tree[i+1][j+1].C_am_vec[k_d] + (1 - w_d) * tree[i+1][j+1].C_am_vec[k_d - 1] 
            C_u_log_am = w_u_log * tree[i+1][j].C_am_log_vec[k_u_log] + (1 - w_u_log) * tree[i+1][j].C_am_log_vec[k_u_log - 1] 
            C_d_log_am = w_d_log * tree[i+1][j+1].C_am_log_vec[k_d_log] + (1 - w_d_log) * tree[i+1][j+1].C_am_log_vec[k_d_log - 1] 


            tree[i][j].C_am_vec[k] = max( A_ijk - K, (p * C_u_am + (1-p) * C_d_am ) * np.exp( -r * dT) ) 
            tree[i][j].C_am_log_vec[k] = max( A_log_ijk - K, (p * C_u_log_am + (1-p) * C_d_log_am) * np.exp(-r * dT) ) 


print("Asian(Lin):    {}".format(tree[0][0].C_vec[0]))
print("Asian(Log): {}".format(tree[0][0].C_log_vec[0]))
print("Asian(american)(Lin):    {}".format(tree[0][0].C_am_vec[0]))
print("Asian(american(Log): {}".format(tree[0][0].C_am_log_vec[0]))



