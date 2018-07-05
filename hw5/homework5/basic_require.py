"""

Price an arithmetic average call with the following payoff using the binomial tree model.
Payofft = max(Save,t − K, 0),
where Save,t is the arithmetic average of stock prices from the issue date until the current
time point t.

Inputs: St , K, r, q, σ, T − t, M, n, Save,t, passing time, number of simulations, number of repetitions. 
Outputs: Option values for both methods and 95% confidence interval for Monte Carlo simulation

"""
import numpy as np 
from math import * 
from numpy.linalg import *

def Amax(i,j,u,St,Save_t,t,delta_t):
    passing_period = t/delta_t + 1
    d = 1/u
    return ((St * (1 - u ** (i - j + 1) ) / (1 - u) + St * u ** (i - j) * d * (1 - d ** j)) / (1 - d) + Save_t * passing_period - St ) / (i + passing_period)
    #return (St * (1 - np.power(u, i - j + 1) ) / (1 - u) + St * np.power(u, i - j) * d * (1 - np.power(d, j)) / (1 - d) + Save_t * passing_period - St ) / (i + passing_period)

def Amin(i,j,u,St,Save_t,t,delta_t):
    passing_period = t/delta_t + 1
    d = 1/u
    return (St * (1 - np.power(d, j +1 )) / (1 - d) + St * np.power(d, j) * u * (1 - np.power(u, i - j)) / (1 - u) + Save_t * passing_period - St ) / (i + passing_period )

def create_node(stock_price,M,i,j,u,passing_time,delta_t):
    node = []
    node.append(stock_price*(u**(2*j-i)))
    Max_av = Amax(i,j,u,stock_price,Save_t,passing_time,delta_t)
    Min_av = Amin(i,j,u,stock_price,Save_t,passing_time,delta_t)
    delta_av = (Max_av - Min_av)/(M-1)
    temp = []
    for k in range(0,M):
        temp.append(Min_av+k*delta_av)
        temp.append(0)
        node.append(list(temp))
        temp = []
    return node

def average(S0,r,q,sigma,T):
    return (np.log(S0)+(r-q-(sigma**2)/2)*T)

def variance(S0,r,q,sigma,T):
    return (sigma**2)*T

def find_interpolation_object(next_price,temp_nodes):
    q = 1
    if next_price<temp_nodes[q][0]:
        return temp_nodes[q][0],temp_nodes[q+1][0],temp_nodes[q][1],temp_nodes[q+1][1]
    for q in range(1,M):
        if next_price >= temp_nodes[q][0] and next_price<=temp_nodes[q+1][0]:
            return temp_nodes[q][0],temp_nodes[q+1][0],temp_nodes[q][1],temp_nodes[q+1][1]
    return temp_nodes[q+1][0],temp_nodes[q+1][0],temp_nodes[q+1][1],temp_nodes[q+1][1]

def interpolation(x_1,x_2,y_1,y_2,x):
    if abs(x_2-x_1)<0.00000000001:
        return max((y_2+y_1)/2,0)
    else:
        '''
        if y_1+((y_2-y_1)/(x_2-x_1))*(x-x_1)<0:
            print(x_1,x_2,y_1,y_2,x)
        '''
        return max(y_1+((y_2-y_1)/(x_2-x_1))*(x-x_1),0)

def construct_tree(St,K,n,r,q,u,M,T,passing_time):
    delta_t = T/n
    d = 1/u
    neutral_prob = (exp((r-q)*delta_t)-d)/(u-d)
    tree = []
    for i in range(0,n+1):
        layer = []
        nodes_num = i + 1
        for j in range(0,nodes_num):
            layer.append( create_node(St,M,i,j,u,passing_time,delta_t) )
        tree.append(layer)
    #tree put_in_answer
    for i in range(0,n+1):
        for j in range(1,M+1):
#            tree[n][i][j][1] = max(tree[n][i][0]-tree[n][i][j][0],0)
            tree[n][i][j][1] = max( tree[n][i][j][0]-K , 0 )
    #backwardation
    # i = 0 to n-1
    for i in range(0,n):
    # j  = 0 to n-i
        for j in range(0,n-i):
            for k in range(1,M+1):
                A_u = ( tree[n-1-i][j][k][0]*( n-i+int(passing_time/delta_t) +1 ) + tree[n-i][j+1][0] )/(n-i+1+int(passing_time/delta_t) +1 )
                A_d = ( tree[n-1-i][j][k][0]*( n-i+int(passing_time/delta_t)+1 ) + tree[n-i][j][0] ) /(n-i+1+int(passing_time/delta_t) +1)
                x_up_1,x_up_2,y_up_1,y_up_2=find_interpolation_object( A_u/(n-i+1+int(passing_time/delta_t) +1) , tree[n-i][j+1] )
                x_down_1,x_down_2,y_down_1,y_down_2=find_interpolation_object( A_d/(n-i+1+int(passing_time/delta_t)+1) , tree[n-i][j] )
                tree[n-1-i][j][k][1] = max( ( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,A_u ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,A_d) )*exp(-r*delta_t),0)
#                tree[n-1-i][j][k][1] = max( ( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j+1][0])/(n-i+1) ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j][0])/(n-i+1) ))*exp(-r*delta_t) , tree[n-1-i][j][k][0]-K ,0) 

    return tree
    #print(tree)

St = 50
K=50
r = 0.1
q = 0.05
sigma = 0.8
diff_T  = 0.25
M = 100
n = 100
Save_t = 50
passing_time = 0
simulations = 10000
repetitions = 20


delta_t = diff_T/n
u = exp( sigma*(delta_t**0.5) )
d = 1/u
neutral_prob = (exp((r-q)*delta_t)-d)/(u-d)

tree_1 = construct_tree(St,K,n,r,q,u,M,diff_T,passing_time)

'''
t_period = passing_time/delta_t
all_repetitions = []
for i in range(0,repetitions):
    all_simulations = []
    for j in range(0,simulations):
        one_simulations = [St]    
        for k in range(0,n):
            mu = average( one_simulations[-1],r,q,sigma,delta_t )
            std = variance( one_simulations[-1],r,q,sigma,delta_t ) ** 0.5
            one_simulations.append( exp( np.random.normal( mu , std ) ) )
        Save_T = float( sum(one_simulations[1:]) + Save_t * ( t_period + 1 ) ) / ( t_period + len(one_simulations) )
        all_simulations.append( max(Save_T-K,0) )
    all_repetitions.append( np.mean(all_simulations)*exp( -r * diff_T ) )

print(np.mean(all_repetitions))
print(np.std(all_repetitions))
'''
print("Asian Call Price is: {}".format(tree_1[0][0][1][1]))
print(tree_1[100][100])
print(tree_1[100][99])
#print(tree_1[0][0][1][1])
