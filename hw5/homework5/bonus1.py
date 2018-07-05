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
import matplotlib.pyplot as plt
def Amax(i,j,u,St,Save_t,t,delta_t):
    d = 1/u

    return ( St *  ( ( u*(1-(u**j) )/(1-u) ) + ( (u**j) * d * ((1-(d**(i-j)))/(1-d) ) ) ) + Save_t*( (passing_time/delta_t)+1 )  )/ ((passing_time/delta_t)+i+1)
#    else:
#        return  ( St* ( ( u* (1- (u**j) )/(1-u) ) ) + Save_t*( (passing_time/delta_t)+1 ) ) /(i+1+(passing_time/delta_t) )

def Amin(i,j,u,St,Save_t,t,delta_t):
    d = 1/u
    #if j>0:
    return (St* ( ( d*(1-(d**(i-j)) )/(1-d)) + ( (d**(i-j) ) * u * ( (1-(u**j) )/(1-u)))) + Save_t*( (passing_time/delta_t)+1 ) ) / ((passing_time/delta_t)+i+1)
    #else:
    #    return (St* ( ( d* (1-(d**(i-j)) )/(1-d))) + Save_t*( (passing_time/delta_t)+1 ) ) / ( i+1+(passing_time/delta_t) ) 

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

def create_node_log(stock_price,M,i,j,u,passing_time,delta_t):
    node = []
    node.append(stock_price*(u**(2*j-i)))
    Max_av = Amax(i,j,u,stock_price,Save_t,passing_time,delta_t)
    Min_av = Amin(i,j,u,stock_price,Save_t,passing_time,delta_t)
    delta_av = ( np.log(Max_av) - np.log(Min_av) )/(M-1)
    temp = []
    for k in range(0,M):
        temp.append( Min_av*(np.exp(delta_av)**k) )
#        print(np.exp(Min_av+k*delta_av))
        temp.append(0)
        node.append(list(temp))
        temp = []
    return node


def average(S0,r,q,sigma,T):
    return (np.log(S0)+(r-q-(sigma**2)/2)*T)

def variance(S0,r,q,sigma,T):
    return (sigma**2)*T

def find_interpolation_object(next_price,temp_nodes,M):
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
                A_u = ( tree[n-1-i][j][k][0]*( n-i+int(passing_time/delta_t) ) + tree[n-i][j+1][0] )/(n-i+1+int(passing_time/delta_t) )
                A_d = ( tree[n-1-i][j][k][0]*( n-i+int(passing_time/delta_t) ) + tree[n-i][j][0] ) /(n-i+1+int(passing_time/delta_t) )
                x_up_1,x_up_2,y_up_1,y_up_2=find_interpolation_object( A_u/(n-i+1+int(passing_time/delta_t) ) , tree[n-i][j+1] ,M)
                x_down_1,x_down_2,y_down_1,y_down_2=find_interpolation_object( A_d/(n-i+1+int(passing_time/delta_t)) , tree[n-i][j],M )
                tree[n-1-i][j][k][1] = max( ( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,A_u ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,A_d) )*exp(-r*delta_t),0)
#                tree[n-1-i][j][k][1] = max( ( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j+1][0])/(n-i+1) ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j][0])/(n-i+1) ))*exp(-r*delta_t) , tree[n-1-i][j][k][0]-K ,0) 
    return tree
    #print(tree)

def construct_tree_1(St,K,n,r,q,u,M,T,passing_time):
    delta_t = T/n
    d = 1/u
    neutral_prob = (exp((r-q)*delta_t)-d)/(u-d)
    tree = []
    for i in range(0,n+1):
        layer = []
        nodes_num = i + 1
        for j in range(0,nodes_num):
            layer.append( create_node_log(St,M,i,j,u,passing_time,delta_t) )
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
                A_u = ( tree[n-1-i][j][k][0]*( n-i+int(passing_time/delta_t) ) + tree[n-i][j+1][0] )/(n-i+1+int(passing_time/delta_t) )
                A_d = ( tree[n-1-i][j][k][0]*( n-i+int(passing_time/delta_t) ) + tree[n-i][j][0] ) /(n-i+1+int(passing_time/delta_t) )
                x_up_1,x_up_2,y_up_1,y_up_2=find_interpolation_object( A_u/(n-i+1+int(passing_time/delta_t) ) , tree[n-i][j+1],M )
                x_down_1,x_down_2,y_down_1,y_down_2=find_interpolation_object( A_d/(n-i+1+int(passing_time/delta_t)) , tree[n-i][j],M )
                tree[n-1-i][j][k][1] = max( ( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,A_u ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,A_d) )*exp(-r*delta_t),0)
#                tree[n-1-i][j][k][1] = max( ( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j+1][0])/(n-i+1) ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j][0])/(n-i+1) ))*exp(-r*delta_t) , tree[n-1-i][j][k][0]-K ,0) 
    return tree


'''
str_St = input('Enter St: ')
str_K = input('Enter K: ')
str_r = input('Enter r: ')
str_q = input('Enter q: ')
str_sigma = input('Enter sigma: ')
str_diff_T = input('Enter T-t: ')
str_M = input('Enter M: ')
str_n = input('Enter n: ')
str_Save_t = input('Enter Save,t: ')
str_passing_time = input('Enter passing_time: ')
str_simulations = input('Enter simulations: ')
str_repetitions = input('Enter repetitions: ')
'''

'''
St = float(str_St)
K = float(str_K)
r = float(str_r)
q = float(str_q)
sigma = float(str_sigma)
diff_T = float(str_diff_T)
M = float(str_M)
n = float(str_n)
Save_t = float(str_Save_t)
passing_time = int(str_passing_time)
simulations = int(str_simulations)
repetitions = int(str_repetitions)
'''
#50 40 0.1 0.05 0.3 0.5 10 100 60 1 10000 20
St = 50
K=40
r = 0.1
q = 0.05
sigma = 0.3
diff_T  = 0.5
n = 10
Save_t = 60
passing_time = 0.1
simulations = 10000
repetitions = 40



delta_t = diff_T/n
u = exp( sigma*(delta_t**0.5) )
d = 1/u
neutral_prob = (exp((r-q)*delta_t)-d)/(u-d)
#iteration_t = passing_time/delta_t
#iteration_t = int(iteration_t)

#true_price = construct_tree(St,K,100,r,q,u,M,diff_T,passing_time)
#tree_1 = construct_tree_1(St,K,n,r,q,u,M,diff_T,passing_time)
#tree_2 = construct_tree(St,K,n,r,q,u,M,diff_T,passing_time)

error_1 = []
error_2 = []
M = [50,100,150,200,250,300,350,400,450,500]
for ele in M:
    print(ele)
    tree_2 = construct_tree_1(St,K,n,r,q,u,ele,diff_T,passing_time)
    tree_1 = construct_tree(St,K,n,r,q,u,ele,diff_T,passing_time)
    print(tree_1[0][0][1][1],tree_2[0][0][1][1])
    error_1.append(tree_1[0][0][1][1])
    error_2.append(tree_2[0][0][1][1])
plt.plot(M,error_1,label='Linearly')
plt.plot(M,error_2,label='Logarithmically')
plt.legend(loc='best')
plt.title("Option Value versus M")
plt.xlabel("M Value")
plt.ylabel("Option Value")
plt.grid(True)

plt.show()


#print("Asian Call Price is: {}".format(tree_1[0][0][1][1]))