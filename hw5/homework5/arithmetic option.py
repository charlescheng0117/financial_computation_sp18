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

def Amax(i,j,u,St):
    d = 1/u
#    return
    if i-j>0:
        return (St*(((1-u**(j+1))/(1-u)) + ((u**j) * d * ((1-d**(i-j))/(1-d)))))/(i+1)
    else:
        return  (St*(((1-u**(j+1))/(1-u))))/(i+1)
def Amin(i,j,u,St):
    d = 1/u
    if j>0:
        return (St* (((1-(d**(i-j+1)))/(1-d)) + ((d**(i-j)) * u * ((1-u**(j))/(1-u)))))/(i+1)
    else:
        return (St* (((1-(d**(i-j+1)))/(1-d))))/(i+1) 

def create_node(stock_price,M,i,j,u):
    node = []
    node.append(round(stock_price*(u**(2*j-i)),3))
    Max_av = round(Amax(i,j,u,stock_price),3)
    Min_av = round(Amin(i,j,u,stock_price),3)
    delta_av = (Max_av - Min_av)/(M-1)
    temp = []
    for k in range(0,M):
        temp.append(round(Min_av+k*delta_av,3))
        temp.append(0)
        node.append(list(temp))
        temp = []
    return node

def average(S0,r,q,sigma,T):
    return (np.log(S0)+(r-q-(sigma**2)/2)*T)

def variance(S0,r,q,sigma,T):
    return (sigma**2)*T
'''
def construct_tree(St,Save_t,K,n,r,u,M):
    tree = []
    root = []
    layer = []
    root.append(St)
    root.append([Save_t 0])
    layer.append(root)
    tree.append(layer)
    for i in range(1,n+1):
        layer = []
        nodes_num = i + 1
        for j in range(0,nodes_num):
            stock_price = (i-j)*np.log(u) - j * np.log(u) + np.log(St)
'''
def find_interpolation_object(next_price,temp_nodes):
    for q in range(1,M):
        if next_price >= temp_nodes[q][0] and next_price<=temp_nodes[q+1][0]:
            break
    return round(temp_nodes[q][0],2),round(temp_nodes[q+1][0],2),round(temp_nodes[q][1],3),round(temp_nodes[q+1][1],3)


def interpolation(x_1,x_2,y_1,y_2,x):
    if x_2-x_1==0:
        return (y_2+y_1)/2
    else:
        return y_1+((y_2-y_1)/(x_2-x_1))*(x-x_1)

def construct_tree(St,K,n,r,q,u,M,T):
    delta_t = T/n
    d = 1/u
    neutral_prob = (exp((r-q)*delta_t)-d)/(u-d)
    tree = []
    for i in range(0,n+1):
        layer = []
        nodes_num = i + 1
        for j in range(0,nodes_num):
            layer.append( create_node(St,M,i,j,u) )
        tree.append(layer)
    #tree put_in_answer
    for i in range(0,n+1):
        for j in range(1,M+1):
#            tree[n][i][j][1] = max(tree[n][i][0]-tree[n][i][j][0],0)
            tree[n][i][j][1] = round(max(tree[n][i][j][0]-K,0),3)
    #backwardation
    # i = 0 to n-1
    for i in range(0,n):
    # j  = 0 to n-i
        for j in range(0,n-i):
            for k in range(1,M+1):
                x_up_1,x_up_2,y_up_1,y_up_2=find_interpolation_object( (tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j+1][0])/(n-i+1),tree[n-i][j+1] )
                x_down_1,x_down_2,y_down_1,y_down_2=find_interpolation_object( (tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j][0])/(n-i+1),tree[n-i][j] )
                
                if i == 15 :
                    if j==2 :
                        print(n-i,j,k)
                        print("average price is: {}".format(tree[n-1-i][j][k][0]))
                        print( "Au is: {}".format((tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j+1][0])/(n-i+1) ))
                        print( "Ad is: {}".format((tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j][0])/(n-i+1) )) 
                        print(x_up_1,x_up_2,y_up_1,y_up_2)
                        print( interpolation(x_up_1,x_up_2,y_up_1,y_up_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j+1][0])/(n-i+1) ))
                        print(x_down_1,x_down_2,y_down_1,y_down_2)
                        print( interpolation(x_down_1,x_down_2,y_down_1,y_down_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j][0])/(n-i+1) ))               
                
                tree[n-1-i][j][k][1] = max( round(( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j+1][0])/(n-i+1) ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j][0])/(n-i+1) ))*exp(-(r-q)*delta_t),3) , tree[n-1-i][j][k][0]-K , 0) 

    return tree
    #print(tree)


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

St = 50
K=50
r = 0.1
q = 0
sigma = 0.4
diff_T  = 1
M = 4
n = 20
Save_t = 60
passing_time = 1
simulations = 1000
repetitions = 40



delta_t = diff_T/n
u = exp( sigma*(delta_t**0.5) )
d = 1/u
neutral_prob = (exp((r-q)*delta_t)-d)/(u-d)
#print(u,d,neutral_prob)
iteration_t = passing_time/delta_t
iteration_t = int(iteration_t)
#print(neutral_prob)
'''
#monte carlo
time_scale = 1/365
sample_point = int( np.ceil( diff_T*365) )
all_repetitions = []
for i in range(0,repetitions):
    all_simulations = []
    for j in range(0,simulations):
        one_simulations = [St]    
        for k in range(0,sample_point):
            mu = average( one_simulations[-1],r,q,sigma,time_scale )
            std = variance( one_simulations[-1],r,q,sigma,time_scale ) ** 0.5
            one_simulations.append( exp( np.random.normal( mu , std ) ) )
        Save_T = float( sum(one_simulations[1:]) + Save_t * ( np.ceil(passing_time * 365) + 1 ) ) / ( np.ceil(passing_time*365) + len(one_simulations) )
        all_simulations.append( max(Save_T-K,0) )
    all_repetitions.append( np.mean(all_simulations)*exp( -r * diff_T ) )

print(np.mean(all_repetitions))
print(np.std(all_repetitions))
'''



#S0*(u**(i-j))*(d**(j))=St ===> S0 = St * u**j * d**(i-j) = St  * u ** ( 2j - i ) ==> 

tree_1 = construct_tree(St,K,n,r,q,u,M,diff_T)
#print(tree_1[4])
#print(tree_1)
print(tree_1[0][0][1][1])