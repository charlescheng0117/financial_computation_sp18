import numpy as np 
from math import * 
from numpy.linalg import *

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

def find_interpolation_object_binary(next_price,temp_nodes):
    q = int((len(temp_nodes)-1)/2)
    temp = list(temp_nodes)
#    print(next_price,temp_nodes[q][0],temp_nodes[q+1][0])
    if next_price<temp_nodes[q][0]:
        return temp_nodes[q][0],temp_nodes[q+1][0],temp_nodes[q][1],temp_nodes[q+1][1]
    elif next_price >= temp_nodes[q][0] and next_price<=temp_nodes[q+1][0]:
        return temp_nodes[q][0],temp_nodes[q+1][0],temp_nodes[q][1],temp_nodes[q+1][1]
    elif next_price < temp[q][0]:
        return find_interpolation_object_binary(next_price,temp[0:q+1])
    else:
        return find_interpolation_object_binary(next_price,temp[q:])

def find_interpolation_object_interpolation(next_price,temp_nodes):
    q = 1
    if next_price<temp_nodes[q][0]:
        return temp_nodes[q][0],temp_nodes[q+1][0],temp_nodes[q][1],temp_nodes[q+1][1]
    else:
        q = int( (next_price-temp_nodes[1][0])/((temp_nodes[M-1][0]-temp_nodes[1][0])/M) )+1
        return temp_nodes[q][0],temp_nodes[q+1][0],temp_nodes[q][1],temp_nodes[q+1][1]
#    return temp_nodes[q+1][0],temp_nodes[q+1][0],temp_nodes[q+1][1],temp_nodes[q+1][1]

def interpolation(x_1,x_2,y_1,y_2,x):
    if abs(x_2-x_1)<0.00000000001:
        return max((y_2+y_1)/2,0)
    else:
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
                x_up_1,x_up_2,y_up_1,y_up_2=find_interpolation_object( A_u/(n-i+1+int(passing_time/delta_t) ) , tree[n-i][j+1] )
                x_down_1,x_down_2,y_down_1,y_down_2=find_interpolation_object( A_d/(n-i+1+int(passing_time/delta_t)) , tree[n-i][j] )
                tree[n-1-i][j][k][1] = max( ( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,A_u ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,A_d) )*exp(-r*delta_t),0)
#                tree[n-1-i][j][k][1] = max( ( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j+1][0])/(n-i+1) ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j][0])/(n-i+1) ))*exp(-r*delta_t) , tree[n-1-i][j][k][0]-K ,0) 
    return tree


def construct_tree_binary(St,K,n,r,q,u,M,T,passing_time):
    delta_t = T/n
    d = 1/u
    neutral_prob = (exp((r-q)*delta_t)-d)/(u-d)
    tree = []
    tstep = int(passing_time/delta_t)
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
                A_u = ( tree[n-1-i][j][k][0]*( n-i+tstep ) + tree[n-i][j+1][0] )/(n-i+1+tstep )
                A_d = ( tree[n-1-i][j][k][0]*( n-i+tstep ) + tree[n-i][j][0] ) /(n-i+1+tstep )
                x_up_1,x_up_2,y_up_1,y_up_2=find_interpolation_object_binary( A_u/(n-i+1+tstep ) , tree[n-i][j+1] )
                x_down_1,x_down_2,y_down_1,y_down_2=find_interpolation_object_binary( A_d/(n-i+1+tstep) , tree[n-i][j] )
                print(x_up_1,x_up_2,y_up_1,y_up_2)
                tree[n-1-i][j][k][1] = max( ( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,A_u ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,A_d) )*exp(-r*delta_t),0)
#                tree[n-1-i][j][k][1] = max( ( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j+1][0])/(n-i+1) ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j][0])/(n-i+1) ))*exp(-r*delta_t) , tree[n-1-i][j][k][0]-K ,0) 
    return tree
def construct_tree_interpolation(St,K,n,r,q,u,M,T,passing_time):
    delta_t = T/n
    d = 1/u
    neutral_prob = (exp((r-q)*delta_t)-d)/(u-d)
    tstep = int(passing_time/delta_t)
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
                A_u = ( tree[n-1-i][j][k][0]*( n-i+tstep ) + tree[n-i][j+1][0] )/(n-i+1+tstep )
                A_d = ( tree[n-1-i][j][k][0]*( n-i+tstep ) + tree[n-i][j][0] ) /(n-i+1+tstep )
                x_up_1,x_up_2,y_up_1,y_up_2=find_interpolation_object_interpolation( A_u/(n-i+1+tstep ) , tree[n-i][j+1] )
                x_down_1,x_down_2,y_down_1,y_down_2=find_interpolation_object_interpolation( A_d/(n-i+1+tstep) , tree[n-i][j] )
                tree[n-1-i][j][k][1] = max( ( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,A_u ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,A_d) )*exp(-r*delta_t),0)
#                tree[n-1-i][j][k][1] = max( ( neutral_prob * interpolation(x_up_1,x_up_2,y_up_1,y_up_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j+1][0])/(n-i+1) ) + (1-neutral_prob)*interpolation(x_down_1,x_down_2,y_down_1,y_down_2,(tree[n-1-i][j][k][0]*(n-i)+tree[n-i][j][0])/(n-i+1) ))*exp(-r*delta_t) , tree[n-1-i][j][k][0]-K ,0) 
    return tree



St = 50
K=40
r = 0.1
q = 0.05
sigma = 0.3
diff_T  = 0.5
M = 100
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

tree_1 = construct_tree_binary(St,K,n,r,q,u,M,diff_T,passing_time)


print("Asian Call Price is: {}".format(tree_1[0][0][1][1]))