import numpy as np
import os
import sys

if len(sys.argv) != 2:
    print("Usage: python3 hw4_bonus1.py <input>")
    sys.exit()

epsilon = 0.01  # compare 2 float
DEBUG = False


class Node(object):

    def __init__(self, St, list_Smax):
        self.St = St
        self.list_Smax = list_Smax
        self.list_put_eu = []
        self.list_put_am = []

    def display(self):
        print(self.St)
        print(self.list_Smax)
        if self.list_put_eu:
            print(self.list_put_eu)

    def compute_put(self, up_child, down_child):
        for s in self.list_Smax:
            u_put_eu, d_put_eu = -1, -1
            backup_u_put_eu = 0 # for case 2: put value for Smax == Stu

            u_put_am, d_put_am = -1, -1
            backup_u_put_am = 0 # for case 2: put value for Smax == Stu
            # up
            for i, uS_max in enumerate(up_child.list_Smax):
                #print("uSt - uS_max = {} * {} - {}: ".format(u, St, uS_max) + str(abs(u * St - uS_max)))
                if (abs(u * self.St - uS_max)) < epsilon:
                    backup_u_put_eu = up_child.list_put_eu[i]
                    backup_u_put_am = up_child.list_put_am[i]
                    #print("uSt, uS_max, backup = {} {} {}".format(u*St, uS_max, backup_u_put))

                if (abs(s - uS_max)) < epsilon:
                    u_put_eu = up_child.list_put_eu[i]
                    u_put_am = up_child.list_put_am[i]
                
            if u_put_eu == -1: # can't find S_max == uS_max
                u_put_eu = backup_u_put_eu
                u_put_am = backup_u_put_am

            # down
            for i, dS_max in enumerate(down_child.list_Smax):
                if abs(s - dS_max) < epsilon:
                    d_put_eu = down_child.list_put_eu[i]
                    d_put_am = down_child.list_put_am[i]
            
            put_eu = ( u_put_eu * p + d_put_eu * (1 - p) ) * np.exp(-r * dT)
            put_am = ( u_put_am * p + d_put_am * (1 - p) ) * np.exp(-r * dT)
            put_am = max(s - self.St, put_am)
            #print("u_put * p + d_put * (1 - p) * e^(-r dT) = {} * {} + {} * {} * e^(-r dT) = {}".format(u_put, p, d_put, (1 - p), put))
            self.list_put_eu.append(put_eu)
            self.list_put_am.append(put_am)
            #input("stop")


def print_line():
    print("----------------------------------------")


in_f = open(sys.argv[1], "r")
line = in_f.readline()
line = line.strip().split(" ")
St, r, q, sigma, t, T, S_max_t, n, num_simulations, num_repetitions = [float(ele) for ele in line]
n = int(n)
num_simulations = int(num_simulations)
num_repetitions = int(num_repetitions)
print("St          = " + str(St))
print("r           = " + str(r))
print("q           = " + str(q))
print("sigma       = " + str(sigma))
print("t           = " + str(t))
print("T           = " + str(T))
print("S_max_t     = " + str(S_max_t))
print("n           = " + str(n))
print("simulations = " + str(num_simulations))
print("repetitions = " + str(num_repetitions))


dT = (T - t) / n
u = np.exp(sigma * np.sqrt(dT))
d = 1/u
p = (np.exp((r-q) * dT) - d)/(u - d)

print("u           = " + str(u))
print("d           = " + str(d))
print("p           = " + str(p))
print("dT          = " + str(dT))

print_line()
print(" Bonus 1")

u_t = 0  # S_max_t / St for each node

tree = []

root = Node(St, [S_max_t])
tree.append([root])

for i in range(1, n + 1):  # period 1, 2, ..., n
    #print("### depth {} ###".format(i))
    sub_t = []
    
    j = 0 # the most top node
    prev_St = tree[i-1][j].St
    cur_St = u * prev_St
    tmp_list_Smax = [max(cur_St, S_max_t) ]
    node_ij = Node(cur_St, tmp_list_Smax)
    if DEBUG:
        node_ij.display()
    sub_t.append(node_ij)

    for j in range(1, i): # j = 1, .., i - 1
        # find St first
        prev_St = tree[i-1][j].St
        cur_St = u * prev_St

        # fine Smax quickly
        tmp_list_Smax = []
        top_Smax = tree[i - j][0].St
        cur_Smax = top_Smax
        for k in range(j):
            tmp_list_Smax.append(cur_Smax)
            cur_Smax = d * cur_Smax
       
        cur_Smax = max(cur_St, S_max_t)
        if not (cur_Smax in tmp_list_Smax):
            tmp_list_Smax.append(cur_Smax)
        
        node_ij = Node(cur_St, tmp_list_Smax)
        if DEBUG:
            node_ij.display()
        sub_t.append(node_ij)

    j = i # the most bottom node
    prev_St = tree[i - 1][j - 1].St
    cur_St = d * prev_St
    tmp_list_Smax = [S_max_t]

    node_ij = Node(cur_St, tmp_list_Smax)
    if DEBUG:
        node_ij.display()
    sub_t.append(node_ij)

    tree.append(sub_t)

print_line()
# get put value for depth n node    
depth = n
sub_t = tree[n]
for node in sub_t:
    list_put_eu = []
    for s in node.list_Smax:
        list_put_eu.append(max(s - node.St, 0))
    list_put_am = list_put_eu[:]

    node.list_put_eu = list_put_eu
    node.list_put_am = list_put_am
    #print(list_put_eu)


for i in range(n - 1, -1, -1): # i = n, n-1, ..., 1
    #print("### depth {} """.format(i))
    pre_sub_t = tree[i + 1]
    cur_sub_t = tree[i]

    for idx, node in enumerate(cur_sub_t):
        up_child = pre_sub_t[idx]
        down_child = pre_sub_t[idx + 1]
        if DEBUG:
            print("up_child is: ")
            up_child.display()
            print("down_child is: ")
            down_child.display()
            print("self is :")
            node.display()
        node.compute_put(up_child, down_child)
        
        if DEBUG:
            print("EU")
            print(node.list_put_eu)
            print("AM")
            print(node.list_put_am)

print("### Answer for Bonus 1 ###")
print("European: {}".format(root.list_put_eu[0]))
print("American: {}".format(root.list_put_am[0]))
