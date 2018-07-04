import numpy as np
import matplotlib.pyplot as plt

linear_f = open("./linear.out", "r")
log_f = open("./log.out", "r")
list_log_option = []

def get_M_option(f):
    list_M = []
    list_option = []
    for line in f:
        line = line.strip().split(" ")
        M, option = line
        M = int(M)
        option = float(option)

        list_M.append(M)
        list_option.append(option)
    return list_M, list_option

list_M, list_linear_option = get_M_option(linear_f)
list_M, list_log_option = get_M_option(log_f)

fig = plt.figure(1)
plt.plot(list_M, list_linear_option)
plt.plot(list_M, list_log_option)
plt.xlabel("M")
plt.ylabel("Option value")
plt.title("Linear vs Log Method")
plt.savefig("bonus1.png")

"""
fig = plt.figure(2)
plt.plot(list_M, list_log_option)
plt.xlabel("M")
plt.ylabel("Option value")
plt.title("Log Method")
plt.savefig("log.png")
"""
