# -*- coding: utf-8 -*-
"""
Ch3_transition.py

Created on Tue Dec  8 09:20:42 2020

Chapter 3.3.4 from the book

    PUBLIC ECONOMICS, by Burkhard Heer, Springer 2019
    
    see pages 84-88

@author: Burkhard Heer
"""

import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10,6)

# parameter values
alpha = 0.36
beta = 0.40
n = 0.1
bigt = 20 # number of transition periods


k = (alpha/n)**(1/(1-alpha))	# steady state value of capital stock 
k0 = 0.7 * k		# initial capital stock 
kfinal = 0.7 * k    # final capital stock
kt = np.zeros(bigt+1)   # time path for capital stock
kt[0] = k0

# function that computes the transition given 
# initial capital stock k0 and guess k1
# input: x, value of capital stock in period 1
# output: k[bigt]-kfinal
def kdyn(x):
    kt[1] = x
    for i in range(1,bigt):
        temp = (1+alpha*kt[i]**(alpha-1)) / (1+n)**2 * (kt[i-1]+kt[i-1]**alpha - (1+n)*kt[i])
        kt[i+1] = (kt[i]+kt[i]**alpha) / (1+n) -temp;
    return kt[bigt] - kfinal

# tests if kdyn is computing the transition for our guess
kguess = 5.35
ksolution = kdyn(kguess)
print("test of initial guess in kdyn: " + str(ksolution))
print("rows in kt: " +str(len(kt)))


x = scipy.optimize.fsolve(kdyn, kguess)
print("solution for kt[1]: " + str(x))

plt.plot(kt)
plt.show()