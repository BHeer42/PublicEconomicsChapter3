# -*- coding: utf-8 -*-
"""
Ch3_olg_dyn1.py

Created on Tue Dec  8 09:20:42 2020

Chapter 3.2 from the book

    PUBLIC ECONOMICS, by Burkhard Heer, Springer 2019
    
    see page 77

@author: Burkhard Heer
"""


import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10,6)

# parameter values
alpha = 0.36
beta = 0.40
n = 0.1
t0 = 20 # number of transition periods
# steady-state capital stock
k = beta / (1+beta) * ((1-alpha) / (1+n) )**(1/(1-alpha))
# initial capital stock
k0 = k / 3

    
def kdyn(x):
    y = beta / (1+beta) * (1-alpha) / (1+n) * x**(alpha)
    return y



kt = []   # initial value of vector
kt.append(k0)   # kt[0] = k0

# iterate over periods i=0,..t0-1
for i in range(int(t0)):
    k1 = kdyn(k0)
    kt.append(k1)
    k0=k1

print("rows in kt: " +str(len(kt)))
plt.plot(kt)
plt.show()