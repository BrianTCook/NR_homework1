from __future__ import division
import numpy as np
import math as m
import matplotlib.pyplot as plt

# Problem 1, Part A

#simple implementation, overflow error for the last two k values
def factorial(x):
    if x == 0:
        return 1
    else:
        return x*factorial(x-1)
    
def Poisson(L, k):
    
    numerator = np.float64(L**k * np.exp(-L))
    
    if factorial(k) < 2**63 - 1:
        
        denominator = np.int64(factorial(k))
    
    else:
        
        denominator = np.float64(factorial(k)) #can't store 40! as a 64-bit integer
        
    return numerator/denominator

Lk_vals = [[1, 0], [5, 10], [3, 21], [2.6, 40]]

print_strs = []

for Lk_pair in Lk_vals:
    L, k = Lk_pair
    print_str = 'for lambda = %.01f and k = %i, P_{lambda}(k) = %.09e'%(L, k, Poisson(L,k))
    print(print_str)
    print_strs.append(print_str)
    
with open('homework1problem1parta_print.txt', 'a') as f:
    for print_str in print_strs:
        f.write("%s \n" % print_str)
