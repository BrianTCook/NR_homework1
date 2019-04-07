# Problem 2, Part G

import numpy as np
import matplotlib.pyplot as plt
import math as m
from random import sample

random_numbers = np.loadtxt('homework1problem1partb_randomnumbers.txt')
a, b, c = np.loadtxt('abc.txt')
Aval = np.loadtxt('Aval.txt')


def mergesort(lst, N_lst):
    
    flag = 0
    while flag == 0:
    
        if len(lst) >= 2:

            if len(lst)%2 == 0:
                midpoint = len(lst)/2
            if len(lst)%2 != 0:
                midpoint = (len(lst) + 1)/2

            midpoint = int(midpoint) #indices must be integers
                
            #step 1, NUR5.pdf pseudocode
            left = lst[:midpoint]
            right = lst[midpoint:]
            
            #step 2
            N_left, N_right = len(left), len(right)
            left = mergesort(left, N_left)
            right = mergesort(right, N_right)

            i = 0 #left index
            j = 0 #right index
            k = 0 #lst index

            #steps 3 and 4
            while i < len(left) and j < len(right):
                if left[i] <= right[j]:
                    lst[k] = left[i]
                    i += 1 #go to next element in left half
                else:
                    lst[k] = right[j]
                    j += 1 #go to next element in right half
                k += 1 #go to next element in list that is the result of merging the two halves


            #in cases of odd lst one of the halves will be iterated over before, in which case we simply need to 
            #include the remaining elements from the longer of the two halves
            while i < len(left):
                lst[k] = left[i]
                i += 1
                k += 1
            while j < len(right):
                lst[k] = right[j]
                j += 1
                k += 1
                
        if len(lst) == N_lst:
        
            flag = 1
            
    return lst
    
'''
from part D
'''

def n_withoutconsts(x):
	return (x/b)**(a-3) * np.exp(-(x/b)**c)

def pxdx(x):
	return Aval * n_withoutconsts(x) * 4 * np.pi * x**2

x_finding = np.linspace(1e-8, 5, 10000)
pxdx_vals = [pxdx(x) for x in x_finding]
max_pxdx = max(pxdx_vals)

random_numbers = list(random_numbers)

def radius_finder(): #method of sampling list of random numbers

	flag = 0
	while flag == 0:

		#rejection sampling
		randnums = sample(random_numbers, 2)

		x = 5*randnums[0] #xmax = 5
		y = max_pxdx*randnums[1] #normalization computed in 1A

		#print('y = %.03f and pxdx(x) = %.03f'%(y, pxdx(x)))

		if y < pxdx(x):
			r = x
			flag = 1

	return r
  
def spherical_coordinates():

	r = radius_finder()

	randnums = sample(random_numbers, 2)
	
	theta = m.acos(1-2*randnums[0])
	phi = 2*np.pi*randnums[1]
	
	return [r, phi, theta]

def creates_satellites(N_satellites):
    
	return [spherical_coordinates() for i in range(N_satellites)]

'''
couldn't think of a way to store the haloes from part e so doing the same calculation twice
'''

haloes = [creates_satellites(100) for i in range(1000)]

all_xs = []
for halo in haloes:
    for satellite in halo:
        r, theta, phi = satellite
        all_xs.append(r)

N_bins_2e = 20

X_bins = np.logspace(-4,np.log10(5),N_bins_2e+1)
Y_bins = [0 for i in range(N_bins_2e)]
Ys = [[] for i in range(N_bins_2e)]

for x in all_xs:
    i = 0
    while i < len(X_bins)-1:
        xi, xf = X_bins[i], X_bins[i+1]

        if x >= xi and x <= xf:
            Y_bins[i] += 1
            Y_bin = Ys[i]
            Y_bin.append(x)
            break
        else:
            i += 1
    
    
max_index = Y_bins.index(max(Y_bins))

list_to_be_sorted = Ys[max_index]
sorted_list = mergesort(list_to_be_sorted, len(list_to_be_sorted))

median_index = int(0.5*len(sorted_list))
sixteen_index = int(0.16*len(sorted_list))
eightyfour_index = int(0.84*len(sorted_list))

print('median of this bin is %.03f'%(sorted_list[median_index]))
print('16th percentile of this bin is %.03f'%(sorted_list[sixteen_index]))
print('84th percentile of this bin is %.03f'%(sorted_list[eightyfour_index]))

with open('homework1problem2partg_print.txt', 'a') as f:
	f.write('median of this bin is %.03f\n'%(sorted_list[median_index]))
	f.write('16th percentile of this bin is %.03f\n'%(sorted_list[sixteen_index]))
	f.write('84th percentile of this bin is %.03f\n'%(sorted_list[eightyfour_index]))

values_to_be_histed = [0 for i in range(len(haloes))]

#only need the max_index
for i in range(len(haloes)):
    halo = haloes[i]
    xs = []
    for satellite in halo:
        r, theta, phi = satellite
        xs.append(r)
        
    for x in xs:

        xi, xf = X_bins[max_index], X_bins[max_index+1]

        if x >= xi and x <= xf:
            values_to_be_histed[i] += 1

mean_values_to_be_histed = sum(values_to_be_histed)/len(values_to_be_histed)

min_val, max_val = min(values_to_be_histed), max(values_to_be_histed)
N_bins = int(max_val - min_val)
histed = [0 for i in range(N_bins+1)]
for val in values_to_be_histed:
    ind = val - min_val
    histed[ind] += 1
    
X_plot = np.linspace(min_val, max_val, len(histed))
Y_plot = [histed[i]/sum(histed) for i in range(len(histed))]

'''
from problem 1, part A
'''

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

Y_Poisson = [Poisson(mean_values_to_be_histed, int(x)) for x in X_plot]

plt.plot(X_plot, Y_plot, label='Data')
plt.plot(X_plot, Y_Poisson, label=r'Poission, $\lambda = %.02f$'%mean_values_to_be_histed)
plt.xlabel(r'$k$', fontsize=14)
plt.ylabel(r'$P_{\lambda}(k)$', fontsize= 14)
plt.legend(loc='best')
plt.savefig('homework1problem2partgfigure1.pdf')

print('The distribution of satellites in the most populated radial bin has been compared to the Poissonian distribution!')
