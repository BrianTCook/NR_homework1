# Problem 2, Part E

import numpy as np
import matplotlib.pyplot as plt
import math as m
from random import sample

random_numbers = np.loadtxt('homework1problem1partb_randomnumbers.txt')
a, b, c = np.loadtxt('abc.txt')
Aval = np.loadtxt('Aval.txt')

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

haloes = [creates_satellites(100) for i in range(1000)]

all_xs = []
for halo in haloes:
    for satellite in halo:
        r, theta, phi = satellite
        all_xs.append(r)
        
def N(x):
    return Aval * n_withoutconsts(x) * 4 * np.pi * x**2

N_bins = 20

#N(x) distribution
X = np.logspace(-4, np.log10(5), 100)
Y =[N(x) for x in X]

#Distribution of satellite radii

X_bins = np.logspace(-4,np.log10(5),N_bins+1)
Y_bins = [0 for i in range(N_bins)]
Ys = [[] for i in range(N_bins)]

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

X_bins_toshow = X_bins[0:len(Y_bins)]
len_all_xs = len(all_xs)
Y_bins_toshow = [8*Y/len_all_xs for Y in Y_bins]

plt.plot(X, Y, 'k', label=r'$N(x)$')
plt.semilogx(X_bins_toshow, Y_bins_toshow, 'r', label='Satellite Distribution*')
plt.xlabel(r'$x$', fontsize=14)
plt.legend(loc='best')
plt.savefig('homework1problem2partefigure1.pdf')

print('Done comparing satellite distribution to N(x)!')
