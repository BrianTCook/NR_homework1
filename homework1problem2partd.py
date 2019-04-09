# Problem 2, Part D

import numpy as np
import math as m #trig functions for sampling theta, phi
from random import sample

random_numbers = np.loadtxt('homework1problem1partb_randomnumbers.txt')
a, b, c = np.loadtxt('abc.txt')
Aval = np.loadtxt('Aval.txt')

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
	
	return [r, theta, phi]

def creates_satellites(N_satellites):
    
	return [spherical_coordinates() for i in range(N_satellites)]

satellites = creates_satellites(100)

i = 1
with open('homework1problem2partd_print.txt', 'a') as f:
	for satellite in satellites:
		print('Satellite %i: (r,theta,phi) = (%.03f, %.03f, %.03f)'%(i, satellite[0], satellite[1], satellite[2]))
		f.write('Satellite %i: (r,theta,phi) = (%.03f, %.03f, %.03f)\n'%(i, satellite[0], satellite[1], satellite[2]))
		i += 1
