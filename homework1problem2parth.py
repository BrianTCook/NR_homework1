# Problem 2, Part H

import numpy as np

a_len, b_len, c_len = 15, 16, 26

a_range = [1.1 + 0.1 * i for i in range(a_len)]
b_range = [0.5 + 0.1 * i for i in range(b_len)]
c_range = [1.5 + 0.1 * i for i in range(c_len)]

xmin_PartA, xmax_PartA = 1e-16, 5
N_sample = 100

Avals = []

#calculates A = A(a,b,c)

def f(a,b,c):

	def n_withoutconsts(x):
		return (x/b)**(a-3) * np.exp(-(x/b)**c)

	def integrand(x):
		return n_withoutconsts(x) * 4*np.pi * x**2

	def simpson(f, x_init, x_final, N_simp): #Simpson's integration rule, \int_{a}^{b} f(x) dx with N sample points

		h = (x_final-x_init)/N_simp

		I = f(a) + f(b)

		odds = [4*f(a + k*h) for k in range(1,N_simp,2)]
		evens = [2*f(a + k*h) for k in range(2,N_simp,2)]
		I += sum(odds) + sum(evens)

		I *= h/3.

		return I

	def A_finder(xmin, xmax):

		denominator = simpson(integrand, xmin, xmax, N_sample)

		return 1/denominator

	return A_finder(1e-16,5)

AValsTable = [[a,b,c,f(a,b,c)] for a in a_range for b in b_range for c in c_range]

with open('homework1problem2parth_table.txt', 'a') as f:
	for item in AValsTable:
		item = tuple(item)
		f.write("%.05f %.05f %.05f %.05f\n" % item)
    
#trilinear interpolation (see trilinear interpolation wikipedia page, section titled "alternative algorithm")
# a0 <= a <= a1, b0 <= b <= b1, c0 <= c <= c1

def A_interpolated(a,a0,a1,b,b0,b1,c,c0,c1):

	#doesn't access the table but does the same calculation
	f000 = f(a0,b0,c0)
	f100 = f(a1,b0,c0)
	f010 = f(a0,b1,c0)
	f110 = f(a1,b1,c0)
	f001 = f(a0,b0,c1)
	f101 = f(a1,b0,c1)
	f011 = f(a0,b1,c1)
	f111 = f(a1,b1,c1)

	#interpolation polynomial coefficients

	x0 = ( (f000*a1*b1*c1)/((a0-a1)*(b0-b1)*(c1-c0))  + (f001*a1*b1*c0)/((a0-a1)*(b0-b1)*(c0-c1))  
	+ (f010*a1*b0*c1)/((a0-a1)*(b0-b1)*(c0-c1))  + (f011*a1*b0*c0)/((a0-a1)*(b0-b1)*(c1-c0))  
	+ (f100*a0*b1*c1)/((a0-a1)*(b0-b1)*(c0-c1))  + (f101*a0*b1*c0)/((a0-a1)*(b0-b1)*(c1-c0))  
	+ (f110*a0*b0*c1)/((a0-a1)*(b0-b1)*(c1-c0))  + (f111*a0*b0*c0)/((a0-a1)*(b0-b1)*(c0-c1)) )

	x1 = ( (f000*b1*c1)/((a0-a1)*(b0-b1)*(c0-c1))  + (f001*b1*c0)/((a0-a1)*(b0-b1)*(c1-c0))  
	+ (f010*b0*c1)/((a0-a1)*(b0-b1)*(c1-c0))  + (f011*b0*c0)/((a0-a1)*(b0-b1)*(c0-c1))  
	+ (f100*b1*c1)/((a0-a1)*(b0-b1)*(c1-c0))  + (f101*b1*c0)/((a0-a1)*(b0-b1)*(c0-c1))  
	+ (f110*b0*c1)/((a0-a1)*(b0-b1)*(c0-c1))  + (f111*b0*c0)/((a0-a1)*(b0-b1)*(c1-c0)) )

	x2 = ( (f000*a1*c1)/((a0-a1)*(b0-b1)*(c0-c1))  + (f001*a1*c0)/((a0-a1)*(b0-b1)*(c1-c0))  
	+ (f010*a1*c1)/((a0-a1)*(b0-b1)*(c1-c0))  + (f011*a1*c0)/((a0-a1)*(b0-b1)*(c0-c1))  
	+ (f100*a0*c1)/((a0-a1)*(b0-b1)*(c1-c0))  + (f101*a0*c0)/((a0-a1)*(b0-b1)*(c0-c1))  
	+ (f110*a0*c1)/((a0-a1)*(b0-b1)*(c0-c1))  + (f111*a0*c0)/((a0-a1)*(b0-b1)*(c1-c0)) )

	x3 = ( (f000*a1*b1)/((a0-a1)*(b0-b1)*(c0-c1))  + (f001*a1*b1)/((a0-a1)*(b0-b1)*(c1-c0))  
	+ (f010*a1*b0)/((a0-a1)*(b0-b1)*(c1-c0))  + (f011*a1*b0)/((a0-a1)*(b0-b1)*(c0-c1))  
	+ (f100*a0*b1)/((a0-a1)*(b0-b1)*(c1-c0))  + (f101*a0*b1)/((a0-a1)*(b0-b1)*(c0-c1))  
	+ (f110*a0*b0)/((a0-a1)*(b0-b1)*(c0-c1))  + (f111*a0*b0)/((a0-a1)*(b0-b1)*(c1-c0)) )

	x4 = ( (f000*c1)/((a0-a1)*(b0-b1)*(c1-c0))  + (f001*c0)/((a0-a1)*(b0-b1)*(c0-c1))  
	+ (f010*c1)/((a0-a1)*(b0-b1)*(c0-c1))  + (f011*c0)/((a0-a1)*(b0-b1)*(c1-c0))  
	+ (f100*c1)/((a0-a1)*(b0-b1)*(c0-c1))  + (f101*c0)/((a0-a1)*(b0-b1)*(c1-c0))  
	+ (f110*c1)/((a0-a1)*(b0-b1)*(c1-c0))  + (f111*c0)/((a0-a1)*(b0-b1)*(c0-c1)) )

	x5 = ( (f000*b1)/((a0-a1)*(b0-b1)*(c1-c0))  + (f001*b1)/((a0-a1)*(b0-b1)*(c0-c1))  
	+ (f010*b0)/((a0-a1)*(b0-b1)*(c0-c1))  + (f011*b0)/((a0-a1)*(b0-b1)*(c1-c0))  
	+ (f100*b1)/((a0-a1)*(b0-b1)*(c0-c1))  + (f101*b1)/((a0-a1)*(b0-b1)*(c1-c0))  
	+ (f110*b0)/((a0-a1)*(b0-b1)*(c1-c0))  + (f111*b0)/((a0-a1)*(b0-b1)*(c0-c1)) )

	x6 = ( (f000*a1)/((a0-a1)*(b0-b1)*(c1-c0))  + (f001*a1)/((a0-a1)*(b0-b1)*(c0-c1))  
	+ (f010*a1)/((a0-a1)*(b0-b1)*(c0-c1))  + (f011*a1)/((a0-a1)*(b0-b1)*(c1-c0))  
	+ (f100*a0)/((a0-a1)*(b0-b1)*(c0-c1))  + (f101*a0)/((a0-a1)*(b0-b1)*(c1-c0))  
	+ (f110*a0)/((a0-a1)*(b0-b1)*(c1-c0))  + (f111*a0)/((a0-a1)*(b0-b1)*(c0-c1)) )

	x7 = ( (f000)/((a0-a1)*(b0-b1)*(c0-c1))  + (f001)/((a0-a1)*(b0-b1)*(c1-c0))  
	+ (f010)/((a0-a1)*(b0-b1)*(c1-c0))  + (f011)/((a0-a1)*(b0-b1)*(c0-c1))  
	+ (f100)/((a0-a1)*(b0-b1)*(c1-c0))  + (f101)/((a0-a1)*(b0-b1)*(c0-c1))  
	+ (f110)/((a0-a1)*(b0-b1)*(c0-c1))  + (f111)/((a0-a1)*(b0-b1)*(c1-c0)) )

	return x0 + x1*a + x2*b + x3*c + x4*a*b + x5*a*c + x6*b*c + x7*a*b*c

print('Table of A(a,b,c) values is computed and trilinear interpolation scheme is written!')
