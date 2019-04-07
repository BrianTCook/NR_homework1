# Problem 2, Part F

import numpy as np

random_numbers = np.loadtxt('homework1problem1partb_randomnumbers.txt')
a, b, c = np.loadtxt('abc.txt')

eps = 1e-3

def RootFinder(x_init, y, yprime):

	deltax = y(x_init)/yprime(x_init) #initialize separation between current and previous root estimate
	x = x_init #initialize x, will converge to an accuracy of eps    
	
	while abs(deltax) > eps:
        
        	deltax = (y(x)/yprime(x)) # deltax = x_{i} - x_{i-1}
        	x -= deltax
        
	return x

def n_withoutconsts(x):
	return (x/b)**(a-3) * np.exp(-(x/b)**c)

def n_withoutconsts_prime(x):
    
	numerator = b**3 * (x/b)**a * np.exp(-(x/b)**c) * (a - c*(x/b)**c - 3)
	denominator = x**4
	return numerator / denominator

def n_withoutconsts_primeprime(x):
    
	#analytic expression of n''(x), needed for natural cubic spline
    
	numerator = b**3 * (x/b)**a * np.exp(-(x/b)**c) * (a**2 - 2*a*c*(x/b)**c - 7*a - \
                        c**2 * (x/b)**c + c**2 * (x/b)**(2*c) + 7*c*(x/b)**c + 12)
	denominator = x**5
    
	return numerator/denominator

def N(x):
	return n_withoutconsts(x) * 4 * np.pi * x**2

def N_prime(x):
    
	return 4 * np.pi * (2*x * n_withoutconsts(x) + x**2 * n_withoutconsts_prime(x))

def N_primeprime(x):
    
	return 4 * np.pi * (2*(x * n_withoutconsts_prime(x) + n_withoutconsts(x)) \
                        + (x**2 * n_withoutconsts_primeprime(x) + 2 * x * n_withoutconsts_prime(x)))

extrema_guesses = [0.2, 1.8]
Max = [RootFinder(x, N_prime, N_primeprime) for x in extrema_guesses]

def func_to_find_roots(x):
	return N(x) - 0.5 * N(Max[0])

def func_to_find_roots_prime(x):
	return N_prime(x)

XStart = [ex + 1e-3 for ex in extrema_guesses]
Roots = [RootFinder(x, func_to_find_roots, func_to_find_roots_prime) for x in XStart]
CleanedRoots = [root for root in Roots if str(root) != 'nan']

print_strs = []

print_str_init = 'The roots of N(x) = N_{max}/2 for the parameters a, b, c are:'
print_strs.append(print_str_init)

for root in CleanedRoots:
	print_strs.append(str(root))

for string in print_strs:
	print(string)

with open('homework1problem2partf_print.txt', 'a') as f:
	for print_str in print_strs:
		f.write("%s \n" % print_str)
