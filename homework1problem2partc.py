# Problem 2, Part C

import numpy as np

random_numbers = np.loadtxt('homework1problem1partb_randomnumbers.txt')
a, b, c = np.loadtxt('abc.txt')

def n_withoutconsts(x):
    return (x/b)**(a-3) * np.exp(-(x/b)**c)

def n_withoutconsts_prime(x):
    
    numerator = b**3 * (x/b)**a * np.exp(-(x/b)**c) * (a - c*(x/b)**c - 3)
    denominator = x**4
    return numerator / denominator

eps = 1e-15 # ~ 10 * machine precision for convergence of Ridder's method

def n_withoutconsts_prime_ridders(x):

    def centraldifference(f, x, h):
        return (f(x+h) - f(x-h))/(2*h)

    h, d = 1e-1, 2. #initial step size and degree to which it will decrease iteratively

    D = [[] for i in range(20)]

    i = 0
    while i < len(D):
        D[i] = centraldifference(n_withoutconsts, x, h)
        h /= d
        i += 1

    flag, i = 0, 1
    while flag == 0:

        diff = D[i] - D[i-1]

        if diff < eps:
            ToKeep = D[i]
            flag = 1

        else:
            i += 1

        if i == len(D):
            print('Did not converge to desired accuracy! Will keep most accurate derivative')
            ToKeep = D[i]
            flag = 1

    return ToKeep

nprime_analytic = n_withoutconsts_prime(b)
nprime_numerical_ridders = n_withoutconsts_prime_ridders(b)
    
print('Analytic value of dn/dx evaluated at x = b is %.12f'%nprime_analytic)
print('Numerical value of dn/dx evaluated at x = b is %.12f'%n_withoutconsts_prime_ridders(b))
print('Fractional error is %.06e'%(abs(nprime_analytic-nprime_numerical_ridders)/abs(nprime_analytic)))

with open('homework1problem2partc_print.txt', 'a') as f:
	f.write('Analytic value of dn/dx evaluated at x = b is %.12f\n'%nprime_analytic)
	f.write('Numerical value of dn/dx evaluated at x = b is %.12f\n'%n_withoutconsts_prime_ridders(b))
	f.write('Fractional error is %.06e\n'%(abs(nprime_analytic-nprime_numerical_ridders)/abs(nprime_analytic)))
