# Problem 2, Part A

import numpy as np

random_numbers = np.loadtxt('homework1problem1partb_randomnumbers.txt')

a = 1.4*random_numbers[9] + 1
b = 1.5*random_numbers[61] + 0.5
c = 2.5*random_numbers[18] + 1.5

abc = [a,b,c]
np.savetxt('abc.txt', abc)

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

N_sample = 100 #number of sample points for integration scheme

def A_finder(xmin, xmax):

    denominator = simpson(integrand, xmin, xmax, N_sample)

    return 1/denominator

xmin_PartA, xmax_PartA = 1e-16, 5.
Aval = A_finder(xmin_PartA, xmax_PartA)
np.savetxt('Aval.txt', [Aval] )
print_str = 'The parameters and associated normalization constant are a = %.03f, b = %.03f, c = %.03f, and A = %.03e.'%(a,b,c,Aval)
print_strs = []
print(print_str)
print_strs.append(print_str)

with open('homework1problem2parta_print.txt', 'a') as f:
    for print_str in print_strs:
        f.write("%s \n" % print_str)
