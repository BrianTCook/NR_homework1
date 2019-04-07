# Problem 2, Part B

import numpy as np
import matplotlib.pyplot as plt

random_numbers = np.loadtxt('homework1problem1partb_randomnumbers.txt')
a, b, c = np.loadtxt('abc.txt')

def n_withoutconsts(x):
    return (x/b)**(a-3) * np.exp(-(x/b)**c)

def n_withoutconsts_primeprime(x):
    
    #analytic expression of n''(x), needed for natural cubic spline
    
    numerator = b**3 * (x/b)**a * np.exp(-(x/b)**c) * (a**2 - 2*a*c*(x/b)**c - 7*a - \
                        c**2 * (x/b)**c + c**2 * (x/b)**(2*c) + 7*c*(x/b)**c + 12)
    denominator = x**5
    
    return numerator/denominator

X = [1e-4, 1e-2, 1e-1, 1, 5]
Y = [n_withoutconsts(x) for x in X]
YPrimePrime = [n_withoutconsts_primeprime(x) for x in X]

plt.figure()
plt.loglog(X, Y, 'ko', label = 'Data')
plt.annotate(r'$A = \langle N_{sat} \rangle = 1$', xy = (0.05, 0.5), xycoords = 'axes fraction', fontsize = 12)

#linear interpolation using Lagrange formula

for i in range(1, len(X)):
    x1, x2 = X[i-1], X[i]
    y1, y2 = Y[i-1], Y[i]
    
    XLinear = np.linspace(x1, x2, 10)
    YLinear = [(xlin/x1)**(np.log10(y2/y1)/np.log10(x2/x1)) * y1 for xlin in XLinear]
    
    if i == 1:
        plt.loglog(XLinear, YLinear, 'r', label='linear interpolation')
    else:
        plt.loglog(XLinear, YLinear, 'r')

#natural cubic spline

def spline(Y, YPrimePrime, xvals, j, x): #evaluated at xvals[j]

    def ASpl(xvals, j, x):
        return (xvals[j+1] - x)/(xvals[j+1] - xvals[j])
    def BSpl(xvals, j, x):
        return 1 - ASpl(xvals, j, x)
    def CSpl(xvals, j, x):
        return (1/6.)*(ASpl(xvals, j, x)**3 - ASpl(xvals, j, x))*(xvals[j+1]-xvals[j])**2
    def DSpl(xvals, j, x):
        return (1/6.)*(BSpl(xvals, j, x)**3 - BSpl(xvals, j, x))*(xvals[j+1]-xvals[j])**2
       
    return ASpl(xvals, j, x)*Y[j] + BSpl(xvals, j, x)*Y[j+1] + \
           CSpl(xvals, j, x)*YPrimePrime[j] + DSpl(xvals, j, x)*YPrimePrime[j+1]
    
XSpline = np.linspace(1e-4, 5, 20) #too many points overemphasizes parts of curve with large second derivative
YSpline = []

X_log10, XSpline_log10 = np.log10(X), np.log10(XSpline)
Y_log10, YPrimePrime_log10 = np.log10(Y), np.log10(YPrimePrime)

for x in XSpline_log10:
    for j in range(len(X_log10) - 1):
        if X_log10[j] <= x and X_log10[j+1] >= x:
            YSpline.append(spline(Y_log10, YPrimePrime_log10, X_log10, j, x))

YSpline = [10**ys for ys in YSpline]

plt.loglog(XSpline, YSpline,'g', label='(natural) cubic spline')
            
plt.legend(loc='best')
plt.xlabel(r'$x$', fontsize = 16)
plt.ylabel(r'$n(x)$', fontsize = 16)
plt.tight_layout()
plt.savefig('homework1problem2partbfigure1.pdf')

print('Done creating the interpolation figure!')
