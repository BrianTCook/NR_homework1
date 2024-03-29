#Problem 3, Part A

import numpy as np
import glob

haloes = glob.glob('satgals_*.txt')

#cleaning the .txt file

data_for_all_haloes = [[] for i in range(len(haloes))]

'''
getting the info from *.txt files
'''

for i in range(len(haloes)):

	halo = haloes[i]

	'''	
	data = []

	for line in open(halo):
		line_stripped = line.strip()

		if '#' not in line_stripped:

			Nsats_or_x = np.float64(line_stripped.split()[0])
			data.append(Nsats_or_x)

	this would be the way to do it if not for being really big .txt files, so the list comprehension will be unreadable
	
	not using given number of haloes just in case there is a bug
	'''

	data = [np.float64(line.strip().split()[0]) for line in open(halo) if '#' not in line.strip()]

	xs = data[1:len(data)]
	Nsats = len(data)-1

	data_for_all_haloes[i] = Nsats, xs

'''
from 2A
'''

N_sample = 500 #number of sample points for integration scheme

def simpson(f, x_init, x_final, N_simp): #Simpson's integration rule, \int_{a}^{b} f(x) dx with N sample points

    h = (x_final-x_init)/N_simp

    I = f(a) + f(b)

    odds = [4*f(a + k*h) for k in range(1,N_simp,2)]
    evens = [2*f(a + k*h) for k in range(2,N_simp,2)]
    I += sum(odds) + sum(evens)

    I *= h/3.
    
    return I

def Aval(a,b,c):

	def n_withoutconsts(x):
	    return (x/b)**(a-3) * np.exp(-(x/b)**c)

	def integrand(x):
	    return n_withoutconsts(x) * 4*np.pi * x**2

	def simpson(f, x_init, x_final, N_simp): #Simpson's integration rule, \int_{a}^{b} f(x) dx with N sample points

	    h = (x_final-x_init)/N_simp

	    I = f(a) + f(b)

	    odds = [4*f(a + k*h) for k in range(1,N_simp,2)]
	    evens = [2*f(a + k*h) for k in range(2,N_simp,2)]
	    I += np.sum(odds) + np.sum(evens)

	    I *= h/3.
	    
	    return I

	def A_finder(xmin, xmax):

	    denominator = simpson(integrand, xmin, xmax, N_sample)

	    return 1/denominator

	return A_finder(1e-12, 5)

'''
need to make N a function of both x and the parameters \mat{p} = (a,b,c)
'''

def N(x, a, b, c):
	return Aval(a,b,c) * (x/b)**3 * np.exp(-(x/b)**c) * 4 * np.pi * x**2

for i in range(len(haloes)):

	N_satellites, xs = data_for_all_haloes[i]

	def neglogL(a, b, c):

		FirstTermSummand = [np.log(N(x,a,b,c)) for x in xs]
		FirstTerm = np.sum(FirstTermSummand)

		# analytic integral is from 0 to infinity, numeral f diverges at zero and integrand is below machine precision for large enough x
		SecondTerm = N_satellites

		return -(FirstTerm + SecondTerm)


	'''
	from NUR6.pdf, this is how the DOWNHILL simplex method would be applied. Need 4 sample points to move around in 3D parameter space
	'''

	avals = [1.1, 1.6, 2.1, 2.5] 
	bvals = [0.5, 1.0, 1.5, 2.0] 
	cvals = [1.5, 2.5, 3.5, 4.0]

	f_simplex_dict = {(a,b,c):neglogL(a,b,c) for a in avals for b in bvals for c in cvals}
	f_simplex_list_sorted = sorted(f_simplex_dict, key = f_simplex_dict.get)
	
	'''
	using a dictionary would allow for manipulation of the values of -ln L(p)
	but the utility of dictionaries is in no small part because the order doesn't matter
	'''

	N_simplex = len(f_simplex_dict)

	fN, f0 = f_simplex_dict_sorted[N_simplex-1], f_simplex_list_sorted[0]
	fracerror = abs(fN - f0)/(abs(fN + f0)/2)
	
	delta = 1e-5

	while fracerror > delta:
		
		a_centroid = sum(avals)/len(avals)
		b_centroid = sum(bvals)/len(bvals)
		c_centroid = sum(cvals)/len(cvals)
				
		atry = 2*a_cenroid - aN
		btry = 2*b_centroid - bN
		ctry = 2*c_centroid - cN

		ftry = neglogL(atry, btry, ctry)

		if f0 <= ftry and ftry < fN:

			aNminusone, bNminusone, cNminusone = atry, btry, ctry
	
		if ftry < f0:

			aexp, bexp, cexp = 2*atry - a_centroid, 2*btry - b_centroid, 2*ctry - c_centroid

			fexp = neglogL(aexp, bexp, cexp)

			if fexp < ftry:

				aNminusone, bNminusone, cNminusone = aexp, bexp, cexp	

			else:

				aNminusone, bNminusone, cNminusone = atry, btry, ctry

		atry, btry, ctry = 0.5*(a_centroid + aN), 0.5*(b_centroid + bN), 0.5*(c_centroid + cN)

		fNminusone = neglogL(aNminusone, bNminusone, cNminusone)

		if ftry < fNminusone:
	
			atry, btry, ctry = aNminusone, bNminusone, cNminusone			
			
		for i in range(1, len(f_simplex_list) - 2):
			
			ai, bi, ci = 0.5*(a_centroid + ai), 0.5*(b_centroid + bi), 0.5*(c_centroid + ci)

		'''	
		would then replace the least accurate key:value pair in the dictionary with the most accurate one
		'''

