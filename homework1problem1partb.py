from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Problem 1, Part B

a_LCG, c_LCG, m_LCG = 1664525, 1013904223, 2**32 #from Numerical Recipes textbook
a1, a2, a3, = 21, 35, 4 #good values for 64-bit XOR-shift method according to lecture

def rng(x):
    newval = (a_LCG*x + c_LCG) % m_LCG

    newval = newval ^ (newval << a1)
    newval = newval ^ (newval >> a2)
    newval = newval ^ (newval << a3)

    return newval
    
seed = 912309568757263 #needs to be an integer
print_strs = []
print_str = 'the seed for this list of random numbers is %i.'%seed
print(print_str)
print_strs.append(print_str)

with open('homework1problem1partb_print.txt', 'a') as f:
    for print_str in print_strs:
        f.write("%s \n" % print_str)

N_total = 1000000
rns = [0 for i in range(N_total)] #need to initialize rns list

rns[0] = seed
for i in range(1,N_total):
    rns[i] = rng(rns[i-1])
    
max_rns = max(rns)
random_numbers = [rns[i]/max_rns for i in range(N_total)]
    
with open('homework1problem1partb_randomnumbers.txt', 'a') as f:
    for rn in random_numbers:
        f.write("%.06f \n" % rn)

N_scatter = 1000
X_scatter = random_numbers[0:N_scatter-1]
Y_scatter = random_numbers[1:N_scatter]

cmap = cm.rainbow(np.linspace(0.0, 1.0, 1000))

plt.figure()
plt.scatter(X_scatter, Y_scatter, c=cmap, s=8)
plt.xlabel('random_numbers[0:N_scatter-1]')
plt.ylabel('random_numbers[1:N_scatter]')
plt.savefig('homework1problem1partbfigure1.pdf')

plt.figure()
plt.hist(random_numbers, bins = 20)
plt.xlabel('random number value')
plt.ylabel('N(generated numbers)')
plt.savefig('homework1problem1partbfigure2.pdf')
