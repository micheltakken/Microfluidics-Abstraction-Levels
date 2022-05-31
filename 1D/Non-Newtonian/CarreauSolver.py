import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative
from math import *
from scipy.integrate import simps

# 100 linearly spaced numbers
resolution = 500
z_axis = np.linspace(0, 0.00005, resolution)
velocity = np.zeros(resolution)
gamma = np.zeros(resolution)
tau = np.zeros(resolution)

mu0 = 1.56e-5
mui = 1.56e-6
a = 0.644
n = 0.392
labda = 0.11
dp = 0.0168
Length = 0.0005 # meter
coeff = dp/Length

def my_newton(f, df, x0, tol):
    # output is an estimation of the root of f
    # using the Newton Raphson method
    # recursive implementation
    if abs(f(x0)) < tol:
        return x0
    else:
        return my_newton(f, df, x0 - f(x0)/df(x0), tol)

def deriv(x):
    return derivative(f,x)

n_prime = (n-1)/a
mu_prime = mu0-mui
l_prime = labda**a
f_prime = lambda x: mui + mu_prime*(1+l_prime*(x**a))**n_prime + mu_prime*(n-1)*(l_prime*x**a)*((1+l_prime*(x**a))**(n_prime-1))

for z in range(z_axis.size - 1):
    z_loc = z_axis[z+1]
    tau[z+1] = (z_loc*coeff)
    f = lambda x: x*(mui+(mu_prime)*(1+l_prime*(x**a))**(n_prime)) - (z_loc*coeff)
    gamma[z+1] = my_newton(f,f_prime, 0.002, 1e-12)
    velocity[z+1] = simps(gamma,z_axis)

gamma = np.flip(gamma)

for z in range(z_axis.size):
    velocity[-z-1] = simps(gamma[1:z+2],z_axis[1:z+2])

print(f'Umax = {simps(gamma,z_axis)}')

total = np.concatenate((np.flip(velocity),velocity),axis=None)

np.savetxt("1D-NonNewtonian.csv", total, delimiter=",")
plt.plot(total,'b')
plt.show()
