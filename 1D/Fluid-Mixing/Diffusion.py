import numpy as np
import matplotlib.pyplot as plt
from math import *

def compute_AMI (data):
    sigma = np.std(data)
    average = np.average(data)
    return (sigma/average)

# 100 linearly spaced numbers
resolution = 1000
y_axis = np.linspace(-0.5,0.5, resolution)
concentration = np.zeros(resolution)
concentration += 0.5

length = 0.000860 # meters
width = 0.000060 # meters
speed = 0.01 # m/s
diffusion = 1e-9# m^2/s
epsilon = 1e-8
coeff = 1
n = 1

x = length/width # Dimensionless x-location

Pe = (speed*width)/diffusion

while coeff > epsilon:
    coeff = exp(-(x*pi*pi*(2*n-1)**2)/Pe)
    product = ((1-cos(pi*(2*n-1)))/(2*n-1))
    for y in range(y_axis.size):
        y_loc = y_axis[y]
        s = sin(pi*(2*n-1)*y_loc)
        concentration[y] += (s*product*coeff)/pi
    n+=1

concentration = 1-concentration

AMI = compute_AMI(concentration)

print(f'The AMI of the 1D approach is: {AMI}')

np.savetxt("1D-Mixing.csv", concentration, delimiter=",")

plt.plot(y_axis,concentration)
plt.show()
