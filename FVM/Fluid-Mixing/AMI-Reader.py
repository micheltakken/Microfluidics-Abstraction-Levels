import csv
import numpy as np

def compute_AMI (data):
    sigma = np.std(data)
    average = np.average(data)
    return (sigma/average)


Data_2D = np.loadtxt('2D/Mixing2DFVM.csv', delimiter=',', skiprows=1, usecols = 3)
Data_3D = np.loadtxt('3D/Mixing3DFVM.csv', delimiter=',', skiprows=1, usecols = 0)

AMI_2D = compute_AMI(Data_2D)
AMI_3D = compute_AMI(Data_3D)

print(f'The AMI of the 2D data is {AMI_2D}')
print(f'The AMI of the 3D data is {AMI_3D}')

