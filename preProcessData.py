
# This file does some preprocessing to the recorded data files. Please run this file after running saveData.py


import numpy as np 

from numpy import genfromtxt

photoDiodeData = genfromtxt('positionPhotoDiodes_4thrun.txt', delimiter=',')

photoDiodeData = np.delete(photoDiodeData , 0 , 1)

photoDiodeData = np.array(photoDiodeData)

zer = np.zeros((1 , np.shape(photoDiodeData)[0]));

photoDiodeData = np.concatenate((photoDiodeData,zer.T), axis = 1)

np.savetxt("photoDiodeData_4thrun.txt", photoDiodeData, delimiter="," , newline = '\n')


print photoDiodeData
print np.shape(photoDiodeData)