import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('path dependence.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    parameter = []
    V = []
    t = []
    
    # loop over the rows
    for row in readCSV:
        parameter.append(row[0])
        V.append(row[1])
        t.append(row[2])
        
    # map to float array
    parameter = list(map(float, parameter))
    V = list(map(float, V))  
    t = list(map(float, t))   


# normalise
norm_V = []
norm_t = []
for i in range(len(V)):
    norm_V.append(V[i]/max(V))
    norm_t.append(t[i]/max(t))


# graph
plt.figure()
plt.plot(parameter, norm_V, label = 'V(S,t=0)')
plt.plot(parameter, norm_t, label = 'Computation time')
plt.xlabel('K')
plt.ylabel('Normalised scale')
plt.legend()
plt.title('Floating strike Asian call with N=200,000')