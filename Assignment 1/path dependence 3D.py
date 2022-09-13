import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


# Read csv file1
with open('path dependence.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    N = []
    K = []
    V = []
    t = []
    
    # loop over the rows
    for row in readCSV:
        N.append(row[0])
        K.append(row[1])
        V.append(row[2])
        t.append(row[3])
        
    # map to float array
    N = list(map(float, N))
    K = list(map(float, K))  
    V = list(map(float, V))        
    t = list(map(float, t))  

# normalise
V_norm = []
t_norm = []
for i in range(len(V)):
    V_norm.append(V[i]/max(V))
    t_norm.append(t[i]/max(t))
    
print(V)
print(V[0]/V[-1])
print(V_norm)


# graph
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(N, K, V)
ax.set_xlabel('N')
ax.set_ylabel('K')
ax.set_zlabel('V(S,t=0)')

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(N, K, t)
ax.set_xlabel('N')
ax.set_ylabel('K')
ax.set_zlabel('t')

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(N, K, V_norm, 'b', label = 'V(S,t=0)')
ax.scatter(N, K, t_norm, 'r', label = 'time')
ax.set_xlabel('N')
ax.set_ylabel('K')
plt.legend()













