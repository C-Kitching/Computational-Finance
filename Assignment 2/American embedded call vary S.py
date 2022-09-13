import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('American embedded call vary S.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    S = []
    V = []
    
    # loop over the rows
    for row in readCSV:
        S.append(row[0])
        V.append(row[1])
        
    # map to float array
    S = list(map(float, S))
    V = list(map(float, V))     

# numerical values    
plt.figure()
plt.plot(S,V, color = 'red', label = "V(S,t=0)")

# paramters
R = 1
F = 50
C = 0.285
alpha = 0.01
r = 0.0114
T = 2
X = 50.5
Cp = 60
kappa = 0.125

# terminal option
V_T = [max(F, R*S[i]) for i in range(len(S))]
plt.plot(S, V_T, label = 'max(F,RS)', linestyle = '--')

V_T2 = [max(Cp, R*S[i]) for i in range(len(S))]
plt.plot(S, V_T2, label = 'max($C_{P}$,RS)', linestyle = '--')

#optimal decision
for i in range(len(V)):
    point_V = []
    point_S = []    
    if(V[i] == Cp or V[i] == R*S[i] or V[i] == F):
      point_V.append(V[i])
      point_S.append(S[i])
      break
plt.scatter(point_V, point_S, label = 'Optimal decision point', color = 'lime', s= 50)
plt.scatter(16, 50, color = 'lime', s= 50)

# graph details
plt.xlabel('S')
plt.ylim(0,150)
plt.xlim(0,150)
plt.ylabel('Value')
plt.title('American sytle convertible bond with embedded call option')
plt.legend()















