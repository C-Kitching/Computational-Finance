import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('American penalty vary S.csv') as csvfile:
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
kappa = 0.125

# terminal option
V_T = [max(F, R*S[i]) for i in range(len(S))]
plt.plot(S, V_T, label = 'max(F,RS)', linestyle = '--')

#optimal decision
opt_V = [V[i] for i in range(len(V)) if V[i] < V_T[i]]
opt_S = [S[i] for i in range(len(S)) if V[i] < V_T[i]]

plt.scatter(opt_S, opt_V, color = 'lime', s = 20, marker = 'x')


# graph details
plt.xlabel('S')
plt.ylabel('Value')
plt.title('American option price')
plt.legend()












