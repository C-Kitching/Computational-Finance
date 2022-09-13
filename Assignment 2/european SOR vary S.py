import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('european SOR vary S.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    S = []
    V = []
    #ana = []
    
    # loop over the rows
    for row in readCSV:
        S.append(row[0])
        V.append(row[1])
        #ana.append(row[2])
        
    # map to float array
    S = list(map(float, S))
    V = list(map(float, V))
    #ana = list(map(float, ana))      

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

# terminal option
V_T = [max(F, R*S[i]) for i in range(len(S))]
plt.plot(S, V_T, label = 'max(F,RS)', linestyle = '--')

# AS+B
ASB = [R*S[i]-2*C/(alpha +r)*np.exp(-alpha*T) for i in range(len(S))]
plt.plot(S, ASB, label = 'V(S,T)', linestyle = '--')

# analytic results
#plt.plot(S, ana, label = 'Analytic', color = 'lime')

# graph details
plt.xlabel('S')
plt.ylim(0, max(V) + 100)
plt.ylabel('Value')
plt.title('Option price varying with S\nThomas solver')
plt.legend()



















