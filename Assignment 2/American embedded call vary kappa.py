import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('American embedded call vary kappa.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    S = []
    V1 = []
    V2 = []
    V3 = []
    
    # loop over the rows
    for row in readCSV:
        S.append(row[0])
        V1.append(row[1])
        V2.append(row[2])
        V3.append(row[3])
        
    # map to float array
    S = list(map(float, S))
    V1 = list(map(float, V1))     
    V2 = list(map(float, V2))
    V3 = list(map(float, V3))  

# numerical values    
plt.figure()
plt.plot(S,V1, label = "$\\beta$ = 0.0625")
plt.plot(S,V2, label = "$\\beta$ = 0.125")
plt.plot(S,V3, label = "$\\beta$ = 0.1875")

# paramters
R = 1
F = 50
C = 0.285
alpha = 0.01
r = 0.0114
T = 2
X = 50.5
Cp = 60

# terminal option
V_T = [max(F, R*S[i]) for i in range(len(S))]
plt.plot(S, V_T, label = 'max(F,RS)', linestyle = '--')

V_T2 = [max(Cp, R*S[i]) for i in range(len(S))]
plt.plot(S, V_T2, label = 'max($C_{P}$,RS)', linestyle = '--')

#theta line
plt.vlines(50.87/2, 0, max(V1), linestyle = '--', label = '$\\theta(t)$')



# graph details
plt.xlabel('S')
plt.xlim(0,70)
plt.ylim(47.5,62.5)
plt.ylabel('V(S,t=0)')
plt.title('American sytle convertible bond with embedded call option')
plt.legend()

plt.figure()
diff = [V3[i] - V1[i] for i in range(len(V1))]
plt.vlines(50.87, 0, max(diff), linestyle = '--', label = '$\\theta(t)$', color = 'red')
plt.plot(S, diff, label = 'Difference')
plt.legend()
plt.xlabel('S')
plt.ylabel('Differrence')
plt.title('Difference in option value with different $\kappa$')







