import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('european thomas vary S data 1.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    S1 = []
    V1 = []
    
    # loop over the rows
    for row in readCSV:
        S1.append(row[0])
        V1.append(row[1])
        
    # map to float array
    S1 = list(map(float, S1))
    V1 = list(map(float, V1))   
    
# Read csv file1
with open('european thomas vary S data 2.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    S2 = []
    V2 = []
    
    # loop over the rows
    for row in readCSV:
        S2.append(row[0])
        V2.append(row[1])
        
    # map to float array
    S2 = list(map(float, S2))
    V2 = list(map(float, V2)) 

# numerical values    
plt.figure()
plt.plot(S1,V1, color = 'blue', label = "$\\beta$ = 1\n$\sigma$ = 0.4")
plt.plot(S2,V2, color = 'red', label = "$\\beta$ = 0.869\n$\sigma$ = 0.668")

# paramters
R = 1
F = 50
C = 0.285
alpha = 0.01
r = 0.0114
T = 2
sigma = 0.668
X = 50.5
kappa = 0.125

sigma1 = 0.4
beta1 = 1
sigma2 = 0.668
beta2 = 0.869

# terminal option
V_T = [max(F, R*S1[i]) for i in range(len(S1))]
plt.plot(S1, V_T, label = 'max(F,RS)', linestyle = '--')

# AS+B
ASB = [R*np.exp(-(kappa+r)*T)*S1[i]+X*R*np.exp(-r*T)*(1-np.exp(-kappa*T))+(C/(alpha+r))*(1-np.exp(-r*T)) for i in range(len(S1))]
plt.plot(S1, ASB, label = '$V(S\\rightarrow\infty,0)$', linestyle = '--')

# graph details
plt.xlabel('S')
plt.ylim(0, max(V1) + 100)
plt.ylabel('Value')
plt.title('Option price varying with S\nThomas solver')
plt.legend()

# difference
plt.figure()
diff = [abs(V1[i] - V2[i]) for i in range(len(V1))]
plt.xlabel('S')
plt.ylabel('Absolute difference')
plt.title('Option value difference for two parameter sets')
plt.plot(S1,diff)








