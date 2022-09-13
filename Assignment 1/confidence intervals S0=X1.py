import csv
import matplotlib.pyplot as plt
import numpy as np

# Read csv file
with open('confidence_intervals_S=X1.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    N = []
    V = []
    upper = []
    lower = []
    
    # loop over the rows
    for row in readCSV:
        N.append(row[0])
        V.append(row[1])
        upper.append(row[2])
        lower.append(row[3])
        
    # map to float array
    N = list(map(float, N))
    V = list(map(float, V)) 
    upper = list(map(float, upper))
    lower = list(map(float, lower)) 


# graph
plt.figure()
#plt.plot(N, V, 'c', label = 'Monte-Carlo approximation', zorder = 1)
plt.plot(N, upper, zorder = 2, label = 'upper 95% confidence interval')
plt.plot(N, lower, zorder = 3, label = 'lower 95% confidence interval')
plt.hlines(-1074.78, min(N), max(N), colors = 'r', label = 'Analytical result', zorder=4)
plt.xlabel('N')
plt.ylabel('V(S,t=0)')
plt.legend()
plt.title('Current portfolio value for $S_{0}=X1$ \nwith 95% confidence intervals for 100 calculations')

# analystic
# X1 = -1074.78
# X2 = -997.365