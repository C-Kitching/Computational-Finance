import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('path dependence with confidence intervals K.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    K = []
    V = []
    upper = []
    lower = []
    
    # loop over the rows
    for row in readCSV:
        K.append(row[0])
        V.append(row[1])
        upper.append(row[2])
        lower.append(row[3])
        
    # map to float array
    K = list(map(float, K))
    V = list(map(float, V))  
    upper = list(map(float, upper))
    lower = list(map(float, lower))  


# graph
plt.figure()
plt.plot(K, V, 'c', label = 'Monte-Carlo approximation', zorder = 1)
plt.plot(K, upper, zorder = 2, label = 'upper 95% confidence interval')
plt.plot(K, lower, zorder = 3, label = 'lower 95% confidence interval')
plt.xlabel('K')
plt.ylabel('V(S,t=0)')
plt.legend()
plt.title('Floating strike Asian call option with $N=200,000$')