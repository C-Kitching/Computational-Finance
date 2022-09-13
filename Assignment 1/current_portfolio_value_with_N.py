import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('current_portfolio_value_with_N_S=X1.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    N = []
    V = []
    
    # loop over the rows
    for row in readCSV:
        N.append(row[0])
        V.append(row[1])
        
    # map to float array
    N = list(map(float, N))
    V = list(map(float, V))        


# graph
plt.figure()
plt.plot(N, V, label = 'Monte-Carlo approximation', zorder = 1)
plt.hlines(-1074.78, min(N), max(N), colors = 'red', label = 'Analytical result', zorder=2)
plt.xlabel('N')
plt.ylabel('V(S,t=0)')
plt.legend()
plt.title('Current portfolio value for $S_{0}=X1$')

# analystic
# X2 = -1074.78
# X1 = -997.365





