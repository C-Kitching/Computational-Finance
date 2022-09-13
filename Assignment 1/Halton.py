import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('Halton.csv') as csvfile:
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



new_N = N[0:50]
exact = -1074.78
log_diff = []
for i in range(len(new_N)):
    log_diff.append(np.log(V[i]-exact))
logN = np.log(new_N)
plt.figure()
plt.scatter(logN, log_diff)



# graph
plt.figure()
plt.plot(N, V, 'c', label = 'Monte-Carlo approximation', zorder = 1)
plt.hlines(-1074.78, min(N), max(N), colors = 'r', label = 'Analytical result', zorder=4)
plt.xlabel('N')
plt.ylabel('V(S,t=0)')
plt.ylim(-1075, -1073.75)
plt.legend()
plt.title('Current portfolio value for $S_{0}=X1$\nusing a Halton Sequence')