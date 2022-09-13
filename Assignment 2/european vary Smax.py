import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('european thomas Smax.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    S = []
    V = []
    t = []
    
    # loop over the rows
    for row in readCSV:
        S.append(row[0])
        V.append(row[1])
        t.append(row[2])
        
    # map to float array
    S = list(map(float, S))
    V = list(map(float, V))
    t = list(map(float, t))      



plt.figure(figsize = (9,5))

plt.subplot(1,2,1)
plt.sca(S, V)
plt.xlabel('$S_{max} \hspace{1} (n*F)$')
plt.ylabel('$V(S,t=0)$')
plt.title('Option value')

plt.subplot(1,2,2)
plt.plot(S, t)
plt.xlabel('$S_{max} \hspace{1} (n*F)$')
plt.ylabel('$t_{comp}$')
plt.title('Computation time')


# common title
plt.suptitle('Effect of varying $S_{max}$')
plt.tight_layout()
plt.subplots_adjust(top=0.88)


