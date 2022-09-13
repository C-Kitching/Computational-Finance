import csv
import matplotlib.pyplot as plt
import numpy as np

# Read csv file
with open('data.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    r = []
    P = []
    V = []
    
    # loop over the rows
    for row in readCSV:
        r.append(row[0])
        P.append(row[1])
        V.append(row[2])
        
    # map to float array
    r = list(map(float, r))
    P = list(map(float, P))
    V = list(map(float, V))        


# graph details
plt.figure()
plt.plot(r, P, label = 'P(r, t=0, T)')
plt.plot(r, V, label = 'V(r, t=0, T)')
plt.xlabel('r')
plt.ylabel('P and V')
plt.title('Comp Finance - Mini Task 2')
plt.legend()
plt.show()
        