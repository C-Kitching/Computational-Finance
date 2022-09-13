import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('American vary step.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    n= []
    #R = []
    #c = []
    V = []
    t = []
    
    # loop over the rows
    for row in readCSV:
        n.append(row[0])
        #R.append(row[1])
        #c.append(row[2])
        V.append(row[1])
        t.append(row[2])
        
        
    # map to float array
    n = list(map(float, n))
    #R = list(map(float, R))  
    #c = list(map(float, c))
    V = list(map(float, V)) 
    t = list(map(float, t))
    

plt.figure(figsize = (9,5))

plt.subplot(1,2,1)
plt.plot(n, V)
plt.xlabel('$i_{max}/j_{max}$')
plt.ylabel('$V(S,t=0)$')
plt.title('Option value')

plt.subplot(1,2,2)
plt.plot(n,t)
plt.hlines(1, min(n), max(n), linestyle = '--', color = 'r')
plt.xlabel('$i_{max}/j_{max}$')
plt.ylabel('$t_{comp}$')
plt.title('Computation time')


# common title
plt.suptitle('Effect of varying $i_{max}/j_{max}$')
plt.tight_layout()
plt.subplots_adjust(top=0.88)

