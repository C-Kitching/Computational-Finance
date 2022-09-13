import csv
import matplotlib.pyplot as plt
import numpy as np

# Read csv file
with open('portfolio_value_at_exercise.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    S = []
    V = []
    
    # loop over the rows
    for row in readCSV:
        S.append(row[0])
        V.append(row[1])
        
    # map to float array
    S = list(map(float, S))
    V = list(map(float, V))        


# graph details
plt.figure()
plt.plot(S, V, label = 'P(S, T)')
plt.xlabel('S')
plt.ylabel('V(S,T)')
plt.title('Portfolio value at exercise')
plt.legend()
plt.show()