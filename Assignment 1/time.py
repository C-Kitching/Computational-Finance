import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('time.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    N = []
    standard_time = []
    antithetic_time = []
    Halton_time = []
    standard_value = []
    antithetic_value = []
    halton_value = []
    
    # loop over the rows
    for row in readCSV:
        N.append(row[0])
        standard_time.append(row[1])
        antithetic_time.append(row[2])
        Halton_time.append(row[3])
        standard_value.append(row[4])
        antithetic_value.append(row[5])
        halton_value.append(row[6])
        
    # map to float array
    N = list(map(float, N))
    standard_time = list(map(float, standard_time))       
    antithetic_time = list(map(float, antithetic_time))    
    Halton_time = list(map(float, Halton_time)) 
    standard_value = list(map(float,standard_value))
    antithetic_value = list(map(float,antithetic_value))
    halton_value = list(map(float,halton_value))

# errors
exact = -1074.78
standard_error = [abs((V-exact)/exact) for V in standard_value]
antithetic_error = [abs((V-exact)/exact) for V in antithetic_value]
halton_error = [abs((V-exact)/exact) for V in halton_value]

# grpah
plt.figure()
plt.plot(N, standard_error, label = 'Standard MC')
plt.plot(N, antithetic_error, label = 'Antithetic variables')
plt.plot(N, halton_error, label = 'Halton sequenece')
plt.xlabel('N')
plt.ylabel('Relative error')
plt.title('Error comparison')
plt.legend()
plt.show()


# graph
plt.figure()
plt.plot(N, standard_time, label = 'Standard MC')
plt.plot(N, antithetic_time, label = 'Antithetic variables')
plt.plot(N, Halton_time, label = 'Halton sequence')
plt.legend()
plt.xlabel('N')
plt.ylabel('t')
plt.title('Computation time comparison')