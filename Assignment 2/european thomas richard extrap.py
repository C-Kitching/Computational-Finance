import csv
import matplotlib.pyplot as plt
import numpy as np


# Read csv file1
with open('european thomas richardson extrap.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    n = []
    R = []
    t = []
    
    # loop over the rows
    for row in readCSV:
        n.append(row[0])
        R.append(row[1])
        t.append(row[2])
        
    # map to float array
    n = list(map(float, n))
    R = list(map(float, R))
    t = list(map(float, t))      

# numerical values    
#plt.figure()
#plt.plot(n,R, color = 'red')
print(R)
# normalise
t_norm = [t[i]/max(t) for i in range(len(t))]
R_abs = [abs(R[i]) for i in range(len(R))]
R_norm = [R[i]/max(R_abs) for i in range(len(R_abs))]

print(R_norm)
print(t_norm)

plt.figure()
plt.plot(np.log(n), R_norm, label = 'Richardson extrapolation')
plt.plot(np.log(n), t_norm, label = 'Computation time')
plt.legend()
plt.xlabel('log(n)')
plt.ylabel('Normalised Scale')
plt.title('Optimal $i_{max}$ and $j_{max}$')


# paramters
R = 1
F = 50
C = 0.285
alpha = 0.01
r = 0.0114
T = 2
X = 50.5




