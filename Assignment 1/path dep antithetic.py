import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns


# Read csv file1
with open('path dep antithetic.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    N = []
    t = []
    
    # loop over the rows
    for row in readCSV:
        N.append(row[0])
        t.append(row[1])
        
    # map to float array
    N = list(map(float, N))
    t = list(map(float, t))
    


ax = sns.regplot(N, t, ci=99)
ax.hlines(10, min(N), max(N), color = 'red', linestyle = '--')
ax.set_xlabel('N')
ax.set_ylabel('Time')
ax.set_title('Computation time')

slope, intercept, r_value, p_value, std_err = stats.linregress(N, t)
print('m = ', slope)
print('c = ', intercept)
print(r_value**2)    

print((10+intercept)/slope)

