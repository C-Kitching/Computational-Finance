import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# Read csv file1
with open('standard MC convergence.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    N = []
    diff = []
    
    # loop over the rows
    for row in readCSV:
        N.append(row[0])
        diff.append(row[1])
        
    # map to float array
    N = list(map(float, N))
    diff = list(map(float, diff))       


N.pop(0)
N.pop(0)
diff.pop(0)
diff.pop(0)

slope, intercept, r_value, p_value, std_err = stats.linregress(N, diff)
y = []
for i in range(len(N)):
    y.append(slope*N[i]+intercept)

print(slope)
print(r_value**2)    
    
    
# graph
plt.figure()
plt.plot(N, y)
plt.scatter(N, diff)
plt.xlabel('log(N)')
plt.ylabel('$log(V_{N}-V_{exact})$')
plt.title('Halton convergence')