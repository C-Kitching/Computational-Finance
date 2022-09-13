import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


# Read csv file1
with open('normal_test.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ',')
    
    # empty arrays for data
    r1 = []
    r2 = []
    r3 = []
    r4 = []
    st1 = []
    st2 = []
    st3 = []
    st4 = []
    
    # loop over the rows
    for row in readCSV:
        r1.append(row[0])
        r2.append(row[1])
        r3.append(row[2])
        r4.append(row[3])
        st1.append(row[0])
        st2.append(row[1])
        st3.append(row[2])
        st4.append(row[3])
        
    # map to float array
    r1 = list(map(float, r1))
    r2 = list(map(float, r2))
    r3 = list(map(float, r3))
    r4 = list(map(float, r4))
    st1 = list(map(float, st1))
    st2 = list(map(float, st2))
    st3 = list(map(float, st3))
    st4 = list(map(float, st4))
    
plt.figure()
pseudo = r1 + r2
halton = r3 + r4

pseudo_JB = stats.jarque_bera(pseudo)
halton_JB = stats.jarque_bera(halton)
p_pseudo = pseudo_JB.pvalue
p_halton = halton_JB.pvalue

print(halton_JB)

plt.figure()
plt.hist(halton, bins = 100, label = 'Halton, p_value = %s' %p_halton, alpha = 0.5)
plt.hist(pseudo, bins = 100, label = 'Pseduo, p_value = %s' %p_pseudo, alpha = 0.5)
plt.xlabel('x')
plt.ylabel('Frequency')
plt.title('Pseudo vs Halton random numbers')
plt.legend()




st_halton = st1 + st2
st_pseudo = st3 + st4
plt.figure()
plt.hist(st_halton, bins = 100, alpha = 0.5, label = 'Halton')
plt.hist(st_pseudo, bins = 100, alpha = 0.5, label = 'Pseudo')
plt.xlabel('x')
plt.ylabel('Frequency')
plt.title('Pseudo vs Halton final stock path values')
plt.legend()




