#plot the probability of X > n/2 as a function of margin (m) and number of voters (n).
#Given population size N  and perceived margin m,  n= min(N, 1/m^2) and X is selected from binomial distribution with parameters n and p  = pA(1) - pB(1).
# we consider p \in [0.5,0.51].

import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
# Define the range of z values
legend_properties = {'weight':'bold'}
p_values = np.linspace(0.50001,0.51, 10000)

def transpose(l1, l2):
    # iterate over list l1 to the length of an item
    for i in range(len(l1[0])):
        # print(i)
        row =[]
        for item in l1:
            # appending to new list with values and index positions
            # i contains index position and item contains values
            row.append(item[i])
        l2.append(row)
    return l2

def computeFc(alpha, c,N):
  s_A = 0.1 + c/2
  s_B = 0.4
  return   c - pow(s_A + s_B, alpha - 0.5)/(pow(np.abs(s_A - s_B), alpha)*np.sqrt(N))

def equilibriumWinProb(c, N):
  s_A = 0.1 + c/2
  s_B = 0.4
  p_A = s_A/(s_A + s_B)
  n_A  = np.ceil(N*(s_A + s_B))
  return 1 - st.binom.cdf(np.ceil(n_A/2), n_A, p_A)

fig, ax = plt.subplots(figsize=(8, 6))
cmin=0.6
cmax=1
c = cmax
N = 10000
itr =0
NArray =[10000,100000,1000000,10000000,100000000]
idxalpha = 0
alphaVal = 17
nVal = 5
alphaArray = [0.5, 0.55, 0.6, 0.65, 0.7,0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 1, 1.025, 1.05, 1.1,1.15, 1.2]
cArray = [ [0] * nVal for _ in range(alphaVal)]
winProb= [ [0] * nVal for _ in range(alphaVal)]
plots= [ [0] * alphaVal for _ in range(nVal)]
for alpha in alphaArray:
  idx = 0
  for N in NArray:
    cmin=0.6
    cmax=1
    c = cmax
    f_c= computeFc(alpha, c,N)
    while np.abs(f_c) > 0.0000001:
      c = (cmax + cmin)/2
      f_c = computeFc(alpha, c,N)
      if f_c > 0:
        cmax = c
      else:
        cmin = c
    cArray[idxalpha][idx] = c
    winProb[idxalpha][idx] = equilibriumWinProb(c, N)
    idx = idx + 1
  idxalpha = idxalpha +1
transpose(winProb, plots)
i=0
while i < nVal:
  plt.plot(alphaArray,plots[i-nVal], marker='o', label=f'N={NArray[i-nVal]}')
  i = i+1
print(plots)
plt.legend()
plt.savefig('samplefigure.png')

plt.show()