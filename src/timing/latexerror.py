import subprocess
import re
import numpy as np
import matplotlib.pyplot as plt

s = 32
arr = []
#number = "([0-9]+.[0-9]+e[+-][0-9]+)"
number = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
la = "\\begin{tabular}{ccccc}\n"
la += " N & L2C (m) & L2C (n) & C2L (m) & C2L (n) \\\\ \n"
la += "\\hline \\n"
for n in range(1, 4):
    N = s * 2**(n+2)
    h0 = subprocess.check_output(f"l2cprec -N{N} -d0 -a1 -n1 -R1", shell=True)
    h1 = subprocess.check_output(f"l2cprec -N{N} -d1 -a1 -n1 -R1", shell=True)
    h0 = np.array(re.findall(number, h0.decode())).astype(float)
    h1 = np.array(re.findall(number, h1.decode())).astype(float)
    la += f"{N:d} & {h0[0]:2.2e} ({h0[1]:2.1f}) & {h0[2]:2.2e} ({h0[3]:2.1f}) & {h1[0]:2.2e} ({h1[1]:2.1f}) & {h1[2]:2.2e} ({h1[3]:2.1f}) \\\\ \n"
    arr.append([N, h0[0], h0[2], h1[0], h1[2]])
la += "\\hline \n"
la += "\\end{tabular}"
print(la)

a = np.array(arr)
plt.loglog(a[:, 0], a[:, 1], 'b', a[:, 0], a[:, 3], 'b:')
plt.loglog(a[:, 0], a[:, 2], 'r', a[:, 0], a[:, 4], 'r:')
plt.show()
