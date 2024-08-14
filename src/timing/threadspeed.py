import subprocess
import re
import numpy as np
import matplotlib.pyplot as plt 

s = 32
arr = []
N = np.array([2**17, 2**20, 2**23])
#number = "([0-9]+.[0-9]+e[+-][0-9]+)"
number = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
la = "\\begin{tabular}{ccccc}\n"
la += " N & 1 & 2 & 4 & 1 & 2 & 4\\\\ \n"
la += "\\hline \n"

tt = np.zeros((3, 6))
t0 = np.zeros((3, 6))
for i, n in enumerate(N):
    la += f' {n} &'
    for j, t in enumerate([1, 2, 4]):
        h0 = subprocess.check_output(f"l2c -N{n} -d0 -r100 -s{s} -t{t} -v0", shell=True)
        h1 = subprocess.check_output(f"l2c -N{n} -d7 -r100 -s{s} -t{t} -v0", shell=True)
        h0 = np.array(re.findall(number, h0.decode())).astype(float)
        h1 = np.array(re.findall(number, h1.decode())).astype(float)
        tt[i, j] = float(h0[1])
        tt[i, j+3] = float(h1[1])
    t0[i, :3] = tt[i, :3] / tt[i, 0]
    t0[i, 3:] = tt[i, 3:] / tt[i, 3] 
    la += f' {t0[i, 0]} & {t0[i, 1]} & {t0[i, 2]} & {t0[i, 3]} & {t0[i, 4]} & {t0[i, 5]} \\\\ \n'
la += "\\hline \n"
la += "\\end{tabular}"
print(la)
plt.loglog([1, 2, 4], np.transpose(t0[:, :3]), 'k', [1, 2, 4], np.transpose(t0[:, 3:]), 'k:', [1, 2, 4], t0[:, 0]/np.array([1, 2, 4]), 'k--')
plt.show()
