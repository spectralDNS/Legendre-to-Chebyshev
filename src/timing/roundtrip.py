import subprocess
import re
import numpy as np
import matplotlib.pyplot as plt

number = "([0-9]+.[0-9]+e[+-][0-9]+)"

s = 32
arr = []
la = "\\begin{tabular}{cc}\n"
la += " N & r=0 & r=1/2  \\\\ \n"
la += "\\hline \\n"
for n in range(2, 16):
    N = s * 2**(n+2)
    h = subprocess.check_output(f"l2c -N{N} -d2 -l0 -m0.0 -v1 -R1", shell=True)
    g = re.search(number, h.decode())
    t = float(g.groups()[0])
    h1 = subprocess.check_output(f"l2c -N{N} -d2 -l0 -m0.5 -v1 -R1", shell=True)
    g1 = re.search(number, h1.decode())
    t1 = float(g1.groups()[0])
    #h2 = subprocess.check_output(f"l2c -N{N} -d2 -l0 -m0.5 -v1 -R1 -M18", shell=True)
    #g2 = re.search(number, h2.decode())
    #t2 = float(g2.groups()[0])
    #h3 = subprocess.check_output(f"l2c -N{N} -d2 -l1 -m0.5 -v1 -R1 -M18", shell=True)
    #g3 = re.search(number, h3.decode())
    #t3 = float(g3.groups()[0])

    la += f"{N:d} & {t:2.1e} & {t1:2.1e} \\\\ \n"
    arr.append([N, t, t1])
la += "\\hline \n"
la += "\\end{tabular}"
print(la)

a = np.array(arr)
plt.loglog(a[:, 0], a[:, 1], 'k', a[:, 0], a[:, 2], 'k--')
#plt.loglog(a[:, 0], a[:, 3], 'r', a[:, 0], a[:, 4], 'r:')
plt.ylim([1e-16, 1e-10])
plt.legend(['r = 0', 'r = 1/2'])
plt.ylabel('$E_{\\infty}$')
plt.xlabel('N')
plt.savefig('accuracy_roundtrip.png')
plt.show()
