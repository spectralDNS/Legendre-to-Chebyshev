"""
Compute the Chebyshev coefficients of the Lambda function with high precision.
"""
from mpmath import mp
import sympy as sp
import numpy as np
import scipy
import math
import matplotlib.pyplot as plt
np.set_printoptions(precision=16, suppress=False)

mp.dps = 100

q = sp.Rational(1, 4)

def lamxf(x, nx=7):
    y = 1 - 0.015625/(x+0.25)**2
    s = np.sqrt(x+0.25)
    if nx == 2:
        return y/s
    y += 0.0025634765625/(x+0.25)**4
    if nx == 3:
        return y/s
    y -= 0.0012798309326171875/(x+0.25)**6
    if nx == 4:
        return y/s
    y += 0.0013435110449790955/(x+0.25)**8
    if nx == 5:
        return y/s
    y -= 0.0024328966392204165/(x+0.25)**10
    if nx == 6:
        return y/s
    y += 0.006754237533641572/(x+0.25)**12
    return y/s

def taux(x, nx=7):
    y = 1 - 1/64/(x+q)**2
    if nx == 2:
        return y
    y += 21/8192/(x+q)**4
    if nx == 3:
        return y
    y -= 671/524288/(x+q)**6
    if nx == 4:
        return y
    y += 180323/134217728/(x+q)**8
    if nx == 5:
        return y
    y -= 20898423/8589934592/(x+q)**10
    if nx == 6:
        return y
    y += 7426362705/1099511627776/(x+q)**12
    return y

def lamx(x, nx=7):
    return taux(x, nx) / mp.sqrt(x+q)

def lamb(x):
    return mp.gamma(x+mp.mpf('1/2')) / mp.gamma(x+1)

def lambi(x, a):
    return mp.sqrt(2*a/(x+1)) * mp.gamma(2*a/(x+1)+mp.mpf('1/2')) / mp.gamma(2*a/(x+1)+1)

def lambiG(x, a):
    return mp.sqrt(mp.sqrt(2*a**2/(x+1))) * mp.gamma(mp.sqrt(2*a**2/(x+1))+mp.mpf('1/4')) / mp.gamma(mp.sqrt(2*a**2/(x+1))+mp.mpf('3/4'))

def chebquad(N, a, nx=8, fun=lambi):
    xj = lambda k, N: mp.cos(mp.pi*(2*k+1)/(2*N))
    res = []
    for j in range(nx):
        s = 0
        for k in range(N):
            x = xj(k, N)
            s += fun(x, a)*mp.cos(j*mp.pi*(2*k+1)/(2*N))
        res.append(s*2/N)
    res[0] /= 2
    return res

def chebval(x, c):
    if len(c) == 1:
        c0 = c[0]
        c1 = 0
    elif len(c) == 2:
        c0 = c[0]
        c1 = c[1]
    else:
        x2 = 2*x
        c0 = c[-2]
        c1 = c[-1]
        for i in range(3, len(c) + 1):
            tmp = c0
            c0 = c[-i] - c1
            c1 = tmp + c1*x2
    return c0 + c1*x


def lambG(x, a=64, nx=8):
    c = chebquad(10, a, nx, lambiG)
    return [chebval(2*(a/(i+q))**2-1, c) / mp.sqrt(i+q) for i in x]

def lambGf(x, a=64, nx=8):
    c = chebquad(10, a, nx, lambiG)
    return [chebval(2*(a/(float(i)+0.25))**2-1, np.array(c).astype(np.double)) / np.sqrt(float(i)+0.25) for i in x]

def fx(y, a, nx=8):
    e = [lamb(x) for x in y]
    f = lambC(y, a, nx)
    return [abs(f[i] - e[i]) for i in range(len(y))]

def fg(y, a, nx=8):
    e = [lamb(x) for x in y]
    f = lambGf(y, a, nx)
    return [abs((f[i] - e[i])/e[i]) for i in range(len(y))]

dfx = lambda x, nx: float(lamx(x, nx) - lamb(x))
dff = lambda x, nx: float(lamxf(float(x), nx) - lamb(x))

#np.savetxt(f"Lambdavals{M}.dat", np.array(a[:8]).astype(float), fmt="%.16e", delimiter=',', newline=',')
#np.savetxt(f"Lambdavals{M}C.dat", np.array(chebquad(20, M)).astype(float), fmt="%.16e", delimiter=',', newline=',')

F = lambda x, nx: mp.log10(mp.absmax(lamx(x, nx)-lamb(x)))-mp.log10(math.nextafter(lamxf(2000), 1e8)-lamxf(2000))
#F = lambda x, nx: mp.log10(mp.absmax((lamx(x, nx)-lamb(x))/lamb(x)))-mp.log10(2.22e-16)
#F = lambda x, nx: mp.log10(mp.absmax((lamx(x, nx)-lamb(x))))-mp.log10(2.22e-16)

mp.findroot(lambda x: F(x, 4), 10)

y = mp.linspace(64, 2000, 200)
plt.semilogy(y, abs(np.array([dfx(i, 2) for i in y])), 'b',
             y, fx(y, 10, 6), 'r',
             y, fg(y, 64, 3), 'k',
             y, abs(np.array([dff(i, 2) for i in y])), 'm')