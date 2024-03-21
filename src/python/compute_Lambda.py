"""
Compute the Chebyshev coefficients of the Lambda function with high precision.
"""
from mpmath import mp
import sympy as sp
import numpy as np
import scipy
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

def lamx(x, nx=7):
    y = 1 - 1/64/(x+q)**2
    s = mp.sqrt(x+q)
    if nx == 2:
        return y/s
    y += 21/8192/(x+q)**4
    if nx == 3:
        return y/s
    y -= 671/524288/(x+q)**6
    if nx == 4:
        return y/s
    y += 180323/134217728/(x+q)**8
    if nx == 5:
        return y/s
    y -= 20898423/8589934592/(x+q)**10
    if nx == 6:
        return y/s
    y += 7426362705/1099511627776/(x+q)**12
    return y/s

def lambN1(x):
    return np.exp(scipy.special.loggamma(x + 0.5) - scipy.special.loggamma(x+1)) * np.sqrt(x)

def lambN2(x):
    return np.exp(scipy.special.loggamma(x + 0.25) - scipy.special.loggamma(x+0.75)) * np.sqrt(x)

def lambHH(x):
    return mp.gamma(1/x+mp.mpf('1/4')) / mp.gamma(1/x+mp.mpf('3/4')) / mp.sqrt(x)

def lambs(x):
    return mp.gamma(1/(x+q)+mp.mpf('1/2')) / mp.gamma(1/(x+q)+1) / mp.sqrt(x+q)

def lambHH(x):
    return mp.gamma(1/x+mp.mpf('1/4')) / mp.gamma(1/x+mp.mpf('3/4')) / mp.sqrt(x)

def lambH(x):
    return mp.gamma(x+mp.mpf('1/4')) / mp.gamma(x+mpf('3/4'))

def lamb(x):
    return mp.gamma(x+mp.mpf('1/2')) / mp.gamma(x+1)

def lambi(x, a):
    return mp.sqrt(2*a/(x+1)) * mp.gamma(2*a/(x+1)+mp.mpf('1/2')) / mp.gamma(2*a/(x+1)+1)

def lambiH(x, a):
    return mp.sqrt(2*a/(x+1)) * mp.gamma(2*a/(x+1)+mp.mpf('1/4')) / mp.gamma(2*a/(x+1)+mp.mpf('3/4'))

def lambiG(x, a):
    return mp.sqrt(mp.sqrt(2*a**2/(x+1))) * mp.gamma(mp.sqrt(2*a**2/(x+1))+mp.mpf('1/4')) / mp.gamma(mp.sqrt(2*a**2/(x+1))+mp.mpf('3/4'))

def lambiE(x, a):
    return mp.gamma(mp.sqrt(2/(x+1))*a+mp.mpf('1/2')) / mp.gamma(mp.sqrt(2/(x+1))*a+1)

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

def curve(x, a):
    return 10 + a / np.log2(x)

#def lambC(x, a=64, nx=8):
#    c = chebquad(10, a, nx)
#    return chebval(-1+2*a/x, np.array(c).astype(float))/np.sqrt(x)

def lambC(x, a=64, nx=8):
    c = chebquad(10, a, nx)
    return [chebval(-1+2*a/i, c) / mp.sqrt(i) for i in x]

def lambG(x, a=64, nx=8):
    c = chebquad(10, a, nx, lambiG)
    return [chebval(2*(a/(i+q))**2-1, c) / mp.sqrt(i+q) for i in x]

def fx(y, a, nx=8):
    e = [lamb(x) for x in y]
    f = lambC(y, a, nx)
    return [abs(f[i] - e[i]) for i in range(len(y))]

def fg(y, a, nx=8):
    e = [lamb(x) for x in y]
    f = lambG(y, a, nx)
    return [abs(f[i] - e[i]) for i in range(len(y))]

dfx = lambda x, nx: float(lamx(x, nx) - lamb(x))
dff = lambda x, nx: float(lamxf(float(x), nx) - lamb(x))

#np.savetxt(f"Lambdavals{M}.dat", np.array(a[:8]).astype(float), fmt="%.16e", delimiter=',', newline=',')
#np.savetxt(f"Lambdavals{M}C.dat", np.array(chebquad(20, M)).astype(float), fmt="%.16e", delimiter=',', newline=',')

#L = []
#for i in range(64):
#    L.append(lamb(i))

F = lambda x, nx: mp.log10(mp.absmax(lamx(x, nx)-lamb(x)))-mp.log10(mp.mpf('2.22e-16'))
#mp.findroot(lambda x: F(x, 7), 10)

y = mp.linspace(400, 1000, 100)
plt.semilogy(y, abs(np.array([dfx(i, 2) for i in y])), 'b',
             y, fx(y, 10, 6), 'r',
             y, fg(y, 10, 6), 'k',
             y, abs(np.array([dff(i, 2) for i in y])), 'm')