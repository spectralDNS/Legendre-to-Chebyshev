"""
Arbitrary precision Legendre to Chebyshev transforms
"""
from mpmath import mp

mp.dps = 50

#Lambda = lambda z: mp.exp(mp.loggamma(z + mp.mpf('1/2')) - mp.loggamma(z+1))
Lambda = lambda z: mp.gamma(z+mp.mpf('1/2')) / mp.gamma(z+1)

def leg2chebmp(u, b, N):
    a = [Lambda(i) for i in range(N)]
    for n in range(0, N, 2):
        for i in range(0, N-n):
            b[i] += a[n//2]*a[n//2+i]*u[n+i]
    b[0] /= 2
    for i in range(0, N):
        b[i] *= (2/mp.pi)

def cheb2legmp(u, b, N):
    dn = [Lambda(i-1)/(2*i) for i in range(1, N//2)]
    dn.insert(0, 0)
    a = [1/(2*Lambda(i)*i*(i+mp.mpf('1/2'))) for i in range(1, N)]
    a.insert(0, 2/mp.sqrt(mp.pi))
    un = u.copy()
    for i in range(2, N):
        un[i] = u[i]*i
    for n in range(N):
        b[n] = mp.sqrt(mp.pi)*a[n]*un[n]
    for n in range(2, N, 2):
        for i in range(0, N-n):
            b[i] -= dn[n//2]*a[n//2+i]*un[n+i]
    for n in range(N):
        b[n] *= (n+mp.mpf('1/2'))

if __name__ == '__main__':
    N = 10
    #u = mp.matrix([mp.mpf(i) for i in range(N)])
    u = mp.matrix([mp.mpf('1')]*N)
    b = mp.matrix([mp.mpf('0')]*N)
    c = mp.matrix([mp.mpf('0')]*N)
    leg2chebmp(u, b, N)
    cheb2legmp(b, c, N)
