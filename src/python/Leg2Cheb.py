import numpy as np
from l2c import Lambda
from mpi4py_fft import fftw
from numpy.polynomial.chebyshev import chebvander
import numba as nb

Mxy = lambda x, y: Lambda((y-x)/2)*Lambda((y+x)/2)
Lxy = lambda x, y: -1/(x+y+1)/(y-x)*Lambda((y-x-2)/2)*Lambda((y+x-1)/2)

class Leg2Cheb:
    def __init__(self, input_array, M=18, s=36, fun=Mxy):
        self.M = M
        N = input_array.shape[0]
        self.L = L = np.ceil(np.log2(N/s)-2).astype(int)
        self.nb = 2**(np.arange(L+1)+1)-1
        self.h = 2**(L-np.arange(L)-1)
        self.s = s = np.ceil(N/2**(L+2)).astype(int)
        self.N = s*2**(L+2)
        self.a = Lambda(np.arange(self.N, dtype='d'))
        xj = np.cos((np.arange(M)+0.5)*np.pi/M)
        dct = fftw.dctn(np.zeros((M, M)), axes=(0, 1), type=2)
        self.TT = self.conversionmatrix(M)
        self.T, self.A = [], []
        for lv in range(L):
            self.A.append([])
            h = s*self.h[lv]
            self.T.append(self.vandermonde(h, M))
            for block in range(self.nb[lv]):
                i0 = 2*block*s*2**(L-lv)
                j0 = 2*s*2**(L-lv)+i0
                for q in range(2):
                    y = (j0+q*2*h+(xj+1)*h)[None, :]
                    for p in range(q+1):
                        x = (i0+p*2*h+(xj+1)*h)[:, None]
                        A_hat = fun(x, y)
                        A_hat = dct(A_hat)/M**2
                        A_hat[0] /= 2
                        A_hat[:, 0] /= 2
                        self.A[-1].append(A_hat.copy())

    def __call__(self, input_array):
        N = len(input_array)
        input_array = np.pad(input_array, [0, self.N-N])
        output_array = np.zeros_like(input_array)
        for oe in (0, 1):
            self.apply(oe, input_array[oe::2], output_array[oe::2])
        self.direct(output_array, input_array, self.a, 2*self.s,
                    self.nb[self.L], self.N)
        output_array *= (2/np.pi)
        output_array[0] *= 0.5
        return output_array[:N].copy()

    def apply(self, oe, invec, outvec):
        M, L, s, T, TT = self.M, self.L, self.s, self.T, self.TT
        wk, ck = [], []
        for lv in range(self.L):
            wk.append(np.zeros((self.nb[lv], 2, M)))
            ck.append(np.zeros((self.nb[lv], 2, M)))
        wk[L-1][:] = np.dot(invec[2*s:].reshape((-1, s)),
                            T[L-1][oe]).reshape((-1, 2, M))
        for lv in range(L-1, 0, -1):
            wk[lv-1] += (np.dot(wk[lv][1:, 0], TT[0].T) +
                         np.dot(wk[lv][1:, 1], TT[1].T)).reshape((-1, 2, M))
        for lv in range(L-1, -1, -1):
            A = iter(self.A[lv])
            for block in range(self.nb[lv]):
                for q in range(2):
                    for p in range(q+1):
                        ck[lv][block, p] += np.dot(next(A), wk[lv][block, q])
        for lv in range(self.L-1):
            for p in range(2):
                ck[lv+1][:-1, p] += np.dot(ck[lv], TT[p]).reshape((-1, M))
        outvec[:-2*s] = np.dot(ck[L-1].reshape((-1, M)), T[L-1][oe].T).ravel()

    @staticmethod
    @nb.njit(cache=True)
    def direct(oa, ia, a, h, Nb, N):
        for n in range(h//2):
            oa[:(N-2*n)] += a[n]*a[n:(N-n)]*ia[2*n:]
        for block in range(Nb):
            i0, j0 = block*h, h+block*h
            for n in range(0, h, 2):
                a0 = i0+(n+h)//2
                oa[i0:i0+h-n] += a[(n+h)//2]*(a[a0:a0+h-n]*ia[j0+n:j0+h])

    def vandermonde(self, h, M):
        V = chebvander(np.linspace(-1, 1, 2*h+1)[:-1], M-1)
        return np.vstack((V[::2], V[1::2])).reshape((2, h, M))

    def conversionmatrix(self, M):
        T = np.zeros((2, M, M))
        k = np.arange(M)
        X = np.cos((k+0.5)*np.pi/M)[None, :]
        dct = fftw.dctn(np.zeros((M, M)), axes=(1,))
        for q in range(2):
            T[q] = dct(np.cos(k[:, None]*np.arccos((X+2*q-1)/2)))/M
            T[q, :, 0] /= 2
        return T

class Cheb2Leg(Leg2Cheb):
    def __init__(self, input_array, M=20):
        Leg2Cheb.__init__(self, input_array, M=M, fun=Lxy)
        k = np.arange(self.N, dtype='d')
        k[0] = 1
        self.dn = Lambda((k[::2]-2)/2)/k[::2]
        self.a = 1/(2*Lambda(k)*k*(k+0.5))
        self.a[0] = 2/np.sqrt(np.pi)
        self.dn[0] = 0

    @staticmethod
    @nb.njit(cache=True)
    def direct(oa, ia, a, dn, h, Nb, N):
        for n in range(1, h//2):
            oa[:(N-2*n)] -= dn[n]*a[n:(N-n)]*ia[2*n:]
        for block in range(Nb):
            i0, j0 = block*h, h+block*h
            for n in range(0, h, 2):
                a0 = i0+(n+h)//2
                oa[i0:i0+h-n] -= dn[(n+h)//2]*(a[a0:a0+h-n]*ia[j0+n:j0+h])

    def __call__(self, input_array):
        N = len(input_array)
        input_array.resize(self.N, refcheck=False)
        output_array = np.zeros_like(input_array)
        w0 = input_array.copy()
        w0[2:] *= np.arange(2, self.N)
        for oe in (0, 1):
            self.apply(oe, w0[oe::2], output_array[oe::2])
        output_array[:] += np.sqrt(np.pi)*self.a*w0
        self.direct(output_array, w0, self.a, self.dn, 2*self.s,
                    self.nb[self.L], self.N)
        output_array *= (np.arange(self.N)+0.5)
        return output_array[:N].copy()

if __name__ == '__main__':
    from shenfun import legendre
    N = 500
    u = np.ones(N)
    L2C = Leg2Cheb(u)
    C2L = Cheb2Leg(u)
    a = np.zeros(N)
    b = np.zeros(N)
    a = legendre.leg2cheb(u, a, axis=0)
    b = legendre.cheb2leg(a, b, axis=0)
    a1 = L2C(u)
    a2 = C2L(a1)
    error_d = np.linalg.norm(b-u, np.inf)
    error_f = np.linalg.norm(a2-u, np.inf)
    assert error_d < 1e-8 and error_f < 1e-8
    print(f'Errornorm direct: {error_d:2.6e}')
    print(f'Errornorm FMM   : {error_f:2.6e}')

