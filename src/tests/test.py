import numpy as np
import l2c
#import pytest

#@pytest.mark.parametrize('N', (1000, 5000, 10000))
def test_accuracy(N):
    u = np.ones(N)
    b = np.zeros_like(u)
    c = np.zeros_like(u)
    C = l2c.Leg2Cheb(u, b, 36, 2, 500)
    b = C(u, b, 0)
    c = C(b, c, 1)
    assert np.linalg.norm(c-1) < 1e-8
    del C

if __name__ == '__main__':
    test_accuracy(1000)
    test_accuracy(5000)
    test_accuracy(10000)
