# This code was taken from the class notebook FDHighOrder
import numpy

def cosspace(a, b, n=50):
    return (a + b)/2 + (b - a)/2 * (numpy.cos(numpy.linspace(-numpy.pi, 0, n)))


def vander_chebyshev(x, n=None):
    """
    Generates a Vandermond matrix using Chebyshev polynomials.
    """
    if n is None:
        n = len(x)
    T = numpy.ones((len(x), n))
    if n > 1:
        T[:,1] = x
    for k in range(2,n):
        T[:,k] = 2 * x * T[:,k-1] - T[:,k-2]
    return T


def chebeval(z, n=None):
    """
    Build matrices to evaluate the n-term Chebyshev expansion and its 
    derivatives at point(s) z
    
    Returns three matrices, which can be pre-multiplied with a coefficient 
    vector of length n to produce
      Tz:   the Chebyshev expansion of a function evaluated at points in z
      dTz:  first derivative of the Chebyshev expansion evaluated at points z
      ddTz: second derivative of the Chebyshev expansion evaluated at points x
    """
    z = numpy.array(z, ndmin=1)

    if n is None:
        n = len(z)
    Tz = vander_chebyshev(z, n)
    dTz = numpy.zeros_like(Tz)
    dTz[:,1] = 1
    dTz[:,2] = 4*z
    ddTz = numpy.zeros_like(Tz)
    ddTz[:,2] = 4

    for n in range(3,n):
        dTz[:,n]  = n * (2*Tz[:,n-1] + dTz[:,n-2]/(n-2))
        ddTz[:,n] = n * (2*dTz[:,n-1] + ddTz[:,n-2]/(n-2))

    return [Tz, dTz, ddTz]


class exact_tanh:
    def __init__(self, k=1, x0=0):
        self.k = k
        self.x0 = x0
    def u(self, x):
        return numpy.tanh(self.k*(x - self.x0))
    def du(self, x):
        return self.k * numpy.cosh(self.k*(x - self.x0))**(-2)
    def ddu(self, x):
        return -2 * self.k**2 * numpy.tanh(self.k*(x - self.x0)) * numpy.cosh(self.k*(x - self.x0))**(-2)
