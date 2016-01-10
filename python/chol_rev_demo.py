import numpy as np
from numpy.random import randn
from numpy.random import standard_normal as randnt
from numpy.linalg import cholesky as chol

from linalg_rev import chol_rev

# Example setup:
N = 500
A = np.cov(randn(N, 3*N)) # a symmetric +ve definite matrix
L = chol(A)               # lower-triangular Cholesky decomposition
Lbar = np.tril(randnt(L.shape)) # declare what dF/dL is for some function F

# Push the derivative back through the Cholesky, Abar = dF/d(tril(A))
Abar = chol_rev(L, Lbar)

# Perturb the input, look at the perturbation of the output and check
# for consistency.
dA = np.cov(randn(N, 3*N))
eps = 1e-5
dL = (chol(A + (eps/2)*dA) - chol(A - (eps/2)*dA)) / eps
dF1 = np.dot(dL.ravel(), Lbar.ravel())
dF2 = np.dot(dA.ravel(), Abar.ravel())

print('Error: %g' % float(dF1 - dF2))
