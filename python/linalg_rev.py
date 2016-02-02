import numpy as _np

from scipy.linalg import solve_triangular as _solve_t

def _Phi(A):
    """Helper routine: take tril, and halve the diagonal"""
    N = A.shape[0]
    A = _np.tril(A)
    # Could offer in-place version of routine with:
    # A[triu_indices(N,1)] = 0
    # But slower!
    A.ravel()[::N+1] *= 0.5
    return A

def chol_rev_np(L, L_bar, inplace=False):
    """Push derivatives back through Cholesky using analytic expression."""
    # "A_bar = L.T \ _Phi(L' * L_bar) / L"
    A_bar = _solve_t(L, _solve_t(L, _Phi(_np.dot(L.T, L_bar)).T, \
            trans=1, lower=True).T, trans=1, lower=True)
    A_bar = _Phi(A_bar + A_bar.T)
    if inplace:
        # Didn't really do stuff inplace, just being compatible with chol_rev
        L_bar[:] = A_bar
    else:
        return A_bar

# TODO For large matrices, someone could make the pure python+numpy version
# faster by porting the blocked algorithm, as in the pure Matlab code. It will
# still be slower than calling the Fortran code though.

try:
    import larmpack as _larmpack
    def chol_rev(L, L_bar, inplace=False):
        """Push derivatives back through Cholesky using differentiated Fortran."""
        N = L.shape[0]
        assert(L.shape == (N, N))
        assert(L.shape == L_bar.shape)
        info = _np.array(0)
        if inplace:
            _larmpack.dpofrt("L", N, L, L_bar, info)
            assert(info == 0)
        else:
            A_bar = L_bar.copy()
            _larmpack.dpofrt("L", N, L, A_bar, info)
            assert(info == 0)
            return A_bar
except:
    import sys
    sys.stderr.write('WARNING: Fortran not compiled; chol_rev will be slow.\n')
    chol_rev = chol_rev_np

