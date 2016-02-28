"""Python code for pushing derivatives through the Cholesky decomposition

Implements forwards- and reverse-mode update rules from:
    Differentiation of the Cholesky decomposition,
    Iain Murray, February 2016.
    http://arxiv.org/abs/1602.07527
"""

import numpy as np
from numpy import tril
from scipy.linalg import solve_triangular as _solve_triangular

try:
    import larmpack as _larmpack
    FORTRAN_COMPILED = True
except:
    FORTRAN_COMPILED = False

# Some comments and documentation refer to the infix matrix multiplication
# operator "@", available in Python >= 3.5 with NumPy >= 1.10. For this release
# I have used the dot() function in the actual code, to maintain compatibility
# with older Python+NumPy.
#
# I haven't profiled or optimized. The pure Python code was ported more-or-less
# straight from pseudo-code, mainly to check the pseudo-code was correct! For
# really large matrices, time is mainly spent inside BLAS calls. However, things
# could probably be improved for smaller matrices.

def _st(A, b, trans=0):
    """
    solve triangular system "tril(A) @ x = b", returning x
    
    if trans==1, solve "tril(A).T @ x = b" instead.
    """
    if b.size == 0:
        return b
    else:
        return _solve_triangular(A, b, trans=trans, lower=True)

def _Phi(A):
    """Return lower-triangle of matrix and halve the diagonal"""
    A = tril(A)
    A[np.diag_indices_from(A)] *= 0.5
    return A

def _chol_symbolic_fwd(L, Sigma_dot):
    """
    Forwards-mode differentiation through the Cholesky decomposition
    
    This version uses a "one-line" symbolic expression to return L_dot
    where "_dot" means sensitivities in forwards-mode differentiation,
    and Sigma = L @ L.T.
    """
    # invL = inv(L)
    # return L @ _Phi(invL @ Sigma_dot @ invL.T)
    return np.dot(L, _Phi(_st(L, _st(L, Sigma_dot.T).T)))

def _chol_symbolic_rev(L, Lbar):
    """
    Reverse-mode differentiation through the Cholesky decomposition
    
    This version uses a short symbolic expression to return
    tril(Sigma_bar) where "_bar" means sensitivities in reverse-mode
    differentiation, and Sigma = L @ L.T.
    """
    P = _Phi(np.dot(L.T, Lbar))
    #invL = inv(L)
    #return _Phi(invL.T @ (P + P.T) @ invL)
    return _Phi(_st(L, _st(L, (P + P.T), 1).T, 1))

def _level2partition(A, j):
    """Return views into A used by the unblocked algorithms"""
    # diagonal element d is A[j,j]
    # we access [j, j:j+1] to get a view instead of a copy.
    rr = A[j, :j]      # row  
    dd = A[j, j:j+1]   # scalar on diagonal   /      \
    B = A[j+1:, :j]    # Block in corner      | r d  |
    cc = A[j+1:, j]    # column               \ B c  /
    return rr, dd, B, cc

def _chol_unblocked(A, inplace=False):
    """
    Cholesky decomposition, mirroring LAPACK's DPOTF2
    
    Intended to illustrate the algorithm only. Use a Cholesky routine
    from numpy or scipy instead.
    """
    if not inplace:
        A = A.copy()
    for j in range(A.shape[0]):
        rr, dd, B, cc = _level2partition(A, j)
        dd[:] = np.sqrt(dd - np.dot(rr, rr))
        cc[:] = (cc - np.dot(B, rr)) / dd
    return A

def _chol_unblocked_fwd(L, Adot, inplace=False):
    """
    Forwards-mode differentiation through the Cholesky decomposition
    
    Obtain L_dot from Sigma_dot, where "_dot" means sensitivities in
    forwards-mode differentiation, and Sigma = L @ L.T.

    This version uses an unblocked algorithm to update sensitivities
    Adot in place. tril(Adot) should start containing Sigma_dot, and
    will end containing the L_dot. The upper triangular part of Adot
    is untouched, so take tril(Adot) at the end if triu(Adot,1) did
    not start out filled with zeros.

    If inplace=False, a copy of Adot is modified instead of the
    original. The Abar that was modified is returned.
    """
    if not inplace:
        Adot = Adot.copy()
    for j in range(L.shape[0]):
        rr, dd, B, cc = _level2partition(L, j)
        rdot, ddot, Bdot, cdot = _level2partition(Adot, j)
        ddot[:] = (ddot/2 - np.dot(rr, rdot)) / dd
        cdot[:] = (cdot - np.dot(Bdot, rr) - np.dot(B, rdot) - cc*ddot) / dd
    return Adot

def _chol_unblocked_rev(L, Abar, inplace=False):
    """
    Reverse-mode differentiation through the Cholesky decomposition
    
    Obtain tril(Sigma_bar) from L_bar, where "_bar" means sensitivities
    in reverse-mode differentiation, and Sigma = L @ L.T.

    This version uses an unblocked algorithm to update sensitivities
    Abar in place. tril(Abar) should start containing L_bar, and will
    end containing the tril(Sigma_bar). The upper triangular part of
    Adot is untouched, so take tril(Abar) at the end if triu(Abar,1)
    did not start out filled with zeros. Alternatively, (tril(Abar) +
    tril(Abar).T) will give the symmetric, redundant matrix of
    sensitivities.
    
    If inplace=False, a copy of Abar is modified instead of the
    original. The Abar that was modified is returned.
    """
    if not inplace:
        Abar = Abar.copy()
    for j in range(L.shape[0] - 1, -1, -1): # N-1,N-2,...,1,0
        rr, dd, B, cc = _level2partition(L, j)
        rbar, dbar, Bbar, cbar = _level2partition(Abar, j)
        dbar -= np.dot(cc, cbar) / dd
        dbar /= dd # / These two lines could be
        cbar /= dd # \ done in one operation
        rbar -= dbar*rr         # / These two lines could be done
        rbar -= np.dot(cbar, B) # \ with one matrix multiply
        Bbar -= np.dot(cbar[:,None], rr[None,:])
        dbar /= 2
    return Abar

def _chol_unblocked_fortran_rev(L, L_bar, inplace=False):
    """Like _chol_unblocked_rev but calling Fortran code."""
    # Checking important. Fortran code could segfault with bad input!
    N = L.shape[0]
    assert(L.shape == (N, N))
    assert(L.shape == L_bar.shape)
    info = np.array(0)
    if inplace:
        A_bar = L_bar
    else:
        A_bar = L_bar.copy()
    _larmpack.dpo2ft("L", N, L, A_bar, info)
    assert(info == 0)
    return A_bar

def _level3partition(A, j, k):
    """Return views into A used by the blocked algorithms"""
    # Top left corner of diagonal block is [j,j]
    # Block size is NB = (k-j)
    R = A[j:k, :j]     # Row block                     /      \
    D = A[j:k, j:k]    # triangular block on Diagonal  |      |
    B = A[k:, :j]      # Big corner block              | R D  |
    C = A[k:, j:k]     # Column block                  \ B C  /
    return R, D, B, C

def _chol_blocked(A, NB=256, inplace=False):
    """Cholesky decomposition, mirroring LAPACK's DPOTRF
    
    Intended to illustrate the algorithm only. Use a Cholesky routine
    from numpy or scipy instead."""
    if not inplace:
        A = A.copy()
    for j in range(0, A.shape[0], NB):
        k = min(N, j + NB)
        R, D, B, C = _level3partition(A, j, k)
        D -= tril(np.dot(R, R.T))
        _chol_unblocked(D, inplace=True)
        C -= np.dot(B, R.T)
        #C[:] = C @ inv(tril(D)).T
        C[:] = _st(D, C.T).T
    return A

def _chol_blocked_fwd(L, Adot, NB=256, inplace=False):
    """
    Forwards-mode differentiation through the Cholesky decomposition
    
    Obtain L_dot from Sigma_dot, where "_dot" means sensitivities in
    forwards-mode differentiation, and Sigma = L @ L.T.

    This version uses a blocked algorithm to update sensitivities Adot
    in place. tril(Adot) should start containing Sigma_dot, and will
    end containing the L_dot. Take tril() of the answer if
    triu(Adot,1) did not start out filled with zeros. Unlike the
    unblocked routine, if the upper triangular part of Adot started
    with non-zero values, some of these will be overwritten.

    If inplace=False, a copy of Adot is modified instead of the
    original. The Abar that was modified is returned.
    """
    if not inplace:
        Adot = Adot.copy()
    for j in range(0, L.shape[0], NB):
        k = min(N, j + NB)
        R, D, B, C = _level3partition(L, j, k)
        Rdot, Ddot, Bdot, Cdot = _level3partition(Adot, j, k)
        Ddot[:] = tril(Ddot) - tril(np.dot(Rdot, R.T) + np.dot(R, Rdot.T))
        #chol_unblocked_fwd(D, Ddot, inplace=True) # slow in Python
        Ddot[:] = _chol_symbolic_fwd(D, Ddot + tril(Ddot, -1).T)
        Cdot -= (np.dot(Bdot, R.T) + np.dot(B, Rdot.T))
        #Cdot[:] = (Cdot - C@Ddot.T) @ inv(tril(D)).T
        Cdot[:] = _st(D, Cdot.T - np.dot(Ddot, C.T)).T
    return Adot

def _chol_blocked_rev(L, Abar, NB=256, inplace=False):
    """
    Reverse-mode differentiation through the Cholesky decomposition
    
    Obtain tril(Sigma_bar) from L_bar, where "_bar" means sensitivities
    in reverse-mode differentiation, and Sigma = L @ L.T.

    This version uses a blocked algorithm to update sensitivities Abar
    in place. tril(Abar) should start containing L_bar, and will end
    containing the tril(Sigma_bar). Take tril(Abar) at the end if
    triu(Abar,1) did not start out filled with zeros. Alternatively,
    (tril(Abar) + tril(Abar).T) will give the symmetric, redundant
    matrix of sensitivities.
    
    Unlike the unblocked routine, if the upper triangular part of Abar
    started with non-zero values, some of these will be overwritten.

    If inplace=False, a copy of Abar is modified instead of the
    original. The Abar that was modified is returned.
    """
    if not inplace:
        Abar = Abar.copy()
    for k in range(L.shape[0], -1, -NB):
        j = max(0, k - NB)
        R, D, B, C = _level3partition(L, j, k)
        Rbar, Dbar, Bbar, Cbar = _level3partition(Abar, j, k)
        #Cbar[:] = Cbar @ inv(tril(D))
        Cbar[:] = _st(D, Cbar.T, trans=1).T
        Bbar -= np.dot(Cbar, R)
        Dbar[:] = tril(Dbar) - tril(np.dot(Cbar.T, C))
        #chol_unblocked_rev(D, Dbar, inplace=True) # slow in Python
        Dbar[:] = _chol_symbolic_rev(D, Dbar)
        Rbar -= (np.dot(Cbar.T, B) + np.dot(Dbar + Dbar.T, R))
    return Abar

def _chol_blocked_fortran_rev(L, L_bar, NB=None, inplace=False):
    """Like chol_blocked_rev but calling Fortran code. NB is ignored."""
    # Checking important. Fortran code could segfault with bad input!
    N = L.shape[0]
    assert(L.shape == (N, N))
    assert(L.shape == L_bar.shape)
    info = np.array(0)
    if inplace:
        A_bar = L_bar
    else:
        A_bar = L_bar.copy()
    _larmpack.dpofrt("L", N, L, A_bar, info)
    assert(info == 0)
    return A_bar

# The public functions that will be listed in the module documentation:
def chol_fwd(L, Adot, NB=256, inplace=False):
    """
    Forwards-mode differentiation through the Cholesky decomposition
    
    Obtain L_dot = dL/dx from Sigma_dot = dSigma/dx, where Sigma = L @ L.T.

    The input tril(Adot) should contain Sigma_dot. The lower-triangle of
    the answer will contain L_dot. Take tril() of the answer if input
    triu(Adot,1) was not filled with zeros.
    
    If inplace=True the answer is placed into the lower triangle of Adot
    as well as being returned. If inplace=False, a copy of Adot is modified
    and returned, and the original array is not touched.

    The optimal block-size NB depends on the machine and size of L, but
    should not change the answer (beyond the usual round-off errors).
    """
    return _chol_blocked_fwd(L, Adot, NB, inplace)
if FORTRAN_COMPILED:
    _chol_rev = _chol_blocked_fortran_rev
else:
    _chol_rev = _chol_blocked_rev
def chol_rev(L, Abar, NB=256, inplace=False):
    """
    Reverse-mode differentiation through the Cholesky decomposition
    
    Obtain Sigma_bar = df/dSigma from L_bar = df/dL, where Sigma = L @ L.T.

    The input tril(Abar) should contain L_bar. The lower-triangle of
    the answer will contain Sigma_bar. Take tril() of the answer if input
    triu(Abar,1) was not filled with zeros.
    
    If inplace=True the answer is placed into the lower triangle of Abar
    as well as being returned. If inplace=False, a copy of Abar is modified
    and returned, and the original array is not touched.

    The optimal block-size NB depends on the machine and size of L, but
    should not change the answer (beyond the usual round-off errors).
    If the Fortran version is compiled, any user-specified NB will be
    ignored, and the block size will be chosen by LAPACK.
    """
    return _chol_rev(L, Abar, NB, inplace)
 

# Testing code follows

def _trace_dot(A, B):
    """_trace_dot(A, B) = trace(A @ B) = A.ravel() @ B. ravel()"""
    return np.dot(A.ravel(), B.ravel())

def _testme(N):
    """Exercise each function using NxN matrices"""
    import scipy as sp
    from time import time
    if N > 1:
        Sigma = np.cov(sp.randn(N, 2*N))
        Sigma_dot = np.cov(sp.randn(N, 2*N))
    elif N == 1:
        Sigma = np.array([[sp.rand()]])
        Sigma_dot = np.array([[sp.rand()]])
    else:
        assert(False)
    tic = time()
    L = np.linalg.cholesky(Sigma)
    toc = time() - tic
    print('Running np.linalg.cholesky:')
    print('   Time taken: %0.4f s' % toc)
    tic = time()
    L_ub = tril(_chol_unblocked(Sigma))
    toc = time() - tic
    print('Unblocked chol works: %r'
            % np.all(np.isclose(L, L_ub)))
    print('   Time taken: %0.4f s' % toc)
    tic = time()
    L_bl = tril(_chol_blocked(Sigma))
    toc = time() - tic
    print('Blocked chol works: %r'
            % np.all(np.isclose(L, L_bl)))
    print('   Time taken: %0.4f s' % toc)
    tic = time()
    Ldot = _chol_symbolic_fwd(L, Sigma_dot)
    toc = time() - tic
    hh = 1e-5 # finite-difference step-size
    L2 = np.linalg.cholesky(Sigma + Sigma_dot*hh/2)
    L1 = np.linalg.cholesky(Sigma - Sigma_dot*hh/2)
    Ldot_fd = (L2 - L1) / hh
    print('Symbolic chol_fwd works: %r'
            % np.all(np.isclose(Ldot, Ldot_fd)))
    print('   Time taken: %0.4f s' % toc)
    tic = time()
    Ldot_ub = tril(_chol_unblocked_fwd(L, Sigma_dot))
    toc = time() - tic
    print('Unblocked chol_fwd works: %r'
            % np.all(np.isclose(Ldot, Ldot_ub)))
    print('   Time taken: %0.4f s' % toc)
    tic = time()
    Ldot_bl = tril(_chol_blocked_fwd(L, Sigma_dot))
    toc = time() - tic
    print('Blocked chol_fwd works: %r'
            % np.all(np.isclose(Ldot, Ldot_bl)))
    print('   Time taken: %0.4f s' % toc)
    Lbar = tril(sp.randn(N, N))
    tic = time()
    Sigma_bar = _chol_symbolic_rev(L, Lbar)
    toc = time() - tic
    Delta1 = _trace_dot(Lbar, Ldot)
    Delta2 = _trace_dot(Sigma_bar, Sigma_dot)
    print('Symbolic chol_rev works: %r'
            % np.all(np.isclose(Delta1, Delta2)))
    print('   Time taken: %0.4f s' % toc)
    tic = time()
    Sigma_bar_ub = _chol_unblocked_rev(L, Lbar)
    toc = time() - tic
    Delta3 = _trace_dot(Sigma_bar_ub, Sigma_dot)
    print('Unblocked chol_rev works: %r'
            % np.all(np.isclose(Delta1, Delta3)))
    print('   Time taken: %0.4f s' % toc)
    tic = time()
    Sigma_bar_bl = _chol_blocked_rev(L, Lbar)
    toc = time() - tic
    Delta4 = _trace_dot(Sigma_bar_bl, Sigma_dot)
    print('Blocked chol_rev works: %r'
            % np.all(np.isclose(Delta1, Delta4)))
    print('   Time taken: %0.4f s' % toc)
    if FORTRAN_COMPILED:
        tic = time()
        Sigma_bar_f = chol_rev(L, Lbar)
        toc = time() - tic
        Delta5 = _trace_dot(Sigma_bar_f, Sigma_dot)
        print('Fortran chol_rev works: %r'
                % np.all(np.isclose(Delta1, Delta5)))
        print('   Time taken: %0.4f s' % toc)
        tic = time()
        Sigma_bar_fub = _chol_unblocked_fortran_rev(L, Lbar)
        toc = time() - tic
        Delta6 = _trace_dot(Sigma_bar_fub, Sigma_dot)
        print('Fortran unblocked chol_rev works: %r'
                % np.all(np.isclose(Delta1, Delta6)))
        print('   Time taken: %0.4f s' % toc)
    else:
        print('Fortran chol_rev not compiled.')

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    else:
        N = 500
    _testme(N)

