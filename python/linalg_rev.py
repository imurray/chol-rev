import larmpack as _larmpack
import numpy as _np

def chol_rev(L, Lbar, inplace=False):
    N = L.shape[0]
    assert(L.shape == (N, N))
    assert(L.shape == Lbar.shape)
    info = _np.array(0)
    if inplace:
        _larmpack.dpofrt("L", N, L, Lbar, info)
        assert(info == 0)
    else:
        Abar = Lbar.copy()
        _larmpack.dpofrt("L", N, L, Abar, info)
        assert(info == 0)
        return Abar

