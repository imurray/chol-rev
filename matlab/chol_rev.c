#include "mex.h"
#include <string.h>

#if __STDC_VERSION__  >= 199901L
/*  C99 says stdint.h header exists.
    Hopefully will work on future systems for some time. */
#   include <stdint.h>
#else
/*  Otherwise the following should work on most current systems. */
#   ifndef int64_t
#       ifdef _WIN32
#           define int64_t LONGLONG
#       else
#           define int64_t long long
#       endif
#   endif
#   ifndef int32_t
#       define int32_t int
#   endif
#endif

#ifdef BLAS64INT
#   define INTEGER int64_t
#else
#   define INTEGER int32_t
#endif

/* fix for systems known not to do fortran name-mangling */
#if defined(_WIN32) || defined(__hpux)
#   define dpofrt_ dpofrt
#endif

extern void dpofrt_(
        const char *uplo, const INTEGER *n, double *a, const INTEGER *lda,
        double *c, const INTEGER *ldc, const INTEGER *info);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double* X;
    INTEGER nn, info;

    if (nrhs != 2 || nlhs > 1)
        mexErrMsgTxt("Usage: X_bar = chol_rev_rf(L, L_bar)");
    nn = mxGetN(prhs[0]);
    if ((mxGetM(prhs[0]) != nn) || (mxGetNumberOfDimensions(prhs[0]) != 2))
        mexErrMsgTxt("Error: L must be a square matrix");
    if ((mxGetN(prhs[1]) != nn) || (mxGetM(prhs[1]) != nn) || (mxGetNumberOfDimensions(prhs[1]) != 2))
        mexErrMsgTxt("Error: L_bar must be a square matrix matching L");

    plhs[0] = mxCreateDoubleMatrix(nn, nn, mxREAL);
    X = mxGetPr(plhs[0]);

    if (nn == 0)
        return;

    memcpy(X, mxGetPr(prhs[1]), nn*nn*sizeof(double));
    dpofrt_("L", &nn, mxGetPr(prhs[0]), &nn, X, &nn, &info);
    if (info != 0)
        mexErrMsgTxt("Error: unable to push sensitivity back through Cholesky decomposition.");
}

