#!/usr/bin/env python

import os, os.path

cmd = 'f2py -c -m larmpack larmpack.pyf'
cmd += ' ' + os.path.join('..', 'dpofrt.f')
cmd += ' ' + os.path.join('..', 'dpo2ft.f')
val = os.system(cmd + ' --link-lapack_opt')
if val != 0:
    # Anaconda's commercial "accelerate" numpy+scipy package uses Intel's MKL,
    # but doesn't ship with a complete LAPACK library. However, its "BLAS"
    # libraries have all of the LAPACK routines we need in them, so try
    # asking just for those if compilation failed.
    val = os.system(cmd + ' --link-blas_opt')

