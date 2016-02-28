#!/usr/bin/env python

import subprocess, os.path

# In Anaconda with python3, f2py3 was available and f2py was not.
# On Ubuntu 14.04, I saw that f2py2.7 and f2py3.4 were available.
# Try to find most specific f2py version:
import sys, distutils.spawn
full_ver = "%d.%d" % sys.version_info[:2] # e.g., "2.7" or "3.5"
f2py = distutils.spawn.find_executable('f2py' + full_ver)
if not f2py:
    ver = str(sys.version_info[0]) # e.g., "2" or "3"
    f2py = distutils.spawn.find_executable('f2py' + ver)
    if not f2py:
        f2py = 'f2py'

cmd = [f2py, '-c', '-m', 'larmpack', 'larmpack.pyf',
        os.path.join('..', 'dpofrt.f'), os.path.join('..', 'dpo2ft.f'),
        '--link-lapack_opt']
val = subprocess.call(cmd)
if val != 0:
    # Anaconda can be configured to use Intel's commercial MKL library,
    # but doesn't ship with a complete LAPACK library. However, its "BLAS"
    # libraries have all of the LAPACK routines we need in them, so try
    # asking just for those if compilation failed.
    cmd[-1] = '--link-blas_opt'
    val = subprocess.call(cmd)

try:
    import larmpack
    print('\nIt appears the Fortran-based module is now available.\n')
except:
    print("\nFortran-based module isn't available." + \
            "Compilation must have failed.\n")

