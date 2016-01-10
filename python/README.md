Python Cholesky Reverse
=======================

This is a first attempt at providing a routine to push derivatives through
a Cholesky decomposition in Python. It uses the Fortran routines from the
parent directory; see the README.md there for more information.

On my machine I can run:
```
    python compile.py
    python chol_rev_demo.py
```
A small error is reported (e.g., `1e-8` or smaller), indicating the routine
is working properly.

You need f2py (which comes with numpy these days), and the ability to
compile Fortran and C, including BLAS and LAPACK routines as used by your
install of numpy+scipy. (Although see the caveat in compile.py -- things
sometimes manage to work if you don't have complete LAPACK libraries.)


Status and open issues
----------------------

See the parent directory for more information about the routines.

I haven't tested things much in Python yet, beyond the quick demo. I also
haven't done anything to check the performance overhead of the f2py
wrapping. It would probably be better to write code that deals directly
with the C-style matrices that NumPy seems to prefer.

I've only built it on my own Linux machine with the system's Python and
with Anaconda. I'm not yet experienced with distributing Python projects.

Unlike the Matlab wrapper, I currently don't do anything to try to check
whether BLAS and LAPACK use 64-bit integers, or set the compiler options
appropriately. When (in future?) NumPy/SciPy uses BLAS/LAPACK with 64-bit
integers, we should detect that and give f2py the `*_ilp64.f` routines
(unless it sets compiler options like `-fdefault-integer-8` or `-i8` in
this situation?).

