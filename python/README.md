Python Cholesky Differentiation
===============================

This directory contains Python routines to push derivatives through the
Cholesky decomposition. See the README.md in the parent directory for a
little more information, and the [note on arXiv](http://arxiv.org/abs/1602.07527)
for a full explanation. The parent directory also links to some
Matlab/Octave code, which contains a small Gaussian process demo.

A simple demo can be run with:
```
    python chol_rev_demo.py
```
A small error should be reported (e.g., `1e-8` or smaller), indicating that
the reverse-mode differentiation routine is working properly.

The routines described in the [note on arXiv](http://arxiv.org/abs/1602.07527)
are available in `chol_diff.py`. You can run
```
    python chol_diff.py
```
to test these routines, and get some simplistic but illustrative
timings. You can try different matrix sizes with an optional argument:
```
    python chol_diff.py 4000
```
The overhead of using Python becomes negligible for large matrices, which
is also where the advantage of blocked versions is largest. For small
matrices it may be better to use compiled code.


Linking with Fortran code
-------------------------

You can compile a Fortran version of the reverse-mode differentiation
routine with
```
    python compile.py
```
and try the demos above again. See the parent README.md for more
information about the Fortran code. The Fortran version isn't always
faster! For large matrices the algorithm matters most: the unblocked
Fortran version is much slower than the blocked Python+Numpy version.

You need f2py (which comes with NumPy these days), python development
headers (the python-dev or python3-dev package in Debian/Ubuntu), and the
ability to compile Fortran and C, including BLAS and LAPACK routines as
used by your install of NumPy+SciPy. (Although see the caveat in compile.py
-- things sometimes manage to work if you don't have complete LAPACK
libraries.)

I haven't tested things much in Python yet, beyond the quick demo. I
also haven't done anything to check the performance overhead of the f2py
wrapping. It would probably be better to write code that deals directly
with the C-style matrices that NumPy seems to prefer.

I've only built it on my own Linux machine with the system's Python and
with Anaconda. My experience with distributing Python projects is limited.

Unlike the Matlab wrapper, I don't yet do anything to try to check
whether BLAS and LAPACK use 64-bit integers, or set the compiler options
appropriately. When (in future?) NumPy/SciPy uses BLAS/LAPACK with
64-bit integers, we should detect that and give f2py the `*_ilp64.f`
routines (unless it sets compiler options like `-fdefault-integer-8` or
`-i8` in this situation automatically?).

