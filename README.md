Cholesky Reverse
================

This mini project provides a routine that let's you differentiate
expressions involving Cholesky decompositions, using 'reverse-mode
differentiation' or 'backpropagation'.

The LAPACK FORTRAN 77 code `dpotrf.f` implements the Cholesky decomposition
using blocked level-3 BLAS routines so it's fast. This implementation is
widely used, such as by Matlab/Octave, R, and SciPy. The file `dpofrt.f` in
this repository is a new companion routine, which takes derivatives with
respect to a Cholesky decomposition from `dpotrf.f` and replaces them with
derivatives with respect to elements of the original positive definite
input matrix.

If you wish to use BLAS and LAPACK binaries that use 64-bit indexes, make sure
to use the relevant compiler flag (possibly `-fdefault-integer-8` or `-i8`), or
use the derived file `dpotrf_ilp64.f`. See `make_ilp64.sh` for more details.
It's important to check the size of integers used by BLAS and LAPACK: they are
traditionally 32-bit even on 64-bit machines, but builds with 64-bit integers
are becoming more common.


How it can be used
------------------

The matlab and python sub-directories demonstrate how to compile this
routine for Matlab/Octave, call it, and check consistency. The matlab
directory has a toy Gaussian process demo, to show how the new routines
could be used within a larger gradient computation. There is also a slow
Matlab version of the Fortran code, which may be easier to read than the
Fortran for anyone wishing to port the algorithm.

What follows is a high level description of how this routine could be used
in general:

Assume you wish to differentiate a scalar function `f(chol(A(W)))` with
respect to an array of input values `W`. `A` is positive-definite matrix
that depends on the inputs, and `L = chol(A)` is its lower-diagonal
Cholesky decomposition. If you can compute derivatives `L_bar = df/dL`, the
code provided here will convert them into `A_bar = df/dA`. These
derivatives should in turn be passed on to reverse-mode/backpropagation
code that can compute `W_bar = df/dW`, the final required result.

(NB don't compute the 4-dimensional array with elements `dA_{ij}/dW_{kl}`
-- use backpropagation to compute `df/dW_{kl}` from `A_bar = df/dA` without
producing large intermediate objects! A lot of Gaussian process code in the
wild does create large intermediate objects, or has inefficient loops that
forward propagate some of the gradient computation. The demo in the Matlab
directory includes a comparison of the different approaches.)

We're only providing the core primitive for pushing derivatives through the
Cholesky decomposition. You'll have to do the rest of the differentiation
by hand, or provide this routine as a useful primitive to an automatic
differentiation package.

A great review of how to differentiate linear algebra is:

[An extended collection of matrix derivative results for forward and
reverse mode automatic differentiation, Mike B Giles (2008)](https://people.maths.ox.ac.uk/gilesm/files/NA-08-01.pdf)

Giles provides pseudo-code and simple Matlab code for pushing derivatives
through the Cholesky decomposition. However, this algorithm doesn't use
block matrix operations like LAPACK routines do, so is comparatively slow.


This code is deliberately archaic
---------------------------------

LAPACK sticks to FORTRAN 77 and 6-letter function names. For the moment,
we've been masochistic and done the same. The last two or three letters of
a LAPACK routine represent the operation, and we have reversed these
letters to indicate a reverse operation.

The idea is that lots of projects already use LAPACK. If our routines are
just like those, in principle it should be possible to build them with
existing tool-chains for many use cases. However, many people use binary
releases of LAPACK, and working out how to compile new routines against
BLAS and LAPACK is not easy for everyone. In the long term, hopefully the
importance of performant reverse mode primitives will become more widely
appreciated. Then routines like this one might ship pre-compiled with
binary distributions of projects like Octave and SciPy.


Status and open issues
----------------------

I've thrown this up in the hope it will be useful. It is not heavily tested
and there are some known rough edges. I should probably make these github
issues...

I've only built it on my own linux machine. Tweaks to the Matlab/Octave and
Python demos that help them work on other platforms would be appreciated.

Currently the routine only lets you specify a lower-triangular Cholesky
decomposition. Fixing that to allow upper-triangular matrices would be a
quick job.

The unblocked routine `dpo2ft.f` is neat in that it only touches one
triangle of the input/output matrix of derivatives. The blocked routine
`dpofrt.f` has a couple of messy parts, where it creates unneeded values in
the other triangle of the matrix, and then overwrites them with zeros at
the end. This code seems needlessly inefficient, but I didn't immediately
see how to make it better, given the standard BLAS routines that exist.
However, I'm sure these parts could be improved.

Is there any demand to allow the Cholesky decomposition to be in one
triangle and the gradients to be in another. That is, to have two `UPLO`
arguments one for `A` and one for `C`? If the messy details around the
diagonal (previous paragraph) could be resolved, maybe the Cholesky
decomposition and gradients could optionally be packed into one `Nx(N+1)`
array, saving memory. Would it be worth the hassle?

In the end I think we should have "LARMPACK", a library containing the
reverse-mode functions for everything in LAPACK. This project is an
experiment, starting with one LAPACK routine, and seeing how that goes
before working on the rest. I don't think differentiating the rest of the
routines should be hard (in principle it can be done automatically,
although some manual work can make the results neater).

