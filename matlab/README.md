Matlab/Octave Cholesky Reverse
==============================

This is a first attempt at providing a routine to push derivatives through
a Cholesky decomposition in Matlab/Octave. It uses the Fortran routines
from the parent directory; see the README.md there for more information.

On my machine I can type within Octave or Matlab:
```
    compile
    chol_rev_demo
```
A small error is reported (e.g., `1e-8` or smaller), indicating the routine
is working properly.

`compile.m` builds a `chol_rev` mex routine for Matlab or Octave. Only
tested on Linux. While I've tried to make it more general, small tweaks are
probably still required for other platforms. If you don't have a Fortran
compiler, it's not going to work!

There are also two pure m-file implementations: `chol_rev_unblocked.m` and
`chol_rev.m` a blocked algorithm, which requires `chol_rev_unblocked.m` as
a helper routine. If you don't care about speed, you can simply use
`chol_rev.m`. This m-file version was actually written first, and the
Fortran code was ported from it. Similarly the m-file code might be a good
starting point for ports to other languages.

A Gaussian process demo `chol_rev_gp_demo.m` shows that this routine can
give competitive (sometimes faster) times when used within a larger
gradient computation. In other applications, for example if wanting to run
[Hamiltonian Monte Carlo](http://arxiv.org/abs/1206.1901) on
[a parameterization of a model using a Cholesky
decomposition](http://homepages.inf.ed.ac.uk/imurray2/pub/10hypers/), it's
hard to see how to implement the method at all without a `chol_rev` routine.


32 bit vs 64 bit indexing
-------------------------

On 64-bit platforms, MATLAB has used BLAS and LAPACK with 64-bit integers
from version 7.8 (R2009a). Octave can be built with 64-bit integers using
the configure option `--enable-64`, although it normally isn't, even on
64-bit platforms (circa 2016). Running mex libraries built for the
wrong-sized integers will lead to crashes. The compile script attempts to
do the right thing. Please be careful if compiling yourself, or modifying
`compile.m`.

