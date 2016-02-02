Matlab/Octave Cholesky Reverse
==============================

This is a first attempt at providing a Matlab/Octave routine to push
derivatives backwards through the Cholesky decomposition. See the README.md
in the parent directory for more information.

On my machine I can type within Octave or Matlab:
```
    compile
    chol_rev_demo
```
A small error is reported (e.g., `1e-8` or smaller), indicating the routine
is working properly.

If you omit the compilation step, it should still work, but you'll be
warned that the core `chol_rev` routine will probably be slower.

The fast compiled version uses the Fortran routines from the parent
directory; see the README.md there both for more information about the
Fortran code and for the general idea. I've only tested compilation on
Linux. While I've tried to make it more general, small tweaks are probably
still required for other platforms. If you don't have a Fortran compiler,
it's not going to work!

The two pure m-file routines `chol_rev_unblocked.m` and
`chol_rev_blocked.m` were written as prototypes for the Fortran routines in
the parent directory. These may be also be a good starting point for ports
to other languages. They resulted from differentiating the blocked and
unblocked LAPACK routines for the Cholesky decomposition.

The `chol_rev.m` routine is the same as `chol_rev_blocked.m`, except with
the call to the unblocked routine replaced with an analytically derived
expression that's faster to evaluate in Matlab or Octave.

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

