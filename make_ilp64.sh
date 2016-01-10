#!/bin/sh

set -e

# This script created the *_ilp64.f files. Although they're derived files, I've
# committed them to version control anyway to simplify things for end users.
# The rest of this comment explains the need for these files.
#
# Some builds of BLAS/LAPACK use 64-bit integer indexes (originally not the
# case, even on 64-bit machines). To link with those, we must use 64-bit
# integers too.
#
# One option is to set appropriate compiler flags:
#     gfortran -fdefault-integer-8
#     ifort -i8
#     pgf77 -i8
#     crayftn -s integer64
#     ... ?
# But it seems a pain to detect the compiler and set these flags correctly in
# all situations. Moreover, there isn't always an appropriate compiler flag. For
# example, the NAG compiler seems to only have an option that doubles the
# integer and floating point precision at the same time, sometimes creating
# 128-bit floats.
#
# An alternative is to create versions of the routines that use 8-byte integers.
# The first step is to replace the types of variables and parameters. That's
# what the sed below dow. The types of all arguments to functions must also be
# clear. That's why I have used the typed parameter "ONE_I" instead of "1" in
# places. Without that, this script wouldn't be enough, and the code would
# probably crash.
#
# For the moment, I'm sticking to the old crufty style of the LAPACK routines
# that I'm following, and not using modern Fortran. INTEGER*8 isn't in the
# Fortran standard, and might not necessarily always mean 64-bits, but seems to
# work widely. Maybe in the future I should just shift to using a recent Fortran
# standard however.

for a in dpo2ft.f dpofrt.f ; do
    sed 's/INTEGER  \? \?/INTEGER*8 /' "$a" > "${a%.f}"_ilp64.f
done

