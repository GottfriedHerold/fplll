/**
  Include this file to use the k-tuple GaussSieve for arbitrary template parameters.

  This file just wraps around Sieve_templ.h
*/

/**
The macros
SIEVE_GAUSS_SINGLE_THREADED
and
SIEVE_GAUSS_MULTI_THREADED
determine whether compile for single-threaded or multi-threaded mode or both.
If either macro is defined, we support the corresponding mode.
If neither is set, we fall back to a default (which is "both")
*/

#ifndef SIEVE_GAUSS_H
#define SIEVE_GAUSS_H

// temporary, for as long as we do not support multithreaded variant.
#define SIEVE_GAUSS_SINGLE_THREADED

#define INCLUDE_SIEVE_TEMPL_FROM_SIEVE_GAUSS_H
#include "Sieve_templ.h"
#undef INCLUDE_SIEVE_TEMPL_FROM_SIEVE_GAUSS_H

#endif  // main include guard for SieveGauss.h
