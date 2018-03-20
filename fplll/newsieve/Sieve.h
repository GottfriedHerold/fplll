// TODO: Documentation.

#ifndef SIEVE_H
#define SIEVE_H

// temporary, for as long as we do not support multithreaded variant.
#define SIEVE_GAUSS_SINGLE_THREADED

#define INCLUDE_SIEVE_TEMPL_FROM_SIEVE_H
#include "Sieve_templ.h"
#undef INCLUDE_SIEVE_TEMPL_FROM_SIEVE_H

#endif  // main include guard for Sieve.h
