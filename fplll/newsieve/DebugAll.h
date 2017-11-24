/**
  This file collects and specifies all debugging preprocessor macros that
  enable checks.
*/

#ifndef DEBUG_ALL_H
#define DEBUG_ALL_H

// clang-format kills the indenting for #defines
// clang-format off

// This causes generic implementations of lattice points to produce verbos output.
// Used to track which functions get called.
// #define DEBUG_SIEVE_LP_TRACEGENERIC

// enables all silent checks
#define DEBUG_SIEVE_SILENT_ALL

// verbos output for Lazy
#define DEBUG_SIEVE_LAZY_ALL


#ifdef DEBUG_SIEVE_SILENT_ALL
// tests that LP classes were initialized prior to using.
  #define DEBUG_SIEVE_LP_INIT

// tests that operations on LP have matching dimensions
  #define DEBUG_SIEVE_LP_MATCHDIM

// tests that the RNG engine is called with valid parameters
  #define DEBUG_SIEVE_MTPRNG

// Tests that accesses to the mu and g-matrices only access entries that are potentially non-zero.
// While code may reasonably access such entries that are always zero, we give an error, because
// it's more likely to be a bug rather than intended.
  #define DEBUG_SIEVE_LOWERTRIANGULAR_MUG

// Tests whether some lazy evaluation actually actually loses some levels in the process.
// This may be caused by missing features of an approximation.
//  #define DEBUG_SIEVE_ALL_APPROX_AT_SAME_LEVEL

#endif  // DEBUG_SIEVE_SILENT_ALL

// enables modules which use a pointer to Sieve<...> to be used standalone.
// Because this actually changes semantics by disabling some features,
// do NOT use this outside of testing, because it actually breaks things.

#ifdef DEBUG_SIEVE_STANDALONE_MODULES_ALL

  #define DEBUG_SIEVE_STANDALONE_SAMPLER

  #define DEBUG_SIEVE_STANDALONE_QUEUE

#endif



#ifdef DEBUG_SIEVE_LAZY_ALL
  #define DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
  #define DEBUG_SIEVE_LAZY_TRACE_EVALS
  #define DEBUG_SIEVE_LAZY_TRACE_APPROX
#endif


// Outputs late-initializations of modules from GaussSieve
//#define DEBUG_SIEVE_INITIALIZATIONS

#ifdef DEBUG_SIEVE_INITIALIZATIONS
  #define DEBUG_SIEVE_TRACEINITIATLIZATIONS(STR) std::cout << STR << std::endl << std::flush;
#else
  #define DEBUG_SIEVE_TRACEINITIATLIZATIONS(STR)
#endif


// verbosely trace calls of generic functions for LPs.
// Disabled by default:


#ifdef DEBUG_SIEVE_LP_TRACEGENERIC
  #define DEBUG_TRACEGENERIC(STR) std::cout << STR << std::endl << std::flush;
#else  // DEBUG_SIEVE_LP_TRACEGENERIC
  #define DEBUG_TRACEGENERIC(STR)
#endif  // DEBUG_SIEVE_LP_TRACEGENERIC

// clang-format on

#endif
