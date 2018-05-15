// This file intentionally has no (traditional) include guards!

// clang-format off
// turned off for this file, which is only preprocessor stuff
// because we want to have proper indentation for (nested) preprocessor #if's in this file.
// Ideally, we would want to use different clang-settings for this file.
// Unfortunately, there is no "clang-format option value" comment in the file that overrides
// clang-format's options on a per-file basis; I hope such a feature will be added at some point.

/**
  This file contains some preprocessor magic to create the files
  SieveGauss.h
  Sieve.h
  Sieve.cpp
  (each of the above 3 files essentially just consists of an #include for this file with different
  preprocessor symbols set. The one that #included this file is dubbed "caller")
  This setup is to avoid a complicated Makefile setup and to avoid having to keep the three files
  consistent.
*/

/**
  SieveGauss.h is the #include-file a user should use who want to compile the ktuplesieve for the
  user's choice of template parameters (and, probably, *only* for this choice)

  Sieve.cpp contains explicit instantiations of the ktuplesieve for some specific (default) choices
  of SieveTraits template parameters. Compiling this will give a library that a user can use without
  having to recompile the ktuplesieve, but the user is restricted to the set of choices of template
  parameters that were set in Sieve.cpp

  Sieve.h is the #include-file a user should use for the precompiled library (made from Sieve.cpp).
  This works by using "extern template" for the same choices as in Sieve.cpp.
  Of course, the user has to link against the precompiled ktuplesieve library that was made from
  Sieve.cpp.
*/

/**
  Depending on the caller (i.e. the file that #included us), we set one of the following three
  macro constants:
  INCLUDE_SIEVE_TEMPL_FROM_SIEVE_GAUSS_H
  INCLUDE_SIEVE_TEMPL_FROM_SIEVE_H
  INCLUDE_SIEVE_TEMPL_FROM_SIEVE_CPP
*/

// Note that the caller .h files all have include guards.
// Furthermore, they do not (even indirectly) #include each other and Sieve.cpp includes neither of
// them:
// This would cause trouble with above macros getting set; in the case of Sieve.cpp, including
// Sieve.h would be outright wrong, because Sieve.h has extern template statements.
// The non-standard include guard used here is to raise an error if this is violated.

#ifdef SIEVE_TEMPL_H
  #error Including Sieve_templ.h twice
#endif
#define SIEVE_TEMPL_H

// Ensure that exactly one of the caller macros is set.
// We use #if defined(...) because there is no elifdef
#if defined(INCLUDE_SIEVE_TEMPL_FROM_SIEVE_GAUSS_H)
  #ifdef INCLUDE_SIEVE_TEMPL_FROM_SIEVE_H
    #error Mutually exclusive INCLUDE_SIEVE_TEMPL_* macros: GAUSS_H and SIEVE_H
  #endif
  #ifdef INCLUDE_SIEVE_TEMPL_FROM_SIEVE_CPP
    #error Mutually exclusive INCLUDE_SIEVE_TEMPL_* macros: GAUSS_H and SIEVE_CPP
  #endif
#elif defined(INCLUDE_SIEVE_TEMPL_FROM_SIEVE_H)
  #ifdef INCLUDE_SIEVE_TEMPL_FROM_SIEVE_CPP
    #error Mutually exclusive INCLUDE_SIEVE_TEMPL_* macros: SIEVE_H and SIEVE_CPP
  #endif
#elif !defined(INCLUDE_SIEVE_TEMPL_FROM_SIEVE_CPP)
  #error Sieve_templ.h must only be included by one of SieveGauss.h, Sieve.h, Sieve.cpp
#endif

/**
  The macros
  SIEVE_GAUSS_SINGLE_THREADED
  and
  SIEVE_GAUSS_MULTI_THREADED
  determine whether compile for single-threaded or multi-threaded mode or both.
  If either macro is defined, we support the corresponding mode.
  If neither is set, we fall back to a default (which is "both")
*/

// temporary:
#ifdef SIEVE_GAUSS_MULTI_THREADED
  #error Multithreaded version not supported yet
#endif
#ifndef SIEVE_GAUSS_SINGLE_THREADED
  #error Need to set single-threaded for now
#endif

// If neither of the above macros is set, we set a default value (compile both).
// We also set a default value for the template argument:
// SIEVE_GAUSS_DEFAULT_THREADED is set to true or false and determines the default template parameter
// (which is multi-threaded iff we compiled it)

#if !defined (SIEVE_GAUSS_SINGLE_THREADED) && !defined (SIEVE_GAUSS_MULTI_THREADED)
  #define SIEVE_GAUSS_SINGLE_THREADED
  #define SIEVE_GAUSS_MULTI_THREADED
  #define SIEVE_GAUSS_DEFAULT_THREADED true
#elif defined(SIEVE_GAUSS_MULTI_THREADED)
  #define SIEVE_GAUSS_DEFAULT_THREADED true
#elif defined(SIEVE_GAUSS_SINGLE_THREADED)
  #define SIEVE_GAUSS_DEFAULT_THREADED false
#endif

// We forward declare some templates and set default parameters.

#include "DefaultIncludes.h"
#include "Typedefs.h"

// we now include the "main (internal) header file SieveJoint.h for the Gauss Sieve.
// This file may be included twice with different values
// Note that this header file requires GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED to be set to either true or false.
// SieveJoint.h and SieveJoint_impl.h are designed in a way that allows them to be included twice
// with a different value of this macro.

#ifdef GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
  #error Macro already defined. This must never happen.
#endif

#ifdef SIEVE_GAUSS_SINGLE_THREADED
  #define GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED false
  #include "SieveJoint.h"
  #undef GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
#endif // SIEVE_GAUSS_SINGLE_THREADED

#ifdef SIEVE_GAUSS_MULTI_THREADED
  #define GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED true
  #include "SieveJoint.h"
  #undef GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
#endif // SIEVE_GAUSS_MULTI_THREADED

#include "Sampler.h"
#include "TerminationConditions.h"
#include "DefaultTermConds.h"

namespace GaussSieve
{
template<class SieveTraits, bool MT> class Sieve;

template<class CoefficientType>
using TupleSieve = Sieve<DefaultSieveTraits<CoefficientType, SIEVE_GAUSS_DEFAULT_THREADED, -1>, SIEVE_GAUSS_DEFAULT_THREADED>;

#ifdef SIEVE_GAUSS_SINGLE_THREADED
template<class CoefficientType>
using SieveST = Sieve<DefaultSieveTraits<CoefficientType, false, -1> , false>;
#endif
#ifdef SIEVE_GAUSS_MULTI_THREADED
template<class SieveTraits>
using SieveMT = Sieve<DefaultSieveTraits,CoefficientType, true, -1>, true>;
#endif

}  // end namespace GaussSieve

// decativated in Sieve.h, activated in the others
#if defined(INCLUDE_SIEVE_TEMPL_FROM_SIEVE_CPP) || defined(INCLUDE_SIEVE_TEMPL_FROM_SIEVE_GAUSS_H)

  /**
    Single-threaded-only implementation files
  */
  #ifdef SIEVE_GAUSS_SINGLE_THREADED
    #define GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED false
    #include "SieveJoint_impl.h"
    #include "SieveST_impl.h"
    #include "SieveST2_impl.h"
    #include "SieveST3_impl.h"
    #undef GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
  #endif // SIEVE_GAUSS_SINGLE_THREADED

  /**
    Multi-threaded-only implementation files
  */

  #ifdef SIEVE_GAUSS_MULTI_THREADED
    #define GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED true
    #include "SieveJoint_impl.h"
    #include "SieveMT.cpp"
    #undef GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
  #endif // SIEVE_GAUSS_MULTI_THREADED

  /**
    Implementation files for either case:
  */

  #include "DefaultTermConds_impl.h"
  #include "GaussQueue_impl.h"
  #include "GPVSampler_impl.h"
  #include "Sampler_impl.h"
  #include "UniformSampler_impl.h"
  #include "GPVSamplerExtended_impl.h"
  #include "GPVSamplerCVP_impl.h"

#endif // of block that is conditional on caller!=Sieve.h

// Sieve.cpp has explicit instantiation.
// Sieve.h has extern template.
// SieveGauss.h has nothing (instantiate as needed)

#if defined(INCLUDE_SIEVE_TEMPL_FROM_SIEVE_CPP) || defined(INCLUDE_SIEVE_TEMPL_FROM_SIEVE_H)

  #ifdef INCLUDE_SIEVE_TEMPL_FROM_SIEVE_CPP
    #define EXTERN_IN_H
  #else
    #define EXTERN_IN_H extern
  #endif

  namespace GaussSieve
  {

  #ifdef SIEVE_GAUSS_SINGLE_THREADED
    EXTERN_IN_H template class Sieve<DefaultSieveTraits<long,      false, -1>, false>;
    // EXTERN_IN_H template class Sieve<DefaultSieveTraits<double,    false, -1>, false>;
    EXTERN_IN_H template class Sieve<DefaultSieveTraits<mpz_class, false, -1>, false>;

    // testing only

    EXTERN_IN_H template class Sieve<DefaultSieveTraits<int32_t, false, -1, fplll::ZZ_mat<mpz_t>>,false>;

    EXTERN_IN_H template class Sampler<DefaultSieveTraits<long,      false, -1>, false>;
    //  EXTERN_IN_H template class Sampler<DefaultSieveTraits<double,    false, -1>, false>;
    EXTERN_IN_H template class Sampler<DefaultSieveTraits<mpz_class, false, -1>, false>;

    EXTERN_IN_H template class TerminationCondition<DefaultSieveTraits<long,      false, -1>, false>;
    //  EXTERN_IN_H template class TerminationCondition<DefaultSieveTraits<double,    false, -1>, false>;
    EXTERN_IN_H template class TerminationCondition<DefaultSieveTraits<mpz_class, false, -1>, false>;
  #endif

  #ifdef SIEVE_GAUSS_MULTI_THREADED
    EXTERN_IN_H template class Sieve<DefaultSieveTraits<long,      true, -1>, false>;
    //  EXTERN_IN_H template class Sieve<DefaultSieveTraits<double,    true, -1>, false>;
    EXTERN_IN_H template class Sieve<DefaultSieveTraits<mpz_class, true, -1>, false>;

    EXTERN_IN_H template class Sampler<DefaultSieveTraits<long,      true, -1>, true>;
    //  EXTERN_IN_H template class Sampler<DefaultSieveTraits<double,    true, -1>, true>;
    EXTERN_IN_H template class Sampler<DefaultSieveTraits<mpz_class, true, -1>, true>;

    EXTERN_IN_H template class TerminationCondition<DefaultSieveTraits<long,      true, -1>, true>;
    //  EXTERN_IN_H template class TerminationCondition<DefaultSieveTraits<double,    true, -1>, true>;
    EXTERN_IN_H template class TerminationCondition<DefaultSieveTraits<mpz_class, true, -1>, true>;
  #endif

  }  // end namespace GaussSieve

#endif // of block that is only there for Sieve.h and Sieve.cpp

// clang-format on
