//This is the header file for the Gauss Sieve (both single-threaded and multi-threaded version).
//The actual classes are available as SieveST<ZT> = Sieve<ZT,false> resp. SieveMT = Sieve<ZT,true>
//Include SieveST.h resp. SieveMT.h if you want only the single-threaded resp. multi-threaded versions.
//Don't include more than one of SieveGauss.h, SieveST.h, SieveMT.h

//Note:
//With C++ lacking a good "static_if" that works with templates (as of now), we use some pre-processor shenanigans to avoid code duplication:
//The "true" header file for the single-threaded Gauss Sieve is SieveJoint.h (with GAUSS_SIEVE_SINGLE_THREADED set to true)
//The "true" header file for the multi -threaded Gauss Sieve is SieveJoint.h (with GAUSS_SIEVE_MULTI_THREADED set to true)
//(setting both/neither gives an error)
//More complicated functions go to SieveJoint.cpp (where we use preprocess macros to change details)
//If the implementations are completely seperate, we use SieveST resp. SieveMT (there should be a comment in SieveJoint in this case).

#ifndef SIEVE_GAUSS_BOTH_THREADED_H
#define SIEVE_GAUSS_BOTH_THREADED_H

#define SIEVE_GAUSS_DEFAULT_THREADED false

#ifdef GAUSS_SIEVE_MULTI_THREADED
#undef GAUSS_SIEVE_MULTI_THREADED
#endif

#define GAUSS_SIEVE_SINGLE_THREADED
#define GAUSS_SIEVE_IS_MULTI_THREADED false
#include "SieveJoint.h"
#include "SieveJoint.cpp"
#include "SieveST.cpp"
#undef  GAUSS_SIEVE_SINGLE_THREADED
#define GAUSS_SIEVE_MULTI_THREADED
#undef  GAUSS_SIEVE_IS_MULTI_THREADED
#define GAUSS_SIEVE_IS_MULTI_THREADED true
#include "SieveJoint.h"
#include "SieveJoint.cpp"
#include "SieveMT.cpp"
#undef  GAUSS_SIEVE_MULTI_THREADED

#endif
