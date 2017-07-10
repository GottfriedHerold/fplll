#define SIEVE_GAUSS_SINGLE_THREADED //only single-threaded for now.

//SIEVE_GAUSS_SINGLE_THREADED : if defined by user, means that we request the single-threaded implementation
//SIEVE_GAUSS_MULTI_THREADED  : if defined by user, means that we request the multi -threaded implementation
//If neither is set, we default to both.


#if !defined (SIEVE_GAUSS_SINGLE_THREADED) && !defined (SIEVE_GAUSS_MULTI_THREADED)
#define SIEVE_GAUSS_SINGLE_THREADED
#define SIEVE_GAUSS_MULTI_THREADED
#define SIEVE_GAUSS_DEFAULT_THREADED true
#elif defined(SIEVE_GAUSS_MULTI_THREADED)
#define SIEVE_GAUSS_DEFAULT_THREADED true
#elif defined(SIEVE_GAUSS_SINGLE_THREADED)
#define SIEVE_GAUSS_DEFAULT_THREADED false
#endif

#ifdef  GAUSS_SIEVE_IS_MULTI_THREADED
#undef  GAUSS_SIEVE_IS_MULTI_THREADED
#endif

#ifdef SIEVE_GAUSS_SINGLE_THREADED //if Single-Threaded implementation is requested:
#define GAUSS_SIEVE_IS_MULTI_THREADED false
//#include "SieveJoint.h"
#include "SieveJoint.cpp"
#include "SieveST.cpp"
#undef GAUSS_SIEVE_IS_MULTI_THREADED
#endif // SIEVE_GAUSS_SINGLE_THREADED

#ifdef SIEVE_GAUSS_MULTI_THREADED //If Multi-Threaded implementation is requested:
#define GAUSS_SIEVE_IS_MULTI_THREADED true
//#include "SieveJoint.h"
#include "SieveJoint.cpp"
#include "SieveMT.cpp"
#undef GAUSS_SIEVE_IS_MULTI_THREADED
#endif // SIEVE_GAUSS_MULTI_THREADED
