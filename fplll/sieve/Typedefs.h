#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include "LatticePointsNew.h"

//various typedef declarations that control the types used by our classes.
//This typedef defines the return type that the Sampler should have.

namespace GaussSieve{
template<class ET,bool MT, int nfixed> using GaussSampler_ReturnType = CompressedPoint<ET,MT,nfixed>;

// unfortunately, trigonometric functions to compute pi don't have constexpr variants on all compilers we want to support, so we just define pi directly
long double constexpr   pi_long     = 3.14159265358979323846264338327950288419716939937510L;
double constexpr        pi_double   = 3.14159265358979323846264338327950288419716939937510;
long double constexpr   pi          = 3.14159265358979323846264338327950288419716939937510L;

};


#endif
