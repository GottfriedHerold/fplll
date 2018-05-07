#ifndef TEST_GSO_LATTICE_POINT_H
#define TEST_GSO_LATTICE_POINT_H

#include <type_traits>
#include "../GSOLatticePoint.h"

#include "fplll/defs.h"
#include "fplll/nr/nr.h"
#include "vector"
#include "gmpxx.h"
#include <iostream>
#include "TestLatticePoint.h"

bool test_GSO_LP()
{
  using GaussSieve::GSOLatticePoint;
  using GaussSieve::MaybeFixed;

  //typedef GaussSieve::GSOLatticePoint<long, -1> LPvar;
  
  //GaussSieve::StaticInitializer<LPvar> init1 (MaybeFixed<-1>{10});
  //GaussSieve::StaticInitializer<LPvar> init2 (MaybeFixed<-1>{10});
  //GaussSieve::StaticInitializer<LPfix> init3 (MaybeFixed<10>{10});
  //GaussSieve::StaticInitializer<LPGMP> init4 (MaybeFixed<10>{10});

  return true;
};

#endif
