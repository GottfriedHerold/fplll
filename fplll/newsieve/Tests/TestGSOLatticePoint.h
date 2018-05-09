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
#include "../LatticeBases.h"

bool test_GSO_LP(bool output = true)
{
  using GaussSieve::GSOLatticePoint;
  using GaussSieve::MaybeFixed;
  using GSOType = fplll::MatGSO<fplll::Z_NR<long>, fplll::FP_NR<double>>;



  //typedef GaussSieve::GSOLatticePoint<fplll::FP_NR<double>,fplll::FP_NR<double>, GSOType, -1> LPvar;

  //GaussSieve::StaticInitializer<LPvar> init1 (MaybeFixed<-1>{10});
  //GaussSieve::StaticInitializer<LPvar> init2 (MaybeFixed<-1>{10});
  //GaussSieve::StaticInitializer<LPfix> init3 (MaybeFixed<10>{10});
  //GaussSieve::StaticInitializer<LPGMP> init4 (MaybeFixed<10>{10});

  return true;
};

#endif
