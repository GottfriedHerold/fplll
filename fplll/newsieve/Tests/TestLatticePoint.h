#error currently unused

/**
  (templated) unit test for general lattice point classes.
  Used as subroutine for specific tests
*/

#ifndef TEST_GENERAL_LATTICE_POINT
#define TEST_GENERAL_LATTICE_POINT

// InitArg const &init_arg is the argument to the static_initializer
template<class LatticePoint, class InitArg>
bool test_general_lattice_point(InitArg const &init_arg, unsigned int dim)
{
  static_assert(GaussSieve::IsALatticePoint<LatticePoint>::value,"");
  GaussSieve::StaticInitializer<LatticePoint> init1{init_arg};
  GaussSieve::StaticInitializer<LatticePoint> init2{init_arg}; // double-init should work
  return true;
}

#endif
