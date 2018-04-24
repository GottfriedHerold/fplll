/**
  (templated) unit test for general lattice point classes.
  Used as subroutine for specific tests
*/

#ifndef TEST_GENERAL_LATTICE_POINT
#define TEST_GENERAL_LATTICE_POINT

#include <sstream>

// InitArg const &init_arg is the argument to the static_initializer
template<class LatticePoint, int nfixed, class InitArg>
bool test_general_lattice_point(InitArg const &init_arg, GaussSieve::MaybeFixed<nfixed> const  = {})
{
  static constexpr unsigned int dim = nfixed;
  static_assert(dim > 0, "Invalid dimension for test");
  static_assert(dim > 2, "Too small dimension for test"); // some accidential equalities would screw up tests
  static_assert(dim < 1000, "Too large dimension for test"); // fear of overflows
  static_assert(GaussSieve::IsALatticePoint<LatticePoint>::value,"Not a lattice Point");
  GaussSieve::StaticInitializer<LatticePoint> init1{init_arg};
  GaussSieve::StaticInitializer<LatticePoint> init2{init_arg}; // double-init should work

  std::stringstream channel;
  channel << init1;
  channel >> init1;
  channel.clear();

  using GaussSieve::MaybeFixed;

  LatticePoint X1;
  LatticePoint X2{dim};
  LatticePoint X3(GaussSieve::MaybeFixed<-1>{dim});
  LatticePoint Y1, Y2, Y3;

  std::vector<long> vec;
  vec.reserve(dim);
  std::array<long,dim> arr;
  long carr[dim];
  long carr2[dim];

  for(unsigned int i = 0; i < dim; ++i)
  {
  vec[i] = i;
  arr[i] = i*i;
  carr[i]=10*i;
  carr2[i]=dim - 1;
  }

  X1 = GaussSieve::make_from_any_vector<LatticePoint>(vec, MaybeFixed<-1>(dim));
  X2 = GaussSieve::make_from_any_vector<LatticePoint>(vec, dim);
  X3 = GaussSieve::make_from_any_vector<LatticePoint>(vec, dim);
  assert(X1==X2);
  assert(X1==X3);

  Y1 = GaussSieve::make_from_any_vector<LatticePoint>(vec,  MaybeFixed<dim>{});
  Y2 = GaussSieve::make_from_any_vector<LatticePoint>(arr,  MaybeFixed<dim>());
  Y3 = GaussSieve::make_from_any_vector<LatticePoint>(carr, MaybeFixed<dim>());
  assert(Y1==X1);
  assert(Y1!=Y2);
  assert(Y1!=Y3);
  assert(X1 <= X2);
  assert(X2 >= X1);
  assert(Y1 < Y2);
  assert(Y2 > Y1);
  assert(!(X1<X2));
  static constexpr long real_norm2  = (dim * (dim-1) * (2*dim - 1)) / 6; // norm^2 of vec
  assert(X1.get_norm2() == real_norm2 );
  assert(Y3.get_norm2() == 100* real_norm2);
  X2 = X1.make_copy(); // does not change anything.
  X3 = GaussSieve::make_from_any_vector<LatticePoint>(carr2, MaybeFixed<-1>(dim));
  assert(X3.get_norm2() == dim * (dim-1)* (dim-1));
  X2 = X3 - X1; // so X2 is now X1 in reverse.
  assert(X2.get_norm2() == real_norm2);
  assert(X2+ X1 == X3);
  assert( compute_sc_product(X1,X3) == ((dim-1)*(dim-1)*dim) / 2 );
  assert( X1.get_dim() == dim );
//  assert( X1.get_internal_rep_size() == 10);
  assert( Y1.get_dim() == dim);
//  assert( Y1.get_internal_rep_size() == 10);
  assert(X1.is_zero() == false);
  X1 = X1 - X1;
  assert(X1.is_zero() == true);
  X2.fill_with_zero();
  assert(X2.is_zero() == true);

  return true;
}

#endif
