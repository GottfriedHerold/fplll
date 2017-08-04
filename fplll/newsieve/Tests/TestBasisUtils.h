#ifndef TEST_BASIS_UTILS_H
#define TEST_BASIS_UTILS_H

#include "../LatticeBases.h"
#include "fplll.h"
#include "../Typedefs.h"
#include "gmpxx.h"

bool test_basis_utils()
{
  int constexpr dim = 20;
  fplll::ZZ_mat<mpz_t> B;
  B.resize(dim, dim);
  srand (1);

  //generates GM lattice
  B.gen_qary_prime(1, 10*dim);

  fplll::lll_reduction(B, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA, fplll::LM_WRAPPER);

  std::cout << "Generating random basis: Result is" << std::endl;
  std::cout << B << std::endl << std::flush;

  using Traits = GaussSieve::DefaultSieveTraits<mpz_class, false, -1>;



  GaussSieve::SieveLatticeBasis<Traits,false> sieve_basis(B);

  double maxbistar2 = sieve_basis.get_maxbistar2();

  std::cout << "maxbistar = " << maxbistar2 << std::endl;

  std::vector<std::vector<double>> mu = sieve_basis.get_mu_matrix();
  std::cout << "Mu-matrix:" << std::endl;
  for(uint_fast16_t i=0; i < dim ; ++i)
  {
    for(uint_fast16_t j=0;j<dim;++j)
    {
    std::cout << mu[i][j] << " ";
    }
    std::cout << std::endl;
  }

  return true;
}


#endif
