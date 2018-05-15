#ifndef TEST_LATTICE_BASIS_H
#define TEST_LATTICE_BASIS_H

#include "../DefaultIncludes.h"
#include "../LatticeBases.h"

// explicitly instantiate:


template class GaussSieve::LatticeBasis<long, long, 20, 10>;
template class GaussSieve::LatticeBasis<long, long, -1, 10>;
template class GaussSieve::LatticeBasis<long, long, -1, -1>;
template class GaussSieve::LatticeBasis<long, long, 20, -1>;
//GaussSieve::LatticeBasis<long, long

bool test_lattice_basis(bool output = true)
{
  static constexpr unsigned int RANK = 10;
  static constexpr unsigned int DIM  = 20;

  using LB  = GaussSieve::LatticeBasis<long, long, DIM, RANK>;
  using LB2 = GaussSieve::LatticeBasis<long>;
  using LongDIM = long[DIM];
  using LongRANK = long[RANK];
  LongDIM basis_vecs[RANK];
  LongRANK gram_vecs[RANK];
  std::vector<std::vector<double>> mu_vecs;
  mu_vecs.clear();
  for(unsigned int i=0; i < RANK; ++i)
  {
    std::vector<double> v;
    v.reserve(RANK);
    for (unsigned int j = 0; j < RANK; ++j)
    {
      v.push_back(double(i)/double(j+1)); // just some arbitrary meaningless values
    }
    mu_vecs.push_back(v);
  }
  for(unsigned int i=0; i < RANK; ++i)
  {
    for(unsigned int j = 0; j < RANK; ++j)
    {
      gram_vecs[i][j] = (i+1) * (j+1);
    }
    for(unsigned int j= 0; j < DIM; ++j)
    {
      basis_vecs[i][j] = 1000 * (i+1) + j;
    }
  }
  LB latb1{DIM, RANK, basis_vecs, mu_vecs, gram_vecs};
  std::cout << latb1 << std::endl;

  return true;
}

#endif // include guards
