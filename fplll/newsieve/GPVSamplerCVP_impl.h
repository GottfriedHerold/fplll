#ifndef GPV_SAMPLER_CVP_IMPL_H
#define GPV_SAMPLER_CVP_IMPL_H

#include "DefaultIncludes.h"
#include "GPVSamplerCVP.h"
#include "LatticeBases.h"
#include "Sampler.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr.h"
#include "fplll/nr/nr_Z.inl"
#include "fplll/svpcvp.h"

#include "gmp.h"
//#include <string> // in DefaultIncludes.h
//#include <sstream>

namespace GaussSieve
{
template <class SieveTraits, bool MT, class Engine, class Sseq>
void GPVSamplerCVP<SieveTraits, MT, Engine, Sseq>::custom_init(
    SieveLatticeBasis<SieveTraits, MT> const &input_basis)
{
  assert(!initialized);
#ifndef DEBUG_SIEVE_STANDALONE_SAMPLER
  assert(sieveptr != nullptr);
#else
  assert(sieveptr == nullptr);
#endif

  dim          = input_basis.ambient_dimension;
  lattice_rank = input_basis.lattice_rank;
  mu_matrix    = make_array_or_vector<mystd::decay_t<decltype(mu_matrix) >>(Parent::basis_ptr->get_mu_matrix(),lattice_rank, lattice_rank);
//  mu_matrix    = input_basis.get_mu_matrix();

  //assert(start_babai < lattice_rank);  // use strictly less to prevent always outputting 0-vector
  assert(lattice_rank <= dim);

  s2pi.resize(lattice_rank);
  maxdeviations.resize(lattice_rank);
  basis.resize(lattice_rank);
  basis_for_cvp.resize(lattice_rank, dim);

  auto const maxbistar2 = input_basis.get_maxbistar2();

  double const st_dev = maxbistar2 * 15.2;

  for (uint_fast16_t i = 0; i < lattice_rank; ++i)
  {

    double const maxdev_nonsc =
        st_dev / convert_to_double(input_basis.get_g(i, i));  // g_(i,i) = ||b^*_i||^2

    s2pi[i]          = maxdev_nonsc / GaussSieve::pi;
    maxdeviations[i] = sqrt(maxdev_nonsc) * cutoff;

    basis[i] = input_basis.get_basis_vector(i).make_copy();

    //std::cout << basis[i] << std::endl;

    //TODO: make mpz_import work
    for (uint_fast16_t j = 0; j < dim; ++j)
    {
      //std::cout << basis[i][j] << " ";
      mpz_t tmp;
      mpz_init(tmp);
      std::stringstream ss;
      ss<<basis[i][j];
      //std::cout<< ss.str() << " ";
      mpz_init_set_str(tmp, ss.str().c_str(),10);
      //mpz_import(tmp, 1, -1, sizeof(basis[i][j]), 0, 0, &(basis[i][j])); // TO BE TESTED
      //gmp_printf ("tmp is ", tmp);
      basis_for_cvp[i][j] = tmp;
    }
    //std::cout <<"basis_for_cvp[i]: " << basis_for_cvp[i] << std::endl;

  }

  using RetType = typename SieveTraits::GaussSampler_ReturnType;

  if (static_init_plainpoint != nullptr)
  {
    assert(false);
  }
  if (static_init_rettype != nullptr)
  {
    assert(false);
  }

  static_init_rettype    = new StaticInitializer<RetType>(MaybeFixed<SieveTraits::get_nfixed>{dim});
  static_init_plainpoint = new StaticInitializer<typename SieveTraits::PlainPoint>(
      MaybeFixed<SieveTraits::get_nfixed>{dim});

  initialized = true;
  std::cout <<"initialization of GPVSamplerCVP is finshed" << std::endl;
}

// TODO: Not using thread / Engine / MTPRNG correctly
template <class SieveTraits, bool MT, class Engine, class Sseq>
typename SieveTraits::GaussSampler_ReturnType
GPVSamplerCVP<SieveTraits, MT, Engine, Sseq>::sample(int const thread)
{
  assert(initialized);
#ifdef DEBUG_SIEVE_STANDALONE_SAMPLER
  assert(sieveptr == nullptr);
#else
  assert(sieveptr != nullptr);
#endif
  typename SieveTraits::PlainPoint vec;
  vec.fill_with_zero();

  std::vector<double> shifts(lattice_rank, 0.0);
  std::vector<fplll::Z_NR<mpz_t>> cvp_sol;
  cvp_sol.resize(dim);

  while (vec.is_zero())
  {
#ifdef PROGRESSIVE
    uint_fast16_t i = this->get_progressive_rank();
#else
    uint_fast16_t i = lattice_rank;
#endif
    // run GPV sampling

    //start_babai = 0; // TODO: make it work for any start_babai
    while (i > start_babai)
    {
      --i;

      /* Does the same as below but with a faster for-loop.
         However, requires to store the coos-vector.

         Gotti: actually, the main difference is that
         the above version allows to operate on a narrow data type for most of the time.
         The extra storage requirement does not matter if allocation is done only once.
       *
      for (uint_fast16_t j = lattice_rank-1; j > i; --j)  // adjust shifts
      {
        shifts[i] -= coos[j] * (mu_matrix[j][i]);
      }
       */


      //double const newc = sample_fp_gaussian<double, Engine>(s2pi[i], shifts[i], engine.rnd(thread));

      long const newcoeff = sample_z_gaussian_VMD<long, Engine>(
          s2pi[i], shifts[i], engine.rnd(thread), maxdeviations[i]);  // coefficient of b_j in vec.

      vec += basis[i] * newcoeff;

      for (uint_fast16_t j = 0; j < i; ++j)  // adjust shifts
      {
        shifts[j] -= newcoeff * (mu_matrix[i][j]);
      }
    }

    std::vector<fplll::Z_NR<mpz_t>> target_vec;
    target_vec.resize(dim);

    //std::cout<<"vec before integer shift: " << vec << std::endl;

    // generate random shift to move vec from the lattice
    // convert vec to target_vec with Z_NR<mpz_t>-entries suirable for the fplll cvp-oracle

    //TODO: make mpz_import work
    for (uint_fast16_t j = dim-15; j < dim; ++j)
    {
        int coo_shift = sample_uniform2<Engine>(-200, 200,engine.rnd(thread));
        vec[j]+= coo_shift;
        //mpz_t tmp;
        //mpz_import(tmp, 1, -1, sizeof vec[i], 0, 0, &vec[i]);
        //target_vec[i] = tmp;
    }
    for (uint_fast16_t j = 0; j < dim; ++j)
    {
        mpz_t tmp;
        mpz_init(tmp);
        std::stringstream ss;
        ss<<vec[j];
        //std::cout<< ss.str() << " ";
        mpz_init_set_str(tmp, ss.str().c_str(),10);
        target_vec[j] = tmp;
    }
    //std::cout<<"vec after integer shift: " << vec << std::endl;
    //std::cout << "target_vec: " << target_vec << std::endl;

    std::vector<fplll::Z_NR<mpz_t>> sol_coord;
    sol_coord.resize(lattice_rank);


    // alternative to Babai: calling cvp-solver from fplll
    //int cvp_stat = fplll::closest_vector(fplll::ZZ_mat<mpz_t> &basis_for_cvp, const vector<fplll::Z_NR<mpz_t>> &int_target,
    //          vector<fplll::Z_NR<mpz_t>> &sol_coord, int method = fplll::CVPM_FAST, int flags = fplll::CVP_DEFAULT);

    int cvp_stat = fplll::closest_vector(basis_for_cvp,target_vec,sol_coord,fplll::CVPM_FAST,fplll::CVP_DEFAULT);

    //std::cout  << "cvp_stat = " << cvp_stat <<std::endl;
    std::cout << "cvp returns: " << sol_coord << std::endl;

    assert(cvp_stat==0 && "error in cvp!");

    for (uint_fast16_t i = 0; i < dim; ++i)
    {
      cvp_sol[i] = 0;
      for (uint_fast16_t j = 0; j < lattice_rank; ++j)
      {
        cvp_sol[i].addmul(sol_coord[j],basis_for_cvp[j][i]);
      }

    }

    //std::cout << "cvp_sol = " << cvp_sol << std::endl;
  }

  // TODO : Fix conversion here.
  typename SieveTraits::GaussSampler_ReturnType ret;
  ret = make_from_znr_vector<typename SieveTraits::GaussSampler_ReturnType>(cvp_sol, dim);
  //std::cout <<"ret = " << ret << std::endl;
  return ret;
}


template <class SieveTraits, bool MT, class Engine, class Sseq>
inline std::ostream &GPVSamplerCVP<SieveTraits, MT, Engine, Sseq>::dump_to_stream(std::ostream &os) const
{
  os << "cutoff" << cutoff <<'\n';
  os << "StartBabai" << start_babai;
  return os;
}

template <class SieveTraits, bool MT, class Engine, class Sseq>
inline std::istream &GPVSamplerCVP<SieveTraits, MT, Engine, Sseq>::read_from_stream(std::istream &is)
{
  if(initialized)
  {
    initialized = false;
    delete static_init_plainpoint;
    static_init_plainpoint = nullptr;
    delete static_init_rettype;
    static_init_rettype = nullptr;
  }
  if (!string_consume(is, "cutoff")) throw bad_dumpread("Could not read GPVSamplerCVP");
  is >> cutoff;
  if (!string_consume(is, "StartBabai")) throw bad_dumpread("Could not read GPVSamplerCVP");
  is >> start_babai;
  sieveptr = nullptr;
  return is;
}


}  // end namespace GaussSieve

#endif
