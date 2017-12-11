// clang-format off

#define USE_REGULAR_QUEUE // only regular queue is implemented for now
                          // For large dimensions priority queue might be faster

#define PROGRESSIVE


/*
  This provides an implementation of k-tuple Gauss sieving
  for k=2 and k=3.
*/

#include "fplll.h"
#include <thread>
#include <chrono>
#include "SieveGauss_main.h"
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <stdio.h>
#include <unistd.h>


using namespace GaussSieve;

/**
 * help function
 */

static void main_usage(char *myself)
{
  std::cout << "Usage: " << myself << " [options]\n"
       << "List of options:\n"
       << "  -k [2|3|]\n"
       << "     2- or 3-sieve;\n"
       << "  -d nnn\n"
       << "     dimension;\n"
       << "  -f filename\n"
       << "     Filename with input basis\n"
       << "  -t nnn\n"
       << "     Targeted norm^2=nnn\n"
       << "  -b nnn\n"
       << "     BKZ preprocessing of blocksize=nnn\n"
       << "  -v\n"
       << "     Verbose mode\n";
  exit(0);
}

/**
 * main function
 */

int main(int argc, char **argv)
{
  char *target_norm_string = NULL;
  char* input_file_name = NULL;
  bool flag_file = false;
  int opt;
  int dim;
  fplll::Z_NR<mpz_t> target_norm;
  target_norm = 0;
  mpz_class target_norm_conv;
  int k=2;
  int beta=0; //beta = 0 means we run LLL

  if (argc < 2)
  {
      std::cout << "Please provide the dimension." << std::endl;
      return -1;
  }

  while ((opt = getopt(argc, argv, "d:t:k:f:b:")) != -1) {
      switch (opt) {
          case 'd':
              dim = atoi(optarg);
              break;
          case 'k':
              k=atoi(optarg);
          case 't':
              target_norm_string = optarg;
              break;
          case 'f':
              input_file_name = optarg;
              flag_file = true;
              break;
          case 'b':
              beta = atoi(optarg);
              break;
          break;
          case 'h':
            main_usage(argv[0]);
            return -1;
          case '?':
            main_usage(argv[0]);
            return -1;
          case ':':
            main_usage(argv[0]);
            return -1;
      }
  }
  
  if (dim==0)
  {
    std::cout << "Please, provide the dimension" << std::endl;
    return -1;
  }

  if (target_norm_string!=NULL)
  {
      target_norm_conv = mpz_class(target_norm_string);
  }
  
  // ZZ_mat is an integer row-oriented matrix. See /nr/matrix.h
  fplll::ZZ_mat<mpz_t> B;
  B.resize(dim, dim);

  if (flag_file) {
      std::ifstream input_file(input_file_name);
      if (input_file.is_open()) {
          std::cout << "reading B from file ..." << std::endl;
          input_file >> B;
          input_file.close();
      }
  }
  else {
      srand (1);
      //generates GM lattice
      B.gen_qary_prime(1, 10*dim);
  }

  std::cout << "beta = " << beta << std::endl;

  /* preprocessing of basis */
  clock_t stime = clock();
  if (beta > 0)
      fplll::bkz_reduction(B, beta, fplll::BKZ_DEFAULT, fplll::FT_DEFAULT, 0);
  else
      fplll::lll_reduction(B, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA, fplll::LM_WRAPPER);

  clock_t etime = clock();
  double secs   = (etime - stime) / (double)CLOCKS_PER_SEC;

  if (beta > 0)
      std::cout << "# [info] BKZ took time " << secs << " s" << std::endl;
  else
      std::cout << "# [info] LLL took time " << secs << " s" << std::endl;


  #ifndef USE_REGULAR_QUEUE
      std::cout << "Use Priority Queue" << std::endl;
  #else
      std::cout << "Use Standard Queue" << std::endl;
  #endif
  
  auto start = std::chrono::high_resolution_clock::now();
  
  bool constexpr multithreaded = false;

  // Define all the types, consts for the sieve
  // template params are <entry type for sieving, single/multi-threaded, is_dim_fixed, entry type of reduced B>
  // here we assume that the entries of reduced B fit in long or int32_t
  
  using Traits = GaussSieve::DefaultSieveTraits<int32_t, false, -1, fplll::ZZ_mat< mpz_t > >;
  //using Traits = GaussSieve::DefaultSieveTraits<long, multithreaded, -1, fplll::ZZ_mat< mpz_t >>;
  
  //instantiate the Sieve class with the basis, k-number of tuples, and verbosity
	Sieve<Traits, multithreaded> test_sieve (B, k, 0);
  
  TerminationCondition<Traits,multithreaded> * termcond;
  
  if (target_norm_conv!=0)
  {
    termcond = new LengthTerminationCondition<Traits, multithreaded> (ConvertMaybeMPZ<long>::convert_to_inttype(target_norm_conv));
  }
  else
  {
    termcond = new MinkowskiTerminationCondition<Traits, multithreaded>;
  }
  
	test_sieve.set_termination_condition(termcond);


  test_sieve.run();
  std::cout << "sv is " << std::endl;
  test_sieve.print_status();


  auto finish = std::chrono::high_resolution_clock::now();
  auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
  std::cout << " Time taken: " << microseconds.count()/1000000.0 << "sec" << std::endl;
  delete termcond;

  return 1;
}

//clang-format on
