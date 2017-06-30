//#define USE_REGULAR_QUEUE //comment out if you use priority-queue

#include "sieve_main.h"
#include "fplll.h"
//#include "LatticePoint.h"
//#include "PointList.h"
#include <thread>
#include <chrono>
#include "SieveGauss.h"
//#include "LatticePoint2.h"
#include <iostream>
#include <fstream>
#include "Sampler.h"
#include <random>
//#include "FilteredList.h"

using namespace fplll;

// Z_NR can be either long, double, or mpz_t

template <class ZT> void test_run_sieve(int dim, std::ofstream &ofs)
{
    //TerminationConditions< Z_NR<ZT> > term_cond; //sets target-legnth as Mink. bound


    //lll_reduction(B, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);
    //Sieve<Z_NR<ZT>, false> Test_Queue (B);
    //Test_Queue.run_2_sieve();

    ZZ_mat<ZT> BTest;
    BTest.resize(dim, dim);
    //BTest.gen_trg(1.1);
    srand (1);
    BTest.gen_qary_prime(1, 10*dim);

    if (dim >= 60)
        bkz_reduction(BTest, 8, BKZ_DEFAULT, FT_DEFAULT, 0);
    else
        lll_reduction(BTest, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);

    Sieve<Z_NR<ZT>, false> Test_Queue (BTest);

    ofs << "sieve is run on B[0]" << BTest[0] << endl;

    auto start = std::chrono::high_resolution_clock::now();
    Test_Queue.run();
    auto finish = std::chrono::high_resolution_clock::now();
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);

    Test_Queue.print_status(-1,ofs);


    ofs << " Time taken: " << microseconds.count()/1000000.0 << "sec" << endl;
}

template <class Z> void sample_gaussians(int number, double s, double center, double cutoff)
{
    std::mt19937_64 engine;
    for (int i=0; i<number;++i)
    {
        cout << GaussSieve::sample_z_gaussian<Z, std::mt19937_64>(s,center,engine, cutoff) << endl;
    }
}



int main(int argc, char **argv)
{
sample_gaussians<long>(50, 10.0, 0.3, 4.0);
sample_gaussians<long>(50, 0.0000001, 0.48, 4.0); //should still be fast.


ZZ_mat<mpz_t> B,u,u_inv;
B.resize(10, 10);
Matrix<FP_NR<double > > r,mu;
Matrix<Z_NR<mpz_t> > g;
    //generates a lower-triangular matrix B; the argument determines (in a complicated way) the bit-size of entries
    //B.gen_trg(1.1);

srand (1);
    //generates GM lattice
B.gen_qary_prime(1, 15);


auto pGSO = new MatGSO<Z_NR<mpz_t>, FP_NR<double>>(B, u, u_inv, 1);

  pGSO->update_gso();
  mu = pGSO->get_mu_matrix();
  r  = pGSO->get_r_matrix();
  g  = pGSO->get_g_matrix();

cout << B << endl;
cout << mu<< endl;
cout << r << endl;
cout << g << endl;


//        //int dim[] = {52, 54, 56, 58, 60, 62, 64};
//        int dim = 62;
//        int length = 7;
//
//    	#ifdef USE_REGULAR_QUEUE
//        std::ofstream ofs("test_sieve_PQ_dim" +to_string(dim) + ".txt");
//		ofs << "WITH PRIORITY QUEUE" << endl;
//	#else
//		std::ofstream ofs("test_sieve_dim"+to_string(dim) + ".txt");
//		ofs << "WITH STANDARD QUEUE" << endl;
//	#endif
//	for (int i=0; i<length; i++) {
//
//
//		ofs << "start sieve on lattice of dim = " << dim[i] << endl;
//        	test_run_sieve<mpz_t>(dim[i], ofs);
//		ofs << "----------------------------------------" << endl;
//	}

//    for (int i=0; i<1; i++)
//    {
//        ofs << "start sieve on lattice of dim =  " << dim << endl;
//        test_run_sieve<mpz_t>(dim, ofs);
//        ofs << "----------------------------------------" << endl;
//    }
//   ofs.close();
  return 1;
}
