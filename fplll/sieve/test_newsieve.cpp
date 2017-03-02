#define USE_REGULAR_QUEUE //comment out if you use priority-queue


#include "sieve_main.h"
#include "fplll.h"
#include "LatticePoint.h"
#include "PointList.h"
#include <thread>
#include <chrono>
#include "SieveGauss.h"
#include "LatticePoint2.h"
#include <iostream>
#include <fstream>

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
    Test_Queue.run_2_sieve();
    auto finish = std::chrono::high_resolution_clock::now();
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);

    Test_Queue.print_status(-1,ofs);
    
    
    ofs << " Time taken: " << microseconds.count()/1000000.0 << "sec" << endl;
}

int main(int argc, char **argv)
{
        //int dim[] = {52, 54, 56, 58, 60, 62, 64};
        int dim = 62;
        int length = 7;

    	#ifdef USE_REGULAR_QUEUE 
        std::ofstream ofs("test_sieve_PQ_dim" +to_string(dim) + ".txt");
		ofs << "WITH PRIORITY QUEUE" << endl;
	#else 
		std::ofstream ofs("test_sieve_dim"+to_string(dim) + ".txt");
		ofs << "WITH STANDARD QUEUE" << endl;
	#endif
//	for (int i=0; i<length; i++) {
//        	
//       		
//		ofs << "start sieve on lattice of dim = " << dim[i] << endl;
//        	test_run_sieve<mpz_t>(dim[i], ofs); 
//		ofs << "----------------------------------------" << endl;
//	}
    
    for (int i=0; i<1; i++)
    {
        ofs << "start sieve on lattice of dim =  " << dim << endl;
        test_run_sieve<mpz_t>(dim, ofs);
        ofs << "----------------------------------------" << endl;
    }
    

   ofs.close();
  return 1;
}
