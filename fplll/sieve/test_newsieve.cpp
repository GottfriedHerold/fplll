
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
    BTest.gen_trg(1.1);
    lll_reduction(BTest, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);
    Sieve<Z_NR<ZT>, false> Test_Queue (BTest);
    
    auto start = std::chrono::high_resolution_clock::now();
    Test_Queue.run_2_sieve();
    auto finish = std::chrono::high_resolution_clock::now();
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);

    Test_Queue.print_status(-1,ofs);
    
    
    cout << " Time taken: " << microseconds.count()/1000000.0 << "sec" << endl;
}

int main(int argc, char **argv)
{
    bool longflag = false; //to be read_in as input
    int dim[] = {60, 62};
    int length = 2;

    if (longflag)
    {
        //for (int i=0; i<length; i++)
        //    test_run_sieve<long>(dim[i]);
    }

    else
    {
        std::ofstream ofs("test_sive.txt");
        ZZ_mat<mpz_t> B;
        B.resize(dim[0], dim[0]);
        B.gen_trg(1.1);
        test_run_sieve<mpz_t>(dim[0], ofs);
    }

    //auto start = std::chrono::high_resolution_clock::now();
//    auto finish = std::chrono::high_resolution_clock::now();
//    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
//    cout << " Time taken: " << microseconds.count()/1000000.0 << "sec" << endl;
  return 1;
}
