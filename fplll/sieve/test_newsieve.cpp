
#include "sieve_main.h"
#include "fplll.h"
//#include "LatticePoint.h"
//#include "PointList.h"
#include <thread>
#include <chrono>
#include "SieveGauss.h"
#include <iostream>
#include <fstream>

using namespace fplll;

// Z_NR can be either long, double, or mpz_t

template <class ZT>
void test_run_sieve(int dim)
{
    TerminationConditions< Z_NR<ZT> > term_cond; //sets target-legnth as Mink. bound
    ZZ_mat<ZT> B;
    B.resize(dim, dim);
    B.gen_trg(1.1);
    lll_reduction(B, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);
    Sieve<Z_NR<ZT>, false> Test_Queue (B);
    //Test_Queue.run_2_sieve(term_cond); //file-name to output the results of tests
    {
    std::ofstream ofs("FILENAME");
    Test_Queue.print_status(-1,ofs);
    }
}

int main(int argc, char **argv)
{

    cout << "Hello " << endl;

    bool longflag = false; //to be read_in as input
    int dim[] = {25, 30, 35, 40, 42, 44, 46, 48, 50, 52, 54, 56};
    int length = 12;

    if (longflag)
    {
        for (int i=0; i<length; i++)
            test_run_sieve<long>(dim[i]);
    }

    else
    {
        test_run_sieve<mpz_t>(dim[0]);
    }

    //auto start = std::chrono::high_resolution_clock::now();
//    auto finish = std::chrono::high_resolution_clock::now();
//    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
//    cout << " Time taken: " << microseconds.count()/1000000.0 << "sec" << endl;
  return 1;
}
