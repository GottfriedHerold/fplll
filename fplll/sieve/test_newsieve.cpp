
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
    
    
    //this is ugly
    ZZ_mat<ZT> Empty_mat;
    MatGSO<Z_NR<ZT>, FP_NR<double>> BGSO(B, Empty_mat, Empty_mat, 0);
    bool upd = BGSO.update_gso();
    
    // returns det(B)^{2/dim}
    FP_NR<double> det = BGSO.get_root_det (1, B.get_rows());
    
    
    
    cout << "det = " << det << endl;
    //cout << "intdet = " << intdet << endl;
    
    FP_NR<double> MinkBound = sqrt(B.get_rows() * det); // \lambda_1(B) <= \sqrt(n * det(B)^{2/n} )
    cout << "MindBound = " << MinkBound << endl;
    Z_NR<ZT> Minkowski;
    Minkowski.set_f(MinkBound);
    
    //Test_Queue.run_2_sieve(Minkowski); //file-name to output the results of tests
    {
        //std::ofstream ofs("FILENAME");
        //Test_Queue.print_status(-1,ofs);
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
