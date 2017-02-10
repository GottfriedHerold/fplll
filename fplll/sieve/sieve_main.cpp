/*
  This provides an implementation of Gauss sieving, including using
  tuples of vectors in fplll framework. The Gauss Sieve code is
  based on Panagiotis Voulgaris's implementation of the Gauss sieve.
*/
#include "sieve_main.h"
#include "fplll.h"
#include "LatticePoint.h"
#include "PointList.h"
#include <thread>
#include <chrono>
#include "SieveGauss.h"

using namespace fplll;

/**
 * help function
 */
//static void main_usage(char *myself)
//{
//  cout << "Usage: " << myself << " [options]\n"
//       << "List of options:\n"
//       << "  -a [2|3|4]\n"
//       << "     2- or 3- or 4-sieve;\n"
//       << "  -f filename\n"
//       << "     Input filename\n"
//       << "  -r nnn\n"
//       << "     Generate a random instance of dimension nnn\n"
//       << "  -t nnn\n"
//       << "     Targeted norm^2=nnn\n"
//       << "  -s nnn\n"
//       << "     Using seed=nnn\n"
//       << "  -b nnn\n"
//       << "     BKZ preprocessing of blocksize=nnn\n"
//       << "  -v\n"
//       << "     Verbose mode\n";
//  exit(0);
//}

/**
 * run sieve
 */

#if 0
template <class ZT> int main_run_sieve(ZZ_mat<ZT> B, Z_NR<ZT> goal_norm, int alg, int ver, int seed)
{
  /*
   GaussSieve<ZT, FP_NR<double>> gsieve(B, alg, ver, seed);
  gsieve.set_goal_norm2(goal_norm);
  if (gsieve.alg == 3)
    gsieve.run_3sieve();
  else if (gsieve.alg == 4)
    gsieve.run_4sieve();
  else
    gsieve.run_2sieve();
   */
  return 0;
}
#endif

/**
 * main function
 */


template<class DT>
void ListTester(ListMultiThreaded<DT> * const Z, GarbageBin<DT> * const gb, int id,int verbose)
{
for (int i=0; i <25; ++i)
{
    int count =0;
    int insertions=0;
    int deletions =0;
    for(MTListIterator<DT> it = Z->begin(); !it.is_end(); ++it)
    {
        ++count;
        DT tmp = *it;
        if(( tmp * tmp + i + 33*id) % 39 == 11)
        {
        Z->unlink(it,*gb);
        if (verbose>=2) cout << "Deleted";
        ++deletions;
        }
        if(( tmp * tmp + i + 18*id) %57 ==13  )
        {
        Z->insert(it, 1000000*id + 100* i + ((tmp * tmp )%93) );
        if (verbose>=2) cout << "Inserted";
        ++insertions;
        }
    }
    if (verbose>=1)
    {
    cout << "Thread " << id <<", iteration " << i << "count: " <<count <<" deletes: "<<deletions << "inserts: "<< insertions<<endl;
    }
    MTListIterator<DT> it = Z->end();
    Z->insert(it, 1000000*id + 100*i + 1);
}
return;
}

template<class ZT> void call_sieve (ZZ_mat<ZT> B, int sieve_verbosity)
{

	TerminationConditions<ZT> term_cond();
	//Sieve<ZT, false> Test_2Sieve (B, term_cond, 0, 0, sieve_verbosity);
	//Sieve<ZT, false> Test_2Sieve;
	//Test_2Sieve.run_2_sieve();
	
	
} 


int main(int argc, char **argv)
{
  //char *input_file_name = NULL;
  //char *goal_norm_s     = NULL;
  //bool flag_verbose = false, flag_file = false;
  //int option, alg, dim = 10, seed = 0, bs = 0;
    cout<< "Hello" << endl;


    //PointListSingleThreaded<long int> X;
//    PointListMultiThreaded<long int> Y;
    //printLatticePoint(p);
    //ListMTNode<int> Z;

    //LatticePoint<long int> p = LatticePoint<long int> (10, 2);
    LatticePoint<Z_NR <long> > p1 (10,3);
    LatticePoint<Z_NR <long> > p2 (10,1);
    Z_NR<long> res;
    sc_product(res,p1,p2);
    //Z_NR<long> X (0);
    //p1.printLatticePoint();
//    LatticePoint<Z_NR <long> > sum = p1+p2;
//    sum.printLatticePoint();
//    LatticePoint<Z_NR <long> > sub = sum - p2;
//    sub.printLatticePoint();
//    cout << res << endl;
    
    // ZZ_mat is an integer row-oriented matrix. See /nr/matrix.h
    ZZ_mat<mpz_t> B;
    int dim = 10;
    B.resize(dim, dim);
    
    //generates a lower-triangular matrix B; the argument determines the bit-size of entries
    B.gen_trg(1.1);
    
    
    // KleinSampler<ZT, F> is templated by two classes; returns NumVect<Z_NR<ZT> of dim = B.NumCols()
    KleinSampler<mpz_t, FP_NR<double>> *Sampler = new KleinSampler<mpz_t, FP_NR<double>>(B, 0, 0);
    NumVect<Z_NR<mpz_t> > sample(dim);
    sample = Sampler->sample();
    cout << sample << endl;

    call_sieve(B, 1);
    
    
    int num_threads=4;
    ListMultiThreaded<int> Z;
    GarbageBin<int>* GarbageBins = new GarbageBin<int> [num_threads];
    std::vector <std::thread> threads;
    //ListTester<int>(&Z,GarbageBins[2],2);
    auto start = std::chrono::high_resolution_clock::now();
    cout<< "Starting threads." << endl;
    for(int i=0; i <num_threads; ++i)
    {
    std::thread tmp(ListTester<int>, &Z, GarbageBins+i,i,0);
    tmp.join();
        //threads.push_back(std::thread(ListTester<int>, &Z, GarbageBins+i,i,0));
    //for(auto &th: threads) th.join();
    }

    cout << "Waiting for them to finish" <<endl;
    //for(auto &th: threads) th.join();
    cout << "Finished" << endl;
    auto finish = std::chrono::high_resolution_clock::now();
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
    cout << " Time taken: " << microseconds.count() << endl;
    delete[] GarbageBins;
#if 0
#if 0
  dot_time = 0;
  dot_num = 0;
  count_bad = 0;
#endif

  alg = 2;

  /* parse */
  if (argc == 1)
  {
    main_usage(argv[0]);
    return -1;
  }
  while ((option = getopt(argc, argv, "a:f:r:t:s:b:v")) != -1)
  {
    switch (option)
    {
    case 'a':
      alg = atoi(optarg);
      if (alg != 2 && alg != 3 && alg != 4)
      {
        cout << " Error, only support 2-, 3- and 4-sieve" << endl;
        exit(1);
      }
    case 'f':
      input_file_name = optarg;
      flag_file       = true;
      break;
    case 'r':
      dim       = atoi(optarg);
      flag_file = false;
      break;
    case 's':
      seed = atoi(optarg);
      break;
    case 'b':
      bs = atoi(optarg);
      break;
    case 'v':
      flag_verbose = true;
      break;
    case 't':
      // ngoal_norm = atol(optarg);
      cout << optarg << endl;
      goal_norm_s = optarg;
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

  /* set lattice */
  ZZ_mat<mpz_t> B;
  if (flag_file)
  {
    ifstream input_file(input_file_name);
    if (input_file.is_open())
    {
      input_file >> B;
      input_file.close();
    }
    else
    {
      cin >> B;
    }
    if (flag_verbose)
    {
      cout << "# [info] reading lattice of dimension " << B.get_rows() << "x" << B.get_cols()
           << endl;
    }
  }
  else
  {
    if (flag_verbose)
    {
      cout << "# [info] generating random lattice of dimension " << dim << endl;
    }
    srand(time(NULL));
    B.resize(dim, dim);
    B.gen_trg(1.1);
  }

  /* set targeted norm */
  Z_NR<mpz_t> goal_norm, max;
  if (goal_norm_s != NULL)
  {
    goal_norm.set_str(goal_norm_s);
  }
  if (goal_norm < 0)
    goal_norm = 0;
  if (flag_verbose)
    cout << "# [info] goal norm^2 is " << goal_norm << endl;

  /* preprocessing of basis */
  clock_t stime = clock();
  if (bs > 0)
    bkz_reduction(B, bs, BKZ_DEFAULT, FT_DEFAULT, 0);
  else
    lll_reduction(B, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);

  clock_t etime = clock();
  double secs   = (etime - stime) / (double)CLOCKS_PER_SEC;
  if (flag_verbose)
  {
    if (bs > 0)
      cout << "# [info] BKZ took time " << secs << " s" << endl;
    else
      cout << "# [info] LLL took time " << secs << " s" << endl;
  }
  // cout << B << endl;

  /* decide integer type */
  stime = clock();
  max   = B.get_max();

#if 1
  if (max < std::numeric_limits<int>::max())
  {
    long goal_norm_l = abs(goal_norm.get_si());
    Z_NR<long> goal_norm_lt;
    goal_norm_lt = goal_norm_l;
    ZZ_mat<long> B2(B.get_rows(), B.get_cols());
    for (int i = 0; i < B.get_rows(); i++)
      for (int j = 0; j < B.get_cols(); j++)
        B2(i, j) = B(i, j).get_si();
    main_run_sieve<long>(B2, goal_norm_lt, alg, flag_verbose, seed);
  }
  else
#endif
 //   main_run_sieve<mpz_t>(B, goal_norm, alg, flag_verbose, seed); //commented out for test

  etime = clock();
  secs  = (etime - stime) / (double)CLOCKS_PER_SEC;
  if (flag_verbose)
  {
    cout << "# [info] sieve took time " << secs << " s" << endl;
/* dot product time */
#if 0
    cout << "# [info] dot_time " << dot_time << endl;
    cout << "# [info] dot_num " << dot_num << endl;
    cout << "# [info] dot_time/dot_number " << (double) dot_time/dot_num << endl;
#endif
  }
#endif
  return 1;
}
