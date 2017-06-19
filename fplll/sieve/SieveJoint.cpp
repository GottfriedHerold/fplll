#undef DO_INCLUDE_SIEVE_JOINT_CPP
#ifndef GAUSS_SIEVE_IS_MULTI_THREADED
#error wrong usage of SieveJoint.cpp --1
#endif

#if GAUSS_SIEVE_IS_MULTI_THREADED == false

#if !defined(SIEVE_GAUSS_SINGLE_THREADED)
#error wrong usage of SieveJoint.cpp -- 2
#endif

#ifndef SIEVE_JOINT_CPP_ST
#define SIEVE_JOINT_CPP_ST
#define DO_INCLUDE_SIEVE_JOINT_CPP
#endif

#elif GAUSS_SIEVE_IS_MULTI_THREADED == true
#if !defined(SIEVE_GAUSS_MULTI_THREADED)
#error wrong usage of SieveJoint.cpp -- 3
#endif

#ifndef SIEVE_JOINT_CPP_MT
#define SIEVE_JOINT_CPP_MT
#define DO_INCLUDE_SIEVE_JOINT_CPP
#endif

#endif

#ifdef DO_INCLUDE_SIEVE_JOINT_CPP

//actual code starts here.
//Be aware that code may be read twice.

#ifndef SIEVE_JOINT_CPP //code in this block only read once.



//End of things included only once.
#endif // SIEVE_JOINT_CPP

template<class ET, int nfixed> //ET : underlying entries of the vectors. Should be a Z_NR<foo> - type.
#if GAUSS_SIEVE_IS_MULTI_THREADED==true
Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>::Sieve(LatticeBasisType B, unsigned int k, unsigned int num_threads, TermCondType const termcond, unsigned int verbosity_, int seed_sampler):
#else
Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>::Sieve(LatticeBasisType B, unsigned int k, TermCondType const termcond, unsigned int verbosity_, int seed_sampler):
#endif
    main_list(),
    main_queue(this),
//    filtered_list(),
    original_basis(B),
    lattice_rank(B.get_rows()),
    ambient_dimension(B.get_cols()), //Note : this means that rows of B form the basis.
    multi_threaded_wanted(GAUSS_SIEVE_IS_MULTI_THREADED),
    #if GAUSS_SIEVE_IS_MULTI_THREADED == true
    num_threads_wanted(num_threads),
    #endif // GAUSS_SIEVE_IS_MULTI_THREADED
    sieve_k(k),
    verbosity(verbosity_),
    term_cond(termcond),
    sieve_status(SieveStatus::sieve_status_init),
    shortest_vector_found(), //TODO : initialize to meaningful value, e.g. any vector of B.
    number_of_collisions(0),
    number_of_points_sampled(0),
    number_of_points_constructed(0),
    current_list_size(0),
    number_of_scprods(0),
    number_of_exact_scprods(0),
    number_of_mispredictions(0)
    #if GAUSS_SIEVE_IS_MULTI_THREADED==true
    ,garbage_bins(nullptr)
    #endif // GAUSS_SIEVE_IS_MULTI_THREADED

{
    #if GAUSS_SIEVE_IS_MULTI_THREADED==true
    if(num_threads_wanted==0) //0 means we take a meaningful default, which is given by thread::hardware_concurrency
    num_threads_wanted = std::max(std::thread::hardware_concurrency(),static_cast<unsigned int>(1)); //Note: hardware_concurrency might return 0 for "unknown".
    #endif
    //unsigned int n = lattice_rank;
    //auto it = main_list.before_begin();
    //assert(main_list.empty()); We don't have a function to check that yet...
    if (verbosity>=2) {cout <<"Initializing list with original basis..." << endl;}
    auto it = main_list.cbegin();
    for (unsigned int i=0; i<lattice_rank; ++i)
    {
        ExactLatticePoint<ET,nfixed> * new_basis_vector = new ExactLatticePoint<ET,nfixed> ( conv_matrixrow_to_lattice_point<ET,nfixed> (original_basis[i]));
        main_list.insert_before(it,  static_cast<CompressedPoint<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed> >(new_basis_vector) );
    }
    current_list_size+=lattice_rank;
//    #if GAUSS_SIEVE_IS_MULTI_THREADED == false
    if(verbosity>=2)    {cout << "Sorting ...";}
    //main_list.sort();
    if(verbosity>=2)    {cout << "is finished." << endl;}
    shortest_vector_found = main_list.cbegin().dereference_exactly_r();
//    #endif // GAUSS_SIEVE_IS_MULTI_THREADED
    //TODO : enable sorting for multithreaded case.
    #if GAUSS_SIEVE_IS_MULTI_THREADED==true
    garbage_bins = new GarbageBin<typename MainListType::DataType>[num_threads_wanted]; //maybe init later.
    #endif
    main_queue.sampler->init(this);
};

template<class ET,int nfixed>
Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>::~Sieve()
{
    #if GAUSS_SIEVE_IS_MULTI_THREADED==true
    delete[] garbage_bins;
    #endif
};

template<class ET, int nfixed>
bool Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>::check_if_done()
{
    return (term_cond->check(this) != 0)?true:false;
}

//template<class ET>
//void Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::run() //runs the sieve specified by the parameters.
//{
//    assert(sieve_k == 2); //for now
//    sieve_status =SieveStatus::sieve_status_running;
//    check_if_done(); //this updates the Minkowski condition, if needed.
//    run_2_sieve();
//    sieve_status = SieveStatus::sieve_status_finished;
//}

#define SIEVE_JOINT_CPP
#endif
