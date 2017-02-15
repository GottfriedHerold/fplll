#undef DO_INCLUDE_SIEVE_JOINT_CPP
#ifndef GAUSS_SIEVE_IS_MULTI_THREADED
#error wrong usage of SieveJoint.cpp -- 1
#endif

#if GAUSS_SIEVE_IS_MULTI_THREADED == false

#if !defined(GAUSS_SIEVE_SINGLE_THREADED) || defined(GAUSS_SIEVE_MULTI_THREADED)
#error wrong usage of SieveJoint.cpp -- 2
#endif

#ifndef SIEVE_JOINT_CPP_ST
#define SIEVE_JOINT_CPP_ST
#define DO_INCLUDE_SIEVE_JOINT_CPP
#endif

#elif GAUSS_SIEVE_IS_MULTI_THREADED == true
#if defined(GAUSS_SIEVE_SINGLE_THREADED) || !defined(GAUSS_SIEVE_MULTI_THREADED)
#error wrong usage of SieveJoint.cpp -- 3
#endif

#ifndef SIEVE_JOINT_CPP_MT
#define SIEVE_JOINT_CPP_MT
#define DO_INCLUDE_SIEVE_JOINT_CPP
#endif

#endif

#ifdef DO_INCLUDE_SIEVE_JOINT_CPP

//actual code starts here.
//Be aware that this code may be read twice.

//template<class ET, GAUSS_SIEVE_IS_MULTI_THREADED>
//void Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::SieveIteration (LatticePoint<ET> &p);

//template<class ET>
//void Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::run_2_sieve()
//{
//    using SieveT = Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>;
//
//    //want to put my basis-vectors into the main_list --Make a separate function for that
//
//    unsigned int n = lattice_rank;
//    auto it = main_list.before_begin();
//    for (unsigned int i=0; i<n; i++) {
//        LatticePoint<ET> p (  conv_matrixrow_to_lattice_point (original_basis[i]) );
//        it = main_list.emplace_after(it, p);
//        //it = main_list.insert_after(it, p);
//    }
//    for ( LatticePoint<ET> & x : main_list) cout << x << endl;
//
//    /* can do main_list.sort here, but I assume original_basis is preporcessed
//
//     */
//
//    int i=0;
//    int MaxIteration = 10;
//
//    LatticePoint<ET> p;
//    NumVect<ET> sample;
//
//    while(i < MaxIteration) // TerminationCondition Here
//    {
//        if (main_queue.empty())
//        {
//            sample = sampler -> sample();
//            p = conv_sample_to_lattice_point(sample);
//        }
//        else
//        {
//            p = main_queue.top(); //why does this return contst?
//            main_queue.pop();
//        }
//
//        SieveIteration(p);
//
//        ++i;
//    }
//
//
//}

#ifndef SIEVE_JOINT_CPP

//template<class ET, GAUSS_SIEVE_IS_MULTI_THREADED>
//void Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::SieveIteration (LatticePoint<ET> &p)
//{
//    //auto it = main_list.before_begin();
//    
//    bool loop = true;
//    
//    while (loop) //while p keeps changing
//    {
//        loop = false;
//        for (auto it1 = main_list.begin(); it1!=main_list.end(); ++it1) {
//            if (p.norm2 < (*it1).norm2) {
//                break;
//            }
//            
//            
//        }
//    }
//    //cout << "running SieveIteration " << endl;
//};

// p1 > p2; reduce p1
template<class ET>
bool check2red (LatticePoint<ET> &p1, const LatticePoint<ET> &p2)
{
    cout << "inside check2red" << endl;
    return true;
}


#endif

#define SIEVE_JOINT_CPP
#endif
