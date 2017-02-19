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

    //cout << "before the reduction: ";
    //p1.printLatticePoint();
    assert(p1.norm2 >= p2.norm2);
    ET sc_prod, abs_2scprod, scalar;
    sc_product(sc_prod, p1, p2);
    abs_2scprod.mul_ui(sc_prod,2);
    abs_2scprod.abs(abs_2scprod);

    // check if |2 * <p1, p2>| <= |p2|^2. If yes, no reduction
    if (abs_2scprod <= p2.norm2)
        return false;

    // compute the (integer) multiple for p1: mult = round(<p1, p2> / |p2|^2)
    FP_NR<double> mult, tmp; //may be can use another type
    mult.set_z(sc_prod); //conversions
    tmp.set_z(p2.norm2);


    mult.div(mult, tmp);
    mult.rnd(mult);
    scalar.set_f(mult); //converts mult to the type suitable for mult_const;


    LatticePoint<ET> res(p2);
    scalar_mult(res, scalar);
    p1 = p1 - res;

    //cout << endl << "after the reduction: ";
    //p1.printLatticePoint();
    return true;
}

//helper function for reading in from streams. Gobbles up str from the stream (and whitespace before/after).
//If str is not on the stream, outputs an error.
//Note that whitespace inside str is OK, as long as it is not at the beginning or end.

bool string_consume(istream &is, std::string const & str, bool elim_ws, bool verbose)
{
    unsigned int len = str.length();
    char *buf = new char[len+1];
    buf[len] = 0; //for error message.
    if (elim_ws)
    {
        is >> std::ws;
    }
    is.read(buf,len);
    if(is.gcount() != len)
    {
        if(verbose)
        {
            cerr << "Failure reading header: Expected to read" << str << endl;
            cerr << "Read only "<<is.gcount() << "bytes. String read was" << buf<<endl;
        }
        return false;
    }
    if(elim_ws)
    {
        is >> std::ws;
    }
    if(str.compare(0,len,buf,len)!=0)
    {
        if(verbose)
        {
            cerr << "Failure reading header: Expected to read" << str << endl;
            cerr << "Read instead:" << buf << endl;
        }
        return false;
    }
    return true;
}

#endif

template<class ET> //ET : underlying entries of the vectors. Should be a Z_NR<foo> - type. Consider making argument template itself.
Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::Sieve(LatticeBasisType B, unsigned int k, TerminationConditions<ET> termcond, unsigned int verbosity_, int seed_sampler):  //move to cpp //TODO:MT
    main_list(),
    main_queue(this),
    original_basis(B),
    lattice_rank(B.get_rows()),
    ambient_dimension(B.get_cols()), //Note : this means that rows of B form the basis.
    multi_threaded_wanted(GAUSS_SIEVE_IS_MULTI_THREADED),
    sieve_k(k),
    sampler(nullptr),
    verbosity(verbosity_),
    term_cond(termcond),
    sieve_status(SieveStatus::sieve_status_init),
    shortest_vector_found(), //TODO : initialize to meaningful value, e.g. any vector of B.
    number_of_collisions(0),
    number_of_points_sampled(0),
    number_of_points_constructed(0),
    current_list_size(0)
{
    sampler = new KleinSampler<typename ET::underlying_data_type, FP_NR<double>>(B, verbosity, seed_sampler);
    //unsigned int n = lattice_rank;
    //auto it = main_list.before_begin();
    //assert(main_list.empty()); We don't have a function to check that yet...
    auto it = main_list.cbegin();
    for (unsigned int i=0; i<lattice_rank; ++i)
    {
        //LatticePoint<ET> p (  conv_matrixrow_to_lattice_point (original_basis[i]) );
        main_list.insert_before(it,  conv_matrixrow_to_lattice_point (original_basis[i])  );
        //it = main_list.insert_after(it, p);
    }
    current_list_size+=lattice_rank;
    main_list.sort();



//TODO : initialize term_condition to some meaningful default.
};

template<class ET>
bool Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::check_if_done()
{
    if(term_cond.do_we_use_default_condition())
    {
        //compute GSO for original_basis
        ZZ_mat<mpz_t> Empty_mat;
        MatGSO<ET, FP_NR<double>> BGSO(original_basis, Empty_mat, Empty_mat, 0);
        bool upd = BGSO.update_gso();
        
        // returns det(B)^{2/dim}
        FP_NR<double> det = BGSO.get_root_det (1, original_basis.get_rows());
        
        //lambda_1^2 = n * det(B)^{2/n}
        FP_NR<double> MinkBound_double = original_basis.get_rows() * det;
        ET Minkowski;
        Minkowski.set_f(MinkBound_double);
        
        cout << "set Mink. bound to: " << Minkowski << endl;
        
        term_cond.set_target_length(Minkowski);
        //compute Minkwoski
        //ET det;

        //FT MatGSO< ZT, FT >::get_root_det in gso.cpp
        //term_cond(set_target_length(Minkowski))
    }
    if(term_cond.do_we_check_length())
    {
        
    }
    //...
    //return ...
}

#define SIEVE_JOINT_CPP
#endif
