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


template<class ET>
bool GaussSieve::check2red (LatticePoint<ET> &p1, const LatticePoint<ET> &p2)
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

// separate chec2Red and perform2Red

//if true, sclara is the multiple s.t. we reduce p1 = p1-sclar * p2;
template<class ET>
bool GaussSieve::check2red_new (const LatticePoint<ET> &p1, const LatticePoint<ET> &p2, ET scalar)
{
    return false;
}


//helper function for reading in from streams. Gobbles up str from the stream (and whitespace before/after).
//If str is not on the stream, outputs an error.
//Note that whitespace inside str is OK, as long as it is not at the beginning or end.

bool GaussSieve::string_consume(istream &is, std::string const & str, bool elim_ws, bool verbose)
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

Z_NR<mpz_t> GaussSieve::compute_mink_bound(ZZ_mat<mpz_t> const & basis)
{
    assert(basis.get_rows() == basis.get_cols()); //Note : Alg might work even otherwise. This assertion failure is just a reminder that this needs to be checked.
    //compute Gram-Schmidt-Orthogonalization.
    ZZ_mat<mpz_t> Empty_mat;
    ZZ_mat<mpz_t> basis2 = basis; //need to copy, as BGSO is not const-specified...
    MatGSO<Z_NR<mpz_t>, FP_NR<double>> BGSO(basis2, Empty_mat, Empty_mat, 0);
    //cout << "before the update... " << endl;
    BGSO.update_gso();
    //cout << "after the update... " << endl;

    FP_NR<double> entry;

    //for (int i=0; i<basis.get_rows(); i++)
    //{
// 	for (int j=0; j<basis.get_rows(); j++)
//        	//cout << BGSO.get_gram(entry, j, j) << endl;
//		cout << (BGSO.get_r(entry, j, j)) << "  " << log(BGSO.get_r(entry, j, j)) << endl;
//        cout << endl;
//   }

    // returns det(B)^{2/dim}

    FP_NR<double> root_det2 = BGSO.get_root_det (1, basis.get_rows());
    FP_NR<double> log_det2 = BGSO.get_log_det (1, basis.get_rows());
    //cout << "root_det2: " << root_det2 << endl;
    //cout << "log_det2: " << log_det2 << endl;

    //lambda_1^2 = n * det(B)^{2/n}

    FP_NR<double> MinkBound_double = 0.114 * root_det2 * static_cast<double> (basis.get_rows() ); //technically, we need to multiply by Hermite's constant in dim n here. We are at least missing a constant factor here.
    //DUE TO [KL79], the best know multiple (for the squared norm) whould be (pi* exp(1)*2^{2*0.099} ~ 0.102).

    //cout << "after MinkBound_double is assigned... " << endl;
    Z_NR<mpz_t> Minkowski;
    cout << "MinkBound_double: " << MinkBound_double << endl;
    Minkowski.set_f(MinkBound_double); // the problem is here
cout << "after MinkBound is set... " << endl;
    return Minkowski;
}

//End of things included only once.
#endif // SIEVE_JOINT_CPP

template<class ET> //ET : underlying entries of the vectors. Should be a Z_NR<foo> - type. Consider making argument template itself.
Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::Sieve(LatticeBasisType B, unsigned int k, TerminationConditions<ET> termcond, unsigned int verbosity_, int seed_sampler):  //move to cpp //TODO:MT
    main_list(),
    main_queue(this),
    original_basis(B),
    lattice_rank(B.get_rows()),
    ambient_dimension(B.get_cols()), //Note : this means that rows of B form the basis.
    multi_threaded_wanted(GAUSS_SIEVE_IS_MULTI_THREADED),
    #if GAUSS_SIEVE_IS_MULTI_THREADED == true
    num_threads_wanted(1),
    #endif // GAUSS_SIEVE_IS_MULTI_THREADED
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
    //D1(0),D2(0),D3(0),D4(0),D5(0),D6(0),D7(0),D8(0),D9(0),D10(0)

{
    sampler = new KleinSampler<typename ET::underlying_data_type, FP_NR<double>>(B, verbosity, seed_sampler);
    //unsigned int n = lattice_rank;
    //auto it = main_list.before_begin();
    //assert(main_list.empty()); We don't have a function to check that yet...
    cout <<"Initializing list with original basis..." << endl << std::flush;
    auto it = main_list.cbegin();
    for (unsigned int i=0; i<lattice_rank; ++i)
    {
        main_list.insert_before(it,  conv_matrixrow_to_lattice_point (original_basis[i])  );
    }
    current_list_size+=lattice_rank;
    cout << "Sorting ..." << std::flush;
    main_list.sort();
    cout << "is finished." << endl << std::flush;


//TODO : initialize term_condition to some meaningful default.
};

template<class ET>
bool Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::check_if_done()
{
    if(term_cond.do_we_use_default_condition())
    {
        //assert(false); //DOES NOT WORK, since compute_mink_bound does not work.
        cout << original_basis.get_cols();
        ET Minkowski = GaussSieve::compute_mink_bound(original_basis);
        if (verbosity>=1) cout << "set Mink. bound to: " << Minkowski << endl;
        term_cond.set_target_length(Minkowski);

        //FT MatGSO< ZT, FT >::get_root_det in gso.cpp
        //term_cond(set_target_length(Minkowski))
    }
    if(term_cond.do_we_check_length())
    {
            if (get_best_length2() <= term_cond.get_target_length()) return true; //TODO : Use current_best or somesuch.
    }
    return false;
}

template<class ET>
void Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::run() //runs the sieve specified by the parameters.
{
    assert(sieve_k == 2); //for now
    sieve_status =SieveStatus::sieve_status_running;
    run_2_sieve();
    sieve_status = SieveStatus::sieve_status_finished;
}



#define SIEVE_JOINT_CPP
#endif
