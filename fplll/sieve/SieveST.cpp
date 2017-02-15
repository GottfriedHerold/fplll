#ifndef SIEVE_GAUSS_SINGLE_THREADED_CPP
#define SIEVE_GAUSS_SINGLE_THREADED_CPP

#if !defined(GAUSS_SIEVE_SINGLE_THREADED) || defined(GAUSS_SIEVE_MULTI_THREADED) || GAUSS_SIEVE_IS_MULTI_THREADED == true
#error SieveST.cpp included with wrong macro settings
#endif


template<class ET>
void Sieve<ET,false>::run_2_sieve()
{
    //using SieveT = Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>;
    
    //want to put my basis-vectors into the main_list --Make a separate function for that
    
    unsigned int n = lattice_rank;
    auto it = main_list.before_begin();
    for (unsigned int i=0; i<n; i++) {
        LatticePoint<ET> p (  conv_matrixrow_to_lattice_point (original_basis[i]) );
        it = main_list.emplace_after(it, p);
        //it = main_list.insert_after(it, p);
    }
    for ( LatticePoint<ET> & x : main_list) cout << x << endl;
    
    /* can do main_list.sort here, but I assume original_basis is preporcessed
     
     */
    
    int i=0;
    int MaxIteration = 10;
    
    LatticePoint<ET> p;
    NumVect<ET> sample;
    
    while(i < MaxIteration) // TerminationCondition Here
    {
        if (main_queue.empty())
        {
            sample = sampler -> sample();
            p = conv_sample_to_lattice_point(sample);
        }
        else
        {
            p = main_queue.top();
            main_queue.pop();
        }
        
        SieveIteration(p);
        
        ++i;
    }
    
    
}


template<class ET>
void Sieve<ET,false>::SieveIteration (LatticePoint<ET> &p)
{
    //auto it = main_list.before_begin();

    bool loop = true;

    while (loop) //while p keeps changing
    {
        loop = false;
        for (auto it1 = main_list.begin(); it1!=main_list.end(); ++it1) {
            if (p.norm2 < (*it1).norm2) {
                break;
            }
            bool check  = check2red(p, *it1);
            
        }
    }
    //cout << "running SieveIteration " << endl;
};



#endif
