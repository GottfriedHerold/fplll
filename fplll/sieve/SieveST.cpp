#ifndef SIEVE_GAUSS_SINGLE_THREADED_CPP
#define SIEVE_GAUSS_SINGLE_THREADED_CPP

#if !defined(SIEVE_GAUSS_SINGLE_THREADED) || GAUSS_SIEVE_IS_MULTI_THREADED == true
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
    //for ( LatticePoint<ET> & x : main_list) cout << x << endl;

    /* can do main_list.sort here, but I assume original_basis is preporcessed

     */

    int i=0;
    int MaxIteration = 30;

    LatticePoint<ET> p;
    NumVect<ET> sample;

    while(i < MaxIteration) // TerminationCondition Here
    {
        if (main_queue.empty())
        {
            sample = sampler -> sample();
            p = conv_sample_to_lattice_point(sample);
            cout << "sampled p: ";
            p.printLatticePoint();
            cout<< endl;
        }
        else
        {
            p = main_queue.top();
            main_queue.pop();
            cout << "popped p: ";
            p.printLatticePoint();
            cout<< endl;
        }


	
        
	SieveIteration(p);

        ++i;
    }


}


template<class ET>
void Sieve<ET,false>::SieveIteration (LatticePoint<ET> &p)
{
    auto it1= main_list.before_begin();
    auto prev = main_list.before_begin();
    bool loop = true;
    

    while (loop) //while p keeps changing
    {
        loop = false;
        for (it1 = main_list.begin(); it1!=main_list.end(); ++it1) {
            if (p.norm2 < (*it1).norm2) {
                break;
            }
            if(check2red(p, *it1)) //p was changed
                loop = true;

        }
    }

    if (p.norm2 == 0)
	{
		//increase the number of collisions
		number_of_collisions++;
		return;
	}
    
    
    //insert p into main_list; 
//     cout << "----- before the insertion ---- " << endl;
//     for (auto it2 = main_list.begin(); it2!=main_list.end(); ++it2) {
//    	(*it2).printLatticePoint();
//    }
    it1 = main_list.emplace_after (it1, p);
//
//    cout << "----- after the insertion ---- " << endl;
//     for (auto it2 = main_list.begin(); it2!=main_list.end(); ++it2) {
//    	(*it2).printLatticePoint();
//    }

    
    prev = it1;
    if(it1!=main_list.end()) ++it1;
    
    while (it1 !=main_list.end()) {
		if (check2red(*it1, p)) //*it was changed, remove it from the list, put
    		{
			cout << "v was found" <<  endl;
                
            if ((*it1).norm2 == 0)
            {
                number_of_collisions++;
                return;
            }
                
			main_queue.emplace(*it1);
			it1 = main_list.erase_after(prev); //+it1 is done
			
		}
		else  
		{
			prev = it1;
			++it1;
		}
		
    }
   
    /* print for debugging */
    //for (it1 = main_list.begin(); it1!=main_list.end(); ++it1) {
    //	(*it1).printLatticePoint();
    //}

};

#endif
