#ifndef SIEVE_GAUSS_MULTI_THREADED_CPP
#define SIEVE_GAUSS_MULTI_THREADED_CPP

#if !defined(SIEVE_GAUSS_MULTI_THREADED) || GAUSS_SIEVE_IS_MULTI_THREADED == false
#error SieveMT.cpp included with wrong macro settings
#endif

template<class ET>
void Sieve<ET,true>::set_num_threads(unsigned int t)
{
    assert(t>0);
    assert(!check_whether_sieve_is_running()); //can only change number of threads in suspended state.
    delete[] garbage_bins;
    num_threads_wanted = t;
    garbage_bins = new GarbageBin<typename MainListType::DataType>[num_threads_wanted];
}

template<class ET>
void Sieve<ET,true>::run_2_sieve()
{
    std::vector<std::thread> threads;
    for (unsigned int i=0;i<num_threads_wanted;++i )
    {
        threads.emplace_back( &Sieve<ET,true>::sieve_2_thread, this, i);
    }
    for (unsigned int i=0;i<num_threads_wanted;++i)
    {
        threads.at(i).join();
    }
}

template<class ET>
void Sieve<ET,true>::sieve_2_thread(int const thread_id)
{
    cout << "Starting thread " << thread_id << endl;
    int const n = get_ambient_dimension();
    while(!check_if_done() )
    {
        auto p=main_queue.true_pop();
        if(p.norm2==0) continue;
        ApproxLatticePoint<ET,false> pApprox (p);

        bool loop = true;

        typename MainListType::Iterator it_comparison_flip=main_list.cend(); //used to store the point where the list elements become larger than p.

        while (loop) //while p keeps changing
        {
            loop = false;
            for (auto it = main_list.cbegin(); it!=main_list.cend(); ++it)
            {
                if (p.norm2 < it.get_true_norm2())
                {
                    it_comparison_flip = it;
                    break;
                }
                bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );
                if(!predict) continue;
                else
                {
                    auto other_point = it.get_exact_point();
                    if(GaussSieve::check2red(p, other_point ) ) //p was changed
                    {
                        if(p.norm2!=0)  pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);
                        loop = true;
                        break;
                    }
                }
            }
        }
        //p no longer changes. it_comparison_flip is iterator to first (shortest) element in the list that is longer than p.
        if (p.norm2 == 0)
        {
            //increase the number of collisions
            number_of_collisions++;
            continue; // refers to outer while(!check_if_done()) - loop
        }
        //insert p into main_list;
        main_list.insert_before(it_comparison_flip,p);
        ++current_list_size;
        for(auto it = it_comparison_flip; it!=main_list.cend(); )   //go through rest of the list.
                                                                    //We know that every element current_list_point=*it is at least as long as p, so we reduce x using p.
        {
            bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,pApprox.get_approx_norm2(),2* pApprox.get_length_exponent()-1,n   );
            if(!predict){++it;continue;} //if prediction is bad, don't even bother to reduce.
            LatticePoint<ET> current_list_point = it.get_exact_point();
            if (GaussSieve::check2red(current_list_point, p)) //We can reduce *it.
            {
			//if (!predict) cerr << "Misprediction 2" << endl;
			//cout << "v was found" <<  endl;

                if (current_list_point.norm2 == 0)
                {
                    number_of_collisions++;
                    break;
                }

                main_queue.push(current_list_point);
                main_list.unlink(it,garbage_bins[thread_id]); //does not invalidate it.
                ++it;
                --current_list_size;

            }
            else
            {
                ++it;
            }
        }
    }
}


template<class ET>
void Sieve<ET,true>::SieveIteration2 (LatticePoint<ET> &p) //note : Queue might output Approx ?
{
    if (p.norm2==0) return; //cerr << "Trying to reduce 0"; //TODO: Ensure sampler does not output 0 (currently, it happens).
    ApproxLatticePoint<ET,false> pApprox (p);
    //simplified the code, because main_list supports deleting AT pos and inserting BEFORE pos now. -- Gotti
    int const n = get_ambient_dimension();
    bool loop = true;

    typename MainListType::Iterator it_comparison_flip=main_list.cend(); //used to store the point where the list elements become larger than p.

    while (loop) //while p keeps changing
    {
        loop = false;
        for (auto it = main_list.cbegin(); it!=main_list.cend(); ++it)
        {
            if (p.norm2 < it.get_true_norm2())
            {
                it_comparison_flip = it;
                break;
            }
	    //Diagnosis(this, pApprox, p, *it, *(it.access_details()) ,n);
            //if(!LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-2,n   ) ) continue;
	    bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );
	    if(!predict) continue;
            if(GaussSieve::check2red(p, *(it.access_details()) ) ) //p was changed
            {
		        if(p.norm2!=0)  pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);
                loop = true;
                break;
            }
        }
    }

    //p no longer changes. it_comparison_flip is iterator to first (shortest) element in the list that is longer than p.
    //If no such element exists, it_comparison_flip refers to after-the-end.

    if (p.norm2 == 0)
	{
		//increase the number of collisions
		number_of_collisions++;
		return;
	}


    //insert p into main_list;
    main_list.insert_before(it_comparison_flip,p);
    ++current_list_size;
    for(auto it = it_comparison_flip; it!=main_list.cend(); ) //go through rest of the list.
                                                              //We know that every element current_list_point=*it is at least as long as p, so we reduce x using p.
    {

        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,pApprox.get_approx_norm2(),2* pApprox.get_length_exponent()-1,n   );
        if(!predict){++it;continue;} //if prediction is bad, don't even bother to reduce.
        LatticePoint<ET> current_list_point = it.get_exact_point();
        if (GaussSieve::check2red(current_list_point, p)) //We can reduce *it.
		{
			//if (!predict) cerr << "Misprediction 2" << endl;
			//cout << "v was found" <<  endl;

            if (current_list_point.norm2 == 0)
            {
                number_of_collisions++;
                break;
            }

			main_queue.push(current_list_point);
			it = main_list.erase(it); //this also moves it forward
			--current_list_size;

		}
		else
		{
		//	prev = it1;
			++it;
		}

    }

    /* print for debugging */
    //for (it1 = main_list.begin(); it1!=main_list.end(); ++it1) {
    //	(*it1).printLatticePoint();
    //}

};


#endif
