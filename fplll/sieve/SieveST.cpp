#ifndef SIEVE_GAUSS_SINGLE_THREADED_CPP
#define SIEVE_GAUSS_SINGLE_THREADED_CPP

#if !defined(SIEVE_GAUSS_SINGLE_THREADED) || GAUSS_SIEVE_IS_MULTI_THREADED == true
#error SieveST.cpp included with wrong macro settings
#endif

template<class ET>
void Sieve<ET,false>::run_2_sieve()
{
    //using SieveT = Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>;

    //want to put my basis-vectors into the main_list
    // -- This is now done by the constructor of Sieve -- Gotti


    //for ( LatticePoint<ET> & x : main_list) cout << x.norm2 << endl;

    int i=0;
    //int MaxIteration = 8000;

    LatticePoint<ET> p;
    //NumVect<ET> sample;

    check_if_done(); //sets up default conditions if not already set. We ignore the return value.
    //ET target_norm = term_cond.get_target_length();
    while (!check_if_done() )
    //while(main_list.cbegin()->norm2 > target_norm)
    {
        p=main_queue.true_pop();

        SieveIteration2(p);

        ++i;
        if (i % 200 == 0) {
            //print_status();
            //cout << "# of collisions: " << number_of_collisions << endl;
            cout << "norm2 of the so far shortest vector: " << main_list.cbegin()->norm2 << endl;

        }
    }

    cout << "sv is " << endl;
    main_list.cbegin()->printLatticePoint();
    print_status();

}


template<class ET>
void Sieve<ET,false>::SieveIteration2 (LatticePoint<ET> &p)
{

    //simplified the code, because main_list supports deleting AT pos and inserting BEFORE pos now. -- Gotti

    bool loop = true;

    typename MainListType::Iterator it_comparison_flip=main_list.cend(); //used to store the point where the list elements become larger than p.

    while (loop) //while p keeps changing
    {
        loop = false;
        for (auto it = main_list.cbegin(); it!=main_list.cend(); ++it)
        {
            if (p.norm2 < (*it).norm2)
            {
                it_comparison_flip = it;
                break;
            }
            if(GaussSieve::check2red(p, *it)) //p was changed
            {
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
        auto current_list_point = *it;
		if (GaussSieve::check2red(current_list_point, p)) //We can reduce *it.
			{
			//cout << "v was found" <<  endl;

            if (current_list_point.norm2 == 0)
            {
                number_of_collisions++;
                break;
            }

			main_queue.push(current_list_point);
			it = main_list.erase(it); //this also increases the iterator it
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
