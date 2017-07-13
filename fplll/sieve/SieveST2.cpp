/* DO NOT INCLUDE THIS FILE DIRECTLY
*/

template<class ET, int nfixed> void Sieve<ET,false,nfixed>::sieve_2_iteration (FastAccess_Point &p) //note : Queue might output Approx ?
{
    int const n = get_ambient_dimension();
    
    /*
    if (p.get_norm2() == 0) return;
    bool loop = true;
    
    typename MainListType::Iterator it_comparison_flip=main_list.cend(); 
    */
    
    
    /*
    if (p.access_exact_point_r().norm2==0) return; //cerr << "Trying to reduce 0"; //TODO: Ensure sampler does not output 0 (currently, it happens).
    ApproximateLatticePoint<ET,nfixed> p_approx (p.access_approximation_r(), n); //makes a copy.
    ET p_exact_norm = p.get_exact_norm2();
    ExactLatticePoint<ET,nfixed> p_exact = p.get_exact_point();

    bool loop = true;

    ET scalar; //reduction multiple output by check2red_new

    typename MainListType::Iterator it_comparison_flip=main_list.cend(); //used to store the point where the list elements become larger than p.
    
    while (loop) //while p keeps changing
    {
        loop = false;
        for (auto it = main_list.cbegin(); it!=main_list.cend(); ++it)
        {
            if (p_exact_norm < it.get_true_norm2())
            {
                it_comparison_flip = it;
                break;
            }
        ++number_of_scprods;
        bool const predict =GaussSieve::compare_abs_approx_scalar_product(p_approx,*it,it->get_norm2_mantissa(), it->get_norm2_exponent() -1,n); //the -1 acts a a constant factor 1/2, related to our target scalar product.
//        bool predict = LatticeApproximations::Compare_Sc_Prod(p_approx,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );
	    if(!predict) continue;

            ++number_of_exact_scprods;
            if ( GaussSieve::check2red_exact<ET,nfixed>(p_exact, it.dereference_exactly_r(), scalar) )
            {
                p = GaussSieve::perform2red_exact_to_compressed<ET,false,nfixed>(p_exact,it.dereference_exactly_r(),scalar);
                p_approx.replace_by(p.access_approximation_r(), n);
                p_exact  = p.get_exact_point();
                p_exact_norm = p_exact.norm2;
                //p = GaussSieve::perform2red(p, *(it.access_details()), scalar);
            //update the approximation of f
                //if (p.norm2!=0) p_approx = static_cast< ApproxLatticePoint<ET,false> >(p);
                loop = true;
                break;
            }
            else
            {
                ++number_of_mispredictions;
            }
        }
    }

    //p no longer changes. it_comparison_flip is iterator to first (shortest) element in the list that is longer than p.
    //If no such element exists, it_comparison_flip refers to after-the-end.

    if (p_exact_norm == 0) //essentially means that p was already inside the list.
	{
        //cout << "p has norm 2" << endl;
		number_of_collisions++;
		return;
	}

    //insert p into main_list;
    main_list.insert_before(it_comparison_flip,p.deep_copy_compressed_point()); //FIXME: This should cause an error without .deep_copy, but doesn't. Something is wrong!
    ++current_list_size;
    if(update_shortest_vector_found(p_exact))
    {
        if(verbosity>=2)
        {
            cout << "New shortest vector found. Norm2 = " << get_best_length2() << endl;
        }
    }

    for(auto it = it_comparison_flip; it!=main_list.cend(); ) //go through rest of the list.
                                                              //We know that every element current_list_point=*it is at least as long as p, so we reduce x using p.
    {
        ++number_of_scprods;
        bool const predict = GaussSieve::compare_abs_approx_scalar_product(p_approx,*it,p_approx.get_norm2_mantissa(),p_approx.get_norm2_exponent()-1,n);
        if(!predict){++it;continue;} //if prediction is bad, don't even bother to reduce.
        //We believe we can reduce *it
        ExactLatticePoint<ET,nfixed> current_list_point = it.dereference_exactly_r();
        ++number_of_exact_scprods;
        if (GaussSieve::check2red_exact(current_list_point, p_exact, scalar)) //We can reduce *it.
		{
			//create new list point
			CompressedPoint<ET,false,nfixed> reduced_point = GaussSieve::perform2red_exact_to_compressed<ET,false,nfixed>(current_list_point, p_exact, scalar);
			//if (!predict) cerr << "Misprediction 2" << endl;
			//cout << "v was found" <<  endl;

            if (reduced_point.get_exact_norm2() == 0) //Note : this cannot happen unless the list contains a non-trivial multiple of p (because the collision would have triggered before).
            {
                number_of_collisions++;
                ++it;
                continue; //was a break. Changed to ++it;continue; -- Gotti
            }

			main_queue.push(std::move(reduced_point));
			it = main_list.erase(it); //this also moves it forward
			--current_list_size;
		}
		else
		{
            ++number_of_mispredictions;
		//	prev = it1;
			++it;
		}
    }

    /* print for debugging */
    //for (it1 = main_list.begin(); it1!=main_list.end(); ++it1) {
    //	(*it1).printLatticePoint();
    //}
}
