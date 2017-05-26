#ifndef SIEVE_GAUSS_SINGLE_THREADED_CPP
#define SIEVE_GAUSS_SINGLE_THREADED_CPP

#if !defined(SIEVE_GAUSS_SINGLE_THREADED) || GAUSS_SIEVE_IS_MULTI_THREADED == true
#error SieveST.cpp included with wrong macro settings
#endif

class next_block;
using namespace LatticeApproximations; // to be able to use ApproxTypeNorm2 to store inner-produces scaled by length

class main_list;
template<class ET>
bool Sieve<ET,false>::update_shortest_vector_found(LPType const & newvector)
{
    if(newvector.norm2 < shortest_vector_found.norm2)
    {
        shortest_vector_found = newvector;
        return true;
    }
    return false;
}

template<class ET>
typename Sieve<ET,false>::LPType Sieve<ET,false>::get_shortest_vector_found()
{
    return shortest_vector_found;
}

template<class ET>
ET Sieve<ET,false>::get_best_length2()
{
    return shortest_vector_found.norm2;
}


//to avoid errs from SieveMT
//template<class ET>
//void Sieve<ET,false>::run_2_sieve()
//{
//    run_sieve (2);
//}

//we should re-use this function for k-sieve, but call different SieveInteration
template<class ET>
void Sieve<ET,false>::run()
{
    //using SieveT = Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>;

    //for ( LatticePoint<ET> & x : main_list) cout << x.norm2 << endl;
    
    cout << "the shortest vector in the input basis has norm2 = " << (main_list.cbegin())->get_details().norm2 << endl;
    sieve_status =SieveStatus::sieve_status_running;
    int i=0;
    //int MaxIteration = 8000;
    term_cond->init(this); //initialisation of termination conditions.
    LatticePoint<ET> p;
    cout << "start " << sieve_k << "-Sieve" << endl;

    //ET target_norm = term_cond.get_target_length();
    while (!check_if_done() )
//	while (i<2)
    //while(main_list.cbegin()->norm2 > target_norm)
    {
        p=main_queue.true_pop();
        if (sieve_k==2)
            SieveIteration2(p);
        else if (sieve_k==3)
        {
            
            //SieveIteration3New_Pointer(p); //DOES NOT WORK
            //SieveIteration3New(p);
            SieveIteration3(p);
        }


        //cout << i <<  " list size" << current_list_size << " Queue: " << main_queue.size() << endl << flush;
        ++i;
//        if (i % 500 == 0) {
//            print_status();
//            cout << "# of collisions: " << number_of_collisions << endl;
//            cout << "norm2 of the so far shortest vector: " << get_best_length2() << endl;
//
//        }
    }
    sieve_status = SieveStatus::sieve_status_finished;

    //Diagnostic and use of the information moved up to the caller.

    /*
    cout << "sv is " << endl;
    main_list.cbegin().access_details()->printLatticePoint();
    print_status();
    */
}

template<class ET>
void Sieve<ET,false>::SieveIteration2 (LatticePoint<ET> &p) //note : Queue might output Approx ?
{
    if (p.norm2==0) return; //cerr << "Trying to reduce 0"; //TODO: Ensure sampler does not output 0 (currently, it happens).
    ApproxLatticePoint<ET,false> pApprox (p);
    //simplified the code, because main_list supports deleting AT pos and inserting BEFORE pos now. -- Gotti
    int const n = get_ambient_dimension();
    bool loop = true;

    ET scalar; //reduction multiple output by check2red_new

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
        ++number_of_scprods;
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );
	    if(!predict) continue;

	/* OLD IMPLEMENTATION: check2red both check and reduces p
	    if(GaussSieve::check2red(p, *(it.access_details()) ) ) //p was changed
            {
		        if(p.norm2!=0)  pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);
		cout << "p = ";
		p.printLatticePoint();
                loop = true;
                break;
            }
	*/
            ++number_of_exact_scprods;
            if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
            {
                p = GaussSieve::perform2red(p, *(it.access_details()), scalar);
            //update the approximation of f
                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);
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

    if (p.norm2 == 0) //essentially means that p was already inside the list.
	{
        //cout << "p has norm 2" << endl;
		number_of_collisions++;
		return;
	}


    //insert p into main_list;
    main_list.insert_before(it_comparison_flip,p);
    ++current_list_size;
    if(update_shortest_vector_found(p))
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
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,pApprox.get_approx_norm2(),2* pApprox.get_length_exponent()-1,n   );
        if(!predict){++it;continue;} //if prediction is bad, don't even bother to reduce.
        LatticePoint<ET> current_list_point = it.get_exact_point();
        ++number_of_exact_scprods;
        if (GaussSieve::check2red_new(current_list_point, p, scalar)) //We can reduce *it.
		{
			//create new list point
			LatticePoint<ET> reduced = GaussSieve::perform2red(p, current_list_point, scalar);
			//if (!predict) cerr << "Misprediction 2" << endl;
			//cout << "v was found" <<  endl;

            if (reduced.norm2 == 0) //Note : this cannot happen unless the list contains a non-trivial multiple of p (because the collision would have triggered before).
            {
                number_of_collisions++;
                ++it;
                continue; //was a break. Changed to ++it;continue; -- Gotti
            }

			main_queue.push(reduced);
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


template<class ET>
void Sieve<ET,false>::SieveIteration3 (LatticePoint<ET> &p)
{
    if (p.norm2==0) return;
    ApproxLatticePoint<ET,false> pApprox (p);
    
    int const n = get_ambient_dimension();
    typename MainListType::Iterator it_comparison_flip=main_list.cend();
    
    //if abs( <p, x1> / (|p||x1|) ) < px1bound, do not put it in the filtered_list
    float  px1bound = 0.33; // TO ADJUST
    
    //'promissingness'
    float x1x2=.33;

    // length_factor determines the difference (multiplicatively) between the legnth of the i-th and (i+1)-st blocks
    float length_factor = 2.0; //TO ADJUST
    
    //number of blocks
    int NumOfBlocks = 0;
    
    auto it = main_list.cbegin();
    
    typename MainListType::Iterator first_element_of_the_block = main_list.cbegin();
    LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_the_current_block =  first_element_of_the_block->get_approx_norm2();
    LatticeApproximations::ApproxTypeNorm2 max_length_of_the_current_block = floor(length_factor * assumed_norm_of_the_current_block + 1);
    
    
    //BlockDivisionType BlockPointers;
    //ApproxLatticePoint<ET,false> first_block_element =*it;
    //BlockPointers[NumOfBlocks] = &first_block_element;
    //BlockPointers[NumOfBlocks] = *it;
    
    LengthDivisionType BlockDivision;
    BlockDivision[NumOfBlocks] = it->get_approx_norm2();
    
    
    // the same number of appendices as the number of blocks
    std::array<AppendixType, 100> appendices;
    
    FilteredListType filtered_list;
    FilterDivisionType last_elements;
    FilterNumOfElems num_of_elements;
    num_of_elements.fill(0);
    
    while (it!=main_list.cend())
    {
        
        //cout << "1: consider a list element of approx norm = " << it->get_approx_norm2() << endl;
        //check if p is still of max. length
        if (p.norm2 < it.get_true_norm2())
        {
            //cout << "reached the position to insert p" << endl;
            it_comparison_flip = it;
            break;
        }
        
        
         //--------------------------------2-red-------------------------------
         
        //check 2-red
        ++number_of_scprods;
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );
        ET scalar;    
       
        if(predict){

                
        ++number_of_exact_scprods;
        
         // preform 2-reduction
        if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
            {
                p = GaussSieve::perform2red(p, *(it.access_details()), scalar);
                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);
                
                //put p back into the queue and break
                if (p.norm2!=0)
                    main_queue.push(p);
                //cout << "put p of size " << p.norm2 << " into the queue " << endl;
                //cout << "pApprox has norm2 " << pApprox.get_approx_norm2() <<  endl;
                //cout << "2-reduced p, break all the loops" << endl;
                return;
                //inner_loop = false;
                //break;
            }
            else
                ++number_of_mispredictions;
        }
        
        
        //--------------------------------3-red-------------------------------
        
        
        //check if we reached the next block; if yes, re-compute max_length_of_the_current_block
        if (it->get_approx_norm2() > max_length_of_the_current_block)
        {
                assumed_norm_of_the_current_block = it->get_approx_norm2();
                max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
                NumOfBlocks++;
                BlockDivision[NumOfBlocks] = assumed_norm_of_the_current_block;
                
        }
        
        
        ApproxTypeNorm2 true_inner_product_px1 = compute_sc_prod(pApprox.get_approx(), it->get_approx(), n);
        
        // scaling factor for <px1>: ||p|| * ||x_1||
        //float scale = (float)(pow(pApprox.get_approx_norm2(), 0.5)) * (float)(pow (it->get_approx_norm2(), 0.5));
        float scale = sqrt ( (float)pApprox.get_approx_norm2() * (float)(it->get_approx_norm2()) );
    
        //cout  << "true_inner_product_px1 / scale " << (float)true_inner_product_px1 / scale << endl;
        
        //if true_inner_product_px1 is large enough, put it in filtered_list
        if (abs((float)true_inner_product_px1 / scale)>px1bound)
        {
            
            typename FilteredListType::iterator it_filter =  filtered_list.begin();
            typename FilteredListType::iterator it_comparison_flip_filter = filtered_list.end();
            
            int appendixCounter = 0;
            
            
            //loop over all elements x2 from the filtered list to 3-reduc (p, x1, x2)
            while ( it_filter != filtered_list.cend())
            {
    
                // OR take the length from BlockDivision[appendixCounter]
                LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_x1 = (it_filter->getApproxVector()).get_approx_norm2();

            

                //in res_upper, we store the upper bound on the inner-product of px2 for x2 from filtered_list
                float res_upper = 0.0;
                
                //Compute_one_third_bound(assumed_norm_of_x1, assumed_norm_of_the_current_block, true_inner_product_px1, x1x2, res_upper);
                
                //check if we even need to iterate over the current block of filtered_list
                if (it_filter->get_sc_prod() > res_upper)
                {
                
                    //check if the top element of the relevant appendix is larger than res_upper
                    bool insertion = false;
                    if (!appendices[appendixCounter].empty())
                    {
                        if ( (appendices[appendixCounter].top().getApproxVector()).get_approx_norm2() > res_upper )
                            insertion = true;
                    }
                
                
                    ApproxTypeNorm2 true_inner_product_x1x2;
                    
                    
                }
                   
                it_filter = last_elements[appendixCounter];
                appendixCounter++;
                it_filter++;
        
            } //end of the loop over filtered_list
            
            //once insert (the insertion happend only into the last block), check if we need to update last_elements[appendixCounter]
            
            bool px1_sign;
            if (true_inner_product_px1>0)
                px1_sign = true;
            else
                px1_sign = false;
                
            FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(*it, abs(true_inner_product_px1), px1_sign);
            filtered_list.insert(it_comparison_flip_filter,new_filtered_point);
           
            
        
            if (num_of_elements[appendixCounter] == 0)
                last_elements[appendixCounter] = it_filter;
            else if ( (last_elements[appendixCounter]->getApproxVector()).get_approx_norm2() < it->get_approx_norm2() )
                
                
            num_of_elements[appendixCounter]++;
            
            
            
            
        }
        
    }
        
        
    
    
}








































































template<class ET>
void Sieve<ET,false>::SieveIteration3New (LatticePoint<ET> &p)
{
    
    if (p.norm2==0) return;
    ApproxLatticePoint<ET,false> pApprox (p);
    
    //cout << "------------" << endl;
    //cout << "Run iteration on p of approx norm " << pApprox.get_approx_norm2() << endl;
    
    int const n = get_ambient_dimension();
    
    ET scalar;
    
    typename MainListType::Iterator it_comparison_flip=main_list.cend();
    
    //if abs( <p, x1> / (|p||x1|) ) < px1bound, do not put it in the filtered_list
    float  px1bound = 0.29; // TO ADJUST
    
    //'promissingness'
    float x1x2=.33;

    // length_factor determines the difference (mult) between the legnth of the i-th and (i+1)-st blocks
    float length_factor = 2.0; //TO ADJUST
    //float length_factor =10.0; //to debug assume we have only 1 block
    
    
    // store target inner-products
    //float px1=.0;
    //float px2=.0;
    
    //number of blocks
    int NumOfBlocks = 1;
    
    auto it = main_list.cbegin();
    
    typename MainListType::Iterator first_element_of_the_block = main_list.cbegin();
    LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_the_current_block =  first_element_of_the_block->get_approx_norm2();
    LatticeApproximations::ApproxTypeNorm2 max_length_of_the_current_block = floor(length_factor * assumed_norm_of_the_current_block + 1);
    
    
    // main loop over main_list
    // this while-loop runs until the p is longer
    
    
    
    while (it!=main_list.cend())
    {
        
        //cout << "1: consider a list element of approx norm = " << it->get_approx_norm2() << endl;
        //check if p is still of max. length
        if (p.norm2 < it.get_true_norm2())
        {
            //cout << "reached the position to insert p" << endl;
            it_comparison_flip = it;
            break;
        }
        
        
         //--------------------------------2-red-------------------------------
         
        //check 2-red
        ++number_of_scprods;
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );
            
       
        if(predict){

                
        ++number_of_exact_scprods;
        
         // preform 2-reduction
        if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
            {
                p = GaussSieve::perform2red(p, *(it.access_details()), scalar);
                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);
                
                //put p back into the queue and break
                if (p.norm2!=0)
                    main_queue.push(p);
                //cout << "put p of size " << p.norm2 << " into the queue " << endl;
                //cout << "pApprox has norm2 " << pApprox.get_approx_norm2() <<  endl;
                //cout << "2-reduced p, break all the loops" << endl;
                return;
                //inner_loop = false;
                //break;
            }
            else
                ++number_of_mispredictions;
        }
        
        
        //--------------------------------3-red-------------------------------
        
        
        
        //check if we reached the next block; if yes, re-compute max_length_of_the_current_block
        if (it->get_approx_norm2() > max_length_of_the_current_block)
        {
                assumed_norm_of_the_current_block = it->get_approx_norm2();
                max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
                NumOfBlocks++;
            
                //cout << "enter the next block" << endl;
                //cout << "assumed_norm_of_current_block = " << assumed_norm_of_the_current_block << endl;
                //cout <<"max_length_of_current_block = " << max_length_of_the_current_block << endl;
                //cout << "filtered_list is:" << endl;
            
                //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;
            
                //cout << endl;
            
            
        }
        
        ApproxTypeNorm2 true_inner_product_px1 = compute_sc_prod(pApprox.get_approx(), it->get_approx(), n);
        
        // scaling factor for <px1>: ||p|| * ||x_1||
        //float scale = (float)(pow(pApprox.get_approx_norm2(), 0.5)) * (float)(pow (it->get_approx_norm2(), 0.5));
        float scale = sqrt ( (float)pApprox.get_approx_norm2() * (float)(it->get_approx_norm2()) );
    
        //cout  << "true_inner_product_px1 / scale " << (float)true_inner_product_px1 / scale << endl;
        
        //if true_inner_product_px1 is large enough, put it in filtered_list
        if (abs((float)true_inner_product_px1 / scale)>px1bound)
        {
            
            //loop over all elements x2 from the filtered list to 3-reduce (p, x1, x2)
            auto it_filter = filtered_list2.cbegin();
            auto next_block = filtered_list2.cend();
            
            //bool filter_loop = true;
            while ( it_filter != filtered_list2.cend())
            {
                LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_x1 = get<0>(it_filter->first);
                
                
                pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2> new_key_pair_bound = make_pair(assumed_norm_of_x1+1, 0);
                
                //typename FilteredListType2:: iterator next_block = filtered_list2.lower_bound(new_key_pair_bound );
                next_block = filtered_list2.lower_bound(new_key_pair_bound );
                
                //returns garbadge if there is only one block in the filtered list
                //LatticeApproximations::ApproxTypeNorm2 norm_of_the_next_block = get<0>(next_block->first);
                //cout << "in filter: " << "block-len = " <<  assumed_norm_of_x1 << " next_block = " << norm_of_the_next_block << endl;
                
                
                //in res_upper, we store the upper bound on the inner-product of px2 for x2 from filtered_list
                float res_upper = 0.0;
                
                // as the second argument pass the value assumed_norm_of_the_current_block
                //Compute_px1_bound(assumed_norm_of_x1, assumed_norm_of_the_current_block, true_inner_product_px1, x1x2, res_upper)
                //the result is stored in res_upper
                Compute_px1_bound(assumed_norm_of_x1, it->get_approx_norm2(), true_inner_product_px1, x1x2, res_upper);
                
                //cout <<  " res_upper = " << (LatticeApproximations::ApproxTypeNorm2)res_upper << endl;
                
                typename FilteredListType2:: iterator itup;
                auto it_tmp = it_filter;
                
                // if the bounds are too large, consider the next length-block
                // otherwise find itup s.t. all inner-products after itup are smaller (i.e. worse) than res_upper. Iterate up until itup;
                // it_filter should point to the element with the smallest inner-product within this block

                
                new_key_pair_bound = make_pair(assumed_norm_of_x1, (LatticeApproximations::ApproxTypeNorm2)res_upper);
                
                //find the iterator s.t. all elements < it have <px1> <= res_upper
                itup = filtered_list2.upper_bound(new_key_pair_bound );
                
                
                //if (it_tmp->first == itup->first )
                //{
                    //cout << "itup = " << get<1>(itup->first) << endl;
                    
                    //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                    //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;
                    
                    //cout << "itup is equal to it_filter, no iteration" << endl;
                    
                //}
                //else
                //{
                    //debugging
                    
                    //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                    //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;
                        
                    //cout << "itup = " << get<1>(itup->first) << endl;
                    //cout << "it_filter = " << get<1>(it_filter->first) << endl;
                    ApproxTypeNorm2 true_inner_product_x1x2;
                    
                    //while (it_filter != itup && it_filter!=filtered_list2.cend())
                    while(it_filter !=itup && it_filter !=next_block)
                    {
                        
                        //cout << "it_filter = " << get<1>(it_filter->first) << endl;
                        //retrieve x2 from it_filters and compare x1x2
                        //predict = LatticeApproximations::Compare_Sc_Prod_3red((it_filter->second).getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);
                        predict = LatticeApproximations::Compare_Sc_Prod_3red((it_filter->second).getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);

                       // cout << " true_inner_product_x1x2 = " << true_inner_product_x1x2 << endl;

                        if (predict)
                        {
                            //cout << "REDUCTION" << endl;
                            
                            //sgn1 and sgn2 store the signs s.t. || p + sgn1*x1 + sgn*x2 || < || p ||
                            // these values are computed in check3red_signs
                            int sgn1 = 0;
                            int sgn2 = 0;
                            
                            
                            //in fistered_list we store abs(<p, x1>); retrieve the true sign of x1 stored in filtered_list
                            int sgn_px1 = 1;
                            int sgn_px2 = 1 ;
                            int sgn_x1x2 = 1;
                            
                            
                            if (true_inner_product_px1 < 0 ) sgn_px1 = -1;
                            if (!(it_filter->second).get_sign()) sgn_px2 = -1;
                            if (true_inner_product_x1x2 < 0 ) sgn_x1x2 = -1;
                            
                            // check
                            //cout << "true_inner_product_px1 " << true_inner_product_px1 << " sgn_px1 " << sgn_px1 << endl;
                            //cout << "true_inner_product_px2 " << (it_filter->second).get_sign() << " sgn_px2 " << sgn_px2 << endl;
                            //cout << "true_inner_product_x1x2 " << true_inner_product_x1x2 << " sgn_x1x2 " << sgn_x1x2 << endl;
                            
                            
                            //check if reduction
                            if (GaussSieve::check3red_signs(p, *(it.access_details()), ((it_filter->second).getApproxVector()).get_details(), sgn_px1, sgn_px2, sgn_x1x2, sgn1, sgn2))
                            {
                               
                                //perfrom actual reduction
                                p = GaussSieve::perform3red(p, *(it.access_details()), ((it_filter->second).getApproxVector()).get_details(), sgn1, sgn2);
                                
                                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);
                                
                                if(p.norm2!=0)
                                    main_queue.push(p);
                                
                                //if(get_best_length2() == 3135573)
                                //{
                                    //cout << "reduced p: RE-START with p on norm2 = " << p.norm2 << endl;
                                    
                                //}
                                //cout << "reduced p: RE-START with p on norm2 = " << p.norm2 << endl;
                                //assert(false);
                                
                            
                                filtered_list2.clear();
                                return;
                
                            }
                        
                        }
                        
                        ++it_filter;
                    }   //while -loop over one length-block
                    
                    //cout << "reached it_up" << endl;
                    
                //} //else-condition
                
                //go the next length_block
                it_filter = next_block;
                
                //cout << "go to length " << get<0>(it_filter->first) << endl;
                
            } // end of while it_filter - loop
            
            
            //cout << "filter while-loop is finished "<< endl;
            
            //insert *it into filtered_list2
            
            bool px1_sign;
            if (true_inner_product_px1>0)
                px1_sign = true;
            else
                px1_sign = false;
                
            FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(*it, abs(true_inner_product_px1), px1_sign);
            //ApproxLatticePoint<ET> point = *it;
            //FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(&point, abs(true_inner_product_px1), px1_sign);

            pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2> new_key_pair;
            new_key_pair = make_pair(assumed_norm_of_the_current_block, abs(true_inner_product_px1));
            
            //cout << "about to insert into filtered_list" << endl;
            //TODO: provide hint to emplace
            filtered_list2.emplace_hint(next_block, new_key_pair, new_filtered_point);
            
            
            //cout << "input x with px = " << true_inner_product_px1  << " and px1_sign = " << px1_sign << endl;
            
            //if(get_best_length2() == 3135573)
            //{
            //    cout <<"filtered_list.size = " << filtered_list2.size() << endl;
            //}
            //if (filtered_list2.size() > 22)
            //    assert(false);
            
        }
        
        else
        {
            //cout << "1: scaled px1 = " <<  (float)true_inner_product_px1 / scale << " is smaller than " << px1bound << endl;
        }
    
        
            
        ++it;
    }
    
    //if(get_best_length2() == 3135573)
   //{
      //  cout << "INSERT p of norm " << p.norm2 << endl;
    //}
    //cout << "INSERT p of norm " << p.norm2 << endl;
    //cout << &it_comparison_flip << endl;
    
    //INSERT p
    main_list.insert_before(it_comparison_flip,p);
    ++current_list_size;
    //cout << "list_size = " <<current_list_size << endl;
    if(update_shortest_vector_found(p))
    {
        if(verbosity>=2)
        {
            cout << "New shortest vector found. Norm2 = " << get_best_length2() << endl;
        }
    }
    
    
    
    //filtered_list2.clear();
    
    //keep the old filtered_list2
    //start iteration from it_comparison_flip
    
    it =it_comparison_flip; //points to the next after p list-element

    
    //cout << "start the lower part of the list" << endl;
    
    assumed_norm_of_the_current_block = it->get_approx_norm2();
    max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
    
    bool reduced_x1 = false; //if true, do not increase ++it and do not include into filtered_list
    
    //now we are reducing *it
    while (it!=main_list.cend())
        
    {
        reduced_x1 = false;
        //cout << "2: consider a list element of approx norm = " << it->get_approx_norm2() << endl;
        //assert(false);
        
        
        //--------------------------------2-red-------------------------------
        ++number_of_scprods;
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );
            
        // preform 2-reduction
        if(predict){

                
        ++number_of_exact_scprods;
        if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
        {
                LatticePoint<ET> reduced = GaussSieve::perform2red(p, *(it.access_details()), scalar);
                
                //put reduced back into the queue
                if (reduced.norm2!=0)
                    main_queue.push(reduced);
                    
                    
                //cout << "put x1 of size " << reduced.norm2 << " into the queue " << endl;
                
                it = main_list.erase(it); //also makes ++it
                --current_list_size;
                reduced_x1 = true;
        }
        else
                ++number_of_mispredictions;
        }
        
        
        //--------------------------------3-red-------------------------------
        
        //check if we reached the next block; if yes, re-compute max_length_of_the_current_block
        if (it->get_approx_norm2() > max_length_of_the_current_block)
        {
                assumed_norm_of_the_current_block = it->get_approx_norm2();
                max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
                NumOfBlocks++;
                //cout << "enter the next block" << endl;
                //cout << "assumed_norm_of_current_block = " << assumed_norm_of_the_current_block << endl;
                //cout <<"max_length_of_current_block = " << max_length_of_the_current_block << endl;
                //cout << "filtered_list is:" << endl;
            
                //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;
            
                //cout << endl;
            
            
        }
        
        ApproxTypeNorm2 true_inner_product_px1 = compute_sc_prod(pApprox.get_approx(), it->get_approx(), n);
        float scale = (float)(pow(pApprox.get_approx_norm2(), 0.5)) * (float)(pow (it->get_approx_norm2(), 0.5));
        
        
        if (abs((float)true_inner_product_px1 / scale)>px1bound)
        {
                auto it_filter = filtered_list2.cbegin();
                
                auto next_block = filtered_list2.cend();
            
                while ( it_filter != filtered_list2.cend() && !reduced_x1 )
                {
                    LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_x1 = get<0>(it_filter->first);
                    pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2> new_key_pair_bound = make_pair(assumed_norm_of_x1+1, 0);
                    
                    next_block = filtered_list2.lower_bound(new_key_pair_bound );
                    
                    //returns garbadge if there is only one block in the filtered list
                    //LatticeApproximations::ApproxTypeNorm2 norm_of_the_next_block = get<0>(next_block->first);
                    //cout << "in filter: " << "block-len = " <<  assumed_norm_of_x1 << " next_block = " << norm_of_the_next_block << endl;
                    
                    
                    float res_upper = 0.0;
                    
                    // as the second argument pass the value assumed_norm_of_the_current_block
                    //Compute_px1_bound(assumed_norm_of_x1, assumed_norm_of_the_current_block, true_inner_product_px1, x1x2, res_upper)
                    Compute_px1_bound(pApprox.get_approx_norm2(), assumed_norm_of_x1, true_inner_product_px1, x1x2, res_upper);
                    //cout <<  " res_upper = " << (LatticeApproximations::ApproxTypeNorm2)res_upper << endl;
                    
                    typename FilteredListType2:: iterator itup;
                    auto it_tmp = it_filter;
                    
                    new_key_pair_bound = make_pair(assumed_norm_of_x1, (LatticeApproximations::ApproxTypeNorm2)res_upper);
                    itup = filtered_list2.upper_bound(new_key_pair_bound );
                    
                    if (it_tmp->first == itup->first )
                    {
                        //cout << "itup = " << get<1>(itup->first) << endl;
                        
                        //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                        //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;
                        
                        //cout << "itup is equal to it_filter, no iteration" << endl;
                        
                    }
                    else
                    {
                        //if (filtered_list2.size()>72)
                        //{
                        //    for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                        //        cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;
                            
                        //    assert(false);
                        //}
                        //cout << "itup = " << get<1>(itup->first) << endl;
                        //cout << "it_filter = " << get<1>(it_filter->first) << endl;
                        ApproxTypeNorm2 true_inner_product_x1x2;
                        
                        while (it_filter != itup && it_filter!=next_block)
                        {
                            predict = LatticeApproximations::Compare_Sc_Prod_3red((it_filter->second).getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);
                            
                            if (predict)
                            {
                                //cout << "REDUCTION FOR x1" << endl;
                                
                                int sgn1 = 0;
                                int sgn2 = 0;
                                
                                int sgn_px1 = 1;
                                int sgn_px2 = 1 ;
                                int sgn_x1x2 = 1;
                                
                                
                                if (true_inner_product_px1 < 0 ) sgn_px1 = -1;
                                if (!(it_filter->second).get_sign()) sgn_px2 = -1;
                                if (true_inner_product_x1x2 < 0 ) sgn_x1x2 = -1;
                                
                                // check
                                //cout << "true_inner_product_px1 " << true_inner_product_px1 << " sgn_px1 " << sgn_px1 << endl;
                                //cout << "true_inner_product_px2 " << (it_filter->second).get_sign() << " sgn_px2 " << sgn_px2 << endl;
                                //cout << "true_inner_product_x1x2 " << true_inner_product_x1x2 << " sgn_x1x2 " << sgn_x1x2 << endl;
                                
                                
                                //check if reduction
                                if (GaussSieve::check3red_signs(*(it.access_details()), p, ((it_filter->second).getApproxVector()).get_details(), sgn_px1, sgn_x1x2, sgn_px2, sgn1, sgn2))
                                {
                                    
                                    //perfrom actual reduction
                                    LatticePoint<ET> reduced = GaussSieve::perform3red(*(it.access_details()), p, ((it_filter->second).getApproxVector()).get_details(), sgn1, sgn2);
                                    
                                    if(reduced.norm2!=0)
                                        main_queue.push(reduced);
                                    
                                    it = main_list.erase(it); //also makes ++it
                                    --current_list_size;
                                    
                                    
                                    //break both while-loops and set the flag not to insert *it into filtered_list
                                    reduced_x1 = true;
                                    break;
                                    
                                    //cout << "reduced x1:  now  x1 is if norm2 = " << p.norm2 << endl;
                                    //assert(false);
                                    
                                }
                            
                            } //if(predict)
                            
                            ++it_filter;
                        } //while -loop over one length-block
                    }//else
                    
                    it_filter = next_block;
                    //cout << "go to length " << get<0>(it_filter->first) << endl;

                }
            
                //cout << "filter while-loop is finished "<< endl;
                if (!reduced_x1)
                {
                    bool px1_sign;
                    if (true_inner_product_px1>0)
                        px1_sign = true;
                    else
                        px1_sign = false;
            
                    FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(*it, abs(true_inner_product_px1), px1_sign);
                    //ApproxLatticePoint<ET> point = *it;
                    //FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(&point, abs(true_inner_product_px1), px1_sign);
                    pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2> new_key_pair;
                    new_key_pair = make_pair(assumed_norm_of_the_current_block, abs(true_inner_product_px1));
                    
                    //if(get_best_length2() == 3135573)
                    //{
                    //    cout << get<0>(new_key_pair) << " " << get<1>(new_key_pair) << endl;
                        //cout << "about to insert into filtered_list from lower-part" << endl;
                        
                    //}
                    
                    filtered_list2.emplace_hint(next_block, new_key_pair, new_filtered_point);
                    //cout << "filtered_list2.size" << filtered_list2.size() << endl;

                    
                }
        
        } //if
        
        else
        {
            //if(get_best_length2() == 3135573)
            //{
                //cout << "scaled px1 = " <<  (float)true_inner_product_px1 / scale << " is smaller than " << px1bound << endl;
            //}
        }
        
        //otherwise ++it was already performed during main_list.erase();
        if (!reduced_x1)
            ++it;
    }
    //cout << "finish the main while-loop, clean filtered_list" << endl;
    filtered_list2.clear();
    
}


template<class ET>
void Sieve<ET,false>::SieveIteration3New_Pointer (LatticePoint<ET> &p)
{
    if (p.norm2==0) return;
    ApproxLatticePoint<ET,false> pApprox (p);
    
    //cout << "------------" << endl;
    //cout << "Run iteration on p of approx norm " << pApprox.get_approx_norm2() << endl;
    
    int const n = get_ambient_dimension();
    
    ET scalar;
    
    typename MainListType::Iterator it_comparison_flip=main_list.cend();
    
    //if abs( <p, x1> / (|p||x1|) ) < px1bound, do not put it in the filtered_list
    float  px1bound = 0.29; // TO ADJUST
    
    //'promissingness'
    float x1x2=.33;

    // length_factor determines the difference (mult) between the legnth of the i-th and (i+1)-st blocks
    float length_factor = 2.0; //TO ADJUST
    //float length_factor =10.0; //to debug assume we have only 1 block
    
    
    // store target inner-products
    //float px1=.0;
    //float px2=.0;
    
    //number of blocks
    int NumOfBlocks = 1;
    
    auto it = main_list.cbegin();
    
    typename MainListType::Iterator first_element_of_the_block = main_list.cbegin();
    LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_the_current_block =  first_element_of_the_block->get_approx_norm2();
    LatticeApproximations::ApproxTypeNorm2 max_length_of_the_current_block = floor(length_factor * assumed_norm_of_the_current_block + 1);
    
    
    // main loop over main_list
    // this while-loop runs until the p is longer
    
    while (it!=main_list.cend())
    {
        
        cout << "1: consider a list element of approx norm = " << it->get_approx_norm2() << endl;
        //check if p is still of max. length
        if (p.norm2 < it.get_true_norm2())
        {
            //cout << "reached the position to insert p" << endl;
            it_comparison_flip = it;
            break;
        }
        
        
         //--------------------------------2-red-------------------------------
         
        //check 2-red
        ++number_of_scprods;
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );
            
       
        if(predict){

                
        ++number_of_exact_scprods;
        
         // preform 2-reduction
        if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
            {
                p = GaussSieve::perform2red(p, *(it.access_details()), scalar);
                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);
                
                //put p back into the queue and break
                if (p.norm2!=0)
                    main_queue.push(p);
                //cout << "put p of size " << p.norm2 << " into the queue " << endl;
                //cout << "pApprox has norm2 " << pApprox.get_approx_norm2() <<  endl;
                cout << "2-reduced p, break all the loops" << endl;
                return;
                //inner_loop = false;
                //break;
            }
            else
                ++number_of_mispredictions;
        }
        
        
        //--------------------------------3-red-------------------------------
        
        //check if we reached the next block; if yes, re-compute max_length_of_the_current_block
        if (it->get_approx_norm2() > max_length_of_the_current_block)
        {
                assumed_norm_of_the_current_block = it->get_approx_norm2();
                max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
                NumOfBlocks++;
            
                //cout << "enter the next block" << endl;
                //cout << "assumed_norm_of_current_block = " << assumed_norm_of_the_current_block << endl;
                //cout <<"max_length_of_current_block = " << max_length_of_the_current_block << endl;
                //cout << "filtered_list is:" << endl;
            
                //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;
            
                //cout << endl;
            
            
        }
        
        ApproxTypeNorm2 true_inner_product_px1 = compute_sc_prod(pApprox.get_approx(), it->get_approx(), n);
        
        // scaling factor for <px1>: ||p|| * ||x_1||
        float scale = (float)(pow(pApprox.get_approx_norm2(), 0.5)) * (float)(pow (it->get_approx_norm2(), 0.5));
        
        //cout  << "true_inner_product_px1 / scale " << (float)true_inner_product_px1 / scale << endl;
        
        //if true_inner_product_px1 is large enough, put it in filtered_list
        if (abs((float)true_inner_product_px1 / scale)>px1bound)
        {
            
            //loop over all elements x2 from the filtered list to 3-reduce (p, x1, x2)
            auto it_filter = filtered_listp.cbegin();
            auto next_block = filtered_listp.cend();
            
            //bool filter_loop = true;
            while ( it_filter != filtered_listp.cend())
            {
                LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_x1 = get<0>(it_filter->first);
                
                
                pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2> new_key_pair_bound = make_pair(assumed_norm_of_x1+1, 0);
                
                //typename FilteredListType2:: iterator next_block = filtered_list2.lower_bound(new_key_pair_bound );
                next_block = filtered_listp.lower_bound(new_key_pair_bound );
                
                //returns garbadge if there is only one block in the filtered list
                LatticeApproximations::ApproxTypeNorm2 norm_of_the_next_block = get<0>(next_block->first);
                cout << "in filter: " << "block-len = " <<  assumed_norm_of_x1 << " next_block = " << norm_of_the_next_block << endl;
                
                
                //in res_upper, we store the upper bound on the inner-product of px2 for x2 from filtered_list
                float res_upper = 0.0;
                
                // as the second argument pass the value assumed_norm_of_the_current_block
                //Compute_px1_bound(assumed_norm_of_x1, assumed_norm_of_the_current_block, true_inner_product_px1, x1x2, res_upper)
                //the result is stored in res_upper
                Compute_px1_bound(assumed_norm_of_x1, it->get_approx_norm2(), true_inner_product_px1, x1x2, res_upper);
                
                cout <<  " res_upper = " << (LatticeApproximations::ApproxTypeNorm2)res_upper << endl;
                
                typename FilteredListTypeP:: iterator itup;
                auto it_tmp = it_filter;
                
                // if the bounds are too large, consider the next length-block
                // otherwise find itup s.t. all inner-products after itup are smaller (i.e. worse) than res_upper. Iterate up until itup;
                // it_filter should point to the element with the smallest inner-product within this block

                
                new_key_pair_bound = make_pair(assumed_norm_of_x1, (LatticeApproximations::ApproxTypeNorm2)res_upper);
                
                //find the iterator s.t. all elements < it have <px1> <= res_upper
                itup = filtered_listp.upper_bound(new_key_pair_bound );
                
                
                if (it_tmp->first == itup->first )
                {
                    //cout << "itup = " << get<1>(itup->first) << endl;
                    
                    //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                    //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;
                    
                    cout << "itup is equal to it_filter, no iteration" << endl;
                    
                }
                else
                {
                    //debugging
                    
                    //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                    //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;
                        
                    //cout << "itup = " << get<1>(itup->first) << endl;
                    //cout << "it_filter = " << get<1>(it_filter->first) << endl;
                    ApproxTypeNorm2 true_inner_product_x1x2;
                    
                    //while (it_filter != itup && it_filter!=filtered_list2.cend())
                    while(it_filter !=itup && it_filter !=next_block)
                    {
                        
                        cout << "it_filter = " << get<0>(it_filter->first) << endl;
                        cout << "filtered_listp.size() = " << filtered_listp.size() << endl;
                        cout << "*it_filter contains vector = " << (*it_filter->second).getApproxVector() << endl;
                        //retrieve x2 from it_filters and compare x1x2
                        predict = LatticeApproximations::Compare_Sc_Prod_3red((*it_filter->second).getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);
                
                        cout << " true_inner_product_x1x2 = " << true_inner_product_x1x2 << endl;

                        if (predict)
                        {
                            //cout << "REDUCTION" << endl;
                            
                            //sgn1 and sgn2 store the signs s.t. || p + sgn1*x1 + sgn*x2 || < || p ||
                            // these values are computed in check3red_signs
                            int sgn1 = 0;
                            int sgn2 = 0;
                            
                            
                            //in fistered_list we store abs(<p, x1>); retrieve the true sign of x1 stored in filtered_list
                            int sgn_px1 = 1;
                            int sgn_px2 = 1 ;
                            int sgn_x1x2 = 1;
                            
                            
                            if (true_inner_product_px1 < 0 ) sgn_px1 = -1;
                            if (!(*it_filter->second).get_sign()) sgn_px2 = -1;
                            if (true_inner_product_x1x2 < 0 ) sgn_x1x2 = -1;
                            
                            // check
                            //cout << "true_inner_product_px1 " << true_inner_product_px1 << " sgn_px1 " << sgn_px1 << endl;
                            //cout << "true_inner_product_px2 " << (it_filter->second).get_sign() << " sgn_px2 " << sgn_px2 << endl;
                            //cout << "true_inner_product_x1x2 " << true_inner_product_x1x2 << " sgn_x1x2 " << sgn_x1x2 << endl;
                            
                            
                            //check if reduction
                            if (GaussSieve::check3red_signs(p, *(it.access_details()), ((*it_filter->second).getApproxVector()).get_details(), sgn_px1, sgn_px2, sgn_x1x2, sgn1, sgn2))
                            {
                               
                                //perfrom actual reduction
                                p = GaussSieve::perform3red(p, *(it.access_details()), ((*it_filter->second).getApproxVector()).get_details(), sgn1, sgn2);
                                
                                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);
                                
                                if(p.norm2!=0)
                                    main_queue.push(p);
                                
                                //if(get_best_length2() == 3135573)
                                //{
                                    //cout << "reduced p: RE-START with p on norm2 = " << p.norm2 << endl;
                                    
                                //}
                                //cout << "reduced p: RE-START with p on norm2 = " << p.norm2 << endl;
                                //assert(false);
                                
                            
                                filtered_listp.clear();
                                return;
                
                            }
                        
                        }
                        
                        ++it_filter;
                    }   //while -loop over one length-block
                    
                    cout << "reached it_up" << endl;
                    
                } //else-condition
                
                //go the next length_block
                it_filter = next_block;
                
                //cout << "go to length " << get<0>(it_filter->first) << endl;
                
            } // end of while it_filter - loop
            
            
            //cout << "filter while-loop is finished "<< endl;
            
            //insert *it into filtered_list2
            
            bool px1_sign;
            if (true_inner_product_px1>0)
                px1_sign = true;
            else
                px1_sign = false;
                
            FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(*it, abs(true_inner_product_px1), px1_sign);
            pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2> new_key_pair;
            new_key_pair = make_pair(assumed_norm_of_the_current_block, abs(true_inner_product_px1));
            
            cout << "about to insert into filtered_list" << endl;
            //TODO: provide hint to emplace
            filtered_listp.emplace_hint(next_block, new_key_pair, &new_filtered_point);
            
            
            //cout << "input x with px = " << true_inner_product_px1  << " and px1_sign = " << px1_sign << endl;
            
            //if(get_best_length2() == 3135573)
            //{
            //    cout <<"filtered_list.size = " << filtered_list2.size() << endl;
            //}
            //if (filtered_list2.size() > 22)
            //    assert(false);
            
        }
        
        else
        {
            cout << "1: scaled px1 = " <<  (float)true_inner_product_px1 / scale << " is smaller than " << px1bound << endl;
        }
    
        
            
        ++it;
    }
    
    //if(get_best_length2() == 3135573)
   //{
      //  cout << "INSERT p of norm " << p.norm2 << endl;
    //}
    //cout << "INSERT p of norm " << p.norm2 << endl;
    //cout << &it_comparison_flip << endl;
    
    //INSERT p
    main_list.insert_before(it_comparison_flip,p);
    ++current_list_size;
    //cout << "list_size = " <<current_list_size << endl;
    if(update_shortest_vector_found(p))
    {
        if(verbosity>=2)
        {
            cout << "New shortest vector found. Norm2 = " << get_best_length2() << endl;
        }
    }
   
    
    it =it_comparison_flip; //points to the next after p list-element

    
    //cout << "start the lower part of the list" << endl;
    
    assumed_norm_of_the_current_block = it->get_approx_norm2();
    max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
    
    bool reduced_x1 = false; //if true, do not increase ++it and do not include into filtered_list
    
    //now we are reducing *it
    while (it!=main_list.cend())
        
    {
        reduced_x1 = false;
        //cout << "2: consider a list element of approx norm = " << it->get_approx_norm2() << endl;
        //assert(false);
        
        
        //--------------------------------2-red-------------------------------
        ++number_of_scprods;
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );
            
        // preform 2-reduction
        if(predict){

                
        ++number_of_exact_scprods;
        if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
        {
                LatticePoint<ET> reduced = GaussSieve::perform2red(p, *(it.access_details()), scalar);
                
                //put reduced back into the queue
                if (reduced.norm2!=0)
                    main_queue.push(reduced);
                    
                    
                //cout << "put x1 of size " << reduced.norm2 << " into the queue " << endl;
                
                it = main_list.erase(it); //also makes ++it
                --current_list_size;
                reduced_x1 = true;
        }
        else
                ++number_of_mispredictions;
        }
        
        
        //--------------------------------3-red-------------------------------
        
        //check if we reached the next block; if yes, re-compute max_length_of_the_current_block
        if (it->get_approx_norm2() > max_length_of_the_current_block)
        {
                assumed_norm_of_the_current_block = it->get_approx_norm2();
                max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
                NumOfBlocks++;
                //cout << "enter the next block" << endl;
                //cout << "assumed_norm_of_current_block = " << assumed_norm_of_the_current_block << endl;
                //cout <<"max_length_of_current_block = " << max_length_of_the_current_block << endl;
                //cout << "filtered_list is:" << endl;
            
                //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;
            
                //cout << endl;
            
            
        }
        
        ApproxTypeNorm2 true_inner_product_px1 = compute_sc_prod(pApprox.get_approx(), it->get_approx(), n);
        float scale = (float)(pow(pApprox.get_approx_norm2(), 0.5)) * (float)(pow (it->get_approx_norm2(), 0.5));
        
        
        if (abs((float)true_inner_product_px1 / scale)>px1bound)
        {
                auto it_filter = filtered_listp.cbegin();
                
                auto next_block = filtered_listp.cend();
            
                while ( it_filter != filtered_listp.cend() && !reduced_x1 )
                {
                    LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_x1 = get<0>(it_filter->first);
                    pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2> new_key_pair_bound = make_pair(assumed_norm_of_x1+1, 0);
                    
                    next_block = filtered_listp.lower_bound(new_key_pair_bound );
                    
                    //returns garbadge if there is only one block in the filtered list
                    //LatticeApproximations::ApproxTypeNorm2 norm_of_the_next_block = get<0>(next_block->first);
                    //cout << "in filter: " << "block-len = " <<  assumed_norm_of_x1 << " next_block = " << norm_of_the_next_block << endl;
                    
                    
                    float res_upper = 0.0;
                    
                    // as the second argument pass the value assumed_norm_of_the_current_block
                    //Compute_px1_bound(assumed_norm_of_x1, assumed_norm_of_the_current_block, true_inner_product_px1, x1x2, res_upper)
                    Compute_px1_bound(pApprox.get_approx_norm2(), assumed_norm_of_x1, true_inner_product_px1, x1x2, res_upper);
                    //cout <<  " res_upper = " << (LatticeApproximations::ApproxTypeNorm2)res_upper << endl;
                    
                    typename FilteredListTypeP:: iterator itup;
                    auto it_tmp = it_filter;
                    
                    new_key_pair_bound = make_pair(assumed_norm_of_x1, (LatticeApproximations::ApproxTypeNorm2)res_upper);
                    itup = filtered_listp.upper_bound(new_key_pair_bound );
                    
                    if (it_tmp->first == itup->first )
                    {
                        //cout << "itup = " << get<1>(itup->first) << endl;
                        
                        //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                        //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;
                        
                        //cout << "itup is equal to it_filter, no iteration" << endl;
                        
                    }
                    else
                    {
                        //cout << "itup = " << get<1>(itup->first) << endl;
                        //cout << "it_filter = " << get<1>(it_filter->first) << endl;
                        ApproxTypeNorm2 true_inner_product_x1x2;
                        
                        while (it_filter != itup && it_filter!=next_block)
                        {
                            predict = LatticeApproximations::Compare_Sc_Prod_3red((*it_filter->second).getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);
                            
                            if (predict)
                            {
                                //cout << "REDUCTION FOR x1" << endl;
                                
                                int sgn1 = 0;
                                int sgn2 = 0;
                                
                                int sgn_px1 = 1;
                                int sgn_px2 = 1 ;
                                int sgn_x1x2 = 1;
                                
                                
                                if (true_inner_product_px1 < 0 ) sgn_px1 = -1;
                                if (!(*it_filter->second).get_sign()) sgn_px2 = -1;
                                if (true_inner_product_x1x2 < 0 ) sgn_x1x2 = -1;
                                
                                // check
                                //cout << "true_inner_product_px1 " << true_inner_product_px1 << " sgn_px1 " << sgn_px1 << endl;
                                //cout << "true_inner_product_px2 " << (it_filter->second).get_sign() << " sgn_px2 " << sgn_px2 << endl;
                                //cout << "true_inner_product_x1x2 " << true_inner_product_x1x2 << " sgn_x1x2 " << sgn_x1x2 << endl;
                                
                                
                                //check if reduction
                                if (GaussSieve::check3red_signs(*(it.access_details()), p, ((*it_filter->second).getApproxVector()).get_details(), sgn_px1, sgn_x1x2, sgn_px2, sgn1, sgn2))
                                {
                                    
                                    //perfrom actual reduction
                                    LatticePoint<ET> reduced = GaussSieve::perform3red(*(it.access_details()), p, ((*it_filter->second).getApproxVector()).get_details(), sgn1, sgn2);
                                    
                                    if(reduced.norm2!=0)
                                        main_queue.push(reduced);
                                    
                                    it = main_list.erase(it); //also makes ++it
                                    --current_list_size;
                                    
                                    
                                    //break both while-loops and set the flag not to insert *it into filtered_list
                                    reduced_x1 = true;
                                    break;
                                    
                                    //cout << "reduced x1:  now  x1 is if norm2 = " << p.norm2 << endl;
                                    //assert(false);
                                    
                                }
                            
                            } //if(predict)
                            
                            ++it_filter;
                        } //while -loop over one length-block
                    }//else
                    
                    it_filter = next_block;
                    //cout << "go to length " << get<0>(it_filter->first) << endl;

                }
            
                //cout << "filter while-loop is finished "<< endl;
                if (!reduced_x1)
                {
                    bool px1_sign;
                    if (true_inner_product_px1>0)
                        px1_sign = true;
                    else
                        px1_sign = false;
            
                    FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(*it, abs(true_inner_product_px1), px1_sign);
                    pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2> new_key_pair;
                    new_key_pair = make_pair(assumed_norm_of_the_current_block, abs(true_inner_product_px1));
                    
                    //if(get_best_length2() == 3135573)
                    //{
                    //    cout << get<0>(new_key_pair) << " " << get<1>(new_key_pair) << endl;
                        //cout << "about to insert into filtered_list from lower-part" << endl;
                        
                    //}
                    
                    filtered_listp.emplace_hint(next_block, new_key_pair, &new_filtered_point);
                    //cout << "filtered_list2.size" << filtered_list2.size() << endl;

                    
                }
        
        } //if
        
        else
        {
            //if(get_best_length2() == 3135573)
            //{
                //cout << "scaled px1 = " <<  (float)true_inner_product_px1 / scale << " is smaller than " << px1bound << endl;
            //}
        }
        
        //otherwise ++it was already performed during main_list.erase();
        if (!reduced_x1)
            ++it;
    }
    //cout << "finish the main while-loop, clean filtered_list" << endl;
    filtered_listp.clear();
}
    




//currently unused diagnostic code.
/*
template<class ET>
void PredictionDiagnosis (Sieve<ET,false> * gs, ApproxLatticePoint<ET,false> const & v1, LatticePoint<ET> const &d1, ApproxLatticePoint<ET,false> const &v2, LatticePoint<ET> const &d2, int dim);
template<class ET>
void PredictionDiagnosis (Sieve<ET,false> * gs, ApproxLatticePoint<ET,false> const & v1, LatticePoint<ET> const &d1, ApproxLatticePoint<ET,false> const &v2, LatticePoint<ET> const &d2, int dim)
{
	static int count =0;
	LatticePoint<ET> c1 = d1;
	LatticePoint<ET> c2 = d2;
	bool actual_red = GaussSieve::check2red(c1,c2);
	//cout << (actual_red?"Red:yes" : "Red:no");
	bool predict_red = LatticeApproximations::Compare_Sc_Prod(v1,v2,v2.get_approx_norm2(),2*v2.get_length_exponent() -2,dim);
	//cout << (predict_red?"Predict:yes" : "Predict:no");

	ET sc_prod, abs_scprod, scalar;
	sc_product(sc_prod, d1, d2);
    	abs_scprod.mul_ui(sc_prod,2);
    	abs_scprod.abs(abs_scprod);
	ET n1true = d1.get_norm2();

	int32_t approxSP = abs(LatticeApproximations::compute_sc_prod(v1.get_approx(),v2.get_approx(),dim));
	int approxExp1 = v1.get_length_exponent();
	int approxExp2 = v2.get_length_exponent();
	int approxExpScP = approxExp1 + approxExp2;
	int32_t n1approx = v1.get_approx_norm2();
	int n1exp = 2*v1.get_length_exponent();

	ET approxSP_real;
	ET n1_real;
	LatticePoint<ET> approxv1_real(dim);
	LatticePoint<ET> approxv2_real(dim);
	approxSP_real = static_cast<long>(approxSP);
	n1_real = static_cast<long>(n1approx);
	for(int i=0;i<dim;++i) approxv1_real[i] = static_cast<long>( (v1.get_approx()) [i]);
	for(int i=0;i<dim;++i) approxv2_real[i] = static_cast<long>( (v2.get_approx()) [i]);
	//stupid:
	for(int i=0;i < approxExp1 ;    ++i) approxv1_real = approxv1_real + approxv1_real;
	for(int i=0;i < approxExp2 ;    ++i) approxv2_real = approxv2_real + approxv2_real;
	for(int i=0;i < approxExpScP ;    ++i) approxSP_real.mul_si(approxSP_real,2);
	for(int i=0;i < n1exp;    ++i) n1_real.mul_si(n1_real,2);

	if(actual_red == true && predict_red ==false)
{
	//misprediction.
	cout << "Misprediction" << endl;
	cout << "v1 =" << v1 << endl;
	cout << "meaning of approx1 = " << approxv1_real << endl;
	cout << "v2 =" << v2 << endl;
	cout << "meaning of approx2 = " << approxv2_real << endl;
	cout << "true absscalar product= " << abs_scprod << endl;
	cout << "approx abssc product = " << approxSP << endl;
	cout << "meaning " << approxSP_real << endl;
	cout << "sqNorm1 = " << n1true << endl;
	cout << "Approx Norm1 = " << n1approx << endl;
	cout << "meaning " << n1_real << endl;

}
else if(count % 100 == 80)
	{
	cout <<"Prediction: ";
	cout << (actual_red?"Red:yes" : "Red: no") << " , ";
	cout << (predict_red?"Predict:yes" : "Predict: no") << endl;
	cout << "v1 =" << v1 << endl;
	cout << "meaning of approx1 = " << approxv1_real << endl;
	cout << "v2 =" << v2 << endl;
	cout << "meaning of approx2 = " << approxv2_real << endl;
	cout << "true absscalar product= " << abs_scprod << endl;
	cout << "approx abssc product = " << approxSP << endl;
	cout << "meaning " << approxSP_real << endl;
	cout << "sqNorm1 = " << n1true << endl;
	cout << "Approx Norm1 = " << n1approx << endl;
	cout << "meaning " << n1_real << endl;
	}
	++count;

	//cout << endl;
}
*/


#endif
