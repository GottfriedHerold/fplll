#ifndef SIEVE_GAUSS_SINGLE_THREADED_CPP
#define SIEVE_GAUSS_SINGLE_THREADED_CPP

#if !defined(SIEVE_GAUSS_SINGLE_THREADED) || GAUSS_SIEVE_IS_MULTI_THREADED == true
#error SieveST.cpp included with wrong macro settings
#endif

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

            SieveIteration3(p);


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

    ET scalar;

    //typename MainListType::Iterator it_comparison_flip=main_list.cend();

    float length_factor = 1.15; //to be verified

    typename MainListType::Iterator first_element_of_LHS_block = main_list.cbegin();
    typename MainListType::Iterator first_element_of_RHS_block = main_list.cbegin();
    LatticeApproximations::ApproxTypeNorm2 pApproxNorm =  pApprox.get_approx_norm2();

    //the largest norm
    auto longest_v = main_list.cend();
    LatticeApproximations::ApproxTypeNorm2 largest_norm = (--longest_v)->get_approx_norm2(); //Attention! Using operation-- here.

    bool outer_loop = true;

    // store target inner-products for the currently considered pair of blocks
    float px1=.0;
    float px2=.0;
    float x1x2=.0;

    int NumOfBlocks = 0;

    //assume we have to lists LHS and RHS. Outer loop is over the blocks of the LHS list. For each block (which defines the length of x1), a corresponding filtering is applied to the RHS list. This is a 'standard' loop. The overhead is the number of blocks.
    
    // loop over the blocks in LHS
    while (outer_loop) {

	NumOfBlocks++;
        cout << "------------" << endl;
        cout << " NumOfBlocks= " <<  NumOfBlocks << endl;

        LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_current_block_LHS =  first_element_of_LHS_block->get_approx_norm2();
        LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_current_block_RHS =  assumed_norm_of_current_block_LHS;

        //max length of elements in the current block

        LatticeApproximations::ApproxTypeNorm2 max_length_of_current_block = floor(length_factor * assumed_norm_of_current_block_LHS + 1);

        //check if max_length_of_current_block is less then the longest vector in the list. If yes, reset

        if(max_length_of_current_block > largest_norm) //this also means we are in the last block on the LHS. TODO: use it for the while-condition
        {
	    cout << "this will be the last LHS block" << endl;
            max_length_of_current_block = largest_norm;
            outer_loop = false;
        }
        
        cout << "assumed_norm_of_current_block_LHS = " << assumed_norm_of_current_block_LHS << endl;
        cout <<"max_length_of_current_block = " << max_length_of_current_block << endl;

        // compute targets <p, x1>, <p,x2> for the block (i,i); p has the largest norm

        LatticeApproximations::Determine_Sc_Prod(pApproxNorm, assumed_norm_of_current_block_LHS, assumed_norm_of_current_block_RHS, x1x2, px1,px2);
        

        //iterate over the RHS list starting from the first_element_of_a_block
        
        auto it = first_element_of_LHS_block;
        
        bool inner_loop = true;
        while (inner_loop || it!=main_list.cend())
        {
            if (p.norm2 < it.get_true_norm2())
            {
		cout << "reached the position to insert p" << endl;
                outer_loop = false;
                break;
            }

            // check if x is still in the assumed block)
            if (it->get_approx_norm2() > max_length_of_current_block)
            {
                assumed_norm_of_current_block_RHS = it->get_approx_norm2();
                max_length_of_current_block =floor(length_factor * assumed_norm_of_current_block_RHS + 1);
                
                cout << "enter the next block on the RHS " << endl;
		cout << "assumed_norm_of_current_block_RHS = " << assumed_norm_of_current_block_RHS << endl;
        	cout <<"max_length_of_current_block = " << max_length_of_current_block << endl;

                //recompute target inner-products

                LatticeApproximations::Determine_Sc_Prod(pApproxNorm, assumed_norm_of_current_block_LHS, assumed_norm_of_current_block_RHS, x1x2, px1,px2);

		cout <<"pApproxNorm " << pApproxNorm << endl;
		cout <<" assumed_norm_of_current_block_LHS " << assumed_norm_of_current_block_LHS << endl;
		cout << "assumed_norm_of_current_block_RHS " << assumed_norm_of_current_block_RHS << endl;
		cout << "x1x2 = "<< x1x2 << " px1 = " << px1 << " px2 = " << px2 << endl;

                //for the first if-cond, change first_element_of_LHH_block to 'it' for the next run of the outer-loop
                if (first_element_of_LHS_block == first_element_of_RHS_block )
		{
		    cout << "LHS block first element is changed from " <<first_element_of_LHS_block->get_approx_norm2() << " to " << it->get_approx_norm2() << endl;
                    first_element_of_LHS_block = it;
		}

		//assert(false);

            }

            // check if 2-red is possible
            ++number_of_scprods;
            bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );
            
	    // preform 2-reduction
	    if(predict){
           	 ++number_of_exact_scprods;
            	if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
            	{
                	p = GaussSieve::perform2red(p, *(it.access_details()), scalar);

                	//put p back into the queue and break
                	if (p.norm2!=0)
                    		main_queue.push(p);

                	cout << "2-reduced p, break the loops" << endl;
                	outer_loop = false;
                	break;
            	}
            	else
                	++number_of_mispredictions;
	    }

            //detC = LatticeApproximations::detConf(px1, px2, x1x2);
            //if(detC > threshold) //decide if we look inside the current block at all


            // check if <p, it> is close to px1

            float true_inner_product_px1 = .0;
            predict = LatticeApproximations::Compare_Sc_Prod_3red(pApprox, *it, n, px1, true_inner_product_px1);
            
            cout << "px1 = " << px1;
            cout << " true_inner_product_px1 = " << true_inner_product_px1 << endl;

	    if(abs(true_inner_product_px1) >.5)
		{
			cout << "missed 2-red" << endl;
			assert(false);
		}

            if(!predict) 
	    {
			++it;
			continue;
	    }

            //loop over all elements x2 from the filtered list to 3-reduce (p, x1, x2)

            for (auto it_filter = filtered_list.cbegin(); it_filter != filtered_list.cend(); ++it_filter)
            {

                //now check if <x1, x2> are close to the target x1x2 -- do we need to check it?
                float true_inner_product_x1x2 = .0;
                predict = LatticeApproximations::Compare_Sc_Prod_3red(it_filter->getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);

		cout << "x1x2 = " << px1;
            	cout << " true_inner_product_px1 = " << true_inner_product_x1x2 << endl;
		assert(false);
                if(!predict) continue;
                
                //from the true inner-products we need only the signs to perform 3 reduction. Assume all the approximations gave the correct signs
                
                float true_inner_product_px2 = it_filter->get_sc_prod();
                
                cout << "true_inner_product_px2 = " << true_inner_product_px2 << endl;
                assert(false);
                
                //if (GaussSieve::check3_red(p, *(it.access_details()), *(it_filter->getApproxVector()).access_details(), ))
                {
                    //p = GaussSieve::perform3red();
                    
                    if(p.norm2!=0)
                        main_queue.push(p);
                    
                    inner_loop = false;
                    outer_loop = false;
                    break;
                }

		
			
			
            }

            //add 'it' to filtered_list
            FilteredPoint<ET, float> new_filtered_point(*it, true_inner_product_px1);
            filtered_list.emplace_back(new_filtered_point);
            cout <<"filtered_list.size = " << filtered_list.size() << endl;
            
            ++it;

        } // end of the inner while-loop

        // filtered list is new for a new block from the LHS

        filtered_list.clear();

    } // end of the outer while-loop

};

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
