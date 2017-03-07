#ifndef SIEVE_GAUSS_SINGLE_THREADED_CPP
#define SIEVE_GAUSS_SINGLE_THREADED_CPP

#if !defined(SIEVE_GAUSS_SINGLE_THREADED) || GAUSS_SIEVE_IS_MULTI_THREADED == true
#error SieveST.cpp included with wrong macro settings
#endif


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

    int i=0;
    //int MaxIteration = 8000;

    LatticePoint<ET> p;
    //NumVect<ET> sample;

    //check_if_done(); //sets up default conditions if not already set. We ignore the return value. -- Moved to run(), Gotti
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
        //if (i % 500 == 0) {
        //    print_status();
            //cout << "# of collisions: " << number_of_collisions << endl;
            //cout << "norm2 of the so far shortest vector: " << get_best_length2() << endl;

        //}
    }

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

    ET scalar; //reduction multiple output by check3red_new

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
	
            if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
            {
		p = GaussSieve::perform2red(p, *(it.access_details()), scalar);

                //update the approximation of f
                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);

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
        if (GaussSieve::check2red_new(current_list_point, p, scalar)) //We can reduce *it.
		{
			//create new list point
			LatticePoint<ET> reduced = GaussSieve::perform2red(p, current_list_point, scalar);
			//if (!predict) cerr << "Misprediction 2" << endl;
			//cout << "v was found" <<  endl;

            if (reduced.norm2 == 0)
            {
                number_of_collisions++;
                break;
            }

			main_queue.push(reduced);
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

template<class ET>
void Sieve<ET,false>::SieveIteration3 (LatticePoint<ET> &p)
{
        
}

//template<class ET>
//void Sieve<ET,false>::run_3_sieve()
//{
//    
//    int i=0; //# of iterations
//    LatticePoint<ET> p;
//    while (! check_if_done()) {
//        
//        p=main_queue.true_pop();
//        
//        
//    }
//    
//}

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
