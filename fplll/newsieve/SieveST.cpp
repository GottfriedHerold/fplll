#ifndef SIEVE_GAUSS_SINGLE_THREADED_CPP
#define SIEVE_GAUSS_SINGLE_THREADED_CPP

#if !defined(SIEVE_GAUSS_SINGLE_THREADED) || GAUSS_SIEVE_IS_MULTI_THREADED == true
#error SieveST.cpp included with wrong macro settings
#endif

namespace GaussSieve{

//We may always assumed that SieveJoint.cpp is prepended before this file


//class next_block;
//using namespace LatticeApproximations; // to be able to use ApproxTypeNorm2 to store inner-produces scaled by length

//class main_list;

template<class SieveTraits>
bool Sieve<SieveTraits,false>::update_shortest_vector_found(FastAccess_Point const & newvector)
{
    if(newvector < (*shortest_vector_found) )
    {
        delete shortest_vector_found;
        shortest_vector_found = new FastAccess_Point (newvector.make_copy());
        return true;
    }
    return false;
}

template<class SieveTraits>
typename Sieve<SieveTraits,false>::FastAccess_Point const & Sieve<SieveTraits,false>::get_shortest_vector_found()
{
    return *shortest_vector_found;
}

template<class SieveTraits>
typename Sieve<SieveTraits,false>::EntryType Sieve<SieveTraits,false>::get_best_length2()
{
    return shortest_vector_found->get_norm2();
}

template<class SieveTraits> void Sieve<SieveTraits,false>::run()
{
    if (verbosity >=2) std::cout << "the shortest vector in the input basis has norm2 = " << get_shortest_vector_found().get_norm2() << std::endl;
    //int MaxIteration = 8000;
    if (term_cond==nullptr)
    {
        if (verbosity >=1) std::cerr << "No termination condition set. Aborting.";
        return;
    }
    if(sieve_k <=1)
    {
        if(verbosity >=1) std::cerr << "sieve_k <=1. Aborting.";
        return;
    }
    sieve_status=SieveStatus::sieve_status_running;
    term_cond->init(this); //initialisation of termination conditions.
    if (verbosity >=2) std::cout << "start " << sieve_k << "-Sieve" << std::endl;

    //dispatch

    switch (sieve_k)
    {
//      case 2: std::cerr << "2-sieve currently deactivated" << std::endl;
        case 2: run_2_sieve(); break;
        case 3: run_3_sieve(); break;
        //default:run_k_sieve(); break;
        default: assert(false);
    }
    sieve_status = SieveStatus::sieve_status_finished;

     //Diagnostic and use of the information moved up to the caller.

    /*
    cout << "sv is " << endl;
    main_list.cbegin().access_details()->printLatticePoint();
    print_status();
    */
}

template<class SieveTraits> void Sieve<SieveTraits,false>::run_2_sieve()
{
//    typename SieveTraits::GaussQueue_ReturnType p;
    int i=0;

    std::cout << "start 2-sieve " << std::endl;

    while (!check_if_done() )
    {
//        p = main_queue.true_pop();

        //convert here???

//        GaussSieve::FastAccess_Point<ET, false, nfixed> p_converted (std::move(p));

        typename SieveTraits::FastAccess_Point p = main_queue.true_pop(); // may need conversion.

//        Sieve<ET,false,nfixed>::sieve_2_iteration(p_converted);
//        std::cout << p << std::endl << std::flush;
        
        //sieve_2_iteration(p);
        hash_sieve_2_iteration(p);
        ++i;
        if (( i % 1000 == 0) && (verbosity >=2))
        {
            std::cout << "[" << i << "]"  << "  |L|=" << current_list_size  << " |Q|=" << main_queue.size() << " #samples = " << number_of_points_sampled << " |sv|= " <<  get_best_length2() << std::endl << std::flush;
        }
    }
}


template<class SieveTraits> void Sieve<SieveTraits,false>::run_3_sieve()
{
    int i=0;
    
    std::cout << "start 3-sieve " << std::endl;
    while (!check_if_done() )
    {
        typename SieveTraits::FastAccess_Point p = main_queue.true_pop();

        sieve_3_iteration(p);
        //sieve_3_iteration_test(p);
        ++i;
        if (( i % 1000 == 0) && (verbosity >=2))
        {
             std::cout << "[" << i << "]"  << "  |L|=" << current_list_size  << " |Q|=" << main_queue.size() << " #samples = " << number_of_points_sampled << " |FL|= " << filtered_list_size << " |sv|= " <<  get_best_length2() <<  std::endl;
        }
    }
}

/*
template<class ET>
void Sieve<ET,false>::run_k_sieve()
{
    LatticePoint<ET> p;
    int i=0;
    while (!check_if_done() )
    {
        p=main_queue.true_pop();
        sieve_k_iteration(p);
        ++i;
        if (( i % 1000 == 0) && (verbosity >=2))
        {
            cout << "[" << i << "]"  << "  |L|=" << current_list_size  << " |Q|=" << main_queue.size() << " #samples = " << number_of_points_sampled << " |sv|= " <<  get_best_length2() << endl;
        }
    }
}
*/


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

}

#include "HyperplaneLSH.h"
#include "SieveST2.cpp"
#include "SieveST3.cpp"
//#include "SieveSTk.cpp"

#endif

