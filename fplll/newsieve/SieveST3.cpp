#define filtered_list_size_max  1000

namespace GaussSieve{


    //The function checks if ||x1 + scalar* x2|| < ||x1||
    // the value <x1,x2> is provided as input
    // scalar is modified
template<class SieveTraits, class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
bool check_2red_with_scprod (typename SieveTraits::FastAccess_Point const &x1,
                             typename SieveTraits::FastAccess_Point const &x2,
                             typename SieveTraits::EntryType const &x1x2, Integer & scalar)
{
  using std::abs;
  using std::round;

  typename SieveTraits::EntryType const abs_2scprod = abs(x1x2 * 2);

    // check if |2 * <x1, x2>| <= |x2|^2. If yes, no reduction
  if (abs_2scprod <= x2.get_norm2())
  {
    return false;
  }

    //compute the multiple mult s.t. res = x1 \pm mult* x2;
  double const mult = convert_to_double( x1x2 ) / convert_to_double( x2.get_norm2() );
  // TODO: Check over- / underflows.
  scalar =  round (mult);
  return true;


}

    // The function checks if
    //    || x1 \pm x2 \pm x3 || < || x1 ||
    // The first arguments is assumed to have the largest norm
    // The last two arguments are modified. They return the correct signs, i.e.
    // x_new = x1 + sgn1*x2 + x3*sng2
    // p_is_max is true if x1 ==p, in which case x3_X stores <x3,x2>
    // otherwise x3_X stores <x3, x1>
template<class SieveTraits, bool p_is_max>
bool check_3red (typename SieveTraits::FastAccess_Point  const &x1,
                 typename SieveTraits::FastAccess_Point const &x2,
                 typename SieveTraits::FlilteredPointType const &x3,
                 typename SieveTraits::EntryType const &x1x2,
                 typename SieveTraits::EntryType const &x3_X,
                 int &sgn2, int &sgn3)


{
        //retrieve the signs of the 3 inner-products;
        //they should satisfy (<p, x1>*<p, x2>*<x1, x2>) < 0
        //check for all 4 combinations that do not lead to a reduction

  typename SieveTraits::EntryType x3x1;
  typename SieveTraits::EntryType x3x2;
  using std::abs;

  if (x1x2 * x3_X * x3.get_sc_prod() >0)
  {
    return false;
  }


  if (!p_is_max) {
    x3x1 = x3_X;
    x3x2 = x3.get_sc_prod();
  }
  else
  {
    x3x1 = x3.get_sc_prod();
    x3x2 = x3_X;
  }


  if (x3x1 <0 && x1x2<0 && x3x2<0 &&
      x2.get_norm2() + (x3.get_point()).get_norm2() <
      2 * ( abs(x3x1 + x1x2 + x3x2) ) )
  {

    sgn2 = 1;
    sgn3 = 1;
        //std::cout << "sgns: case 1" << std::endl;
    return true;
  }

    /* bool f = ((x ^ y) < 0); // true iff x and y have opposite signs*/

    //if (x3x1 <0 && !((x1x2^x3x2) <0) &&
  if (x3x1 <0 && x1x2>0 && x3x2 >0 &&
      x2.get_norm2() + (x3.get_point()).get_norm2() <
      2 * ( -x3x1 + x1x2 + x3x2 ) )
  {

    sgn2 = -1;
    sgn3 = 1;
        //std::cout << "sgns: case 2" << std::endl;
    return true;

  }

  if (x3x1 >0 && x1x2<0 && x3x2>0 &&
      x2.get_norm2() + (x3.get_point()).get_norm2() < 2 * ( x3x1 - x1x2 + x3x2 ) )
  {
    sgn2 = 1;
    sgn3 = -1;
    //std::cout << "sgns: case 3" << std::endl;
    return true;
  }

  if (x3x1 >0 && x1x2>0 && x3x2<0 &&
      (x2.get_norm2() + (x3.get_point()).get_norm2() < 2 * (  x3x1 + x1x2 - x3x2 )) )
  {

    sgn2 = -1;
    sgn3 = -1;
        //std::cout << "sgns: case 4" << std::endl;
    return true;
  }

  return false;
}

    /*
     runs 3-reduction
     TODO: input should be GaussQueue_ReturnType (same for sieve_2_iteration)
     */

template<class SieveTraits> void Sieve<SieveTraits,false>::sieve_3_iteration (typename SieveTraits::FastAccess_Point &p)
{
  using std::abs;
  if (p.is_zero() )
  {
    return; //TODO: Ensure sampler does not output 0 (currently, it happens).
  }



    // ! targets are squared

    //double px1_target = 0.1024;
    //double px1_target  = .1111; // TO ADJUST
  double px1_target = 0.123;

  int scalar; //for 2-reduction

  auto it_comparison_flip=main_list.cend(); //to store the point where the list elements become larger than p.

  FilteredListType filtered_list;

  auto it = main_list.cbegin();

  while (it!=main_list.cend())
  {
    if (p  < (*it) )
    {
      it_comparison_flip = it;
      break;
    }

    EntryType sc_prod_px1 = compute_sc_product(p,*it);
    statistics.increment_number_of_scprods_level1();

        //
        //check for 2-reduction
        //
    if ( check_2red_with_scprod<SieveTraits>(p, *it, sc_prod_px1, scalar) )
    {
      assert(scalar!=0); //should not be 0 in any case
      p-= (*it) * scalar;

      if (p.is_zero() )
      {
        statistics.increment_number_of_collisions();
      }
      else
      {
        main_queue.push(std::move(p));
      }
      return;
    }

    statistics.set_filtered_list_size(0);

    filtered_list.reserve(filtered_list_size_max);

    //
    //compare <p, x1> with px1
    //
    //TODO: Change the computations below to smth. faster/better/cleverer/etc.

    double sc_prod_px1_norm = convert_to_double( sc_prod_px1)*convert_to_double(sc_prod_px1 ) /
                ( convert_to_double ( p.get_norm2()) * convert_to_double( (*it).get_norm2() )) ;
    if (abs(sc_prod_px1_norm) > px1_target)
    {
      //This is a fast iteration accodring to the Internet
      //use 'auto &' to take a reference (copy-constructor is deleted)
      for (auto & filtered_list_point: filtered_list)
      {


        EntryType sc_prod_x1x2 = compute_sc_product(*it, filtered_list_point.get_point());
        statistics.increment_number_of_scprods_level2();

        int sgn2, sgn3;

        //check if || p \pm x1 \pm x2 || < || p ||
        // ! check_3red assumes that the first argument has the largest norm
        if ( check_3red<SieveTraits, true> ( p, *it, filtered_list_point, sc_prod_px1, sc_prod_x1x2, sgn2, sgn3) )
        {

                    //TODO:  RETRIEVE ||p|| from the sc_prods

                    //EntryType pnorm_old = p.get_norm2();

          p += (*it)*sgn2 + (filtered_list_point).get_point() * sgn3;
          //p +=  (*it)*sgn2 + *(filtered_list_point).get_point() * sgn3;

          //FOR DEBUGGING

          /*
          if (p.get_norm2() > pnorm_old)
          {
            std::cout << "bug in computing p " << std::endl;
            assert(false);
          }
          */

          if (p.is_zero() )
          {
            statistics.increment_number_of_collisions();
          }
          else
          {
            main_queue.push(std::move(p));
          }
          return;
        }
      }

      //typename SieveTraits::FlilteredPointType new_filtered_point((*it).make_copy(), sc_prod_px1);
      typename SieveTraits::FlilteredPointType new_filtered_point(&(*it), sc_prod_px1);
      filtered_list.push_back(std::move(new_filtered_point));

#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
      SimHash::BitApproxScalarProduct approx_scprod_res = compute_sc_product_bitapprox_fixed(p, *it);
      statistics.red_stat_sim_hash[static_cast<uint_fast32_t>(approx_scprod_res)]++;
#endif

    } // if  | <p, x1> | >=px1_target
    else
    {
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
      SimHash::BitApproxScalarProduct approx_scprod_res = compute_sc_product_bitapprox_fixed(p, *it);
      statistics.no_red_stat_sim_hash[static_cast<uint_fast32_t>(approx_scprod_res)]++;
#endif
    }

    ++it;
  } //while-loop

  main_list.insert_before(it_comparison_flip,p.make_copy());
  statistics.increment_current_list_size();
  //std::cout << "list_size = " <<current_list_size << std::endl;
  if(update_shortest_vector_found(p))
  {
    if(verbosity>=2)
    {
      std::cout << "New shortest vector found. Norm2 = " << get_best_length2() << std::endl;
    }
  }


  //now p is not the largest
  //it_comparison_flip points to the next after p list-element
  it = it_comparison_flip;
  //for(auto it = it_comparison_flip; it!=main_list.cend(); )
  while (it !=main_list.cend())
  {

        //if true, do not put into the filtered_list
        //if true due to 2-reduction, re-compute the sc_product
    bool x1_reduced = false;

    EntryType sc_prod_px1 = compute_sc_product(p,*it);
    statistics.increment_number_of_scprods_level2();

        //
        //check for 2-reduction
        //

    if ( check_2red_with_scprod<SieveTraits>(*it, p, sc_prod_px1, scalar) )
    {
      assert(scalar!=0); //should not be 0 in any case
      typename SieveTraits::FastAccess_Point v_new = (*it) - (p*scalar);


      if (v_new.is_zero() )
      {
        statistics.increment_number_of_collisions();
      }

      main_queue.push(std::move(v_new));

      it = main_list.erase(it);
      statistics.decrement_current_list_size();

      x1_reduced = true;

            //it was increased by the erase; we may already reach the end, so the code below will segfalut without the if-cond below
      if (it == main_list.cend())
      {
        return;
      }
    }
        //
        // 3-rediction
        //
        // Now x1 is the largest      //x1 can be modified during the check_2_red, hence the computed sc_prod there is no longer valid
    if (x1_reduced)
    {
      sc_prod_px1 = compute_sc_product(p,*it);
      statistics.increment_number_of_scprods_level1();
      x1_reduced = false;
    }

    double sc_prod_px1_norm = convert_to_double( sc_prod_px1 )* convert_to_double( sc_prod_px1 )  /
        ( convert_to_double ( p.get_norm2()) * convert_to_double( (*it).get_norm2() )) ;

    if (std::abs(sc_prod_px1_norm) > px1_target)
    {

      for (auto & filtered_list_point: filtered_list)
      {
        EntryType sc_prod_x1x2 = compute_sc_product(*it, (filtered_list_point).get_point());
        statistics.increment_number_of_scprods_level2();

        int  sgn2, sgn3;

                // ! check_3red assumes that the first argument has the largest norm
        if ( check_3red<SieveTraits, false> ( *it, p, filtered_list_point, sc_prod_px1, sc_prod_x1x2, sgn2, sgn3) )
        {
          typename SieveTraits::FastAccess_Point v_new =(*it) + p*sgn2 + filtered_list_point.get_point() * sgn3;
                    //typename SieveTraits::FastAccess_Point v_new =(*it) + p*sgn2 + *(filtered_list_point).get_point() * sgn3;

          if (v_new.is_zero() )
          {
            statistics.increment_number_of_collisions();
          }
                    //FOR DEBUG
                  /*
                    if (v_new.get_norm2() > (*it).get_norm2())
                    {
                        std::cout << "bug in computing v_new" << std::endl;
                        assert(false);
                    }
                  */
          main_queue.push(std::move(v_new));
          it = main_list.erase(it);

          statistics.decrement_current_list_size();

          x1_reduced = true;

          break; //for-loop over the filtered_list
        }
      } //for-loop


      if (!x1_reduced)
      {
                //typename SieveTraits::FlilteredPointType new_filtered_point((*it).make_copy(), sc_prod_px1);
        typename SieveTraits::FlilteredPointType new_filtered_point(&(*it), sc_prod_px1);
        filtered_list.push_back(std::move(new_filtered_point));
      }

#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
      SimHash::BitApproxScalarProduct approx_scprod_res = compute_sc_product_bitapprox_fixed(p, *it);
      statistics.red_stat_sim_hash[static_cast<uint_fast32_t>(approx_scprod_res)]++;
#endif

    }
    else
    {
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
      SimHash::BitApproxScalarProduct approx_scprod_res = compute_sc_product_bitapprox_fixed(p, *it);
      statistics.no_red_stat_sim_hash[static_cast<uint_fast32_t>(approx_scprod_res)]++;
#endif
    }

    if (!x1_reduced)
    {
      ++it;
    }

        /*
        if (filtered_list.size()>0) {
            std::cout << "filtered.size() = " << filtered_list.size() << std::endl;
        }
         */
  } // 'lower' while-loop

  statistics.set_filtered_list_size(filtered_list.size());
  filtered_list.clear();
}


} //namespace GaussSieve
