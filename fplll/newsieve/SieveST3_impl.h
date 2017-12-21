#ifndef SIEVE_ST_3_IMPL_H
#define SIEVE_ST_3_IMPL_H

/*
 MAIN ROUTINES FOR 3-GAUSS SIEVE
 */

namespace GaussSieve
{
#ifndef USE_ORDERED_LIST
template <class SieveTraits>
void Sieve<SieveTraits, false>::sieve_3_iteration_vec()
{
  using std::abs;
  using std::round;
  using std::max;

  typename SieveTraits::FastAccess_Point p = main_queue.true_pop();

  bool is_p_max;

  int scalar = 0;

  static FilteredListType filtered_list;
  filtered_list.reserve(SieveTraits::filtered_list_size_max);
  
  std::cout << "p: " << p.get_norm2() << std::endl;

start_over:

  filtered_list.clear();

  for (auto it_x1 = main_list.cbegin(); it_x1 != main_list.cend();)
  {
    if (!check_simhash_scalar_product<typename SieveTraits::CoordinateSelectionUsed>(
                                                                                   p, it_x1, SieveTraits::threshold_lvls_3sieve_lb_out,
                                                                                   SieveTraits::threshold_lvls_3sieve_ub_out))
    {
      ++it_x1;
      continue;
    }
    std::cout << "it_x1: " << (*it_x1).get_norm2() <<std::endl;
    // CHECK FOR 2-REDUCTION
    LengthType sc_prod_px1 = 0; // annoyingly warns otherwise
    LengthType cond_x1 = 0;
    if(check2red_max_for_3red(p, it_x1, scalar, sc_prod_px1, cond_x1, is_p_max))
    {
      if (is_p_max)
      {
        p.sub_multiply(*it_x1, scalar);
        if (p.is_zero())
        {
          statistics.increment_number_of_collisions();
          return;
        }
        p.update_bitapprox();
        goto start_over;
      }
      else
      {
        auto v_new = main_list.true_pop_point(it_x1); // also does ++it_x1
        v_new.sub_multiply(p, scalar);
        if (v_new.is_zero())  // this only happens if the list contains a non-trivial multiple of p.
        {
          statistics.increment_number_of_collisions();
        }
        else
        {
          main_queue.push(std::move(v_new));
        }
        continue;  // This increments the iterator in the sense that its point to the next element now
      }
    }
    // could not perform 2-reduction, but possibly 3-reduction
    // compute scaled inner-product: <p, x1> / ( ||p||^2 * ||x1||^2)
    // the result should always be positive
    double const sc_prod_px1_normalized =
    convert_to_double(sc_prod_px1) * convert_to_double(sc_prod_px1) /
    (convert_to_double(p.get_norm2()) * convert_to_double(it_x1->get_norm2()));
    
    // If the scalar product is too small, we cannot perform 3-reduction, so we take the next x1
    if (sc_prod_px1_normalized < SieveTraits::x1x2_target)
    {
      ++it_x1;
      continue;  // for loop over it_x1;
    }
    bool const sign_px1 = (sc_prod_px1 > 0);
    
    //for (auto &filtp_x2 : filtered_list)
    typename std::vector<Filtered_Point>::iterator filtp_x2 = filtered_list.begin();
    //auto filtp_x2 = filtered_list.cbegin();
    while (filtp_x2!=filtered_list.end())
    {
      if (!check_simhash_scalar_product<typename SieveTraits::CoordinateSelectionUsed>(
                                                                                       it_x1, (*filtp_x2).sim_hashes, SieveTraits::threshold_lvls_3sieve_lb_inn,
                                                                                       SieveTraits::threshold_lvls_3sieve_ub_inn)
          || (*filtp_x2).delete_flag)
      {
        ++filtp_x2;
        continue; // next filtp_x2
      }
      std::cout << "filtp_x2: " << (*filtp_x2).ptr_to_exact->get_norm2() <<std::endl;
      LengthType sc_prod_x1x2 = ( (*filtp_x2).sign_flip == sign_px1) ?
                                compute_sc_product(*it_x1, *((*filtp_x2).ptr_to_exact))
                                : -compute_sc_product(*it_x1, *((*filtp_x2).ptr_to_exact));
      statistics.increment_number_of_scprods_level2();
      
      //case 1: p is max and can reduce
      
      // actually, cond_x1 = 2*abs(sc_prod_px1) - it_x1->get_norm2(). TODO: check is works
      if (is_p_max && (*filtp_x2).is_p_max &&
          2 * sc_prod_x1x2 < 2*abs(sc_prod_px1) - it_x1->get_norm2() + (*filtp_x2).twice_abs_sc_prod - (*filtp_x2).ptr_to_exact->get_norm2() )
      {
        statistics.increment_number_of_3reds();
        std::cout << "p is max " << std::endl;
        LengthType const debug_test = p.get_norm2();
        if (sign_px1)
        {
          p -= *it_x1;
        }
        else
        {
          p += *it_x1;
        }
        if ((*filtp_x2).sign_flip)
        {
          p -= *((*filtp_x2).ptr_to_exact);
        }
        else
        {
          p += *((*filtp_x2).ptr_to_exact);
        }
        assert(p.get_norm2() <= debug_test);  // make sure we are making progress.
        if (p.is_zero())
        {
          statistics.increment_number_of_collisions();
          return;
        }
        p.update_bitapprox();
        statistics.increment_number_of_3reds();
        goto start_over;
      }
      
      // case 2: x1 is max
      
      // actually, cond_x1 = 2*abs(sc_prod_px1) - p->get_norm2() (Now p<it_x1). TODO: check is works
      if (!is_p_max && it_x1->get_norm2() > (*filtp_x2).ptr_to_exact->get_norm2() &&
          2 * sc_prod_x1x2 < ( 2*abs(sc_prod_px1) - p.get_norm2() ) + ((*filtp_x2).twice_abs_sc_prod - (*filtp_x2).ptr_to_exact->get_norm2() ))
      {
        statistics.increment_number_of_3reds();
        std::cout << "x1 is max " << std::endl;
        LengthType const debug_test = (*it_x1).get_norm2();
        auto v_new = main_list.true_pop_point(it_x1);  // also performs ++it_x1 !
        // Note: sign_px1 says whether we need to change x1 (i.e.
        //       we need to consider p+/- v_new. We instead flip the global sign
        //       and look at v_new +/- p
        if (sign_px1)
        {
          v_new -= p;
        }
        else
        {
          v_new += p;
        }
        // If sign_px1 == true, we need to invert the sign of x2 because of the global sign flip.
        if ((*filtp_x2).sign_flip != sign_px1)
        {
          v_new -= *((*filtp_x2).ptr_to_exact);
        }
        else
        {
          v_new += *((*filtp_x2).ptr_to_exact);
        }
         assert(v_new.get_norm2() < debug_test);  // make sure we are making progress.
        if (v_new.is_zero())
        {
          statistics.increment_number_of_collisions();
        }
        else
        {
          main_queue.push(std::move(v_new));
        }
        // GOTO NEXT X2, do not put x1 into filtered_list
        goto end_of_x1_loop;
      }
      
      // case 3: x2 from filtered_list is max
      // want to compare 2<x1,x2> + 2|<x1, p>| - p.norm2 + 2|<x2, p>| - x1.norm2 < 0
      // cond_x1 = 2|<x1, p>| - min{p.norm2, x1.norm2}
      
      if (it_x1->get_norm2() < (*filtp_x2).ptr_to_exact->get_norm2() &&
          2 * sc_prod_x1x2 < 2*abs(sc_prod_px1) - p.get_norm2() + (*filtp_x2).twice_abs_sc_prod - it_x1->get_norm2() )
      {
        // TODO: CHECK IF IT_X1 is end();
        auto to_make_it_x1_plus = it_x1;
        if (++to_make_it_x1_plus == main_list.cend())
        {
          //SWAP WILL NOT WORK, IGNORE THIS REDUCTION
           goto end_of_x1_loop;
          // ++it_x1;
          // continue;
        }
        
        std::cout << "filtp_x2->delete_flag = " << filtp_x2->delete_flag << std::endl;
        
        filtp_x2->delete_flag = true;
        
        std::cout << "x2 is max " << std::endl;
        LengthType const debug_test = (*filtp_x2).ptr_to_exact->get_norm2();
        std::cout << "debug_test = " << debug_test << std::endl;
        auto v_new = main_list.true_pop_point((*filtp_x2).it_to_main_list);
        std::cout << "v_new :" << v_new.get_norm2() << std::endl << std::flush;
        if (sign_px1)
        {
          v_new -= p;
        }
        else
        {
          v_new += p;
        }
        std::cout << "if test" << std::endl;
        // If sign_px1 == true, we need to invert the sign of x2 because of the global sign flip.
        if ((*filtp_x2).sign_flip != sign_px1)
        {
          v_new -= *((*filtp_x2).ptr_to_exact);
        }
        else
        {
          v_new += *((*filtp_x2).ptr_to_exact);
        }
        assert(v_new.get_norm2() < debug_test);  // make sure we are making progress.
        if (v_new.is_zero())
        {
          statistics.increment_number_of_collisions();
        }
        else
        {
          main_queue.push(std::move(v_new));
        }
        std::cout << "x2 is max finished " << std::endl;
      }
      std::cout << "increasing filtp " << std::endl;
      ++filtp_x2;
      //continue;
      
    } // while-loop over filtered_list
  
    // add x1 into filtered_list in case it was not reduced
    filtered_list.emplace_back(it_x1, sign_px1, is_p_max, cond_x1, 2*abs(sc_prod_px1) );
    ++it_x1;
    end_of_x1_loop:;
  } // loop over main_list

  if (update_shortest_vector_found(p))
  {
    if (verbosity >= 2)
    {
      std::cout << "New shortest vector found. Norm2 = " << get_best_length2() << std::endl;
    }
  }

  main_list.emplace_back(std::move(p));

}
#endif

/*
  main 3-sieve iteration
  for the input point p, we check if ||p +/- x1 +/- x2|| < max {||p||, ||x1||, ||x2||}
  for all x1, x2 from main_list. In case, max is ||p||, p is modified and the iteration
  is re-started with this new p (via goto). If max is x1 or x2, the max is modified,
  removed from main_list and pushed into main_queue.
  This routines also checks for collisions (0-vector after any of the modifications).

*/


#ifdef USE_ORDERED_LIST
template <class SieveTraits>
void Sieve<SieveTraits, false>::sieve_3_iteration()
{
  using std::abs;
  using std::round;
  typename SieveTraits::FastAccess_Point p = main_queue.true_pop();
  // TODO: may be make  filtered_list a member of the Sieve class

  static FilteredListType filtered_list;
  filtered_list.reserve(SieveTraits::filtered_list_size_max);

start_over:
  // it_comparison_flip stores the point where the list elements become larger than p
  auto it_comparison_flip = main_list.cend();
  filtered_list.clear();

  // ||p|| >= ||x1|| >= ||x2||
  for (auto it_x1 = main_list.cbegin(); it_x1 != main_list.cend(); ++it_x1)
  {
    if (p < (*it_x1))  // TODO: use approx_norm2_p (may be)
    {
      it_comparison_flip = it_x1;
      break;  // we proceed to the case where p is no longer the largest point of the triple
    }

    // if the sim_hash - scalar product between p and it_x1 is bad, don't bother with this x1:
    if (!check_simhash_scalar_product<typename SieveTraits::CoordinateSelectionUsed>(
            p, it_x1, SieveTraits::threshold_lvls_3sieve_lb_out,
            SieveTraits::threshold_lvls_3sieve_ub_out))
    {
      continue;
    }

    LengthType const sc_prod_px1 = compute_sc_product(p, *it_x1);
    statistics.increment_number_of_scprods_level1();

    // check for 2-reduction. We already computed the (exact) scalar product anyway.
    LengthType const abs_sc_prod_px1 = abs(sc_prod_px1);

    // cond_x1 = 2*|<p,x_1>| - ||x_1||^2.
    // If >0, we can perform 2-reduction. It is useful to keep this term, since we will re-use it
    // in the exact check for 3-reduction
    LengthType const cond_x1 = 2 * abs_sc_prod_px1 - it_x1->get_norm2();

    if (cond_x1 > 0)  // We can perform 2-reduction, changing p:
    {
      statistics.increment_number_of_2reds();
      double const mult = convert_to_double(sc_prod_px1) / convert_to_double(it_x1->get_norm2());
      int const scalar  = round(mult);
      p.sub_multiply(*it_x1, scalar);
      if (p.is_zero())  // might move this to after start_over.
      {
        statistics.increment_number_of_collisions();
        return;
      }
      // we start the current iteration all over:
      // This should be faster than main_queue.push(std::move(p)); return;
      statistics.increment_number_of_2reds();
      p.update_bitapprox();
      goto start_over;
    }

    // could not perform 2-reduction, but possibly 3-reduction
    // compute scaled inner-product: <p, x1> / ( ||p||^2 * ||x1||^2)
    // the result should always be positive
    double const sc_prod_px1_normalized =
        convert_to_double(sc_prod_px1) * convert_to_double(sc_prod_px1) /
        (convert_to_double(p.get_norm2()) * convert_to_double(it_x1->get_norm2()));

    // If the scalar product is too small, we cannot perform 3-reduction, so we take the next x1
    if (sc_prod_px1_normalized < SieveTraits::x1x2_target)
    {
      continue;  // for loop over it_x1;
    }

    // From here : x1 is a candidate for 3-reduction and will eventually be put into filtered_list.
    //             To avoid checking the triple (p, x1, x1), we only append to filtered_list after
    //             we iterate over candidates for x2.
    bool const sign_px1 = (sc_prod_px1 > 0);
    for (auto &filtp_x2 : filtered_list)  // Note that we know ||p|| >= ||*it_x1|| >= ||x2||
    {

      // (approximate) check if x1 and x2 are close enough to participate in 3-reduction
      // we can remove it and perform the exact check only, however, it seems to speed-up the search
      if (!check_simhash_scalar_product<typename SieveTraits::CoordinateSelectionUsed>(
              it_x1, filtp_x2.sim_hashes, SieveTraits::threshold_lvls_3sieve_lb_inn,
              SieveTraits::threshold_lvls_3sieve_ub_inn))
      {
        continue;
      }
      LengthType sc_prod_x1x2 = (filtp_x2.sign_flip == sign_px1)
                                    ? compute_sc_product(*it_x1, *(filtp_x2.ptr_to_exact))
                                    : -compute_sc_product(*it_x1, *(filtp_x2.ptr_to_exact));
      statistics.increment_number_of_scprods_level2();
      //      cond_x1 == 2|<p,x_1>| - ||x_1||^2.
      // filp_x2.cond == 2|<p,x_2>| - ||x_2||^2, which we computed and stored in a previous
      // iteration. The condition is equivalent to ||p+/-x_1 +/- x_2||^2 < ||p||^2,
      // where we have a minus sign @x1 iff sign_px1 == true
      // and                        @x2 iff filtp_x2.sign_flip == true
      if (2 * sc_prod_x1x2 < cond_x1 + filtp_x2.cond)  // perform 3-reduction:
      {
        statistics.increment_number_of_3reds();
        // LengthType const debug_test = p.get_norm2();
        if (sign_px1)
        {
          p -= *it_x1;
        }
        else
        {
          p += *it_x1;
        }
        if (filtp_x2.sign_flip)
        {
          p -= *(filtp_x2.ptr_to_exact);
        }
        else
        {
          p += *(filtp_x2.ptr_to_exact);
        }
        // assert(p.get_norm2() <= debug_test);  // make sure we are making progress.
        if (p.is_zero())
        {
          statistics.increment_number_of_collisions();
          return;
        }
        p.update_bitapprox();
        statistics.increment_number_of_3reds();
        goto start_over;
      }
    }
    filtered_list.emplace_back(it_x1, sign_px1, cond_x1);
  }  // end of first part of for-loop where p is the largest of the triple

  /**********

  **********/

  // p no longer changes now. it_comparison_flip is iterator to first (shortest) element in the list
  // that is longer than p. If no such element exists, it_comparison_flip refers to after-the-end.
  for (auto it_x1 = it_comparison_flip; it_x1 != main_list.cend();)  // ++it inside loop body.
  {
    // if <p,x1> is bad, don't bother with x1
    if (!check_simhash_scalar_product<typename SieveTraits::CoordinateSelectionUsed>(
            p, it_x1, SieveTraits::threshold_lvls_3sieve_lb_out,
            SieveTraits::threshold_lvls_3sieve_ub_out))
    {
      ++it_x1;
      continue;
    }
    // x1 is somewhat promising according to the above check
    LengthType const sc_prod_px1 = compute_sc_product(p, *it_x1);
    statistics.increment_number_of_scprods_level1();
    bool const sign_px1 = (sc_prod_px1 > 0);

    // check for 2-reduction. We already computed the (exact) scalar product anyway.
    LengthType twice_abs_sc_prod_px1 = 2 * abs(sc_prod_px1);
    // cond_x1_store = 2*|<p,x_1>| - ||x_1||^2.
    // cond_x1_p     = 2*|<p,x_1>| - ||p||^2.
    LengthType const cond_x1_p = twice_abs_sc_prod_px1 - p.get_norm2();
    if (cond_x1_p > 0)  // In this case, we can perform 2-reduction, changing x1
    {
      statistics.increment_number_of_2reds();
      double const mult = convert_to_double(sc_prod_px1) / convert_to_double(p.get_norm2());
      int const scalar  = round(mult);
      assert(scalar != 0);

      // true_pop_point erases it, sio the place to insert p is incremented to avoid segfaults
      if (it_x1 == it_comparison_flip)
      {
        ++it_comparison_flip;
      }

      auto v_new = main_list.true_pop_point(it_x1);  // also performs ++it_x1 !
      v_new.sub_multiply(p, scalar);
      statistics.increment_number_of_2reds();
      if (v_new.is_zero())
      {
        statistics.increment_number_of_collisions();
      }
      else
      {
        main_queue.push(std::move(v_new));
      }
      continue;  // with next it_x1
    }

    // no 2-reduction possible, consider 3-reductions
    // compute scaled inner-product: <p, x1> / ( ||p||^2 * ||x1||^2)
    // the result should always be positive
    double const sc_prod_px1_normalized =
        convert_to_double(sc_prod_px1) * convert_to_double(sc_prod_px1) /
        (convert_to_double(p.get_norm2()) * convert_to_double(it_x1->get_norm2()));

    // If the scalar product is too small, 3-reduction is unlikely to happen, so we take the next x1
    if (sc_prod_px1_normalized < SieveTraits::x1x2_target)
    {
      ++it_x1;
      continue;  // for loop over it_x1;
    }

    // From here : x1 is a candidate for 3-reduction and will eventually be put into filtered_list.
    //             To avoid checking the triple (p, x1, x1), we only append to filtered_list after
    //             we iterate over candidates for x2.
    for (auto const &filtp_x2 : filtered_list)
    {

      if (!check_simhash_scalar_product<typename SieveTraits::CoordinateSelectionUsed>(
              it_x1, filtp_x2.sim_hashes, SieveTraits::threshold_lvls_3sieve_lb_inn,
              SieveTraits::threshold_lvls_3sieve_ub_inn))
      {
        continue;
      }

      // Note that we know ||p|| < ||x1|| and ||x1|| >= ||x2||, so x1 is the maximum.
      LengthType sc_prod_x1x2 = (filtp_x2.sign_flip == sign_px1)
                                    ? compute_sc_product(*it_x1, *(filtp_x2.ptr_to_exact))
                                    : -compute_sc_product(*it_x1, *(filtp_x2.ptr_to_exact));
      statistics.increment_number_of_scprods_level2();
      // The correct condition here has cond_x1_p = 2|<p,x_1>| - ||p||^2.
      // This differs from the case above, because now x_1 is larger than p.
      if (2 * sc_prod_x1x2 < cond_x1_p + filtp_x2.cond)  // perform 3-reduction on x1
      {
        statistics.increment_number_of_3reds();
        // LengthType const debug_test = it_x1->get_norm2();
        if (it_x1 == it_comparison_flip)
        {
          ++it_comparison_flip;
        }
        auto v_new = main_list.true_pop_point(it_x1);  // also performs ++it_x1 !
        // Note: sign_px1 says whether we need to change x1 (i.e.
        //       we need to consider p+/- v_new. We instead flip the global sign
        //       and look at v_new +/- p
        if (sign_px1)
        {
          v_new -= p;
        }
        else
        {
          v_new += p;
        }
        // If sign_px1 == true, we need to invert the sign of x2 because of the global sign flip.
        if (filtp_x2.sign_flip != sign_px1)
        {
          v_new -= *(filtp_x2.ptr_to_exact);
        }
        else
        {
          v_new += *(filtp_x2.ptr_to_exact);
        }
        // assert(v_new.get_norm2() < debug_test);  // make sure we are making progress.
        if (v_new.is_zero())
        {
          statistics.increment_number_of_collisions();
        }
        else
        {
          main_queue.push(std::move(v_new));
        }
        // we break the inner loop over filtp_x2 and continue the it_x1-loop here
        goto end_of_x1_loop;
      }
      // No 3-reduction for this x2
    }
    // No 3-reduction for for this x1 for any x2
    twice_abs_sc_prod_px1 -= it_x1->get_norm2();
    filtered_list.emplace_back(it_x1, sign_px1, std::move(twice_abs_sc_prod_px1));
    ++it_x1;
  end_of_x1_loop:;
  }  // end of second part of for loop

  // put p into the main_list
  assert(!(p.is_zero()));
  if (update_shortest_vector_found(p))
  {
    if (verbosity >= 2)
    {
      std::cout << "New shortest vector found. Norm2 = " << get_best_length2() << std::endl;
    }
  }

  // could be done between the two for-loops, but it would require making a copy of p
  main_list.insert_before(it_comparison_flip, std::move(p));
}

#endif

}  // namespace GaussSieve


#endif
