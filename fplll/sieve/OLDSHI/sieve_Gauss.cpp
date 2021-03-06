/*
  test small svp challenges
  --------------------------------------
  dim  seed  norm current_status  time
  --------------------------------------
  40    0    1702.46439022965    ok
  50    0    1893.16982862077    ok
  60    0    1943.40088504662    ok
 */


#include "sieve_Gauss.h"
#include "sieve_Gauss_2sieve.cpp"
#include "sieve_Gauss_3sieve.cpp"
#include "sieve_Gauss_4sieve.cpp"

#define REDUCE_TIMING

/**
 * constructor
 */
template<class ZT, class F>
Gauss_sieve<ZT, F>::Gauss_sieve (ZZ_mat<ZT> &B, int alg_arg)
{
  
  /* stats */
  b = B;
  nr = b.getRows();
  nc = b.getCols();
  max_list_size = 0;
  iterations = 0;
  collisions = 0;
  samples = 0;
  verbose = 0;
  goal_sqr_norm = 0;
  mem_lower = pow(2.0, 0.18*nc);
  alg = alg_arg;

  /* sanity check */
  if (alg == 2) {
    cout << "# [info] running 2-sieve" <<endl;
    mult = 0.1;
    add = 200.0;
    iterations_step = 200;
  }
  else if (alg == 3) {
    cout << "# [info] running 3-sieve" <<endl;
    mult = 0.1;
    add = 100.0;
    iterations_step = 50;
  }
  else if (alg == 4) {
    cout << "# [info] running 4-sieve" <<endl;
    mult = 0.1;
    add = 50.0;
    iterations_step = 5;
  }
  else {
    cout << " Error, only support 2-, 3- and 4-sieve" <<endl;
    exit(1);
  }
  
  /* clean up list */
  free_list_queue();

  /* initialize sampler */
  Sampler = new KleinSampler<ZT, F> (b);

  /* initialize list */
  init_list ();

  /* further initialization by randomization */
  //init_list_rand();

  /* done initialization */
  max_list_size = List.size();

  /* output stats */
  cout << "# [info] done initialization, size(List)="
       << List.size() << endl;
  cout << "# [info] done initialization, size(Queue)="
       << Queue.size() << endl;
  cout << "# [info] done initialization, mem_est="
       << mem_lower <<endl;

}


/**
 * deconstructor
 */
template<class ZT, class F>
Gauss_sieve<ZT, F>::~Gauss_sieve ()
{
  free_list_queue();
  free_sampler();
}


/**
 * put matrix vectors to list
 */
template<class ZT, class F>
void Gauss_sieve<ZT, F>::add_mat_list (ZZ_mat<ZT> &B)
{
  Z_NR<ZT> t, current_norm;
  dotProduct(best_sqr_norm, B[0], B[0]);
  ListPoint<ZT>* p;

  for (int i = 0; i < nr; ++i) {
    p = new_listpoint<ZT>(nc);
    MatrixRowToListPoint(B[i], p);

    //cout << "# [info] init: additing point ";
    //cout << p->v << endl;

    if (alg == 3)
      current_norm = update_p_3reduce(p);
    else if (alg == 2)
      current_norm = update_p_2reduce(p);
    else if (alg == 4)
      current_norm = update_p_4reduce(p);
    else {
      cout << " Error, only support 2-, 3- and 4-sieve" <<endl;
      exit(1);
    }
    
    if ((current_norm < best_sqr_norm) && (current_norm > 0)) 
      //if ((current_norm < best_sqr_norm) )
      best_sqr_norm = current_norm;

  }
}


/**
 * init function (used in constructor)
 */
template<class ZT, class F>
void Gauss_sieve<ZT, F>::init_list ()
{
  add_mat_list(b);
}


/**
 * init pool of samples (used later)
 */
template<class ZT, class F>
void Gauss_sieve<ZT, F>::init_list_rand ()
{
  /* after transformation, the size could be large */
  ZZ_mat<mpz_t> NewZ (nr, nc);
  ZZ_mat<ZT> New (nr, nc);
  mpz_t tmp;
  Z_NR<mpz_t> tmpZ;
  mpz_init(tmp);
  long flag = 0;
  FP_NR<double> c, t;
  Z_NR<ZT> x;
  c = 0.0;
  t = 32.0;

  /* init */
  for (int i = 0; i < nr; i++ ) {
    for (int j = 0; j < nc; j++ ) {
      (b[i][j]).get_mpz(tmp);
      NewZ[i][j] = tmp;
    }
  }

  /* randomization */
  for (int i = 0; i < nr; i++ ) {
    for (int k = 0; k < nr; k++ ) {
      if (i != k) {
        x = sampleZ_basic_alt<ZT, FP_NR<double> >(c, t);
        x.get_mpz(tmp);
        tmpZ = tmp;
        (NewZ[i]).addmul(NewZ[k], tmpZ, (NewZ[k]).size());
      }
    }
  }
  
  /* reduce */
  lllReduction(NewZ, 0.99, 0.51, LM_FAST);

  /* set */
  for (int i = 0; i < nr; i++ ) {
    for (int j = 0; j < nc; j++ ) {
      tmpZ = (NewZ[i][j]).getData();
      tmpZ.get_mpz(tmp);
      New[i][j] = tmp;
    }
  }

  /* add to list */
  add_mat_list(New);
  mpz_clear(tmp);
}


/**
 * free list and queue
 */
template<class ZT, class F>
void Gauss_sieve<ZT, F>::free_list_queue ()
{
  /* clean list */
  typename list<ListPoint<ZT>* >::iterator lp_it;
  for (lp_it = List.begin(); lp_it != List.end(); ++lp_it)
    del_listpoint<ZT>(*lp_it);
  List.clear();

  /* clean queue */
  while (!Queue.empty()) {
    del_listpoint<ZT>(Queue.front());
    Queue.pop();
  }
  
  /* clean priority queue */
  while (!Queue_Samples.empty()) {
    del_listpoint<ZT>(Queue_Samples.top());
    Queue_Samples.pop();
  }

}


/**
 * free sampler
 */
template<class ZT, class F>
void Gauss_sieve<ZT, F>::free_sampler ()
{
  delete Sampler;
}



/**
 * set targeted norm^2 
 */
template<class ZT, class F>
void Gauss_sieve<ZT, F>::set_goal_norm2 (Z_NR<ZT> norm) 
{
  goal_sqr_norm = norm;
}


/**
 * set verbose
 */
template<class ZT, class F>
void Gauss_sieve<ZT, F>::set_verbose (bool ver) 
{
  verbose = ver;
}


/**
 * print current info
 */
template<class ZT, class F>
void Gauss_sieve<ZT, F>::print_curr_info () 
{
  if (iterations % iterations_step == 0) {
    cout << "# [info] [" << iterations << "] cols=" << collisions;
    cout << " (" << mult * max_list_size + add << ")";
    cout << " |L|=" << List.size();
    cout << " |Q|=" << Queue.size();
    cout << " |samples|=" << samples;
    cout << " |sv|^2=" << List.front()->norm;
    cout << endl;
    cout << std::flush;
  }
}


/**
 * print final info
 */
template<class ZT, class F>
void Gauss_sieve<ZT, F>::print_final_info () 
{
  long first_size = 0;
  vector<long>::iterator it2 = iters_ls.begin();
  for (typename vector<Z_NR<ZT> >::iterator it1 = iters_norm.begin(); 
       it1 != iters_norm.end(); ++it1, ++it2) {
    if ((*it1) == best_sqr_norm) {
      first_size = (*it2);
      break;
    }
  }
  cout << "# [****] done!" << endl;
  cout << "# [info] [" << iterations << "] cols=" << collisions;
  cout << " (" << mult * max_list_size + add << ")";
  cout << " |L|=" << List.size();
  cout << " |Q|=" << Queue.size();
  cout << " |samples|=" << samples << endl;
  cout << "# [info] max(|L|)=" << max_list_size;
  cout << " log2(max|L|)/n=" << log2(max_list_size)/nc << endl;
  cout << "# [info] true max|L| = " << first_size << endl;;
  cout << "# [info] true log2(max|L|)/n = " << log2(first_size)/nc << endl;
  cout << "# [info] sv is" << endl;
  cout << List.front()->v << endl;
  final_norm.set_z(best_sqr_norm);
  final_norm.sqrt(final_norm, GMP_RNDN);
  cout << "# [info] |sv| = " << final_norm << " (" << best_sqr_norm
       << ")" << endl;
}

template class Gauss_sieve<mpz_t, FP_NR<double> >;
template class Gauss_sieve<long, FP_NR<double> >;
