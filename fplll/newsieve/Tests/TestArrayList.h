#ifndef TEST_ARRAY_LIST_H
#define TEST_ARRAY_LIST_H

#include "../ArrayList.h"

#include <random>
#include <algorithm>
#include <ctime>


bool test_array_list()
{
  static constexpr int LEN = 100000;
  GaussSieve::ArrayList<int, 200000> test_list1;
  std::list<int> stl_list;
  std::vector<int> V;
  V.resize(LEN, 0);
  for(unsigned int i = 0; i< LEN; ++i)
  {
    V[i] = i;
  }
  std::random_device rd;
  std::mt19937 rand_gen(rd());
  std::shuffle(V.begin(), V.end(), rand_gen);
  // V is now a random permutation of 0..LEN - 1

  // dumb - sorting in O(LEN^2), intentionally slow!
  std::clock_t stime = clock();
  for (int i = 0; i < LEN; ++i)
  {
    bool at_end = true;
    for(auto it = stl_list.cbegin(); it != stl_list.cend(); ++it)
    {
      if( (*it) > i)
      {
        stl_list.emplace(it, i);
        at_end = false;
        break;
      }
    }
    if(at_end) { stl_list.emplace(stl_list.cend(),i); }
  }
  std::clock_t etime = clock();
  std::cout << "STL LIST dumb insertion sort took ";
  std::cout << (etime - stime) / (double)CLOCKS_PER_SEC;
  std::cout << " sec." << std::endl;
//  std::cout << "STL LIST" << std::endl;
//  for(auto it = stl_list.cbegin(); it!=stl_list.cend(); ++it)
//  {
//    std::cout << *it << " ";
//  }
//  std::cout << std::endl;

  // same with our arraylist
  stime = clock();
  for( int i=0; i < LEN; ++i)
  {
    bool at_end = true;
    for(auto it = test_list1.cbegin(); !it.is_end(); ++it)
    {
      if( (*it) > i)
      {
        test_list1.emplace(it, i);
        at_end = false;
        break;
      }
    }
    if(at_end)
    {
      auto end_it = test_list1.cend();
      test_list1.emplace_before(end_it, i);
    }
  }
  etime = clock();

  std::cout << "Array List dumb insertion sort took";
  std::cout << (etime - stime) / (double)CLOCKS_PER_SEC;
  std::cout << " sec." << std::endl;

  stime = clock();
  long long acc = 0;
  for(auto it = stl_list.cbegin(); it!=stl_list.cend(); ++it)
  {
    acc += *it;
  }
  etime = clock();
  std::cout << acc << std::endl;
  std::cout << "STL traversal";
  std::cout << (etime - stime) / (double)CLOCKS_PER_SEC;
  std::cout << " sec." << std::endl;

  stime = clock();
  acc = 0;
  for(auto it = test_list1.cbegin(); !it.is_end(); ++it)
  {
    acc += *it;
  }
  etime = clock();
  std::cout << acc << std::endl;
  std::cout << "ARR traversal";
  std::cout << (etime - stime) / (double)CLOCKS_PER_SEC;
  std::cout << " sec." << std::endl;

  stime = clock();
  acc = 0;
  for(auto it = V.cbegin(); it !=V.cend(); ++it)
  {
    acc += *it;
  }
  etime = clock();
  std::cout << acc << std::endl;
  std::cout << "Vector traversal";
  std::cout << (etime - stime) / (double)CLOCKS_PER_SEC;
  std::cout << " sec." << std::endl;



//  for(auto it = test_list1.cbegin(); !it.is_end(); ++it)
//  {
//    std::cout << *it << " ";
//  }
//  std::cout << std::endl;

  return true;

}



#endif
