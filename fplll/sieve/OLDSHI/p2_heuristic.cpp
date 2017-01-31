#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <vector>
#include <numeric>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_on_sphere.hpp>
using namespace std;

const double PI  =3.141592653589793238463;

#define DEBUG 1
#define DEBUG_OUTPUT 0


class CombinationsIndexArray {
  vector<int> index_array;
  int last_index;
public:
  CombinationsIndexArray(int number_of_things_to_choose_from, int number_of_things_to_choose_in_one_combination) {
    last_index = number_of_things_to_choose_from - 1;
    for (int i = 0; i < number_of_things_to_choose_in_one_combination; i++) {
      index_array.push_back(i);
    }
  }
  int operator[](int i) {
    return index_array[i];
  }
  int size() {
    return index_array.size();
  }
  bool advance() {

    int i = index_array.size() - 1;
    if (index_array[i] < last_index) {
      index_array[i]++;
      return true;
    } else {
      while (i > 0 && index_array[i-1] == index_array[i]-1) {
        i--;
      }
      if (i == 0) {
        return false;
      } else {
        index_array[i-1]++;
        while (i < index_array.size()) {
          index_array[i] = index_array[i-1]+1;
          i++;
        }
        return true;
      }
    }
  }
};


template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
  if (a.size() != b.size())
    throw std::domain_error("vector addition must have equal length vectors");
  std::vector<T> result;
  result.reserve(a.size());  
  std::transform(a.begin(), a.end(), b.begin(), 
                 std::back_inserter(result), std::plus<T>());
  return result;
}


template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
  if (a.size() != b.size())
    throw std::domain_error("vector subtraction must have equal length vectors");
  std::vector<T> result;
  result.reserve(a.size());
  std::transform(a.begin(), a.end(), b.begin(),
                 std::back_inserter(result), std::minus<T>());
  return result;
}


/**
 *	Cosine between v1 and v2
 */
template<class T>
double cosine(vector<T> v1, vector<T> v2, int verbose)
{

  if (v1.size() != v2.size())
    throw std::domain_error("vector subtraction must have equal length vectors");

	T lv1 = v1[0] * v1[0];
	T lv2 = v2[0] * v2[0];
	T dot12 = v1[0] * v2[0];

	for (unsigned int i=1; i < v1.size(); i++){
		lv1 += v1[i] * v1[i];
		lv2 += v2[i] * v2[i];
		dot12 += v1[i] * v2[i];
	}

	return double(dot12) / ( sqrt(double(lv1)) * sqrt(double(lv2)) );
}


template<class T>
double angle (vector<T> v1, vector<T> v2, int verbose) {
  double cosa = cosine(v1, v2, verbose);
  return acos (cosa) * 180.0 / PI;
}


/**
 *	The mean of a vector
 */
template<class T>
double mean(vector<T> v1) {
	T sum = v1[0];
	for (unsigned int i=1; i < v1.size(); i++)
		sum += v1[i];
	return double(sum) / double(v1.size());
}


/**
 *	The Covariance
 */
template<class T>
double covariance(vector<T> v1, vector<T> v2){
  if (v1.size() != v2.size())
    throw std::domain_error("vector subtraction must have equal length vectors");

	double mean1 = mean(v1), mean2 = mean(v2);
	double sum = (double(v1[0]) - mean1) * (double(v2[0]) - mean2);
	for (unsigned int i=1; i < v1.size(); i++){
		sum += (double(v1[i]) - mean1) * (double(v2[i]) - mean2);
	}
	return double(sum) / double(v1.size()-1);
}


/**
 *	standard deviation the covariance where both vectors are the same.
 */
template<class T>
double std_dev(vector<T> v1){
	return sqrt(covariance(v1, v1));
}


template<class T>
double distance2 (vector<T> v1, vector<T> v2)
{
	if(v1.size() != v2.size()){
    throw std::domain_error("vector subtraction must have equal length vectors");
	}
	T diff, sum;
	diff = v1[0] - v2[0];
	sum = diff * diff;

	for (unsigned int i=1; i < v1.size(); i++){
		diff = v1[i] - v2[i];
		sum += diff * diff;
    // early abort
    //if (sum > 1)
    //  return 2.0;
	}
	return sqrt(double(sum));
}


template<class T>
double norm(vector<T> v)
{
	T sum1,sum2,sum3,sum4;
  sum1 = 0;
  sum2 = 0;
  sum3 = 0;
  sum4 = 0;
  int n = v.size();
  for(int i = 0; i < n/4*4; i += 4) {
    sum1 = sum1 + v[i]*v[i];
    sum2 = sum2 + v[i+1]*v[i+1];
    sum3 = sum3 + v[i+2]*v[i+2];
    sum4 = sum4 + v[i+3]*v[i+3];
  }
  for(int i = n/4*4; i < n; i++)
    sum1 = sum1 + v[i]*v[i];
  sum1 = sum1+sum2+sum3+sum4;
	return sqrt(double(sum1));
}


void
coutPoint (const vector<double>& v)
{
  for (int i=0; i<v.size();i++){
    cout << v[i] << " ";
  }
  cout << '\n';
}


void check_list_2red (vector< vector<double> > &L)
{

  /* n choose 2 */
  vector<int> index;
  for (int i = 1; i <= L.size(); i++) {
    index.push_back(i);
  }
  CombinationsIndexArray combos(L.size(), 2);
  vector<vector<double> > newpoints(combos.size());
  do  {
    for (int i = 0; i < combos.size(); i++) {
      newpoints[i] = L[index[combos[i]-1]];
    }

    /* fix for higher dimension */
    double dist = distance2(newpoints[0], newpoints[1]);
    double ang = angle(newpoints[0], newpoints[1], 0);

    /* check just half */
    if (ang < 60 || ang > 120) {
      cout << "ERROR CHECK THE CODE " << endl;
      double ang = angle(newpoints[0], newpoints[1], 1);
      cout << "distance " << dist << ", ";
      cout << "angle " << ang << endl;
      coutPoint(newpoints[0]);
      coutPoint(newpoints[1]);
      cout << endl;
    }

#if DEBUG_OUTPUT
    myfile1 << ang << endl;
    myfile2 << dist << endl;
#endif
    
  }
  while (combos.advance());

  cout << "check 2-red OK" << endl;  

}


void check_list_3red (vector< vector<double> > &L)
{

#if DEBUG_OUTPUT
  ofstream myfile1;
  ofstream myfile2;
  myfile1.open ("angle.data");
  myfile2.open ("dist.data");
#endif

  vector<double> tmp1;
  vector<double> tmp2;
  vector<double> tmp3;
  vector<double> tmp4;
  vector<double> tmp5;
  vector<double> tmp6;

  
  /* n choose 2 */
  vector<int> index;
  for (int i = 1; i <= L.size(); i++) {
    index.push_back(i);
  }
  CombinationsIndexArray combos(L.size(), 3);
  vector<vector<double> > newpoints(combos.size());
  
  do  {
    for (int i = 0; i < combos.size(); i++) {
      newpoints[i] = L[index[combos[i]-1]];
    }

    /* check pairwisely */
    double ang = angle(newpoints[0], newpoints[1], 0);
    if (ang < 60 || ang > 120) {
      cout << "ERROR CHECK THE CODE " << endl;exit(1);
    }
    ang = angle(newpoints[0], newpoints[2], 0);
    if (ang < 60 || ang > 120) {
      cout << "ERROR CHECK THE CODE " << endl;exit(1);
    }
    ang = angle(newpoints[1], newpoints[2], 0);
    if (ang < 60 || ang > 120) {
      cout << "ERROR CHECK THE CODE " << endl;exit(1);
    }

    /* check more (easier to use length) */
    tmp1 = newpoints[0] + newpoints[1] + newpoints[2];
    double leng = norm(tmp1);
    tmp1 = newpoints[0] + newpoints[1] - newpoints[2];
    leng = min(leng, norm(tmp1));
    tmp1 = newpoints[0] - newpoints[1] + newpoints[2];
    leng = min(leng, norm(tmp1));
    tmp1 = newpoints[0] - newpoints[1] - newpoints[2];
    leng = min(leng, norm(tmp1));
    double leng2 = norm(newpoints[0]);
    leng2 = max(leng2, norm(newpoints[1]));
    leng2 = max(leng2, norm(newpoints[2]));
    
    if (leng < leng2) {
      cout << "ERROR CHECK THE CODE " << endl;
      coutPoint(newpoints[0]);
      coutPoint(newpoints[1]);
      coutPoint(newpoints[2]);
      exit(1);
    }
    
  }
  while (combos.advance());

  //cout << "check 3-red OK" << endl;
  
#if DEBUG_OUTPUT
  myfile1.close();
  myfile2.close();
#endif
  /*
    for(std::vector<int>::size_type j = 0; j != L.size(); j++) {
    newpoint = L[j];
    for(std::vector<int>::size_type k = 0; k != L.size(); k++) {
    newpoint2 = L[k];
    double dist = euclidean(newpoint, newpoint2);
    double ang = angle(newpoint, newpoint2);
    if (dist > 0) {
    cout << "distance " << dist << ", ";
    cout << "angle " << ang << endl;
    }
    }
    }
  */
}


double cover_sphere_uniform (int n, int verbose)
{

  boost::mt19937 gen;
  boost::uniform_on_sphere<double> sampler(n);
  uint32_t seed;
  struct timeval time2; 
  gettimeofday(&time2,NULL);
  seed = static_cast<unsigned int>((time2.tv_sec * 1000) + (time2.tv_usec / 1000));
  //rand_device >> seed;
  gen.seed(seed);

  /* list of centers */
  vector<vector<double> > L;
  vector<double> tmp;
  vector<double> newpoint;
  vector<double> point2;

  /* print information */
  if (verbose) {
    cout << "# [info] seed = " << seed << endl;
    cout << "# [info] dim = " << n << endl;
    //cout << "# [info] samples = " << k << endl;
    cout << "# [info] bound = " << pow(2, 0.21*n) << endl;
  }
  
  /* main loop */
  int i = 0;
  int flag = 1;
  int cont_cols = 0;
  double dist = 0;
  
  while(1) {

    /* sample new point */
    newpoint = sampler(gen);
    flag = 1;

    /* find near points in L */
    std::vector<vector<double> >::const_iterator it;
    for (it = L.begin(); it != L.end(); ++it) {
      point2 = *it;
      tmp = newpoint + point2;
      dist = norm(tmp);
      if (dist < 1) {
        flag = 0;
        cont_cols ++;
        break;
      }
      tmp = newpoint - point2;
      dist = norm(tmp);
      if (dist < 1) {
        flag = 0;
        cont_cols ++;
        break;
      }
    }

    /* too much collisions */
    //cout << "size " << L.size() << endl;
    //cout << "cols " << cont_cols << endl;    
    if (cont_cols >= L.size()*4 + n)
      break;
    
    /* update if no near points */
    if (flag == 1) {
      //cont_cols = 0;
      L.push_back(newpoint);
    }

    i++;
  }

  /* final stats */
  double ss = L.size();
  if (verbose) {
    cout << "# [info] |L| " << ss << endl;
    cout << "# [info] log2(|L|) / n = " << log(ss)/log(2.0) / n << endl;
  }

  //check_list_2red (L);
  
  return (log(ss)/log(2.0) / n);
}


double cover_sphere_uniform_3red (int n, int verbose)
{

  boost::mt19937 gen;
  boost::uniform_on_sphere<double> sampler(n);
  uint32_t seed;
  struct timeval time2; 
  gettimeofday(&time2,NULL);
  seed = static_cast<unsigned int>((time2.tv_sec * 1000) + (time2.tv_usec / 1000));
  gen.seed(seed);

  /* list of centers */
  vector<vector<double> > L;
  vector<double> tmp;
  vector<double> tmp2;
  vector<double> tmp3;
  vector<double> tmp4;
  vector<double> tmp5;
  vector<double> tmp6;
  vector<double> newpoint;
  vector<double> point2;
  vector<double> point3;

  if (verbose) {
    cout << "# [info] seed = " << seed << endl;
    cout << "# [info] dim = " << n << endl;
    //cout << "# [info] samples = " << k << endl;
    cout << "# [info] bound = " << pow(2, 0.21*n) << endl;
  }
  
  int i = 0;
  int flag = 1;
  int cont_cols = 0;
  double dist = 0;
  
  while(1) {

    newpoint = sampler(gen);
    flag = 1;

    /* check 2-red condition */
    std::vector<vector<double> >::const_iterator it;
    for (it = L.begin(); it != L.end(); ++it) {
      point2 = *it;
      tmp = newpoint + point2;
      dist = norm(tmp);
      if (dist < 1) {
        flag = 0;
        cont_cols ++;
        break;
      }
      tmp = newpoint - point2;
      dist = norm(tmp);
      if (dist < 1) {
        flag = 0;
        cont_cols ++;
        break;
      }
    }

    if (flag == 0)
      continue;

    if (cont_cols >= L.size()*4 + n)
      break;

    /* check 3-red condition */
    if (flag == 1) {
      if (L.size()<=1) {
        L.push_back(newpoint);
        continue;
      }
      std::vector<vector<double> >::const_iterator it1;
      std::vector<vector<double> >::const_iterator it2;
      for (it1 = L.begin(); it1 != L.end(); ++it1) {
        point2 = *it1;
        tmp = newpoint + point2;
        tmp2 = newpoint - point2;
        for (it2 = L.begin(); it2 != L.end(); ++it2) {
          point3 = *it2;

          tmp3 = tmp + point3;
          dist = norm(tmp3);
          if (dist < 1) {
            flag = 0; cont_cols ++;break;
          }
          
          tmp4 = tmp - point3;
          dist = norm(tmp4);
          if (dist < 1) {
            flag = 0; cont_cols ++;break;
          }
          
          tmp5 = tmp2 + point3;
          dist = norm(tmp5);
          if (dist < 1) {
            flag = 0; cont_cols ++;break;
          }
          
          tmp6 = tmp2 - point3;
          dist = norm(tmp6);
          if (dist < 1) {
            flag = 0; cont_cols ++;break;
          }
        }
        if (flag == 0)
          break;
      }

      if (flag == 1) {
        //cont_cols = 0;
        L.push_back(newpoint);
      }
    }
    i++;
  }
  
  double ss = L.size();
  if (verbose) {
    cout << "# [info] |L| " << ss << endl;
    cout << "# [info] log2(|L|) / n = " << log(ss)/log(2.0) / n << endl;
  }

  //check_list_3red (L);
  
  return (log(ss)/log(2.0) / n);
}



int main(int argc, char *argv[])
{
  int updim = 30;
  if (argc > 1) {
    updim = atoi(argv[1]); 
  }

  //double curr = cover_sphere_uniform (updim,1);
  
  ///*
  for (int dim = 8; dim <= updim; dim ++) {
    int times = 200000/dim;
    //int times = 1;
    
    double all = 0.0, curr=0.0;
    for (int i=0; i<times; i++){
      curr = cover_sphere_uniform_3red (dim,0);
      //curr = cover_sphere_uniform (dim,0);
      all += curr;
    }
    //cout << "# Average log2(|L|) / n = " << all / times << endl;
    cout << dim << " " << all / times << endl;
  }
  //*/
  
}


