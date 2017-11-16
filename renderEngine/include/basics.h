#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <iostream>

#include "/Code/Eigen/Eigen"
#include "/Code/Eigen/StdVector"
#include "/Code/Eigen/SparseCore"
#include "Vector3.h"
using namespace std;
using namespace Eigen;
#define pi 3.1415926
#define ASSERT(x) if (!(x)) cout << "error" << endl;

inline double sign(double x){ return x > 0.0 ? 1.0 : -1.0; }

// TODO: shouldn't be necessary, given templates below
inline double clamped(double value, double left, double right)
{
  return value<left ? left : (value > right ? right : value);
}
inline double random()
{
  return -1.0 + (2.0 * (double)rand() / (double)RAND_MAX);
}
inline double random(double min, double max)
{
  return min + (max-min)*((double)rand() / (double)RAND_MAX);
}
template<class T> inline const T& max(const T& left, const T& right)
{	
	return right>left ? right : left;
}
template<class T> inline const T& min(const T& left, const T& right)
{	
	return right<left ? right : left;
}
template<class T> inline const T& clamped(const T& value, const T& left, const T& right)
{	
  return value<left ? left : (value > right ? right : value);
}
inline double sqr(double a)
{
  return a*a;
}
