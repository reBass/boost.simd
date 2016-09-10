
//==============================================================================
//         Copyright 2016        Numscale SAS
//
//          Distributed under the Boost Software License, Version 1.0.
//                 See accompanying file LICENSE.txt or copy at
//                     http://www.boost.org/LICENSE_1_0.txt
//==============================================================================
#include <boost/simd/function/simd/acos.hpp>
#include <boost/simd/constant/mone.hpp>
#include <boost/simd/constant/one.hpp>
#include <boost/simd/constant/eps.hpp>
#include <boost/simd/function/next.hpp>
#include <boost/simd/function/exp2.hpp>
#include <boost/simd/function/ilog2.hpp>
#include <boost/simd/pack.hpp>
#include <exhaustive.hpp>

#include <cmath>
#include <cstdlib>


int main(int argc, char* argv[])
{
  float mini = bs::Mone<float>(); // acos is Nan under
  float maxi = bs::One<float>();  // acos is Nan above
  if(argc >= 2) mini = std::atof(argv[1]);
  if(argc >= 3) maxi = std::atof(argv[2]);
  uint64_t n[100];

 for(int i=0; i <100 ; i++)
 {
   n[i] =0;
 }
  for(float i=mini; i < maxi; i = bs::next(i))
  {
    double di =  double(i)+bs::Eps<double>();
    uint64_t z =  bs::min(uint64_t(99), uint64_t(bs::ilog2(2*bs::ulpdist(bs::accurate_(bs::acos)(di), std::acos(di))+1)));
    n[z]++;
 }
 for(int i=0; i <100 ; i++)
 {
   if(n[i])
     std::cout <<  bs::exp2(float(i))/2 << "  ->" << n[i] << std::endl;
 }

  return 0;
}
