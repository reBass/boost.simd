// -------------------------------------------------------------------------------------------------
//                              Copyright 2016 - NumScale SAS
//
//                   Distributed under the Boost Software License, Version 1.0.
//                        See accompanying file LICENSE.txt or copy at
//                            http://www.boost.org/LICENSE_1_0.txt
// -------------------------------------------------------------------------------------------------

#include <simd_bench.hpp>
#include <boost/simd/function/scalar/acos.hpp>
#include <boost/simd/function/accurate.hpp>

namespace nsb = ns::bench;
namespace bs =  boost::simd;

DEFINE_SCALAR_BENCH(scalar_accurate_acos, bs::accurate_(bs::acos));

DEFINE_BENCH_MAIN()
{
  nsb::for_each<scalar_accurate_acos, NS_BENCH_IEEE_TYPES>(-1, 1);
}
