// -------------------------------------------------------------------------------------------------
//                              Copyright 2016 - NumScale SAS
//
//                   Distributed under the Boost Software License, Version 1.0.
//                        See accompanying file LICENSE.txt or copy at
//                            http://www.boost.org/LICENSE_1_0.txt
// -------------------------------------------------------------------------------------------------

#include <simd_bench.hpp>
#include <boost/simd/function/simd/lookup.hpp>
#include <boost/simd/pack.hpp>

namespace nsb = ns::bench;
namespace bs =  boost::simd;

struct look
{
  template<class T> T operator()(const T & a) const
  {
    using iT =  bd::as_integer_t<T>;
    return bs::lookup(a, iT(1));
  }
};

DEFINE_SIMD_BENCH(simd_lookup, look);

DEFINE_BENCH_MAIN()
{
  nsb::for_each<simd_lookup, NS_BENCH_IEEE_TYPES>(-10, 10);
}
