//==================================================================================================
/**
  Copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
**/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_SCALAR_FUNCTION_COMPARE_NOT_EQUAL_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_SCALAR_FUNCTION_COMPARE_NOT_EQUAL_HPP_INCLUDED

#include <boost/simd/logical.hpp>
#include <boost/simd/detail/dispatch/function/overload.hpp>
#include <boost/config.hpp>

namespace boost { namespace simd { namespace ext
{
  namespace bd = boost::dispatch;
  BOOST_DISPATCH_OVERLOAD ( compare_not_equal_
                          , (typename A0)
                          , bd::cpu_
                          , bd::scalar_< bd::unspecified_<A0> >
                          , bd::scalar_< bd::unspecified_<A0> >
                          )
  {
    BOOST_FORCEINLINE bool operator() ( A0 const& a0, A0 const& a1) const BOOST_NOEXCEPT
    {
      return a0 != a1;
    }
  };
} } }

#endif
