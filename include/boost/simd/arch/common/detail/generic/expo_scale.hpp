//==================================================================================================
/**
  Copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
**/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_DETAIL_GENERIC_EXPO_SCALE_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_DETAIL_GENERIC_EXPO_SCALE_HPP_INCLUDED

#include <boost/simd/detail/constant/maxexponent.hpp>
#include <boost/simd/constant/nbmantissabits.hpp>
#include <boost/simd/function/simd/bitwise_cast.hpp>
#include <boost/simd/function/simd/shift_left.hpp>
#include <boost/simd/detail/dispatch/meta/scalar_of.hpp>

namespace boost { namespace simd { namespace detail
{
  template<typename A0, typename A1>
  BOOST_FORCEINLINE A0 scale(A0 const & y, const A1& k)
  {
    auto ik = k + Maxexponent<A0>();
    ik = shift_left(ik, Nbmantissabits<dispatch::scalar_of_t<A0>>());
    return y*bitwise_cast<A0>(ik);
  }
} } }

#endif
