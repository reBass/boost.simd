//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_DETAIL_GENERIC_EXPO_UTILITIES_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_DETAIL_GENERIC_EXPO_UTILITIES_HPP_INCLUDED
#include <boost/simd/arch/common/detail/tags.hpp>
#include <boost/simd/constant/invlog10_2.hpp>
#include <boost/simd/constant/invlog_2.hpp>
#include <boost/simd/detail/constant/log_2hi.hpp>
#include <boost/simd/detail/constant/log_2lo.hpp>
#include <boost/simd/detail/constant/log10_2hi.hpp>
#include <boost/simd/detail/constant/log10_2lo.hpp>
#include <boost/simd/constant/log_10.hpp>
#include <boost/simd/constant/log_2.hpp>
#include <boost/simd/detail/constant/log_2hi.hpp>
#include <boost/simd/detail/constant/log_2lo.hpp>
#include <boost/simd/detail/constant/maxlog.hpp>
#include <boost/simd/detail/constant/maxlog10.hpp>
#include <boost/simd/detail/constant/maxlog2.hpp>
#include <boost/simd/detail/constant/minlog.hpp>
#include <boost/simd/detail/constant/minlog10.hpp>
#include <boost/simd/detail/constant/minlog2.hpp>
#include <boost/simd/function/fnms.hpp>
#include <boost/simd/function/is_greater_equal.hpp>
#include <boost/simd/function/is_less_equal.hpp>
#include <boost/simd/function/nearbyint.hpp>
#include <boost/simd/function/toint.hpp>

namespace boost { namespace simd
{
  namespace detail
  {
    namespace bd =  boost::dispatch;
    template< typename A0
              , typename Tag
              >
    struct exp_utilities;

    template<class A0> struct exp_utilities <A0,bs::tag::exp_>
    {
      using iA0 = bd::as_integer_t<A0>;
      static BOOST_FORCEINLINE auto isgemaxlog(A0 const& a0) BOOST_NOEXCEPT
        -> decltype(is_greater_equal(a0, Maxlog<A0>()))
      {
        return is_greater_equal(a0, Maxlog<A0>());
      }

      static BOOST_FORCEINLINE auto isleminlog(A0 const& a0) BOOST_NOEXCEPT
      -> decltype(is_less_equal(a0, Minlog<A0>()))
      {
        return is_less_equal(a0, Minlog<A0>());
      }
      static BOOST_FORCEINLINE iA0 reduce(A0 const& a0
                                         , A0& hi, A0& lo, A0& x) BOOST_NOEXCEPT
      {
        A0 k = nearbyint(Invlog_2<A0>()*a0);
        hi = fnms(k, Log_2hi<A0>(), a0); //a0-k*L
        lo = k*Log_2lo<A0>();
        x  = hi-lo;
        return toint(k);
      }
    };

    template<class A0> struct exp_utilities <A0,bs::tag::exp2_>
    {
      using iA0 = bd::as_integer_t<A0>;
      static BOOST_FORCEINLINE auto isgemaxlog(A0 const& a0) BOOST_NOEXCEPT
        -> decltype(is_greater_equal(a0, Maxlog2<A0>()))
      {
        return is_greater_equal(a0, Maxlog2<A0>());
      }

      static BOOST_FORCEINLINE auto isleminlog(A0 const& a0) BOOST_NOEXCEPT
      -> decltype(is_less_equal(a0, Minlog2<A0>()))
      {
        return is_less_equal(a0, Minlog2<A0>());
      }
      static BOOST_FORCEINLINE iA0 reduce(A0 const& a0
                                        , A0&, A0&, A0& x) BOOST_NOEXCEPT
      {
        A0 k = nearbyint(a0);
        x = (a0 - k);
        return toint(k);
      }
    };

    template<class A0> struct exp_utilities <A0,bs::tag::exp10_>
    {
      using iA0 = bd::as_integer_t<A0>;
      static BOOST_FORCEINLINE auto isgemaxlog(A0 const& a0) BOOST_NOEXCEPT
        -> decltype(is_greater_equal(a0, Maxlog10<A0>()))
      {
        return is_greater_equal(a0, Maxlog10<A0>());
      }

      static BOOST_FORCEINLINE auto isleminlog(A0 const& a0) BOOST_NOEXCEPT
      -> decltype(is_less_equal(a0, Minlog10<A0>()))
      {
        return is_less_equal(a0, Minlog10<A0>());
      }
      static BOOST_FORCEINLINE iA0 reduce(A0 const& a0
                                        , A0&, A0&, A0& x) BOOST_NOEXCEPT
      {
        A0 k = nearbyint(Invlog10_2<A0>()*a0);
        x = fnms(k, Log10_2hi<A0>(), a0);
        x -= k*Log10_2lo<A0>();
        return toint(k);
      }
    };

} } }

#endif
