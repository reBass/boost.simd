//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_SCALAR_FUNCTION_EXP10_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_SCALAR_FUNCTION_EXP10_HPP_INCLUDED

#include <boost/simd/detail/dispatch/function/overload.hpp>
#include <boost/simd/arch/common/detail/scalar/exponential.hpp>
#include <boost/simd/function/scalar/modf.hpp>
#include <boost/simd/function/musl.hpp>
#include <boost/simd/function/scalar/exp2.hpp>
#include <boost/simd/arch/common/detail/tags.hpp>
#include <boost/config.hpp>

namespace boost { namespace simd { namespace ext
{
  namespace bd = boost::dispatch;
  namespace bs = boost::simd;

  BOOST_DISPATCH_OVERLOAD ( exp10_
                          , (typename A0)
                          , bd::cpu_
                          , bd::scalar_< bd::floating_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() (A0 a0) const BOOST_NOEXCEPT
    {
      return detail::exponential<A0,bs::tag::exp10_,tag::not_simd_type>::expa(a0);
    }
  };

  BOOST_DISPATCH_OVERLOAD ( exp10_
                          , (typename A0)
                          , bd::cpu_
                          , bs::musl_tag
                          , bd::scalar_< bd::single_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() (const musl_tag &, A0 x) const BOOST_NOEXCEPT
    {
      const float p10[] = {
        1e-7f, 1e-6f, 1e-5f, 1e-4f, 1e-3f, 1e-2f, 1e-1f,
        1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7
      };
      float n, y;
      std::tie(y, n)= bs::modf(x);
      using uiA0 =  bd::as_integer_t<A0>;
      uiA0 u = bitwise_cast<uiA0>(n);
      /* abs(n) < 8 without raising invalid on nan */
      if ((u>>23 & 0xff) < 0x7f+3) {
        if (!y) return p10[(int)n+7];
        y = bs::exp2(3.32192809488736234787031942948939f * y);
        return y * p10[(int)n+7];
      }
      return bs::exp2(3.32192809488736234787031942948939 * x);
    }
  };


} } }


#endif
