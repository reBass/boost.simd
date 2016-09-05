//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_SCALAR_FUNCTION_LOG_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_SCALAR_FUNCTION_LOG_HPP_INCLUDED
#include <boost/simd/function/std.hpp>

#include <boost/simd/arch/common/detail/scalar/logarithm.hpp>
#include <boost/simd/detail/dispatch/function/overload.hpp>
#include <boost/assert.hpp>
#include <boost/config.hpp>
#include <cmath>

namespace boost { namespace simd { namespace ext
{
  namespace bd = boost::dispatch;
  namespace bs = boost::simd;

  BOOST_DISPATCH_OVERLOAD ( log_
                          , (typename A0)
                          , bd::cpu_
                          , bs::std_tag
                          , bd::scalar_< bd::floating_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() (const std_tag &, A0 a0) const BOOST_NOEXCEPT
    {
      return std::log(a0);
    }
  };
  BOOST_DISPATCH_OVERLOAD ( log_
                          , (typename A0)
                          , bd::cpu_
                          , bd::scalar_< bd::floating_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() (A0 a0) const BOOST_NOEXCEPT
    {
      return detail::logarithm<A0,tag::not_simd_type>::log(a0);
    }
  };

  BOOST_DISPATCH_OVERLOAD ( log_
                          , (typename A0)
                          , bd::cpu_
                          , bs::fast_tag
                          , bd::scalar_< bd::single_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() (const fast_tag &, A0 x) const BOOST_NOEXCEPT
    {
      A0 ln2_hi = 6.9313812256e-01, /* 0x3f317180 */
        ln2_lo = 9.0580006145e-06, /* 0x3717f7d1 */
        /* |(log(1+s)-log(1-s))/s - Lg(s)| < 2**-34.24 (~[-4.95e-11, 4.97e-11]). */
        Lg1 =0.66666662693, /*   0xaaaaaa.0p-24*/
        Lg2 =0.40000972152, /*   0xccce13.0p-25*/
        Lg3 =0.28498786688, /*   0x91e9ee.0p-25*/
        Lg4 =0.24279078841; /*   0xf89e26.0p-26*/
      using uiA0 = bd::as_integer_t<A0, unsigned>;
      using iA0 = bd::as_integer_t<A0,   signed>;
      uiA0 ix = bitwise_cast<uiA0>(x);
      iA0 k = 0;
      if (ix < 0x00800000 || ix>>31) {  /* x < 2**-126  */
        if (ix<<1 == 0)
          return Minf<A0>();  /* log(+-0)=-inf */
        if (ix>>31)
          return Nan<A0>(); /* log(-#) = NaN */
        /* subnormal number, scale up x */
        k -= 25;
        x *= 0x1p25f;
        ix = bitwise_cast<iA0>(x);
      } else if (ix >= 0x7f800000) {
        return x;
      } else if (ix == 0x3f800000)
        return 0;
      /* reduce x into [sqrt(2)/2, sqrt(2)] */
      ix += 0x3f800000 - 0x3f3504f3;
      k += bitwise_cast<iA0>(ix>>23) - 0x7f;
      ix = (ix&0x007fffff) + 0x3f3504f3;
      x =  bitwise_cast<A0>(ix);
      A0 f = dec(x);
      A0 s = f/(2.0f + f);
      A0 z = sqr(s);
      A0 w = sqr(z);
      A0 t1= w*fma(w, Lg4, Lg2); //w*(Lg2+w*Lg4);
      A0 t2= z*fma(w, Lg3, Lg1); //z*(Lg1+w*Lg3);
      A0 R = t2 + t1;
      A0 hfsq = Half<A0>()*sqr(f);
      return fma(s, (hfsq+R), k*ln2_lo - hfsq + f + k*ln2_hi);
    }
  };
} } }


#endif
