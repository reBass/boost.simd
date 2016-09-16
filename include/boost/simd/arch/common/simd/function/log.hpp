//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_SIMD_FUNCTION_LOG_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_SIMD_FUNCTION_LOG_HPP_INCLUDED

#include <boost/simd/arch/common/detail/scalar/logarithm.hpp>
#include <boost/simd/arch/common/detail/tags.hpp>
#include <boost/simd/arch/common/detail/simd/logarithm.hpp>
#include <boost/simd/detail/dispatch/function/overload.hpp>
#include <boost/config.hpp>
#include <boost/simd/function/simd/is_lez.hpp>
#include <boost/simd/function/simd/any.hpp>
#include <boost/simd/function/musl.hpp>
#include <boost/simd/function/std.hpp>
#include <boost/simd/detail/constant/log_2hi.hpp>
#include <boost/simd/detail/constant/log_2lo.hpp>

namespace boost { namespace simd { namespace ext
{
  namespace bd = boost::dispatch;
  namespace bs = boost::simd;

  BOOST_DISPATCH_OVERLOAD_IF ( log_
                          , (typename A0, typename X)
                          , (detail::is_native<X>)
                          , bd::cpu_
                          , bs::pack_< bd::floating_<A0>, X>
                          )
  {
    BOOST_FORCEINLINE A0 operator() (A0 const & a0) const BOOST_NOEXCEPT
    {
      return detail::logarithm<A0,tag::simd_type>::log(a0);
    }
  };

  BOOST_DISPATCH_OVERLOAD_IF ( log_
                          , (typename A0, typename X)
                          , (detail::is_native<X>)
                          , bd::cpu_
                          , bs::musl_tag
                          , bs::pack_< bd::single_<A0>, X>
                          )
  {
    BOOST_FORCEINLINE A0 operator() (const musl_tag &, const A0& a0) const BOOST_NOEXCEPT
    {
        /* |(log(1+s)-log(1-s))/s - Lg(s)| < 2**-34.24 (~[-4.95e-11, 4.97e-11]). */
      using uiA0 = bd::as_integer_t<A0, unsigned>;
      using iA0 = bd::as_integer_t<A0,   signed>;
      A0 x =  if_nan_else(is_lez(a0), a0);
      iA0 k(0);
      auto isnez = is_nez(a0);
#ifndef BOOST_SIMD_NO_DENORMALS
      auto test = is_less(a0, Smallestposval<A0>())&&isnez;
      if (any(test))
      {
        k = if_minus(test, k, iA0(25));
        x = if_else(test, x*A0(0x1p25f), x);
      }
#endif
      uiA0 ix = bitwise_cast<uiA0>(x);
      /* reduce x into [sqrt(2)/2, sqrt(2)] */
      ix += 0x3f800000 - 0x3f3504f3;
      k += bitwise_cast<iA0>(ix>>23) - 0x7f;
      ix = (ix&0x007fffff) + 0x3f3504f3;
      x =  bitwise_cast<A0>(ix);
      A0 f = dec(x);
      A0 s = f/(2.0f + f);
      A0 z = sqr(s);
      A0 R =  horn<A0
        , 0x3f2aaaaa  //  0.66666662693 0xaaaaaa.0p-24
        , 0x3eccce13  //  0.40000972152 0xccce13.0p-25
        , 0x3e91e9ee  //  0.28498786688 0x91e9ee.0p-25
        , 0x3e789e26  //  0.24279078841 0xf89e26.0p-26
        >(z)*z;
      A0 hfsq = Half<A0>()*sqr(f);
      A0 dk = tofloat(k);
      A0 r = fma(s, (hfsq+R), dk*Log_2lo<A0>() - hfsq + f + dk*Log_2hi<A0>());
#ifndef BOOST_SIMD_NO_INFINITIES
      return if_else(isnez, if_else(a0 == Inf<A0>(), Inf<A0>(), r), Minf<A0>());
#else
      return if_else(isnez, r, Minf<A0>());
#endif
    }
  };

  BOOST_DISPATCH_OVERLOAD_IF ( log_
                             , (typename A0, typename X)
                             , (detail::is_native<X>)
                             , bd::cpu_
                             , bs::musl_tag
                             , bs::pack_< bd::double_<A0>, X>
                             )
  {
    BOOST_FORCEINLINE A0 operator() (const musl_tag &, const A0& a0) const BOOST_NOEXCEPT
    {
      using uiA0 = bd::as_integer_t<A0, unsigned>;
      using iA0 = bd::as_integer_t<A0,   signed>;
      A0 x =  if_nan_else(is_lez(a0), a0);
      uiA0 hx = bitwise_cast<uiA0>(x) >> 32;
      iA0 k(0);
      auto isnez = is_nez(a0);

#ifndef BOOST_SIMD_NO_DENORMALS
      auto test = is_less(a0, Smallestposval<A0>())&&isnez;
      if (any(test))
      {
        k = if_minus(test, k, iA0(54));
        x = if_else(test, x*A0(0x1p54), x);
      }
#endif
      /* reduce x into [sqrt(2)/2, sqrt(2)] */
      hx += 0x3ff00000 - 0x3fe6a09e;
      k += bitwise_cast<iA0>(hx>>20) - 0x3ff;
      A0 dk = tofloat(k);
      hx = (hx&0x000fffff) + 0x3fe6a09e;
      x = bitwise_cast<A0>(hx<<32 | (bitwise_and(0xffffffffull, bitwise_cast<uiA0>(x))));

      A0 f = dec(x);
      A0 hfsq = Half<A0>()*sqr(f);
      A0 s = f/(2.0f + f);
      A0 z = sqr(s);
//       A0 w = sqr(z);
//       A0 t1 = w*fma(w, fma(w, Lg6, Lg4), Lg2); //w*(Lg2+w*(Lg4+w*Lg6));
//       A0 t2 = z*fma(w, fma(w, fma(w, Lg7, Lg5), Lg3), Lg1); //      z*(Lg1+w*(Lg3+w*(Lg5+w*Lg7)));
//       A0 R = t2 + t1;
      A0 R =  horn<A0
        , 0x3FE5555555555593     //  6.666666666666735130e-01
        , 0x3FD999999997FA04     //  3.999999999940941908e-01
        , 0x3FD2492494229359     //  2.857142874366239149e-01
        , 0x3FCC71C51D8E78AF     //  2.222219843214978396e-01
        , 0x3FC7466496CB03DE     //  1.818357216161805012e-01
        , 0x3FC39A09D078C69F     //  1.531383769920937332e-01
        , 0x3FC2F112DF3E5244     //  1.479819860511658591e-01
        >(z)*z;
      A0 r = fma(s, (hfsq+R), fma(dk, Log_2lo<A0>(), - hfsq + f + dk*Log_2hi<A0>()));
//      A0 r = fma(dk, ln2_hi, fma(s, (hfsq+R), fma(dk, ln2_lo, - hfsq + f)));

#ifndef BOOST_SIMD_NO_INFINITIES
      return if_else(isnez, if_else(a0 == Inf<A0>(), Inf<A0>(), r), Minf<A0>());
#else
      return if_else(isnez, r, Minf<A0>());
#endif
    }
  };

} } }

#endif
  /*
   *   1. Argument Reduction: find k and f such that
   *                      x = 2^k * (1+f),
   *         where  sqrt(2)/2 < 1+f < sqrt(2) .
   *
   *   2. Approximation of log(1+f).
   *      Let s = f/(2+f) ; based on log(1+f) = log(1+s) - log(1-s)
   *               = 2s + 2/3 s**3 + 2/5 s**5 + .....,
   *               = 2s + s*R
   *      We use a special Remez algorithm on [0,0.1716] to generate
   *      a polynomial of degree 14 to approximate R The maximum error
   *      of this polynomial approximation is bounded by 2**-58.45. In
   *      other words,
   *                      2      4      6      8      10      12      14
   *          R(z) ~ Lg1*s +Lg2*s +Lg3*s +Lg4*s +Lg5*s  +Lg6*s  +Lg7*s
   *      (the values of Lg1 to Lg7 are listed in the program)
   *      and
   *          |      2          14          |     -58.45
   *          | Lg1*s +...+Lg7*s    -  R(z) | <= 2
   *          |                             |
   *      Note that 2s = f - s*f = f - hfsq + s*hfsq, where hfsq = f*f/2.
   *      In order to guarantee error in log below 1ulp, we compute log
   *      by
   *              log(1+f) = f - s*(f - R)        (if f is not too large)
   *              log(1+f) = f - (hfsq - s*(hfsq+R)).     (better accuracy)
   *
   *      3. Finally,  log(x) = k*ln2 + log(1+f).
   *                          = k*ln2_hi+(f-(hfsq-(s*(hfsq+R)+k*ln2_lo)))
   *         Here ln2 is split into two floating point number:
   *                      ln2_hi + ln2_lo,
   *         where n*ln2_hi is always exact for |n| < 2000.
   */
