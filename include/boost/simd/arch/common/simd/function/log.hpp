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
                          , bs::fast_tag
                          , bs::pack_< bd::single_<A0>, X>
                          )
  {
    BOOST_FORCEINLINE A0 operator() (const fast_tag &, const A0& a0) const BOOST_NOEXCEPT
    {
      A0 ln2_hi(6.9313812256e-01), /* 0x3f317180 */
        ln2_lo(9.0580006145e-06), /* 0x3717f7d1 */
        /* |(log(1+s)-log(1-s))/s - Lg(s)| < 2**-34.24 (~[-4.95e-11, 4.97e-11]). */
        Lg1(0.66666662693), /*   0xaaaaaa.0p-24*/
        Lg2(0.40000972152), /*   0xccce13.0p-25*/
        Lg3(0.28498786688), /*   0x91e9ee.0p-25*/
        Lg4(0.24279078841); /*   0xf89e26.0p-26*/
      using uiA0 = bd::as_integer_t<A0, unsigned>;
      using iA0 = bd::as_integer_t<A0,   signed>;
      A0 x =  if_nan_else(is_lez(a0), a0);
      iA0 k(0);
      auto isnez = is_nez(a0);
      auto test = is_less(a0, Smallestposval<A0>())&&isnez;
#ifndef BOOST_SIMD_NO_DENORMALS
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
      A0 w = sqr(z);
      A0 t1= w*fma(w, Lg4, Lg2); //w*(Lg2+w*Lg4);
      A0 t2= z*fma(w, Lg3, Lg1); //z*(Lg1+w*Lg3);
      A0 R = t2 + t1;
      A0 hfsq = Half<A0>()*sqr(f);
      A0 dk = tofloat(k);
      A0 r = fma(s, (hfsq+R), dk*ln2_lo - hfsq + f + dk*ln2_hi);
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
                             , bs::fast_tag
                             , bs::pack_< bd::double_<A0>, X>
                             )
  {
    BOOST_FORCEINLINE A0 operator() (const fast_tag &, const A0& a0) const BOOST_NOEXCEPT
    {
      static const double
        ln2_hi(6.93147180369123816490e-01),  /* 3fe62e42 fee00000 */
        ln2_lo(1.90821492927058770002e-10),  /* 3dea39ef 35793c76 */
        Lg1(6.666666666666735130e-01),  /* 3FE55555 55555593 */
        Lg2(3.999999999940941908e-01),  /* 3FD99999 9997FA04 */
        Lg3(2.857142874366239149e-01),  /* 3FD24924 94229359 */
        Lg4(2.222219843214978396e-01),  /* 3FCC71C5 1D8E78AF */
        Lg5(1.818357216161805012e-01),  /* 3FC74664 96CB03DE */
        Lg6(1.531383769920937332e-01),  /* 3FC39A09 D078C69F */
        Lg7(1.479819860511658591e-01);  /* 3FC2F112 DF3E5244 */
      using uiA0 = bd::as_integer_t<A0, unsigned>;
      using iA0 = bd::as_integer_t<A0,   signed>;
      A0 x =  if_nan_else(is_lez(a0), a0);
      uiA0 hx = bitwise_cast<uiA0>(x) >> 32;
      iA0 k(0);
      auto isnez = is_nez(a0);
      auto test = is_less(a0, Smallestposval<A0>())&&isnez;

#ifndef BOOST_SIMD_NO_DENORMALS
      if (any(test))
      {
        k = if_minus(test, k, iA0(54));
        x = if_else(test, x*A0(0x1p54), x);
      }
#endif
      /* reduce x into [sqrt(2)/2, sqrt(2)] */
      hx += 0x3ff00000 - 0x3fe6a09e;
      k += bitwise_cast<iA0>(hx>>20) - 0x3ff;
      hx = (hx&0x000fffff) + 0x3fe6a09e;
      x = bitwise_cast<A0>(hx<<32 | (bitwise_and(0xffffffffull, bitwise_cast<uiA0>(x))));

      A0 f = dec(x);
      A0 s = f/(2.0f + f);
      A0 z = sqr(s);
      A0 w = sqr(z);
      A0 t1 = w*fma(w, fma(w, Lg6, Lg4), Lg2); //w*(Lg2+w*(Lg4+w*Lg6));
      A0 t2 = z*fma(w, fma(w, fma(w, Lg7, Lg5), Lg3), Lg1); //      z*(Lg1+w*(Lg3+w*(Lg5+w*Lg7)));
      A0 R = t2 + t1;
      A0 hfsq = Half<A0>()*sqr(f);
      A0 dk = tofloat(k);
      A0 r = fma(s, (hfsq+R), dk*ln2_lo - hfsq + f + dk*ln2_hi);
#ifndef BOOST_SIMD_NO_INFINITIES
      return if_else(isnez, if_else(a0 == Inf<A0>(), Inf<A0>(), r), Minf<A0>());
#else
      return if_else(isnez, r, Minf<A0>());
#endif
    }
  };

} } }

#endif
