//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_SIMD_FUNCTION_LOG1P_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_SIMD_FUNCTION_LOG1P_HPP_INCLUDED
#include <boost/simd/detail/overload.hpp>

#include <boost/simd/meta/hierarchy/simd.hpp>
#include <boost/simd/function/is_nez.hpp>
#include <boost/simd/function/log.hpp>
#include <boost/simd/function/dec.hpp>
#include <boost/simd/function/inc.hpp>
#include <boost/simd/function/if_plus.hpp>
#include <boost/simd/function/if_zero_else.hpp>
#include <boost/simd/function/is_lez.hpp>

#ifndef BOOST_SIMD_NO_INFINITIES
#include <boost/simd/constant/inf.hpp>
#include <boost/simd/function/if_else.hpp>
#include <boost/simd/function/is_equal.hpp>
#endif

namespace boost { namespace simd { namespace ext
{
  namespace bd = boost::dispatch;

  BOOST_DISPATCH_OVERLOAD_IF ( log1p_
                             , (typename A0, typename X)
                             , (detail::is_native<X>)
                             , bd::cpu_
                             , bs::pack_< bd::floating_<A0>, X>
                          )
  {
    BOOST_FORCEINLINE A0 operator()( const A0& a0) BOOST_NOEXCEPT
    {
      A0 u = inc(a0);
      A0 r = if_plus(is_nez(u),
                    log(u),
                    (a0-dec(u))/u); // cancels errors with IEEE arithmetic
#ifndef BOOST_SIMD_NO_INFINITIES
      r = if_else(is_equal(u, Inf<A0>()),u, r);
#endif
      return r;
    }
  };

  BOOST_DISPATCH_OVERLOAD_IF ( log1p_
                          , (typename A0, typename X)
                          , (detail::is_native<X>)
                          , bd::cpu_
                          , bs::musl_tag
                          , bs::pack_< bd::single_<A0>, X>
                          )
  {
    BOOST_FORCEINLINE A0 operator() (const musl_tag &, const A0& a0) const BOOST_NOEXCEPT
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
      const A0 inca0 =  inc(a0);
      A0 x =  if_nan_else(is_lez(inca0), a0);
      auto isnez = is_nez(inca0);

      A0 uf =  inc(x);
      uiA0 iu = bitwise_cast<uiA0>(uf);
      iu += 0x3f800000 - 0x3f3504f3;
      iA0 k = bitwise_cast<iA0>(iu>>23) - 0x7f;
      /* correction term ~ log(1+x)-log(u), avoid underflow in c/u */
      A0  c = if_else_zero(k < 25, if_else( k >= 2, oneminus(uf-x), x-dec(uf))/uf);

      iu = (iu&0x007fffff) + 0x3f3504f3;
      A0 f =  dec(bitwise_cast<A0>(iu));

      A0 s = f/(2.0f + f);
      A0 z = sqr(s);
      A0 w = sqr(z);
      A0 t1= w*fma(w, Lg4, Lg2); //w*(Lg2+w*Lg4);
      A0 t2= z*fma(w, Lg3, Lg1); //z*(Lg1+w*Lg3);
      A0 R = t2 + t1;
      A0 hfsq = Half<A0>()*sqr(f);
      A0 dk = tofloat(k);
      A0 r =  fma(s, hfsq+R,  fma(dk, ln2_lo, c) - hfsq + f + dk*ln2_hi);
#ifndef BOOST_SIMD_NO_INFINITIES
      return if_else(isnez, if_else(a0 == Inf<A0>(), Inf<A0>(), r), Minf<A0>());
#else
      return if_else(isnez, r, Minf<A0>());
#endif
    }
  };

  BOOST_DISPATCH_OVERLOAD_IF ( log1p_
                             , (typename A0, typename X)
                             , (detail::is_native<X>)
                             , bd::cpu_
                             , bs::musl_tag
                             , bs::pack_< bd::double_<A0>, X>
                             )
  {
    BOOST_FORCEINLINE A0 operator() (const musl_tag &, const A0& a0) const BOOST_NOEXCEPT
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
      const A0 inca0 =  inc(a0);
      A0 x =  if_nan_else(is_lez(inca0), a0);
      auto isnez = is_nez(inca0);

      A0 uf =  inc(x);
      uiA0 hu = bitwise_cast<uiA0>(uf)>>32;
      hu += 0x3ff00000 - 0x3fe6a09e;
      iA0 k = bitwise_cast<iA0>(hu>>20) - 0x3ff;
      /* correction term ~ log(1+x)-log(u), avoid underflow in c/u */
      A0  c = if_else_zero(k < 54, if_else( k >= 2, oneminus(uf-x), x-dec(uf))/uf);

      hu =  (hu&0x000fffff) + 0x3fe6a09e;
      A0 f = bitwise_cast<A0>( bitwise_cast<uiA0>(hu<<32) | (bitwise_and(0xffffffffull, bitwise_cast<uiA0>(f))));
      f = dec(f);

      A0 hfsq = Half<A0>()*sqr(f);
      A0 s = f/(2.0f + f);
      A0 z = sqr(s);
      A0 w = sqr(z);
      A0 t1 = w*fma(w, fma(w, Lg6, Lg4), Lg2);
      A0 t2 = z*fma(w, fma(w, fma(w, Lg7, Lg5), Lg3), Lg1);
      A0 R = t2 + t1;
      A0 dk = tofloat(k);
      A0 r =  fma(s, hfsq+R,  fma(dk, ln2_lo, c) - hfsq + f + dk*ln2_hi);

#ifndef BOOST_SIMD_NO_INFINITIES
      return if_else(isnez, if_else(a0 == Inf<A0>(), Inf<A0>(), r), Minf<A0>());
#else
      return if_else(isnez, r, Minf<A0>());
#endif
    }
  };


} } }

#endif
