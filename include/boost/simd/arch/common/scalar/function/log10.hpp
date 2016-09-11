//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_SCALAR_FUNCTION_LOG10_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_SCALAR_FUNCTION_LOG10_HPP_INCLUDED
#include <boost/simd/function/std.hpp>


#include <boost/simd/arch/common/detail/scalar/logarithm.hpp>
#include <boost/simd/arch/common/detail/tags.hpp>
#include <boost/simd/detail/dispatch/function/overload.hpp>
#include <boost/simd/detail/dispatch/meta/as_integer.hpp>
#include <boost/simd/function/musl.hpp>
#include <boost/config.hpp>
#include <cmath>

namespace boost { namespace simd { namespace ext
{
  namespace bd = boost::dispatch;
  namespace bs = boost::simd;

  BOOST_DISPATCH_OVERLOAD ( log10_
                          , (typename A0)
                          , bd::cpu_
                          , bd::scalar_< bd::int_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() ( A0 a0) const BOOST_NOEXCEPT
    {
      using ui_t = bd::as_integer_t<A0, unsigned>;
      BOOST_ASSERT_MSG( a0 > 0, "log10 is not defined for zero or negative integers." );
      return static_cast<A0>(log10(static_cast<ui_t>(a0)));
    }
  };
  BOOST_DISPATCH_OVERLOAD ( log10_
                          , (typename A0)
                          , bd::cpu_
                          , bd::scalar_< bd::uint8_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() ( A0 a0) const BOOST_NOEXCEPT
    {
      return (a0 >= 100u) ? 2 :
             (a0 >= 10u)  ? 1 : 0;
    }
  };
  BOOST_DISPATCH_OVERLOAD ( log10_
                          , (typename A0)
                          , bd::cpu_
                          , bd::scalar_< bd::uint16_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() ( A0 a0) const BOOST_NOEXCEPT
    {
      return (a0 >= 10000u) ? 4 :
             (a0 >= 1000u)  ? 3 :
             (a0 >= 100u)   ? 2 :
             (a0 >= 10u)    ? 1 : 0;
    }

  };
  BOOST_DISPATCH_OVERLOAD ( log10_
                          , (typename A0)
                          , bd::cpu_
                          , bd::scalar_< bd::uint32_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() ( A0 a0) const BOOST_NOEXCEPT
    {
      return (a0 >= 1000000000u) ? 9 :
             (a0 >= 100000000u) ?  8 :
             (a0 >= 10000000u) ?   7 :
             (a0 >= 1000000u) ?    6 :
             (a0 >= 100000u) ?     5 :
             (a0 >= 10000u) ?      4 :
             (a0 >= 1000u) ?       3 :
             (a0 >= 100u) ?        2 :
             (a0 >= 10u) ?         1 : 0;
    }
  };
  BOOST_DISPATCH_OVERLOAD ( log10_
                          , (typename A0)
                          , bd::cpu_
                          , bd::scalar_< bd::uint64_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() ( A0 a0) const BOOST_NOEXCEPT
    {
      return(a0 >= 10000000000000000000ull) ? 19 :
            (a0 >= 1000000000000000000ull) ?  18 :
            (a0 >= 100000000000000000ull) ?   17 :
            (a0 >= 10000000000000000ull) ?    16 :
            (a0 >= 1000000000000000ull) ?     15 :
            (a0 >= 100000000000000ull) ?      14 :
            (a0 >= 10000000000000ull) ?       13 :
            (a0 >= 1000000000000ull) ?        12 :
            (a0 >= 100000000000ull) ?         11 :
            (a0 >= 10000000000ull) ?          10 :
            (a0 >= 1000000000ull) ?            9 :
            (a0 >= 100000000ull) ?             8 :
            (a0 >= 10000000ull) ?              7 :
            (a0 >= 1000000ull) ?               6 :
            (a0 >= 100000ull) ?                5 :
            (a0 >= 10000ull) ?                 4 :
            (a0 >= 1000ull) ?                  3 :
            (a0 >= 100ull) ?                   2 :
            (a0 >= 10ull) ?                    1 : 0;
    }
  };
  BOOST_DISPATCH_OVERLOAD ( log10_
                          , (typename A0)
                          , bd::cpu_
                          , bs::std_tag
                          , bd::scalar_< bd::floating_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() (const std_tag &, A0 a0) const BOOST_NOEXCEPT
    {
      return std::log10(a0);
    }
  };

  BOOST_DISPATCH_OVERLOAD ( log10_
                          , (typename A0)
                          , bd::cpu_
                          , bd::scalar_< bd::floating_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() ( A0 a0) const BOOST_NOEXCEPT
    {
      return detail::logarithm<A0,tag::not_simd_type>::log10(a0);
    }
  };

  BOOST_DISPATCH_OVERLOAD ( log10_
                          , (typename A0)
                          , bd::cpu_
                          , bs::musl_tag
                          , bd::scalar_< bd::single_<A0> >
                          )
  {
  /* origin: FreeBSD /usr/src/lib/msun/src/e_log10f.c */
  /*
   * ====================================================
   * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
   *
   * Developed at SunPro, a Sun Microsystems, Inc. business.
   * Permission to use, copy, modify, and distribute this
   * software is freely granted, provided that this notice
   * is preserved.
   * ====================================================
   */
    BOOST_FORCEINLINE A0 operator() (const musl_tag &, A0 x) const BOOST_NOEXCEPT
    {
      const A0
        ivln10hi  =  4.3432617188e-01, /* 0x3ede6000 */
        ivln10lo  = -3.1689971365e-05, /* 0xb804ead9 */
        log10_2hi =  3.0102920532e-01, /* 0x3e9a2080 */
        log10_2lo =  7.9034151668e-07, /* 0x355427db */
        //idem log/////////////////////////////////////////////////////////////////////////////////////
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
      //idem log fin/////////////////////////////////////////////////////////////////////////////////////

      A0 hi = f - hfsq;
      hi =  bitwise_and(hi, uiA0(0xfffff000ul));
      A0  lo = fma(s, hfsq+R, f - hi - hfsq);
      A0 dk = k;
      return fma(dk, log10_2lo, fma(lo+hi, ivln10lo, fma(lo, ivln10hi, fma(hi, ivln10hi, dk*log10_2hi))));
    }
  };

   BOOST_DISPATCH_OVERLOAD ( log10_
                          , (typename A0)
                          , bd::cpu_
                          , bs::musl_tag
                          , bd::scalar_< bd::double_<A0> >
                          )
  {
  /* origin: FreeBSD /usr/src/lib/msun/src/e_log10.c */
  /*
   * ====================================================
   * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
   *
   * Developed at SunSoft, a Sun Microsystems, Inc. business.
   * Permission to use, copy, modify, and distribute this
   * software is freely granted, provided that this notice
   * is preserved.
   * ====================================================
   */
    BOOST_FORCEINLINE A0 operator() (const musl_tag &, A0 x) const BOOST_NOEXCEPT
    {
      const A0
        ivln10hi  = 4.34294481878168880939e-01, /* 0x3fdbcb7b, 0x15200000 */
        ivln10lo  = 2.50829467116452752298e-11, /* 0x3dbb9438, 0xca9aadd5 */
        log10_2hi = 3.01029995663611771306e-01, /* 0x3FD34413, 0x509F6000 */
        log10_2lo = 3.69423907715893078616e-13, /* 0x3D59FEF3, 0x11F12B36 */
        //idem log/////////////////////////////////////////////////////////////////////////////////////
        Lg1 = 6.666666666666735130e-01,  /* 3FE55555 55555593 */
        Lg2 = 3.999999999940941908e-01,  /* 3FD99999 9997FA04 */
        Lg3 = 2.857142874366239149e-01,  /* 3FD24924 94229359 */
        Lg4 = 2.222219843214978396e-01,  /* 3FCC71C5 1D8E78AF */
        Lg5 = 1.818357216161805012e-01,  /* 3FC74664 96CB03DE */
        Lg6 = 1.531383769920937332e-01,  /* 3FC39A09 D078C69F */
        Lg7 = 1.479819860511658591e-01;  /* 3FC2F112 DF3E5244 */
      using uiA0 = bd::as_integer_t<A0, unsigned>;
      using iA0 = bd::as_integer_t<A0,   signed>;
      uiA0 hx = bitwise_cast<uiA0>(x) >> 32;
      iA0 k = 0;
      if (hx < 0x00100000 || hx>>31) {
        if(is_eqz(x))
          return Minf<A0>();  /* log(+-0)=-inf */
        if (hx>>31)
          return Nan<A0>(); /* log(-#) = NaN */
        /* subnormal number, scale x up */
        k -= 54;
        x *= 0x1p54;
        hx = bitwise_cast<uiA0>(x) >> 32;
      } else if (hx >= 0x7ff00000) {
        return x;
      } else if (x == One<A0>())
        return Zero<A0>();
      /* reduce x into [sqrt(2)/2, sqrt(2)] */
      hx += 0x3ff00000 - 0x3fe6a09e;
      k += bitwise_cast<iA0>(hx>>20) - 0x3ff;
      hx = (hx&0x000fffff) + 0x3fe6a09e;
      x = bitwise_cast<A0>( (uint64_t)hx<<32 | (bitwise_and(0xffffffffull, bitwise_cast<uiA0>(x))));
      A0 f = dec(x);
      A0 hfsq = Half<A0>()*sqr(f);
      A0 s = f/(2.0f + f);
      A0 z = sqr(s);
      A0 w = sqr(z);
      A0 t1 = w*fma(w, fma(w, Lg6, Lg4), Lg2); //w*(Lg2+w*(Lg4+w*Lg6));
      A0 t2 = z*fma(w, fma(w, fma(w, Lg7, Lg5), Lg3), Lg1); //      z*(Lg1+w*(Lg3+w*(Lg5+w*Lg7)));
      A0 R = t2 + t1;
      //idem log fin/////////////////////////////////////////////////////////////////////////////////////

      A0  hi = f - hfsq;
      hi =  bitwise_and(hi, (Allbits<uiA0>() << 32));

      A0 lo = f - hi - hfsq + s*(hfsq+R);

      /* val_hi+val_lo ~ log10(1+f) + k*log10(2) */
      A0 val_hi = hi*ivln10hi;
      A0 dk = k;
//      A0 y = dk*log10_2hi;
      A0 val_lo = fma(dk, log10_2lo,  fma(lo+hi, ivln10lo,  lo*ivln10hi));

//       /*
//        * Extra precision in for adding y is not strictly needed
//        * since there is no very large cancellation near x = sqrt(2) or
//        * x = 1/sqrt(2), but we do it anyway since it costs little on CPUs
//        * with some parallelism and it reduces the error for many args.
//        */
//       A0 w1 = y + val_hi;
//       val_lo += (y - w1) + val_hi;
//       val_hi = w1;

      return val_lo + val_hi;
    }
  };

} } }


#endif

