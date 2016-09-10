//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_SCALAR_FUNCTION_LOG2_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_SCALAR_FUNCTION_LOG2_HPP_INCLUDED
#include <boost/simd/function/std.hpp>

#include <boost/simd/arch/common/detail/scalar/logarithm.hpp>
#include <boost/simd/function/ilog2.hpp>
#include <boost/simd/meta/is_not_scalar.hpp>
#include <boost/simd/function/musl.hpp>
#include <boost/simd/function/std.hpp>
#include <boost/simd/detail/dispatch/function/overload.hpp>
#include <boost/config.hpp>
#include <cmath>

namespace boost { namespace simd { namespace ext
{
  namespace bd = boost::dispatch;
  namespace bs = boost::simd;

  BOOST_DISPATCH_OVERLOAD ( log2_
                          , (typename A0)
                          , bd::cpu_
                          , bd::scalar_< bd::floating_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() (A0 const& a0) const BOOST_NOEXCEPT
    {
      return detail::logarithm<A0,is_not_scalar_t<A0>>::log2(a0);
    }
  };
  BOOST_DISPATCH_OVERLOAD ( log2_
                          , (typename A0)
                          , bd::cpu_
                          , bd::scalar_< bd::arithmetic_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() (A0 const& a0) const BOOST_NOEXCEPT
    {
      return bs::ilog2(a0);
    }
  };

  BOOST_DISPATCH_OVERLOAD ( log2_
                          , (typename A0)
                          , bd::cpu_
                          , boost::simd::std_tag
                          , bd::scalar_< bd::arithmetic_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() (const std_tag &, A0 a0
                                    ) const BOOST_NOEXCEPT
    {
       return std::log2(a0);
    }
  };

  BOOST_DISPATCH_OVERLOAD ( log2_
                          , (typename A0)
                          , bd::cpu_
                          , bs::musl_tag
                          , bd::scalar_< bd::single_<A0> >
                          )
  {
  /* origin: FreeBSD /usr/src/lib/msun/src/e_log2f.c */
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
  /*
   * Return the base 2 logarithm of x.
   *
   * Reduce x to 2^k (1+f) and calculate r = log(1+f) - f + f*f/2
   * as in log.c, then combine and scale in extra precision:
   *    log2(x) = (f - f*f/2 + r)/log(2) + k
   */
    BOOST_FORCEINLINE A0 operator() (const musl_tag &, A0 x) const BOOST_NOEXCEPT
    {
      const A0
        ivln2hi =  1.4428710938e+00, /* 0x3fb8b000 */
        ivln2lo = -1.7605285393e-04, /* 0xb9389ad4 */
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

      A0  hi = f - hfsq;
      hi =  bitwise_and(hi, uiA0(0xfffff000ul));
      A0  lo = fma(s, hfsq+R, f - hi - hfsq);
      return fma((lo+hi), ivln2lo, lo*ivln2hi + hi*ivln2hi + k);
    }
  };

   BOOST_DISPATCH_OVERLOAD ( log2_
                          , (typename A0)
                          , bd::cpu_
                          , bs::musl_tag
                          , bd::scalar_< bd::double_<A0> >
                          )
  {
  /* origin: FreeBSD /usr/src/lib/msun/src/e_log2.c */
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
  /*
   * Return the base 2 logarithm of x.
   *
   * Reduce x to 2^k (1+f) and calculate r = log(1+f) - f + f*f/2
   * as in log.c, then combine and scale in extra precision:
   *    log2(x) = (f - f*f/2 + r)/log(2) + k
   */
    BOOST_FORCEINLINE A0 operator() (const musl_tag &, A0 x) const BOOST_NOEXCEPT
    {
      const A0
        ivln2hi = 1.44269504072144627571e+00, /* 0x3ff71547, 0x65200000 */
        ivln2lo = 1.67517131648865118353e-10, /* 0x3de705fc, 0x2eefa200 */
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
      /*
       * f-hfsq must (for args near 1) be evaluated in extra precision
       * to avoid a large cancellation when x is near sqrt(2) or 1/sqrt(2).
       * This is fairly efficient since f-hfsq only depends on f, so can
       * be evaluated in parallel with R.  Not combining hfsq with R also
       * keeps R small (though not as small as a true `lo' term would be),
       * so that extra precision is not needed for terms involving R.
       *
       * Compiler bugs involving extra precision used to break Dekker's
       * theorem for spitting f-hfsq as hi+lo, unless double_t was used
       * or the multi-precision calculations were avoided when double_t
       * has extra precision.  These problems are now automatically
       * avoided as a side effect of the optimization of combining the
       * Dekker splitting step with the clear-low-bits step.
       *
       * y must (for args near sqrt(2) and 1/sqrt(2)) be added in extra
       * precision to avoid a very large cancellation when x is very near
       * these values.  Unlike the above cancellations, this problem is
       * specific to base 2.  It is strange that adding +-1 is so much
       * harder than adding +-ln2 or +-log10_2.
       *
       * This uses Dekker's theorem to normalize y+val_hi, so the
       * compiler bugs are back in some configurations, sigh.  And I
       * don't want to used double_t to avoid them, since that gives a
       * pessimization and the support for avoiding the pessimization
       * is not yet available.
       *
       * The multi-precision calculations for the multiplications are
       * routine.
       */

      /* hi+lo = f - hfsq + s*(hfsq+R) ~ log(1+f) */
      A0  hi = f - hfsq;
      hi =  bitwise_and(hi, (Allbits<uiA0>() << 32));
      A0 lo = f - hi - hfsq + s*(hfsq+R);

      A0 val_hi = hi*ivln2hi;
      A0 val_lo = fma(lo+hi, ivln2lo, lo*ivln2hi);

      A0 y = k;
      A0 w1 = y + val_hi;
      val_lo += (y - w1) + val_hi;
      val_hi = w1;
      return val_lo + val_hi;
    }
  };

} } }


#endif
