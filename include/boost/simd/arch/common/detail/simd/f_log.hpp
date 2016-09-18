//==================================================================================================
/*!

  Copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_DETAIL_SIMD_F_LOG_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_DETAIL_SIMD_F_LOG_HPP_INCLUDED

#include <boost/simd/function/any.hpp>
#include <boost/simd/function/dec.hpp>
#include <boost/simd/function/fma.hpp>
#include <boost/simd/function/is_eqz.hpp>
#include <boost/simd/function/if_nan_else.hpp>
#include <boost/simd/function/if_else.hpp>
#include <boost/simd/function/if_else_zero.hpp>
#include <boost/simd/function/if_plus.hpp>
#include <boost/simd/function/is_lez.hpp>
#include <boost/simd/function/is_ltz.hpp>
#include <boost/simd/function/horn.hpp>
#include <boost/simd/function/sqr.hpp>
#include <boost/simd/function/tofloat.hpp>
#include <boost/simd/arch/common/detail/tags.hpp>
#include <boost/simd/constant/mhalf.hpp>
#include <boost/simd/constant/minf.hpp>
#include <boost/simd/detail/constant/log_2hi.hpp>
#include <boost/simd/detail/constant/log_2lo.hpp>
#include <boost/simd/constant/log2_em1.hpp>
#include <boost/simd/detail/constant/log10_ehi.hpp>
#include <boost/simd/detail/constant/log10_elo.hpp>
#include <boost/simd/detail/constant/log10_2hi.hpp>
#include <boost/simd/detail/constant/log10_2lo.hpp>
#include <boost/simd/constant/sqrt_2o_2.hpp>
#include <boost/simd/detail/dispatch/meta/as_integer.hpp>
#include <boost/simd/detail/dispatch/meta/scalar_of.hpp>

#ifndef BOOST_SIMD_NO_NANS
#include <boost/simd/function/is_nan.hpp>
#include <boost/simd/function/logical_or.hpp>
#endif
#ifndef BOOST_SIMD_NO_INFINITIES
#include <boost/simd/constant/inf.hpp>
#include <boost/simd/function/is_equal.hpp>
#endif
#ifndef BOOST_SIMD_NO_DENORMALS
#include <boost/simd/function/abs.hpp>
#include <boost/simd/function/is_less.hpp>
#include <boost/simd/constant/smallestposval.hpp>
#include <boost/simd/constant/twotonmb.hpp>
#include <boost/simd/constant/mlogtwo2nmb.hpp>
#include <boost/simd/constant/mlog2two2nmb.hpp>
#include <boost/simd/constant/mlog10two2nmb.hpp>
#endif

  //////////////////////////////////////////////////////////////////////////////
  // how to compute the various logarithms
  //////////////////////////////////////////////////////////////////////////////
  // The method is mainly taken from the cephes library:
  // first reduce the the data
  // a0 is supposed > 0
  // the input a0 is split into a mantissa and an exponent
  // the mantissa m is between sqrt(0.5) and sqrt(2) and the corresponding exponent is e
  // a0 = m*2^e
  // then the log? calculus is split in two parts (? being nothing: natural logarithm,  2: base 2 logarithm,
  //10 base ten logarithm)
  // as log?(a) = log?(2^e)+log?(m)
  // 1) computing log?(m)
  //   first put x = m-1 (so -0.29 <  x < 0.414)
  //   write log(m)   = log(1+x)   = x + x*x/2 + x*x*x*g(x)
  //   write log2(m)  = log2(1+x)  = C2*log(x),
  //     C2 = log(2)  the multiplication have to be taken seriously as C2 is not exact
  //   write log10(m) = log10(1+x) = C10*log(x),
  //     C10= log(10) the multiplication have to be taken seriously as C10 is not exact
  // then g(x) has to be approximated
  // g is ((log(1+x)/x-1)/x-1/2)/x
  // It is not a good idea to approximate directly log(1+x) instead of g,  because this will lead to bad precision around 1.
  //
  // in this approximation one can choose a best approximation rational function given by remez algorithm.
  // there exist a classical solution which is a polynomial p8 one of degree 8 that gives 0.5ulps everywhere
  // this is what is done in the kernel_t::log impl;
  // Now,  it is possible to choose a rational fraction or a polynomial of lesser degree to approximate g
  // providing faster but less accurate logs.
  // 2) computing log?(2^e)
  // see the explanations relative to each case
  // 3) finalize
  // This is simply treating invalid entries
  // 4) For denormal we use the fact that log(x) =  log?(x*y)-log?(y) and that if y is
  // the constant two2nmb if x is denormal x*y and y are not.
  //////////////////////////////////////////////////////////////////////////////

namespace boost { namespace simd
{
  namespace detail
  {
    namespace bd = boost::dispatch;
    template < class A0 >
    struct logarithm< A0, tag::simd_type, float>
    {

      static BOOST_FORCEINLINE void kernel_log(const A0& a0,
                                    A0& dk,
                                    A0& hfsq,
                                    A0& s,
                                    A0& r,
                                    A0& f) BOOST_NOEXCEPT
    {

//       using uiA0 = bd::as_integer_t<A0, unsigned>;
       using iA0 = bd::as_integer_t<A0,   signed>;
//       A0 x =  if_nan_else(is_lez(a0), a0);
//       iA0 k(0);
//       auto isnez = is_nez(a0);
// #ifndef BOOST_SIMD_NO_DENORMALS
//       auto test = is_less(a0, Smallestposval<A0>())&&isnez;
//       if (any(test))
//       {
//         k = if_minus(test, k, iA0(25));
//         x = if_else(test, x*A0(0x1p25f), x);
//       }
// #endif
//       uiA0 ix = bitwise_cast<uiA0>(x);
//       /* reduce x into [sqrt(2)/2, sqrt(2)] */
//       ix += 0x3f800000 - 0x3f3504f3;
//       k += bitwise_cast<iA0>(ix>>23) - 0x7f;
//       ix = (ix&0x007fffff) + 0x3f3504f3;
//       x =  bitwise_cast<A0>(ix);
//      using i_t = bd::as_integer_t<A0, signed>;
      iA0 k;
      A0 x;
      std::tie(x, k) = fast_(frexp)(a0);
      const iA0 x_lt_sqrthf = if_else_zero(is_greater(Sqrt_2o_2<A0>(), x),Mone<iA0>());
      k += x_lt_sqrthf;
      f = dec(x+bitwise_and(x, x_lt_sqrthf));
      dk = tofloat(k);
      s = f/(Two<A0>()+f);
      A0 z = sqr(s);
      r =  horn<A0
        , 0x3f2aaaaa  //  0.66666662693 0xaaaaaa.0p-24
        , 0x3eccce13  //  0.40000972152 0xccce13.0p-25
        , 0x3e91e9ee  //  0.28498786688 0x91e9ee.0p-25
        , 0x3e789e26  //  0.24279078841 0xf89e26.0p-26
        >(z)*z;
      hfsq = (Half<A0>()* sqr(f));
    }

    static BOOST_FORCEINLINE A0 log(A0 a0) BOOST_NOEXCEPT
    {
#ifndef BOOST_SIMD_NO_DENORMALS
      auto denormal = is_less(bs::abs(a0), Smallestposval<A0>());
      A0 t = A0(0);
      if (any(denormal))
      {
        a0 = if_else(denormal, a0*Twotonmb<A0>(), a0);
        t = if_else_zero(denormal, Mlogtwo2nmb<A0>());
      }
#endif
      A0 dk, hfsq, s, R, f;
      kernel_log(a0, dk, hfsq, s, R, f);
      A0 y =  (dk* Log_2hi<A0>())-
        ((hfsq-(s*(hfsq+R)+(dk*Log_2lo<A0>())))-f);
#ifndef BOOST_SIMD_NO_DENORMALS
      y+= t;
#endif
      return finalize(a0, y);
    }

    static BOOST_FORCEINLINE A0 log2(const A0& a0) BOOST_NOEXCEPT
    {
      A0 dk, hfsq, s, R, f;
      kernel_log(a0, dk, hfsq, s, R, f);
      A0 y = -(hfsq-(s*(hfsq+R))-f)*Invlog_2<A0>()+dk;
      return finalize(a0, y);
    }

    static BOOST_FORCEINLINE A0 log10(const A0& a0) BOOST_NOEXCEPT
    {
      A0 dk, hfsq, s, R, f;
      kernel_log(a0, dk, hfsq, s, R, f);
      A0 y = -(hfsq-(s*(hfsq+R))-f)*Invlog_10<A0>()+dk*Log_2olog_10<A0>();
      return finalize(a0, y);
    }
  private:
    static BOOST_FORCEINLINE A0 finalize(const A0& a0, const A0& y) BOOST_NOEXCEPT
    {
    #ifdef BOOST_SIMD_NO_NANS
      auto test =  is_ltz(a0);
    #else
      auto test =  logical_or(is_ltz(a0), is_nan(a0));
    #endif
      A0 y1 = if_nan_else(test, y);
    #ifndef BOOST_SIMD_NO_INFINITIES
      y1 = if_else(is_equal(a0, Inf<A0>()), a0, y1);
    #endif
      return if_else(is_eqz(a0), Minf<A0>(), y1);
    }
  };
  }
} }

#endif





//       using iA0 = bd::as_integer_t<A0, signed>;
//       using sA0 = bd::scalar_of_t<A0>;
//       static BOOST_FORCEINLINE void kernel_log(const A0& a0,
//                                                A0& fe,
//                                                A0& x,
//                                                A0& x2,
//                                                A0& y) BOOST_NOEXCEPT
//       {
//         iA0 e;
//         std::tie(x, e)= fast_(frexp)(a0);
//         auto xltsqrthf = (x < Sqrt_2o_2<A0>());
//         fe = if_plus(xltsqrthf, tofloat(e), Mone<A0>());
//         x =  dec(if_plus(xltsqrthf, x, x));
//         x2 = sqr(x);
//         // performances informations using this kernel for boost::simd::log
//         // exhaustive and bench tests with g++-4.7 sse4.2 or scalar give:
//         // at most 0.5 ulp  for input in [0, 3.40282e+38]
//         // 2130706656 values computed.
//         // 2127648316 values (99.86%)  within 0.0 ULPs
//         //    3058340 values (0.14%)   within 0.5 ULPs
//         // bench produces  8.9 cycles/value (simd) 34.5 cycles/value (scalar) full computation
//         // bench produces  7.1 cycles/value (simd) 32.2 cycles/value (scalar) with NO_DENORMALS, NO_INVALIDS etc.
//         y =  horn< A0,
//           0x3eaaaaa9, //      3.3333328e-01
//           0xbe800064, //     -2.5000298e-01
//           0x3e4cd0a3, //      2.0001464e-01
//           0xbe2a6aa0, //     -1.6642237e-01
//           0x3e116e80, //      1.4202309e-01
//           0xbe04d6b7, //     -1.2972532e-01
//           0x3e0229f9, //      1.2711324e-01
//           0xbda5dff0  //     -8.0993533e-02
//           >(x)*x*x2;
//       }

//       static BOOST_FORCEINLINE A0 log(const A0& a0) BOOST_NOEXCEPT
//       {
//         A0 z = a0;
// #ifndef BOOST_SIMD_NO_DENORMALS
//         A0 t = Zero<A0>();
//         auto denormal = is_less(bs::abs(z), Smallestposval<A0>());
//         z = if_else(denormal, z*Twotonmb<A0>(), z);
//         t = if_else_zero(denormal, Mlogtwo2nmb<A0>());
// #endif
//         //log(2.0) in double is 6.931471805599453e-01
//         //double(0.693359375f)+double(-0.00021219444f)  is  6.931471805600000e-01 at 1.0e-14 of log(2.0)
//         // let us call Log_2hi 0.693359375f anf Log_2lo -0.00021219444f
//         // We use thi to correct the sum where this could matter a lot
//         // log(a0) = fe*Log_2hi+ (0.5f*x*x +(fe*Log_2lo+y))
//         // These operations are order dependent: the parentheses do matter
//         A0 x, fe, x2, y;
//         kernel_log(z, fe, x, x2, y);
//         y = bs::fma(fe, Log_2lo<A0>(), y);
//         y = bs::fma(Mhalf<A0>(), x2, y);
// #ifdef BOOST_SIMD_NO_DENORMALS
//         return finalize(a0, bs::fma(Log_2hi<A0>(), fe, x+y));
// #else
//         return finalize(a0, bs::fma(Log_2hi<A0>(), fe, x+y+t));
// #endif
//       }

//       static BOOST_FORCEINLINE A0 log2(const A0& a0) BOOST_NOEXCEPT
//       {
//         A0 z =  a0;
// #ifndef BOOST_SIMD_NO_DENORMALS
//         auto denormal = is_less(bs::abs(z), Smallestposval<A0>());
//         z = if_else(denormal, z*Twotonmb<A0>(), z);
//         A0 t = if_else_zero(denormal, Mlog2two2nmb<A0>());
// #endif
//         //here let l2em1 = log2(e)-1, the computation is done as:
//         //log2(a0) = ((l2em1*x+(l2em1*(y+x*x/2)))+(y+x*x/2)))+x+fe for best results
//         // once again the order is very important.
//         A0 x, fe, x2, y;
//         kernel_log(z, fe, x, x2, y);
//         y = bs::fma(Mhalf<A0>(),x2, y);
//         z = bs::fma(x,Log2_em1<A0>(),y*Log2_em1<A0>());
// #ifdef BOOST_SIMD_NO_DENORMALS
//         return finalize(a0, ((z+y)+x)+fe);
// #else
//         return finalize(a0, ((z+y)+x)+fe+t);
// #endif
//       }

//       static BOOST_FORCEINLINE A0 log10(const A0& a0) BOOST_NOEXCEPT
//       {
//         A0 z = a0;
// #ifndef BOOST_SIMD_NO_DENORMALS
//         auto denormal = is_less(bs::abs(z), Smallestposval<A0>());
//         z = if_else(denormal, z*Twotonmb<A0>(), z);
//         A0 t = if_else_zero(denormal, Mlog10two2nmb<A0>());
// #endif
//         // here there are two multiplication: log of fraction by log10(e) and base 2 exponent by log10(2)
//         // and we have to split log10(e) and log10(2) in two parts to get extra precision when needed
//         A0 x, fe, x2, y;
//         kernel_log(z, fe, x, x2, y);
//         y = bs::fma(x2, Mhalf<A0>(), y);
//         z = (x+y)* Log10_elo<A0>();
//         z = bs::fma(y, Log10_ehi<A0>(), z);
//         z = bs::fma( x, Log10_ehi<A0>(), z);
//         z = bs::fma(fe, Log10_2hi<A0>(), z);
// #ifdef BOOST_SIMD_NO_DENORMALS
//         return finalize(a0, bs::fma(fe, Log10_2lo<A0>(), z));
// #else
//         return finalize(a0, bs::fma(fe, Log10_2lo<A0>(), z+t));
// #endif
//       }
//     private:
//       static BOOST_FORCEINLINE A0 finalize(const A0& a0, const A0& y) BOOST_NOEXCEPT
//       {
// #ifdef BOOST_SIMD_NO_NANS
//         auto test = bs::is_ltz(a0);
// #else
//         auto test = bs::logical_or(bs::is_ltz(a0), bs::is_nan(a0));
// #endif
//         A0 y1 = bs::if_nan_else(test, y);
// #ifndef BOOST_SIMD_NO_INFINITIES
//         y1 = if_else(bs::is_equal(a0, bs::Inf<A0>()), a0, y1);
// #endif
//         return if_else(is_eqz(a0), bs::Minf<A0>(), y1);
//       }
//     };
//   }
// } }

// #endif

