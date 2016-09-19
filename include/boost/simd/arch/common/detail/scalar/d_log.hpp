//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_DETAIL_SCALAR_D_LOG_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_DETAIL_SCALAR_D_LOG_HPP_INCLUDED

#ifndef BOOST_SIMD_NO_NANS
#include <boost/simd/function/is_nan.hpp>
#endif
#ifndef BOOST_SIMD_NO_INFINITIES
#include <boost/simd/constant/inf.hpp>
#endif
#include <boost/simd/function/fast.hpp>
#include <boost/simd/function/horn.hpp>
#include <boost/simd/arch/common/detail/tags.hpp>
#include <boost/simd/constant/half.hpp>
#include <boost/simd/constant/invlog_10.hpp>
#include <boost/simd/constant/invlog_2.hpp>
#include <boost/simd/detail/constant/log_2hi.hpp>
#include <boost/simd/detail/constant/log_2lo.hpp>
#include <boost/simd/constant/log_2olog_10.hpp>
#include <boost/simd/constant/minf.hpp>
#include <boost/simd/constant/mone.hpp>
#include <boost/simd/constant/sqrt_2o_2.hpp>
#include <boost/simd/constant/two.hpp>
#include <boost/simd/function/bitwise_and.hpp>
#include <boost/simd/function/divides.hpp>
#include <boost/simd/function/frexp.hpp>
#include <boost/simd/function/genmask.hpp>
#include <boost/simd/function/if_allbits_else.hpp>
#include <boost/simd/function/if_else_zero.hpp>
#include <boost/simd/function/is_equal.hpp>
#include <boost/simd/function/is_eqz.hpp>
#include <boost/simd/function/is_greater.hpp>
#include <boost/simd/function/is_ltz.hpp>
#include <boost/simd/function/logical_or.hpp>
#include <boost/simd/function/minus.hpp>
#include <boost/simd/function/dec.hpp>
#include <boost/simd/function/multiplies.hpp>
#include <boost/simd/function/plus.hpp>
#include <boost/simd/function/if_plus.hpp>
#include <boost/simd/function/sqr.hpp>
#include <boost/simd/function/tofloat.hpp>
#include <boost/simd/function/unary_minus.hpp>
#include <tuple>

namespace boost { namespace simd
{
  namespace detail
  {
    namespace bd = boost::dispatch;

    template < class A0 >
    struct logarithm< A0, tag::not_simd_type, double>
    {
      using uiA0 = bd::as_integer_t<A0, unsigned>;
      using iA0 = bd::as_integer_t<A0,   signed>;

      static BOOST_FORCEINLINE void kernel_log(A0 a0,
                                               iA0 k,
                                               A0& dk,
                                               A0& hfsq,
                                               A0& s,
                                               A0& r,
                                               A0& f) BOOST_NOEXCEPT
      {
        reduce_musl(a0, k, dk, f);
        s = f/(2.0+f);
        A0 z = sqr(s);
        A0 w = sqr(z);
        A0 t1= w*horn<A0,
                      0x3fd999999997fa04ll,
                      0x3fcc71c51d8e78afll,
                      0x3fc39a09d078c69fll
                      > (w);
        A0 t2= z*horn<A0,
                      0x3fe5555555555593ll,
                      0x3fd2492494229359ll,
                      0x3fc7466496cb03dell,
                      0x3fc2f112df3e5244ll
                      > (w);
        r = t2+t1;
        hfsq = 0.5*sqr(f);
      }

//       static BOOST_FORCEINLINE void reduce(const A0& a0, A0& dk, A0& f) BOOST_NOEXCEPT
//       /* reduce x into [sqrt(2)/2, sqrt(2)] */
//       {
//         iA0 k;
//         A0 x;
//         std::tie(x, k) = fast_(frexp)(a0);
//         A0  x_lt_sqrthf = genmask(is_greater(Sqrt_2o_2<A0>(), x));
//         k += bitwise_cast<iA0>(x_lt_sqrthf);
//         f = dec(x+bitwise_and(x, x_lt_sqrthf));
//         dk = tofloat(k);
//        }

      static BOOST_FORCEINLINE void reduce_musl(A0 a0, iA0 k, A0& dk, A0& f) BOOST_NOEXCEPT
      /* reduce x into [sqrt(2)/2, sqrt(2)] */
      {
        using uiA0 = bd::as_integer_t<A0, unsigned>;
        uiA0 hx = bitwise_cast<uiA0>(a0) >> 32;
        hx += 0x3ff00000 - 0x3fe6a09e;
        k += bitwise_cast<iA0>(hx>>20) - 0x3ff;
        hx = (hx&0x000fffff) + 0x3fe6a09e;
        A0 x = bitwise_cast<A0>(hx<<32 | (bitwise_and(0xffffffffull, bitwise_cast<uiA0>(a0))));
        f = dec(x);
        dk = tofloat(k);
       }


      static BOOST_FORCEINLINE bool early_return(A0& a0, A0& res, iA0 & k) BOOST_NOEXCEPT
      {

        uiA0 hx = bitwise_cast<uiA0>(a0) >> 32;
        if (hx < 0x00100000 || hx>>31) {
          if(is_eqz(a0))
          {
            res = Minf<A0>();  /* log(+-0)=-inf */
            return true;
          }
          if (hx>>31)
          {
            res = Nan<A0>(); /* log(-#) = NaN */
            return true;
          }
          /* subnormal number, scale x up */
          k -= 54;
          a0 *= 0x1p54;
          hx = bitwise_cast<uiA0>(a0) >> 32;
          return false;
        }
        else if (hx >= 0x7ff00000)
        {
          res =  a0;
          return true;
        }
        else if (a0 == One<A0>())
        {
          // res = 0.0;
          return true;
        }
        else
          return false;
      }

      static BOOST_FORCEINLINE A0 log(A0 a0) BOOST_NOEXCEPT
      {
        iA0 k = 0;
        A0 res = 0.0;
        if (early_return(a0, res, k)) return res;
        A0 dk, hfsq, s, R, f;
        kernel_log(a0, k, dk, hfsq, s, R, f);
        return  (dk* Log_2hi<A0>())- ((hfsq-(s*(hfsq+R)+(dk*Log_2lo<A0>())))-f);
      }

      static BOOST_FORCEINLINE A0 log2(A0 a0) BOOST_NOEXCEPT
      {
        iA0 k = 0;
        A0 res = 0.0;
        if (early_return(a0, res, k)) return res;
        A0 dk, hfsq, s, R, f;
        kernel_log(a0, dk, hfsq, s, R, f);
        return -(hfsq-(s*(hfsq+R))-f)*Invlog_2<A0>()+dk;
      }

      static BOOST_FORCEINLINE A0 log10(A0 a0) BOOST_NOEXCEPT
      {
        iA0 k = 0;
        A0 res = 0.0;
        if (early_return(a0, res, k)) return res;
        A0 dk, hfsq, s, R, f;
        kernel_log(a0, dk, hfsq, s, R, f);
        return -(hfsq-(s*(hfsq+R))-f)*Invlog_10<A0>()+dk*Log_2olog_10<A0>();
      }
    };

  }
} }

#endif
