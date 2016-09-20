//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_DETAIL_SCALAR_F_LOG_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_DETAIL_SCALAR_F_LOG_HPP_INCLUDED

#ifndef BOOST_SIMD_NO_NANS
#include <boost/simd/function/is_nan.hpp>
#endif
#ifndef BOOST_SIMD_NO_INFINITIES
#include <boost/simd/constant/inf.hpp>
#endif
#ifndef BOOST_SIMD_NO_DENORMALS
#include <boost/simd/constant/mlog10two2nmb.hpp>
#include <boost/simd/constant/mlog2two2nmb.hpp>
#include <boost/simd/constant/mlogtwo2nmb.hpp>
#include <boost/simd/constant/smallestposval.hpp>
#include <boost/simd/constant/twotonmb.hpp>
#include <boost/simd/function/abs.hpp>
#endif
#include <boost/simd/arch/common/detail/tags.hpp>
#include <boost/simd/constant/invlog_2.hpp>
#include <boost/simd/constant/invlog10_2.hpp>
#include <boost/simd/constant/invlog_10.hpp>
#include <boost/simd/constant/log_2olog_10.hpp>
#include <boost/simd/detail/constant/log10_2hi.hpp>
#include <boost/simd/detail/constant/log10_2lo.hpp>
#include <boost/simd/detail/constant/log10_ehi.hpp>
#include <boost/simd/detail/constant/log10_elo.hpp>
#include <boost/simd/constant/log2_em1.hpp>
#include <boost/simd/detail/constant/log_2hi.hpp>
#include <boost/simd/detail/constant/log_2lo.hpp>
#include <boost/simd/constant/sqrt_2o_2.hpp>
#include <boost/simd/constant/mhalf.hpp>
#include <boost/simd/constant/minf.hpp>
#include <boost/simd/constant/nan.hpp>
#include <boost/simd/constant/zero.hpp>

#include <boost/simd/function/dec.hpp>
#include <boost/simd/function/fma.hpp>
#include <boost/simd/function/frexp.hpp>
#include <boost/simd/function/genmask.hpp>
#include <boost/simd/function/horn.hpp>
#include <boost/simd/function/if_plus.hpp>
#include <boost/simd/function/is_eqz.hpp>
#include <boost/simd/function/is_ltz.hpp>
#include <boost/simd/function/sqr.hpp>
#include <boost/simd/function/tofloat.hpp>
#include <boost/simd/detail/dispatch/meta/as_integer.hpp>
#include <boost/simd/detail/dispatch/meta/scalar_of.hpp>
#include <boost/simd/function/regular.hpp>
#include <boost/simd/function/musl.hpp>

namespace boost { namespace simd
{
  namespace detail
  {
    namespace bd = boost::dispatch;
    template < class A0,
               class Style ,
               class Tag = regular_tag,
               class base_A0 = bd::scalar_of_t<A0>
               >
               struct logarithm{};

    template < class A0, class Tag>
    struct logarithm< A0, tag::not_simd_type, Tag, float>
    {
      using iA0 = bd::as_integer_t<A0,   signed>;
      using uiA0 = bd::as_integer_t<A0, unsigned>;

      static BOOST_FORCEINLINE void kernel_log(A0 a0,
                                               iA0  k,
                                               A0& dk,
                                               A0& hfsq,
                                               A0& s,
                                               A0& r,
                                               A0& f) BOOST_NOEXCEPT
      {
        reduce(Tag(), a0, k, dk, f);
        s = f/(2.0f+f);
        A0 z = sqr(s);
        A0 w = sqr(z);
        A0 t1= w*horn<A0, 0x3eccce13, 0x3e789e26>(w);
        A0 t2= z*horn<A0, 0x3f2aaaaa, 0x3e91e9ee>(w);
        r = t2 + t1;
        hfsq = 0.5f*sqr(f);
      }

      //Up to now the regular way is never taken in scalar,  musl is speedier
      static BOOST_FORCEINLINE void reduce(const regular_tag &, const A0& a0, A0& dk, A0& f) BOOST_NOEXCEPT
      // reduce x into [sqrt(2)/2, sqrt(2)]
      {
        iA0 k;
        A0 x;
        std::tie(x, k) = fast_(frexp)(a0);
        A0  x_lt_sqrthf = genmask(is_greater(Sqrt_2o_2<A0>(), x));
        k += bitwise_cast<iA0>(x_lt_sqrthf);
        f = dec(x+bitwise_and(x, x_lt_sqrthf));
        dk = tofloat(k);
      }

      static BOOST_FORCEINLINE void reduce(const musl_tag &, const A0& a0, iA0 k, A0& dk, A0& f) BOOST_NOEXCEPT
      // reduce x into [sqrt(2)/2, sqrt(2)]
      {
        uiA0 ix = bitwise_cast<uiA0>(a0);
        ix += 0x3f800000 - 0x3f3504f3;
        k += bitwise_cast<iA0>(ix>>23) - 0x7f;
        ix = (ix&0x007fffff) + 0x3f3504f3;
        A0 x =  bitwise_cast<A0>(ix);
        f = dec(x);
        dk = tofloat(k);
      }


      static BOOST_FORCEINLINE bool early_return(A0& a0, A0& res, iA0 & k) BOOST_NOEXCEPT
      {

        uiA0 ix = bitwise_cast<uiA0>(a0);
        if (ix < 0x00800000 || ix>>31)
        {
          if (ix<<1 == 0)
          {
            res = Minf<A0>();  /* log(+-0)=-inf */
            return true;
          }
          if (ix>>31)
          {
            res = Nan<A0>(); /* log(-#) = NaN */
            return true;
          }
          /* subnormal number, scale x up */
          k -= 25;
          a0 *= 0x1p25;
          return false;
        }
        else if (ix >= 0x7f800000)
        {
          res =  a0;
          return true;
        }
        else if (ix == 0x3f800000)
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

