//==================================================================================================
/*!

  Copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_DETAIL_SIMD_F_LOG_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_DETAIL_SIMD_F_LOG_HPP_INCLUDED

#include <boost/simd/arch/common/detail/version.hpp>

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
#include <boost/simd/function/genmask.hpp>
#include <boost/simd/function/horn.hpp>
#include <boost/simd/function/sqr.hpp>
#include <boost/simd/function/tofloat.hpp>
#include <boost/simd/arch/common/detail/tags.hpp>
#include <boost/simd/constant/mhalf.hpp>
#include <boost/simd/constant/minf.hpp>
#include <boost/simd/detail/constant/log_2hi.hpp>
#include <boost/simd/detail/constant/log_2lo.hpp>
#include <boost/simd/constant/log_2.hpp>
#include <boost/simd/constant/log2_em1.hpp>
#include <boost/simd/constant/invlog_2.hpp>
#include <boost/simd/constant/invlog10_2.hpp>
#include <boost/simd/constant/invlog_10.hpp>
#include <boost/simd/constant/log_2olog_10.hpp>
#include <boost/simd/detail/constant/log10_ehi.hpp>
#include <boost/simd/detail/constant/log10_elo.hpp>
#include <boost/simd/detail/constant/log10_2hi.hpp>
#include <boost/simd/detail/constant/log10_2lo.hpp>
#include <boost/simd/constant/sqrt_2o_2.hpp>
#include <boost/simd/detail/dispatch/meta/as_integer.hpp>
#include <boost/simd/detail/dispatch/meta/scalar_of.hpp>
#include <boost/core/ignore_unused.hpp>
#include <boost/simd/function/regular.hpp>
#include <boost/simd/function/musl.hpp>

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

namespace boost { namespace simd
{
  namespace detail
  {
    namespace bd = boost::dispatch;
    template < class A0, class Tag >
    struct logarithm< A0, tag::simd_type, Tag, float>
    {
      using iA0 = bd::as_integer_t<A0,   signed>;

      static BOOST_FORCEINLINE void kernel_log(const regular_tag &,
                                               const A0& a0,
                                               A0& dk,
                                               A0& hfsq,
                                               A0& s,
                                               A0& r,
                                               A0& f) BOOST_NOEXCEPT
      {
        // reduce x into [sqrt(2)/2, sqrt(2)]
        iA0 k;
        A0 x;
        std::tie(x, k) = fast_(frexp)(a0);
        A0  x_lt_sqrthf = genmask(is_greater(Sqrt_2o_2<A0>(), x));
        k += bitwise_cast<iA0>(x_lt_sqrthf);
        f = dec(x+bitwise_and(x, x_lt_sqrthf));
        dk = tofloat(k);
        //Compute approximation
        s = f/(Two<A0>()+f);
        A0 z = sqr(s);
        A0 w = sqr(z);
        A0 t1= w*horn<A0, 0x3eccce13, 0x3e789e26>(w);
        A0 t2= z*horn<A0, 0x3f2aaaaa, 0x3e91e9ee>(w);
        r = t2 + t1;
        hfsq = Half<A0>()*sqr(f);
      }

      static BOOST_FORCEINLINE void kernel_log(const musl_tag &,
                                               const A0& a0,
                                               A0& dk,
                                               A0& hfsq,
                                               A0& s,
                                               A0& r,
                                               A0& f) BOOST_NOEXCEPT
      {
        // reduce x into [sqrt(2)/2, sqrt(2)]
        using uiA0 = bd::as_integer_t<A0, unsigned>;
        uiA0 ix = bitwise_cast<uiA0>(a0);
        ix += 0x3f800000 - 0x3f3504f3;
        iA0 k = bitwise_cast<iA0>(ix>>23) - 0x7f;
        ix = (ix&0x007fffff) + 0x3f3504f3;
        A0 x =  bitwise_cast<A0>(ix);
        f = dec(x);
        dk = tofloat(k);
        //Compute approximation
        s = f/(Two<A0>()+f);
        A0 z = sqr(s);
        A0 w = sqr(z);
        A0 t1= w*horn<A0, 0x3eccce13, 0x3e789e26>(w);
        A0 t2= z*horn<A0, 0x3f2aaaaa, 0x3e91e9ee>(w);
        r = t2 + t1;
        hfsq = Half<A0>()*sqr(f);
      }
//        static BOOST_FORCEINLINE void reduce(const regular_tag &, const A0& a0, A0& dk, A0& f) BOOST_NOEXCEPT
//       // reduce x into [sqrt(2)/2, sqrt(2)]
//       {
//         iA0 k;
//         A0 x;
//         std::tie(x, k) = fast_(frexp)(a0);
//         A0  x_lt_sqrthf = genmask(is_greater(Sqrt_2o_2<A0>(), x));
//         k += bitwise_cast<iA0>(x_lt_sqrthf);
//         f = dec(x+bitwise_and(x, x_lt_sqrthf));
//         dk = tofloat(k);
//       }

//       static BOOST_FORCEINLINE void reduce(const musl_tag &, const A0& a0, A0& dk, A0& f) BOOST_NOEXCEPT
//       // reduce x into [sqrt(2)/2, sqrt(2)]
//       {
//         using uiA0 = bd::as_integer_t<A0, unsigned>;
//         uiA0 ix = bitwise_cast<uiA0>(a0);
//         ix += 0x3f800000 - 0x3f3504f3;
//         iA0 k = bitwise_cast<iA0>(ix>>23) - 0x7f;
//         ix = (ix&0x007fffff) + 0x3f3504f3;
//         A0 x =  bitwise_cast<A0>(ix);
//         f = dec(x);
//         dk = tofloat(k);
//       }

#ifndef BOOST_SIMD_NO_DENORMALS
      static BOOST_FORCEINLINE void denormal_prepare(A0& a0,  A0& t, A0 const & fac) BOOST_NOEXCEPT
      {
        auto denormal = is_less(bs::abs(a0), Smallestposval<A0>());
        if (any(denormal))
        {
          a0 = if_else(denormal, a0*Twotonmb<A0>(), a0);
          t = if_else_zero(denormal, fac);
        }
      }

      static BOOST_FORCEINLINE void denormal_finalize(A0& y,const A0& t) BOOST_NOEXCEPT
      {
        y += t;
      }
#endif

      static BOOST_FORCEINLINE A0 log(const A0& a0) BOOST_NOEXCEPT
      {
        VERSION(A0, Tag);
        auto isnez = is_nez(a0);
#ifndef BOOST_SIMD_NO_DENORMALS
        A0 t(0);
        denormal_prepare(a0, t, Mlogtwo2nmb<A0>() );
#endif
        A0 dk, hfsq, s, R, f;
        kernel_log(Tag(), a0, dk, hfsq, s, R, f);
        A0 y =  (dk* Log_2hi<A0>())- ((hfsq-(s*(hfsq+R)+(dk*Log_2lo<A0>())))-f);
#ifndef BOOST_SIMD_NO_DENORMALS
        denormal_finalize(y, t);
#endif
        return finalize(a0, isnez, y);
      }

      static BOOST_FORCEINLINE A0 log2(const A0& a0) BOOST_NOEXCEPT
      {
        auto isnez = is_nez(a0);
#ifndef BOOST_SIMD_NO_DENORMALS
        A0 t(0);
        denormal_prepare(a0, t, Mlog2two2nmb<A0>());
#endif
        A0 dk, hfsq, s, R, f;
        kernel_log(Tag(), a0, dk, hfsq, s, R, f);
        A0 y = -(hfsq-(s*(hfsq+R))-f)*Invlog_2<A0>()+dk;
#ifndef BOOST_SIMD_NO_DENORMALS
        denormal_finalize(y, t);
#endif
        return finalize(a0, isnez, y);
      }

      static BOOST_FORCEINLINE A0 log10(const A0& a0) BOOST_NOEXCEPT
      {
        auto isnez = is_nez(a0);
#ifndef BOOST_SIMD_NO_DENORMALS
        A0 t(0);
        denormal_prepare(a0, t, Mlog10two2nmb<A0>());
#endif
        A0 dk, hfsq, s, R, f;
        kernel_log(Tag(), a0, dk, hfsq, s, R, f);
        A0 y = -(hfsq-(s*(hfsq+R))-f)*Invlog_10<A0>()+dk*Log_2olog_10<A0>();
#ifndef BOOST_SIMD_NO_DENORMALS
        denormal_finalize(y, t);
#endif
        return finalize(a0, isnez, y);
      }

    private:
      template < class LA0>
      static BOOST_FORCEINLINE A0 finalize(const A0 & a0, LA0 const & isnez, const A0& y) BOOST_NOEXCEPT
      {
#ifndef BOOST_SIMD_NO_INFINITIES
        A0 z = if_else(isnez, if_else(a0 == Inf<A0>(), Inf<A0>(), y), Minf<A0>());
#else
        boost::ignore_unused(a0);
        A0 z = if_else(isnez, y, Minf<A0>());
#endif
        return if_nan_else(is_ltz(a0), z);
      }
    };
  }
} }

#endif


