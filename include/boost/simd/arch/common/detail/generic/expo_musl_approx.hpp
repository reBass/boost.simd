//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_DETAIL_GENERIC_EXPO_MUSL_APPROX_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_DETAIL_GENERIC_EXPO_MUSL_APPROX_HPP_INCLUDED
#include <boost/simd/arch/common/detail/tags.hpp>
#include <boost/simd/function/musl.hpp>
#include <boost/simd/constant/two.hpp>
#include <boost/simd/constant/log_2.hpp>
#include <boost/simd/function/regular.hpp>
#include <boost/simd/function/fma.hpp>
#include <boost/simd/function/fnms.hpp>
#include <boost/simd/function/horn.hpp>
#include <boost/simd/function/horn1.hpp>
#include <boost/simd/function/inc.hpp>
#include <boost/simd/function/oneminus.hpp>
#include <boost/simd/function/sqr.hpp>


namespace boost { namespace simd
{
  namespace detail
  {
    namespace bd =  boost::dispatch;

    ////////////////////////////////////////////////////////////////////////////////////
    // approx and finalize for exp. as musl
    template < class A0> struct exp_approx < A0, bs::tag::exp_, double, musl_tag>
    {
    /* origin: FreeBSD /usr/src/lib/msun/src/e_expf.c */
    /*
     * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
     */
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
      static BOOST_FORCEINLINE A0 approx(A0 const& x) BOOST_NOEXCEPT
      {
        A0 const t = sqr(x);
        return fnms(t,
                    horn<A0
                         , 0x3fc555555555553eull
                         , 0xbf66c16c16bebd93ull
                         , 0x3f11566aaf25de2cull
                         , 0xbebbbd41c5d26bf1ull
                         , 0x3e66376972bea4d0ull
                    >(t), x); //x-h*t
      }

      static BOOST_FORCEINLINE A0 finalize(A0 const& x, A0 const& c, A0 const& hi, A0 const& lo) BOOST_NOEXCEPT
      {
        return oneminus(((lo-(x*c)/(Two<A0>()-c))-hi));
      }
    };

    template < class A0> struct exp_approx < A0, bs::tag::exp_, float, musl_tag>
    {
      static BOOST_FORCEINLINE A0 approx(A0 const& x) BOOST_NOEXCEPT
      {
        const A0  P1(1.6666625440e-1f),  /*  0xaaaa8f.0p-26 */
                  P2(-2.7667332906e-3f); /* -0xb55215.0p-32 */
       A0 xx = sqr(x);
        return fnms(xx, fma(xx, P2, P1), x);
      }

      static BOOST_FORCEINLINE A0 finalize( A0 const& x, A0 const& c
                                          , A0 const& hi, A0 const& lo) BOOST_NOEXCEPT
      {
        return  inc(x*c/(2-c)- lo + hi);
      }
    };

  }
} }


#endif
