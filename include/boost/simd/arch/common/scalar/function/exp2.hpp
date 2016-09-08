//==================================================================================================
/**
  Copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
**/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_SCALAR_FUNCTION_EXP2_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_SCALAR_FUNCTION_EXP2_HPP_INCLUDED

#include <boost/simd/function/std.hpp>
#include <boost/simd/detail/dispatch/function/overload.hpp>
#include <boost/simd/arch/common/detail/scalar/exponential.hpp>
#include <boost/simd/arch/common/detail/tags.hpp>
#include <boost/config.hpp>
#include <cmath>

namespace boost { namespace simd { namespace ext
{
  namespace bd = boost::dispatch;
  namespace bs = boost::simd;

  BOOST_DISPATCH_OVERLOAD ( exp2_
                          , (typename A0)
                          , bd::cpu_
                          , bd::scalar_< bd::floating_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() (A0 a0) const BOOST_NOEXCEPT
    {
      return detail::exponential<A0,bs::tag::exp2_,tag::not_simd_type>::expa(a0);
    }
  };

  BOOST_DISPATCH_OVERLOAD ( exp2_
                          , (typename A0)
                          , bd::cpu_
                          , bs::std_tag
                          , bd::scalar_< bd::floating_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() (const std_tag &, A0 a0) const BOOST_NOEXCEPT
    {
      return std::exp2(a0);
    }
  };

  BOOST_DISPATCH_OVERLOAD ( exp2_
                          , (typename A0)
                          , bd::cpu_
                          , bs::musl_tag
                          , bd::scalar_< bd::floating_<A0> >
                          )
  {
    BOOST_FORCEINLINE A0 operator() (const musl_tag &, A0 a0) const BOOST_NOEXCEPT
    {
      return std::exp2(a0);
//       constexpr int TBLSIZE = 16;
//       const float
//         redux = 0x1.8p23f / TBLSIZE,
//         P1    = 0x1.62e430p-1f,
//         P2    = 0x1.ebfbe0p-3f,
//         P3    = 0x1.c6b348p-5f,
//         P4    = 0x1.3b2c9cp-7f;

//       const double exp2ft[TBLSIZE] = {
//         0x1.6a09e667f3bcdp-1,
//         0x1.7a11473eb0187p-1,
//         0x1.8ace5422aa0dbp-1,
//         0x1.9c49182a3f090p-1,
//         0x1.ae89f995ad3adp-1,
//         0x1.c199bdd85529cp-1,
//         0x1.d5818dcfba487p-1,
//         0x1.ea4afa2a490dap-1,
//         0x1.0000000000000p+0,
//         0x1.0b5586cf9890fp+0,
//         0x1.172b83c7d517bp+0,
//         0x1.2387a6e756238p+0,
//         0x1.306fe0a31b715p+0,
//         0x1.3dea64c123422p+0,
//         0x1.4bfdad5362a27p+0,
//         0x1.5ab07dd485429p+0,
//       };
//       return a0;
     }
  };

} } }

#endif
