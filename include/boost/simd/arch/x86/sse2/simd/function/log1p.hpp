//==================================================================================================
/*!
    @file

    @Copyright 2016 Numscale SAS

    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_X86_SSE2_SIMD_FUNCTION_LOG1P_HPP_INCLUDED
#define BOOST_SIMD_ARCH_X86_SSE2_SIMD_FUNCTION_LOG1P_HPP_INCLUDED
#include <boost/simd/detail/overload.hpp>
#include <boost/config.hpp>
#include <boost/simd/function/dec.hpp>
#include <boost/simd/function/inc.hpp>
#include <boost/simd/function/if_else.hpp>
#include <boost/simd/function/if_plus.hpp>
#include <boost/simd/function/is_nez.hpp>
#include <boost/simd/constant/inf.hpp>


namespace boost { namespace simd { namespace ext
{
  namespace bd =  boost::dispatch;
  BOOST_DISPATCH_OVERLOAD ( log1p_
                          , (typename A0)
                          , bs::sse2_
                          , bs::pack_<bd::double_<A0>, bs::sse_>
                         )
  {
    BOOST_FORCEINLINE A0 operator() ( const A0 & a0 ) const BOOST_NOEXCEPT
    {
      //This version is a little less accurate, but stay within 0.5 ulps
      // and is a lot faster in sse2
      A0 u = inc(a0);
      A0 r = if_plus(is_nez(u),
                     log(u),
                     (a0-dec(u))/u); // cancels errors with IEEE arithmetic
#ifndef BOOST_SIMD_NO_INFINITIES
      r = if_else(u == Inf<A0>(),u, r);
#endif
      return r;

    }
  };
} } }

#endif
