//==================================================================================================
/**
  Copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
**/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_SIMD_FUNCTION_EXP_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_SIMD_FUNCTION_EXP_HPP_INCLUDED

#include <boost/simd/detail/overload.hpp>
#include <boost/simd/detail/traits.hpp>
#include <boost/simd/arch/common/detail/simd/exponential.hpp>
#include <cmath>

namespace boost { namespace simd { namespace ext
{
  namespace bd = boost::dispatch;
  namespace bs = boost::simd;

  BOOST_DISPATCH_OVERLOAD_IF( exp_
                            , (typename A0, typename X)
                            , (detail::is_native<X>)
                            , bd::cpu_
                            , bs::pack_< bd::floating_<A0>, X>
                            )
  {
    BOOST_FORCEINLINE A0 operator() (A0 const& a0) const BOOST_NOEXCEPT
    {
      return detail::exponential<A0,bs::tag::exp_,tag::simd_type>::expa(a0);
    }
  };

  BOOST_DISPATCH_OVERLOAD_IF ( exp_
                             , (typename A0, typename X)
                             , (detail::is_native<X>)
                             , bd::cpu_
                             , bs::musl_tag
                             , bs::pack_< bd::floating_<A0>, X>
                          )
  {
    BOOST_FORCEINLINE A0 operator() (const musl_tag &, A0 a0) const BOOST_NOEXCEPT
    {
      using sA0 = bd::scalar_of_t<A0>;
      return detail::exponential<A0,bs::tag::exp_,tag::simd_type,sA0,musl_tag>::expa(a0);
    }
  };

} } }

#endif
