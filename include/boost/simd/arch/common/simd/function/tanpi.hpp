//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_SIMD_FUNCTION_TANPI_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_SIMD_FUNCTION_TANPI_HPP_INCLUDED

#include <boost/simd/arch/common/detail/simd/trigo.hpp>
#include <boost/simd/function/restricted.hpp>
#include <boost/simd/arch/common/detail/tags.hpp>
#include <boost/simd/detail/dispatch/function/overload.hpp>
#include <boost/config.hpp>

namespace boost { namespace simd { namespace ext
{
  namespace bd = boost::dispatch;
  namespace bs = boost::simd;
  BOOST_DISPATCH_OVERLOAD_IF ( tanpi_
                          , (typename A0, typename X)
                          , (detail::is_native<X>)
                          , bd::cpu_
                          , bs::pack_< bd::floating_<A0>, X>
                          )
  {
    BOOST_FORCEINLINE A0 operator() ( A0 const& a0) const BOOST_NOEXCEPT
    {
      return detail::trig_base<A0,tag::pi_tag,tag::simd_type,tag::big_tag>::tana(a0);
    }
  };
  BOOST_DISPATCH_OVERLOAD( tanpi_
                          , (typename A0, typename A1, typename X)
                          , bd::cpu_
                          , bs::pack_< bd::floating_<A0>, X>
                          , bd::scalar_ < bd::unspecified_<A1>>
                          )
  {
    BOOST_FORCEINLINE A0 operator() ( A0 const& a0,  A1 const &) const BOOST_NOEXCEPT
    {
      return detail::trig_base<A0,tag::pi_tag,tag::simd_type,A1>::tana(a0);
    }
  };
  BOOST_DISPATCH_OVERLOAD_IF ( tanpi_
                          , (typename A0, typename X)
                          , (detail::is_native<X>)
                          , bd::cpu_
                          , bs::restricted_tag
                          , bs::pack_< bd::floating_<A0>, X>
                          )
  {
    BOOST_FORCEINLINE A0 operator() (const restricted_tag &,  A0 const& a0) const BOOST_NOEXCEPT
    {
      return detail::trig_base<A0,tag::pi_tag,tag::simd_type,tag::clipped_pio4_tag>::tana(a0);
    }
  };
} } }


#endif
