//==================================================================================================
/**
  Copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
**/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_DETAIL_SIMD_EXPO_BASE_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_DETAIL_SIMD_EXPO_BASE_HPP_INCLUDED

#ifndef BOOST_SIMD_NO_INVALIDS
#include <boost/simd/function/if_allbits_else.hpp>
#include <boost/simd/function/is_nan.hpp>
#endif

#include <boost/simd/function/if_else.hpp>
#include <boost/simd/function/if_zero_else.hpp>
#include <boost/simd/constant/inf.hpp>
#include <boost/simd/arch/common/detail/generic/expo_approx.hpp>
#include <boost/simd/arch/common/detail/generic/expo_utilities.hpp>
#include <boost/simd/arch/common/detail/scalar/expo_scale.hpp>
#include <boost/config.hpp>

namespace boost { namespace simd { namespace detail
{
  template<typename A0, typename Tag, typename Elt, typename Algo>
  struct exponential < A0, Tag, tag::simd_type, Elt, Algo>
  {
    using iA0 = bd::as_integer_t<A0>;
    using approx_t =  exp_approx<A0,Tag,Elt,Algo>;
    using util_t  =  exp_utilities<A0, Tag>;
    // compute exp(ax)
    static BOOST_FORCEINLINE A0 expa(const A0& a0)
    {
      A0 hi, lo, x;
      iA0 k = util_t::reduce(a0, hi, lo, x);
      A0 c = approx_t::approx(x);
      c = approx_t::finalize(x, c, hi, lo);
      c = if_zero_else(util_t::isleminlog(a0), scale(c, k));
      c = if_else(util_t::isgemaxlog(a0), bs::Inf<A0>(), c);
#ifdef BOOST_SIMD_NO_INVALIDS
      return c;
#else
      return if_allbits_else(is_nan(a0), c);
#endif
    }
  };

} } }

#endif
