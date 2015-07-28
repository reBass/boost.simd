//==================================================================================================
/*!
  @file

  Convenience header for Boost.SIMD

  @copyright 2012 - 2015 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)

**/
//==================================================================================================
#ifndef BOOST_SIMD_HPP_INCLUDED
#define BOOST_SIMD_HPP_INCLUDED

#include <boost/simd/arch/tags.hpp>
#include <boost/simd/math.hpp>
#include <boost/simd/config.hpp>

namespace boost { namespace simd
{
  /*!
    @defgroup group-config Configuration
    Configuration options for the library
  **/

  namespace tag
  {
    /*!
      @defgroup group-hierarchy Hierarchies
      Type hierarchies defined by the library
    **/
  }

  namespace concept
  {
    /*!
      @defgroup group-concept Concepts
      Concepts defined by the library
    **/
  }

} }

#endif