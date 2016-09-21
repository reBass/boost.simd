//==================================================================================================
/*!
  @file

  Useful version forward declarations and aliases

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)

**/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_DETAIL_VERSION_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_DETAIL_VERSION_HPP_INCLUDED
#include <boost/simd/detail/dispatch/function/overload.hpp>
#include <boost/config.hpp>
#include <boost/core/demangle.hpp>
#include <type_traits>

namespace boost { namespace simd { namespace detail {

  template<typename T> inline std::string type_id()
  {
    typedef std::is_const<typename std::remove_reference<T>::type>  const_t;
    typedef std::is_lvalue_reference<T>                             lref_t;
    typedef std::is_rvalue_reference<T>                             rref_t;
    std::string s = boost::core::demangle(typeid(T).name());
    s += const_t::value ? " const"  : "";
    s += lref_t::value   ? "&"      : "";
    s += rref_t::value   ? "&&"     : "";
    return s;
  }
  template<typename T> inline std::string type_id( const T& )
  {
    return type_id<T>();
  }

  template < typename T, typename Tag>
  BOOST_FORCEINLINE void version(const T &, const Tag&) BOOST_NOEXCEPT
  {
    static bool ici = true;
     if (ici) {
       std::cout << type_id<T>() << "  " << type_id<Tag>() << std::endl;
       ici = false;
     }
  }
#ifdef BOOST_SIMD_TAG_DEBUG
#define VERSION(T, Tag)  detail::version(T(), Tag())
#else
#define VERSION(T, Tag)
#endif
} } }


#endif
