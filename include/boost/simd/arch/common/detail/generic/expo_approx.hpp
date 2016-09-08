//==================================================================================================
/*!
  @file

  @copyright 2016 NumScale SAS

  Distributed under the Boost Software License, Version 1.0.
  (See accompanying file LICENSE.md or copy at http://boost.org/LICENSE_1_0.txt)
*/
//==================================================================================================
#ifndef BOOST_SIMD_ARCH_COMMON_DETAIL_GENERIC_EXPO_APPROX_HPP_INCLUDED
#define BOOST_SIMD_ARCH_COMMON_DETAIL_GENERIC_EXPO_APPROX_HPP_INCLUDED
#include <boost/simd/arch/common/detail/tags.hpp>
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
    ////////////////////////////////////////////////////////////////////////////////////
    // approx and finalize for exp.
    namespace bd =  boost::dispatch;
    template< typename A0
              , typename Tag
              , typename Base_A0 = bd::scalar_of_t<A0>
              , typename Algo = bs::regular_tag
              >
    struct exp_approx;

    template < class A0> struct exp_approx < A0, bs::tag::exp_, double, regular_tag>
    {
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

    template < class A0> struct exp_approx < A0, bs::tag::exp_, float, regular_tag>
    {
      static BOOST_FORCEINLINE A0 approx(A0 const& x) BOOST_NOEXCEPT
      {
        // remez polynomial of degree 4 on [-0.5 0.5] for (exp(x)-1-x)/sqr(x)
        // tested in range: [-88.3763, 88.3763]
        //2214928753 values (98.98%)  within 0.0 ULPs
        //  22831063 values (1.02%)   within 0.5 ULPs
        //  4.89648 cycles/value (SSE4.2 g++-4.8)
        A0 y = horn<A0,
                    0x3f000000, //  5.0000000e-01
                    0x3e2aa9a5, //  1.6666277e-01
                    0x3d2aa957, //  4.1665401e-02
                    0x3c098d8b, //  8.3955629e-03
                    0x3ab778cf  //  1.3997796e-03
                    >(x);
        return inc(fma(y, sqr(x), x));
      }

      static BOOST_FORCEINLINE A0 finalize( A0 const& , A0 const& c
                                          , A0 const& , A0 const& ) BOOST_NOEXCEPT
      {
        return c;
      }
    };

    ////////////////////////////////////////////////////////////////////////////////////
    // approx and finalize for exp2.
    template < class A0 > struct exp_approx < A0, bs::tag::exp2_, double, regular_tag>
    {

      static BOOST_FORCEINLINE A0 approx(A0 const& x) BOOST_NOEXCEPT
      {
        const A0 t =  sqr(x);
        return fnms(t,
                    horn<A0
                         , 0x3fc555555555553eull
                         , 0xbf66c16c16bebd93ull
                         , 0x3f11566aaf25de2cull
                         , 0xbebbbd41c5d26bf1ull
                         , 0x3e66376972bea4d0ull
                    > (t), x); //x-h*t
      }

      static BOOST_FORCEINLINE A0 finalize(A0 const& x, A0 const& c, A0 const&, A0& ) BOOST_NOEXCEPT
      {
        return oneminus(((-(x*c)/(Two<A0>()-c))-x));
      }
    };

    template < class A0 > struct exp_approx < A0, bs::tag::exp2_, float, regular_tag>
    {
      static BOOST_FORCEINLINE A0 approx(A0 const& x) BOOST_NOEXCEPT
      {
        // remez polynom of degree 4 on [-0.5, 0.5] for (exp2(x)-1-x*log(2))/sqr(x)  tested in range: [-127 127]
        // 2247884800 values computed.
        // 2224606419 values (98.96%)  within 0.0 ULPs
        // 23278381 values (1.04%)   within 0.5 ULPs
        // 3.5 cycles/value  (SSE4.2 g++-4.8)
        A0 y = horn<A0,
                    0x3e75fdf1,  //    2.4022652e-01
                    0x3d6356eb,  //    5.5502813e-02
                    0x3c1d9422,  //    9.6178371e-03
                    0x3ab01218,  //    1.3433127e-03
                    0x3922c8c4   //    1.5524315e-04
                    >(x);
        return inc(fma(y, sqr(x), x*Log_2<A0>()));
      }

      static BOOST_FORCEINLINE A0 finalize( A0 const& , A0 const& c
                                          , A0 const& , A0 const& ) BOOST_NOEXCEPT
      {
        return  c;
      }
    };

    ////////////////////////////////////////////////////////////////////////////////////
    // approx and finalize for exp10.
    template < class A0 > struct exp_approx < A0, bs::tag::exp10_, double, regular_tag>
    {
      static BOOST_FORCEINLINE A0 approx(const A0& x) BOOST_NOEXCEPT
      {
        A0 xx = sqr(x);
        A0 px = x*horn<A0,
                       0x40a2b4798e134a01ull,
                       0x40796b7a050349e4ull,
                       0x40277d9474c55934ull,
                       0x3fa4fd75f3062dd4ull
                       > (xx);
        A0 x2 =  px/(horn1<A0,
                          0x40a03f37650df6e2ull,
                          0x4093e05eefd67782ull,
                          0x405545fdce51ca08ull
                     //   0x3ff0000000000000ull
                          > (xx)-px);
        return inc(x2+x2);
      }

      static BOOST_FORCEINLINE A0 finalize(A0 const&, A0 const& c, A0 const&, A0& ) BOOST_NOEXCEPT
      {
        return c;
      }
    };

    template < class A0 > struct exp_approx < A0, bs::tag::exp10_, float, regular_tag>
    {
      static BOOST_FORCEINLINE A0 approx(A0 const& x) BOOST_NOEXCEPT
      {
        // remez polynom of degree 5 on [-0.5, 0.5]*log10(2) for (exp10(x)-1)/x   tested in range: [-37.9, 38.2308]
        //  2217772528 values computed.
        //  2198853506 values (99.15%)  within 0.0 ULPs
        //  18919022 values (0.85%)   within 0.5 ULPs
        //      5.2 cycles/value  (SSE4.2 g++-4.8)
        return inc(horn<A0,
                        0x40135d8e, //    2.3025851e+00
                        0x4029a926, //    2.6509490e+00
                        0x400237da, //    2.0346589e+00
                        0x3f95eb4c, //    1.1712432e+00
                        0x3f0aacef, //    5.4170126e-01
                        0x3e54dff1  //    2.0788552e-01
                   >(x)*x);
      }

      static BOOST_FORCEINLINE A0 finalize( A0 const& , A0 const& c
                                          , A0 const& , A0 const& ) BOOST_NOEXCEPT
      {
        return  c;
      }
    };
  }
} }


#endif
