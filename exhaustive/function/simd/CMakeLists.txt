## -------------------------------------------------------------------------------------------------
##                              Copyright 2016 - NumScale SAS
##
##                   Distributed under the Boost Software License, Version 1.0.
##                        See accompanying file LICENSE.txt or copy at
##                            http://www.boost.org/LICENSE_1_0.txt
## -------------------------------------------------------------------------------------------------

set(SOURCES
  acos.cpp
  acos.accurate.cpp
  acosh.cpp
  asin.cpp
  atan.cpp
  atanh.cpp
  erf.cpp
  exp.cpp
  exp10.cpp
  exp2.cpp
  expm1.cpp
  fast_acos.cpp
  log.cpp
  log.musl.cpp
  log.regular.cpp
  log1p.cpp
  log1p.musl.cpp
  log10.cpp
  log10.musl.cpp
  log2.cpp
  log2.musl.cpp
  nthroot.cpp
 )

make_exhaustive("function.simd" ${SOURCES})

















