## -------------------------------------------------------------------------------------------------
##                              Copyright 2016 - NumScale SAS
##
##                   Distributed under the Boost Software License, Version 1.0.
##                        See accompanying file LICENSE.txt or copy at
##                            http://www.boost.org/LICENSE_1_0.txt
## -------------------------------------------------------------------------------------------------

set(SOURCES
  exp.cpp
  exp.musl.cpp
  exp10.cpp
  log.cpp
  log.musl.cpp
    )

make_exhaustive("function.scalar" ${SOURCES})
