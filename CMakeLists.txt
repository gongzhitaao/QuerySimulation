CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

project (demo)

OPTION(build_test "Build SimSystem test cases." ON)

SET (CMAKE_CXX_FLAGS "-frounding-math -Wall -std=c++11 -O2")
SET (CMAKE_CXX_FLAGS_DEBUG "-frounding-math -Wall -std=c++11 -g")
SET (CMAKE_CXX_FLAGS_RELEASE "-frounding-math -Wall -std=c++11 -O2")

FIND_PACKAGE(Boost 1.53)

ADD_SUBDIRECTORY (src)

if (build_test)
  enable_testing()
  FIND_PACKAGE(GTest REQUIRED)
  ADD_SUBDIRECTORY(test)
endif()
