CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

project (demo)

OPTION(build_test "Build SimSystem test cases." OFF)

SET (CMAKE_CXX_FLAGS "-frounding-math -Wall -std=c++11 -O3")
SET (CMAKE_CXX_FLAGS_DEBUG "-frounding-math -Wall -std=c++11 -g")
SET (CMAKE_CXX_FLAGS_RELEASE "-frounding-math -Wall -std=c++11 -O3")

FIND_PACKAGE(Boost 1.53)
FIND_PACKAGE(CGAL 1.43)

ADD_SUBDIRECTORY (src)

if (build_test)
  enable_testing()
  FIND_PACKAGE(GTest REQUIRED)
  ADD_SUBDIRECTORY(test)
endif()
