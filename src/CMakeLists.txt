add_library(defs STATIC defs.cpp defs.h)
add_library(walkinggraph STATIC walkinggraph.cpp walkinggraph.h)
add_library(particle STATIC particle.cpp particle.h)
add_library(utils STATIC utils.cpp utils.h)

add_executable(main main.cpp)

target_link_libraries(main defs)
target_link_libraries(main walkinggraph)
target_link_libraries(main particle)
target_link_libraries(main utils)

target_link_libraries(main CGAL CGAL_Core)
target_link_libraries(main gmp)

target_link_libraries(main boost_random)
target_link_libraries(main boost_system)
