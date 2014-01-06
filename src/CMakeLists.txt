add_library(walkinggraph STATIC walkinggraph.cpp walkinggraph.h)
add_library(particle STATIC particle.cpp particle.h defs.h)

add_executable(main main.cpp)
target_link_libraries(main walkinggraph)
target_link_libraries(main particle)
target_link_libraries(main CGAL CGAL_Core)
target_link_libraries(main gmp)
