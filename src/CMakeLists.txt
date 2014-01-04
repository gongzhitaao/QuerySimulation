add_library(particle STATIC particle.cpp particle.h)
add_library(simsystem STATIC simsystem.cpp simsystem.h walkinggraph.h particle.h)
target_link_libraries(simsystem particle)
target_link_libraries(simsystem CGAL CGAL_Core)
target_link_libraries(simsystem gmp)

add_executable(simsys_main main.cpp)
target_link_libraries(simsys_main simsystem)
