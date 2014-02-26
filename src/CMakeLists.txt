add_library(global STATIC global.cpp global.h)
add_library(param STATIC param.cpp param.h)
add_library(walkinggraph STATIC walkinggraph.cpp walkinggraph.h)
add_library(particle STATIC particle.cpp particle.h)
add_library(simulation STATIC simulation.cpp simulation.h)

add_executable(main main.cpp)

target_link_libraries(main simulation)
target_link_libraries(main particle)
target_link_libraries(main walkinggraph)
target_link_libraries(main param)
target_link_libraries(main global)

target_link_libraries(main CGAL CGAL_Core)
target_link_libraries(main gmp)
target_link_libraries(main boost_system)
target_link_libraries(main boost_random)
target_link_libraries(main boost_thread)
