add_library(simsystem STATIC simsystem.cpp simsystem.h walkinggraph.h)
target_link_libraries(simsystem CGAL CGAL_Core)
target_link_libraries(simsystem gmp)

# add_executable(main main.cpp)
