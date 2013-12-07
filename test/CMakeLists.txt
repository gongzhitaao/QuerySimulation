include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(test_simsystem test_simsystem.cpp)
target_link_libraries(test_simsystem gtest gtest_main)
target_link_libraries(test_simsystem pthread)
target_link_libraries(test_simsystem simsystem)

add_test(test-simsystem test_simsystem)
