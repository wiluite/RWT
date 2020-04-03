cmake_minimum_required(VERSION 3.10)
project(RWT VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 17)

set (BOOST_TEST_COMPONENTS unit_test_framework)
find_package (Boost REQUIRED COMPONENTS ${BOOST_TEST_COMPONENTS} )
find_package (Threads REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
if (UNIX)
    add_definitions( -DBOOST_ALL_DYN_LINK )
endif()
add_definitions( -fopenmp)

include(external/external.cmake)

add_executable(sample main.cpp)
target_link_libraries (sample -lgomp Threads::Threads)

add_executable(test_app test.cpp)
target_link_libraries (test_app ${Boost_LIBRARIES} -lgomp  Threads::Threads)
add_test (test_app test_app)
