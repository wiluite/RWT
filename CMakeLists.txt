cmake_minimum_required(VERSION 3.8)
project(RWT VERSION 0.0.0)

set(CMAKE_CXX_STANDARD 17)

set (BOOST_TEST_COMPONENTS unit_test_framework log)
find_package (Boost REQUIRED COMPONENTS ${BOOST_TEST_COMPONENTS} )
include_directories(${Boost_INCLUDE_DIRS})
if (UNIX)
    add_definitions( -DBOOST_ALL_DYN_LINK )
endif()


include(external/external.cmake)

add_executable(sample main.cpp)
target_link_libraries (sample ${Boost_LIBRARIES}  )
