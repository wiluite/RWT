cmake_minimum_required(VERSION 3.10)
project(RWT VERSION 2.0.0)

if ( UNIX AND (CMAKE_CXX_COMPILER MATCHES "clang\\+\\+" OR CMAKE_CXX_COMPILER MATCHES "g\\+\\+" OR CMAKE_CXX_COMPILER MATCHES "c\\+\\+") )
    message ("OS and Compiler OK")
else()
    message (FATAL_ERROR "You should be in UNIX, have GCC or Clang compiler")
endif()

set(CMAKE_CXX_STANDARD 17)
find_package (Threads REQUIRED)
find_package(OpenMP REQUIRED)

set (BOOST_TEST_COMPONENTS unit_test_framework)
find_package (Boost REQUIRED COMPONENTS ${BOOST_TEST_COMPONENTS} )
find_package (Threads REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
add_definitions(-DBOOST_ALL_DYN_LINK)
add_definitions(-fopenmp)
if (CMAKE_CXX_COMPILER MATCHES "clang\\+\\+")
    add_definitions( -stdlib=libc++ ) #-D_LIBCPP_HAS_PARALLEL_ALGORITHMS
endif()

include(external/external.cmake)

if (CMAKE_CXX_COMPILER MATCHES "clang\\+\\+")
    add_executable(sample main.cpp detail/rwt_omp.h)
    target_link_libraries (sample -lomp Threads::Threads -stdlib=libc++)

    add_executable(test_app test.cpp detail/rwt_omp.h)
    target_link_libraries (test_app ${Boost_LIBRARIES} -lomp Threads::Threads -stdlib=libc++)
    add_test (test_app test_app)

else()
    add_executable(sample main.cpp detail/rwt_omp.h)
    target_link_libraries (sample -lgomp Threads::Threads)

    add_executable(test_app test.cpp detail/rwt_omp.h)
    target_link_libraries (test_app ${Boost_LIBRARIES} -lgomp  Threads::Threads)
    add_test (test_app test_app)
endif()

