cmake_minimum_required(VERSION 3.8)
project(RWT VERSION 0.0.0)

include(external/external.cmake)

target_link_libraries(RWT INTERFACE wiluite_RWT_transwarp)

