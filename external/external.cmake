message(STATUS "Installing Transwarp via submodule")
execute_process(COMMAND git submodule update --init -- external/transwarp WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(RWT ${CMAKE_CURRENT_SOURCE_DIR}/external/transwarp/src)




