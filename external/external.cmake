#set(dependency_via_submodule ON)
message(STATUS "Installing Transwarp via submodule")
execute_process(COMMAND git submodule update --init -- external/transwarp WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(wiluite_RWT_transwarp SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/transwarp/src)
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/transwarp EXCLUDE_FROM_ALL)




