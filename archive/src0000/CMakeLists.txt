cmake_minimum_required(VERSION 3.13)
set(CMAKE_EXPORT_COMPILE_COMMANDS "true")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -O3 -lrt -Werror")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -O3 -lrt")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -Wextra -Wpedantic -O3 -lrt")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -Wextra -Wpedantic -O0 -lrt")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE RelWithDebInfo)
#set(CMAKE_BUILD_TYPE Relise)
project(synthesator VERSION "1.1.0" LANGUAGES C CXX)
#get google test
#include(FetchContent)
#FetchContent_Declare(googletest	SOURCE_DIR "${CMAKE_SOURCE_DIR}/googletest/")
#FetchContent_GetProperties(googletest)
#FetchContent_MakeAvailable(googletest)
#googletest_POPULATED
#googletest_SOURCE_DIR
#googletest_BUILD_DIR
#if(NOT googletest_POPULATED)
#	FetchContent_Populate(googletest)
#	add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BUILD_DIR})
#endif()
#include(FindPackageHandleStandardArgs)
#find_path(<<TBD>>_INCLUDE_DIR NAMES etherfabric/ef_vi.h)
#find_library(<<TBD>>_LIBRARY NAMES ciul1)
#find_package_handle_standard_args(<<TBD>> DEFAULT_MSG
#        <<TBD>>_LIBRARY
#        <<TBD>>_INCLUDE_DIR)
#add_library(synthesator_lib SHARED 
add_library(synthesator_lib STATIC 
        src/stdafx.h
        src/stdafx.cpp
        src/getmatrix.cpp
        src/extractor_functions.cpp
        src/factorizator.h
        src/factorizator.cpp
        src/factorizator_functions.cpp
        src/synthesator_functions.cpp
#        src/inmemory_matrix_product.cpp
        src/multiply.cpp
)
#target_include_directories(synthesator_lib
#        PRIVATE
#        ${<<TBD>>_INCLUDE_DIR}
#	lib/capp
#	inc
#)
add_executable(synthesator
	src/extractor.cpp
)
target_link_libraries(synthesator
	synthesator_lib
#        ${<<TBD>>_LIBRARY}
)
install(TARGETS synthesator RUNTIME DESTINATION bin)
#add_executable(synthesator_test src/main_test.cpp)
#target_link_libraries(synthesator_test
#	synthesator_lib
#        ${<<TBD>>_LIBRARY})
#target_link_libraries(synthesator_test synthesator_lib gtest_main gmock_main)

