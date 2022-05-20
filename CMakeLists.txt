cmake_minimum_required(VERSION 3.0.0)
project(iditm3d
    VERSION 0.1.0 
    DESCRIPTION "GSK3D Model"
    LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)

include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

#specify C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#include source files from a directory
aux_source_directory(${CMAKE_SOURCE_DIR}/src SOURCE_FILES_DIR)
add_executable(gsk3d ${SOURCE_FILES_DIR})

#set directory where executable is produced
set_target_properties(gsk3d PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")

#header file directories
target_include_directories(gsk3d PUBLIC
    /home/jtu/mpich/include
    /home/jtu/hdf5/include
    PRIVATE ${CMAKE_SOURCE_DIR}/src
)

#specify additional libraries. full path required
find_library(fdlfootLib
    NAMES libfdlfoot.so
    HINTS ${CMAKE_SOURCE_DIR}/lib)

target_link_libraries(gsk3d PUBLIC
    /home/jtu/hdf5/lib/libhdf5.so
    /home/jtu/hdf5/lib/libhdf5_cpp.so
    ${fdlfootLib}
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(gsk3d PRIVATE OpenMP::OpenMP_CXX)
endif()


#set compile options
target_compile_options(gsk3d PRIVATE
    -Wall
    -Wshadow
    -Wextra
    -Wconversion
    -Wsign-conversion
    -Wuninitialized
    -Wunused
    -fopenmp
    #$<$<CONFIG:DEBUG>:-g3>
    #$<$<CONFIG:DEBUG>:-Og>
    #$<$<CONFIG:RELEASE>:-O3>  
)

# Modify compile flags to change optimization level from O3 to O2
string(REGEX REPLACE "([\\/\\-]O)3" "\\12"
  CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
string(REGEX REPLACE "([\\/\\-]O)3" "\\12"
  CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
string(REGEX REPLACE "([\\/\\-]O)3" "\\12"
  CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL}")
string(REGEX REPLACE "([\\/\\-]O)3" "\\12"
  CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

  #if want to install executable to a different directory
install(TARGETS gsk3d DESTINATION ${CMAKE_BINARY_DIR})