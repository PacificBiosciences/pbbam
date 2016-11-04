########################################################################
# CMake build script for htslib library.
########################################################################

cmake_policy(SET CMP0048 NEW)  # lets us set version in project()
project(htslib VERSION 1.0.0 LANGUAGES CXX C)
cmake_minimum_required(VERSION 3.0)

# project version
set(htslib_VERSION
  "${htslib_VERSION_MAJOR}.${htslib_VERSION_MINOR}.${htslib_VERSION_PATCH}"
)

# build-time options
option(htslib_build_shared "Build htslib as shared library."  OFF)

# determine if we are building a shared lib
if(htslib_build_shared)
    set(PB_LIB_MODE SHARED)
    set(PB_LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
    set(PB_LIB_MODE STATIC)
    set(PB_LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

# main project paths
string(REPLACE "//" "/" HTS_LOCAL_DIR ${CMAKE_CURRENT_LIST_DIR})
set(htslib_RootDir       ${HTS_LOCAL_DIR} CACHE INTERNAL "bla" FORCE)
set(htslib_IncludeDir    ${htslib_RootDir} CACHE INTERNAL "bla" FORCE)
set(htslib_LibDir        ${htslib_RootDir} CACHE INTERNAL "bla" FORCE)

if (NOT ZLIB_INCLUDE_DIRS OR
    NOT ZLIB_LIBRARIES)
    find_package(ZLIB REQUIRED)
endif()

string(REGEX REPLACE "/libz.*" "" ZLIB_LIBRARY_REPLACED ${ZLIB_LIBRARIES})
if(NOT TARGET htslibSrc)
add_custom_target(
    htslibSrc ALL
    "make" CC=${CMAKE_C_COMPILER} AR=${CMAKE_AR} RANLIB=${CMAKE_RANLIB} ZLIB_INC=${ZLIB_INCLUDE_DIRS} ZLIB_DIR=${ZLIB_LIBRARY_REPLACED}
    COMMENT "cd ${htslib_RootDir} && make CC=${CMAKE_C_COMPILER} AR=${CMAKE_AR} RANLIB=${CMAKE_RANLIB} ZLIB_INC=${ZLIB_INCLUDE_DIRS} ZLIB_DIR=${ZLIB_LIBRARY_REPLACED}"
    WORKING_DIRECTORY ${htslib_RootDir}
    VERBATIM
)

add_library(htslib STATIC IMPORTED)
add_dependencies(htslib htslibSrc)
endif()
# target_include_directories(htslib   
#     INTERFACE
#     ${htslib_IncludeDir}
# )



# define symbols for projects that use htslib
set(HTSLIB_INCLUDE_DIRS
    ${htslib_IncludeDir}
    CACHE INTERNAL
    "${PROJECT_NAME}: Include Directories"
    FORCE
)

set(HTSLIB_LIBRARIES
    ${htslib_LibDir}/libhts${PB_LIB_SUFFIX}
    CACHE INTERNAL
    "${PROJECT_NAME}: Libraries"
    FORCE
)

if(APPLE)
    # e.g. libhts.1.dylib
    set(HTSLIB_LIBRARIES_VERSIONED_LINK
        ${htslib_LibDir}/libhts.1${CMAKE_SHARED_LIBRARY_SUFFIX}
        CACHE INTERNAL
        ""
        FORCE
    )
else()
    # e.g. libhts.so.1
    set(HTSLIB_LIBRARIES_VERSIONED_LINK
        ${htslib_LibDir}/libhts${CMAKE_SHARED_LIBRARY_SUFFIX}.1
        CACHE INTERNAL
        ""
        FORCE
    )
endif()
