########################################################################
# CMake build script for PacBioBAM library.
########################################################################

cmake_policy(SET CMP0048 NEW)  # lets us set version in project()
project(PacBioBAM VERSION 0.7.0 LANGUAGES CXX C)
cmake_minimum_required(VERSION 3.0)

# project name & version
set(PacBioBAM_NAME pbbam)
set(PacBioBAM_VERSION
  "${PacBioBAM_VERSION_MAJOR}.${PacBioBAM_VERSION_MINOR}.${PacBioBAM_VERSION_PATCH}"
)

# list build-time options
option(PacBioBAM_build_docs    "Build PacBioBAM's API documentation."                   ON)
option(PacBioBAM_build_tests   "Build PacBioBAM's unit tests."                          ON)
option(PacBioBAM_build_shared  "Build PacBioBAM as shared library as well."             OFF)
option(PacBioBAM_build_tools   "Build PacBioBAM command line utilities (e.g. pbindex)"  ON)
option(PacBioBAM_wrap_csharp   "Build PacBioBAM with SWIG bindings for C#."             OFF)
option(PacBioBAM_wrap_python   "Build PacBioBAM with SWIG bindings for Python."         OFF)
option(PacBioBAM_wrap_r        "Build PacBioBAM with SWIG bindings for R."              OFF)
option(PacBioBAM_use_modbuild  "Build PacBioBAM using Modular Build System."            OFF)
option(PacBioBAM_use_ccache    "Build PacBioBAM using ccache, if available."            ON)
option(PacBioBAM_auto_validate "Build PacBioBAM with auto-validation enabled."          OFF)

# enable ccache, if available 
if(PacBioBAM_use_ccache)
    find_program(CCACHE_FOUND ccache)
    if(CCACHE_FOUND)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
        set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK    ccache)
    endif()
endif()

# Deprecating the "PacBioBAM_build_pbindex" command line option in favor of more
# general "PacBioBAM_build_tools", as we're starting to add new utilities.
#
# That said, I don't want to break current auto tests/builds, so I'm providing a
# warning message so devs are aware.
#
if(DEFINED PacBioBAM_build_pbindex)

    # construct warning message
    set(pbindex_warning "\nDeprecated:\n-DPacBioBAM_build_pbindex\n")
    if (PacBioBAM_build_pbindex)
        set(pbindex_warning "${pbindex_warning} Building as requested,")
    else()
        set(pbindex_warning "${pbindex_warning} Skipping as requested,")
    endif()
    set(pbindex_warning "${pbindex_warning} but support for this option will be removed at some point in the future.\n")
    message(AUTHOR_WARNING "${pbindex_warning} ** Use -DPacBioBAM_build_tools instead. **\n")

    # force PacBioBAM_build_tools option
    set(PacBioBAM_build_tools ${PacBioBAM_build_pbindex} CACHE BOOL
        "Build PacBioBAM with add'l utilities (e.g. pbindex, pbindexdump)." FORCE)
endif()

# enable testing if requested
if(PacBioBAM_build_tests)
    enable_testing()
endif()

# determine if we're generating SWIG bindings
if(PacBioBAM_wrap_csharp OR PacBioBAM_wrap_r OR PacBioBAM_wrap_python)
    set(wrapping_swig TRUE)
else()
    set(wrapping_swig FALSE)
endif()

# determine if we need a shared lib
if(PacBioBAM_build_shared OR wrapping_swig)
    set(BUILD_SHARED_LIBS ON)
    set(htslib_build_shared ON CACHE BOOL "force htslibConfig to export proper library name")
    set(PB_LIB_MODE SHARED)
    set(PB_LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
    set(BUILD_SHARED_LIBS OFF)
    set(PB_LIB_MODE STATIC)
    set(PB_LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()
    
# main project paths
set(PacBioBAM_RootDir       ${PacBioBAM_SOURCE_DIR})
set(PacBioBAM_DocsDir       ${PacBioBAM_RootDir}/docs)
set(PacBioBAM_IncludeDir    ${PacBioBAM_RootDir}/include)
set(PacBioBAM_SourceDir     ${PacBioBAM_RootDir}/src)
set(PacBioBAM_SwigSourceDir ${PacBioBAM_RootDir}/src/swig)
set(PacBioBAM_TestsDir      ${PacBioBAM_RootDir}/tests)
set(PacBioBAM_ToolsDir      ${PacBioBAM_RootDir}/tools)

if(NOT PacBioBAM_OutputDir)
    set(PacBioBAM_OutputDir ${PacBioBAM_RootDir})
else()
    # if SWIG bindings requested
    if(${wrapping_swig})
        message(FATAL_ERROR "SWIG bindings not currently supported in modular build.")
    endif()
endif()

set(PacBioBAM_BinDir        ${PacBioBAM_OutputDir}/bin)
set(PacBioBAM_LibDir        ${PacBioBAM_OutputDir}/lib)
set(PacBioBAM_ThirdPartyDir ${PacBioBAM_RootDir}/third-party)

file(MAKE_DIRECTORY ${PacBioBAM_BinDir})
file(MAKE_DIRECTORY ${PacBioBAM_LibDir})

# use some custom Find*, Use* cmake modules
list(APPEND CMAKE_MODULE_PATH "${PacBioBAM_RootDir}/cmake")

# shared & third-party paths
if(NOT HTSLIB_INCLUDE_DIRS OR
   NOT HTSLIB_LIBRARIES)
    if(HTSLIB_ROOTDIR)
        find_package(htslib
                     PATHS ${HTSLIB_ROOTDIR}/
                     REQUIRED)
    else()
        find_package(htslib
                     PATHS ${PacBioBAM_SOURCE_DIR}/../htslib/
                     REQUIRED)
    endif()
endif()

if(NOT Boost_INCLUDE_DIRS)
    find_package(Boost REQUIRED)
endif()

if (NOT ZLIB_INCLUDE_DIRS OR
    NOT ZLIB_LIBRARIES)
    find_package(ZLIB REQUIRED)
endif()

# shared CXX flags for src & tests
if (MSVC)
    set(PacBioBAM_CXX_FLAGS "/Wall")
else()
    set(PacBioBAM_CXX_FLAGS "-std=c++11 -Wall -Wno-sign-compare")
endif()

# NOTE: -Wno-unused-local-typedefs used to quash clang warnings w/ Boost
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-unused-local-typedefs" HAS_NO_UNUSED_LOCAL_TYPEDEFS)
if(HAS_NO_UNUSED_LOCAL_TYPEDEFS)
    set(PacBioBAM_CXX_FLAGS "${PacBioBAM_CXX_FLAGS} -Wno-unused-local-typedefs")
endif()

if(PacBioBAM_auto_validate)
    add_definitions("-DPBBAM_AUTOVALIDATE=1")
endif()

# For now, keep @rpath out of install names on OS X, as it causes SWIG
# tests to fail.
if(APPLE)
    set(CMAKE_MACOSX_RPATH OFF)
endif()

# Turn on windows-style filepath resolution.
# We need to add this #define early (not just in the C# SWIG wrapper)
if(WIN32 AND PacBioBAM_wrap_csharp)
    add_definitions(-DPBBAM_WIN_FILEPATHS)
endif()

# keep this order (src first, at least)
add_subdirectory(src)

if(PacBioBAM_build_tools)
    add_subdirectory(tools)
endif()

if(PacBioBAM_build_docs)
    add_subdirectory(docs)
endif()

if(PacBioBAM_build_tests)

    if (NOT GTEST_SRC_DIR)
        set(PREBUILT_GTEST_SRC ${PacBioBAM_RootDir}/../../../../prebuilt.tmpout/gtest/gtest_1.7.0/)
        if(EXISTS ${PREBUILT_GTEST_SRC})
            set(GTEST_SRC_DIR ${PREBUILT_GTEST_SRC})
        else()
            set(GTEST_SRC_DIR ../gtest) # keep old fallback behavior for external builds, for now at least
        endif()
    endif()

    add_subdirectory(${GTEST_SRC_DIR} external/gtest/build)
    add_subdirectory(tests)

endif()

