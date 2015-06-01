########################################################################
# CMake build script for PacBioBAM library.
########################################################################

project(PacBioBAM CXX C)
cmake_minimum_required(VERSION 2.8)

# project version
set(PacBioBAM_MAJOR_VERSION 0)
set(PacBioBAM_MINOR_VERSION 0)
set(PacBioBAM_PATCH_VERSION 5)
set(PacBioBAM_VERSION
  "${PacBioBAM_MAJOR_VERSION}.${PacBioBAM_MINOR_VERSION}.${PacBioBAM_PATCH_VERSION}"
)

# list build-time options
option(PacBioBAM_build_docs   "Build PacBioBAM's API documentation."           ON)
option(PacBioBAM_build_tests  "Build PacBioBAM's unit tests."                  ON)
option(PacBioBAM_build_chef   "Build PacBioBAM using chef build."              OFF)
option(PacBioBAM_build_shared "Build PacBioBAM as shared library as well."     OFF)
option(PacBioBAM_wrap_csharp  "Build PacBioBAM with SWIG bindings for C#."     OFF)
option(PacBioBAM_wrap_python  "Build PacBioBAM with SWIG bindings for Python." OFF)
option(PacBioBAM_wrap_r       "Build PacBioBAM with SWIG bindings for R."      OFF)

# --check build-time options --

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
    set(building_shared TRUE)
else()
    set(building_shared FALSE)
endif()

# only allow 'vanilla' pbbam (static lib, no SWIG) if chef build for now
if(PacBioBAM_wrap_chef)
    # if shared lib directly requested
    if(PacBioBAM_build_shared)
        message(FATAL_ERROR "Shared lib not currently supported in chef build.")
    endif()
    # if SWIG bindings requested
    if(${wrapping_swig})
        message(FATAL_ERROR "SWIG bindings not currently supported in chef build.")
    endif()
endif()

# main project paths
set(PacBioBAM_RootDir       ${PacBioBAM_SOURCE_DIR})
set(PacBioBAM_BinDir        ${PacBioBAM_RootDir}/bin)
set(PacBioBAM_DocsDir       ${PacBioBAM_RootDir}/docs)
set(PacBioBAM_IncludeDir    ${PacBioBAM_RootDir}/include)
set(PacBioBAM_SourceDir     ${PacBioBAM_RootDir}/src)
set(PacBioBAM_SwigSourceDir ${PacBioBAM_RootDir}/src/swig)
set(PacBioBAM_TestsDir      ${PacBioBAM_RootDir}/tests)

if(PacBioBAM_build_chef)
    set(PacBioBAM_LibDir        ${PacBioBAM_RootDir}/_output/lib)
else()
    set(PacBioBAM_LibDir        ${PacBioBAM_RootDir}/lib)
    set(PacBioBAM_ThirdPartyDir ${PacBioBAM_RootDir}/third-party)
endif()

file(MAKE_DIRECTORY ${PacBioBAM_LibDir})

# use some custom Find*, Use* cmake modules
list(APPEND CMAKE_MODULE_PATH "${PacBioBAM_RootDir}/cmake")

# shared & third-party paths
if(PacBioBAM_build_chef)
    set(htslib_RootDir    "${PacBioBAM_RootDir}/_output/deplinks/htslib")
    set(htslib_IncludeDir "${htslib_RootDir}/include")
    set(htslib_LibDir     "${htslib_RootDir}/lib")
else()
    set(Htslib_RootDir    ${PacBioBAM_ThirdPartyDir}/htslib)
    set(Htslib_IncludeDir ${Htslib_RootDir})
    set(Htslib_LibDir     ${Htslib_RootDir})
endif()

# shared CXX flags for src & tests
if (MSVC)
    set(PacBioBAM_CXX_FLAGS "/Wall")
else()
    set(PacBioBAM_CXX_FLAGS "-std=c++11 -Wall")
endif()

# NOTE: -Wno-unused-local-typedefs used to quash clang warnings w/ Boost
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-unused-local-typedefs" HAS_NO_UNUSED_LOCAL_TYPEDEFS)
if(HAS_NO_UNUSED_LOCAL_TYPEDEFS)
    set(PacBioBAM_CXX_FLAGS "${PacBioBAM_CXX_FLAGS} -Wno-unused-local-typedefs")
endif()

# keep this order (src first, at least)
add_subdirectory(src)
if(PacBioBAM_build_docs)
    add_subdirectory(docs)
endif()
if(PacBioBAM_build_tests)
    add_subdirectory(tests)
endif()
