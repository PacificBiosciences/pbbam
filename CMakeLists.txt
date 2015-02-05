########################################################################
# CMake build script for PacBioBAM library.
########################################################################

project(PacBioBAM CXX C)
cmake_minimum_required(VERSION 2.8)

# project version
set(PacBioBAM_MAJOR_VERSION 0)
set(PacBioBAM_MINOR_VERSION 0)
set(PacBioBAM_PATCH_VERSION 2)
set(PacBioBAM_VERSION
  "${PacBioBAM_MAJOR_VERSION}.${PacBioBAM_MINOR_VERSION}.${PacBioBAM_PATCH_VERSION}"
)

# build-time options
option(PacBioBAM_build_docs  "Build PacBioBAM's API documentation." ON)
option(PacBioBAM_build_tests "Build PacBioBAM's unit tests." ON)

# main project paths
set(PacBioBAM_RootDir       ${PacBioBAM_SOURCE_DIR})
set(PacBioBAM_BinDir        ${PacBioBAM_RootDir}/bin)
set(PacBioBAM_DocsDir       ${PacBioBAM_RootDir}/docs)
set(PacBioBAM_LibDir        ${PacBioBAM_RootDir}/lib)
set(PacBioBAM_IncludeDir    ${PacBioBAM_RootDir}/include)
set(PacBioBAM_SourceDir     ${PacBioBAM_RootDir}/src)
set(PacBioBAM_TestsDir      ${PacBioBAM_RootDir}/tests)
set(PacBioBAM_ThirdPartyDir ${PacBioBAM_RootDir}/third-party)
file(MAKE_DIRECTORY ${PacBioBAM_LibDir})

# shared & third-party paths
set(Htslib_RootDir    ${PacBioBAM_ThirdPartyDir}/htslib)
set(Htslib_IncludeDir ${Htslib_RootDir})
set(Htslib_LibDir     ${Htslib_RootDir})

# shared CXX flags for src & tests
include(CheckCXXCompilerFlag)
set(PacBioBAM_CXX_FLAGS "-std=c++11 -Wall")
# NOTE: -Wno-unused-local-typedefs used to quash clang warnings w/ Boost
check_cxx_compiler_flag("-Wno-unused-local-typedefs" HAS_NO_UNUSED_LOCAL_TYPEDEFS)
if(HAS_NO_UNUSED_LOCAL_TYPEDEFS)
    set(PacBioBAM_CXX_FLAGS "${PacBioBAM_CXX_FLAGS} -Wno-unused-local-typedefs")
endif()

# keep this order (src first, at least)
add_subdirectory(src)
add_subdirectory(docs)
if(PacBioBAM_build_tests)
    enable_testing()
    add_subdirectory(tests)
endif()
