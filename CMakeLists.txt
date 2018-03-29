########################################################################
# CMake build script for PacBioBAM library.
########################################################################

cmake_policy(SET CMP0048 NEW)  # lets us set version in project()
project(PacBioBAM VERSION 0.17.0 LANGUAGES CXX C)
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
option(PacBioBAM_use_modbuild  "Build PacBioBAM using Modular Build System."            OFF)
option(PacBioBAM_use_ccache    "Build PacBioBAM using ccache, if available."            ON)
option(PacBioBAM_auto_validate "Build PacBioBAM with auto-validation enabled."          OFF)

if(PacBioBAM_build_tests)
    enable_testing()
endif()

set(PacBioBAM_permissive_cigar  OFF)
STRING(TOLOWER "${PBBAM_PERMISSIVE_CIGAR}" PBBAM_PERMISSIVE_CIGAR_LOWER)
if (DEFINED PBBAM_PERMISSIVE_CIGAR)
    set(PacBioBAM_permissive_cigar  ON)
    add_definitions(-DPBBAM_PERMISSIVE_CIGAR)
endif()

# project paths
set(PacBioBAM_RootDir       ${CMAKE_CURRENT_LIST_DIR})
set(PacBioBAM_DocsDir       ${PacBioBAM_RootDir}/docs)
set(PacBioBAM_IncludeDir    ${PacBioBAM_RootDir}/include)
set(PacBioBAM_SourceDir     ${PacBioBAM_RootDir}/src)
set(PacBioBAM_TestsDir      ${PacBioBAM_RootDir}/tests)
set(PacBioBAM_ThirdPartyDir ${PacBioBAM_RootDir}/third-party)
set(PacBioBAM_ToolsDir      ${PacBioBAM_RootDir}/tools)

if(NOT PacBioBAM_OutputDir)
    set(PacBioBAM_OutputDir ${CMAKE_CURRENT_BINARY_DIR})
endif()
set(PacBioBAM_BinDir    ${PacBioBAM_OutputDir}/bin)
set(PacBioBAM_LibDir    ${PacBioBAM_OutputDir}/lib)

set(GeneratedDir ${CMAKE_BINARY_DIR}/generated)
set(GeneratedTestDataDir ${GeneratedDir}/data)
file(MAKE_DIRECTORY ${PacBioBAM_BinDir})
file(MAKE_DIRECTORY ${PacBioBAM_LibDir})
file(MAKE_DIRECTORY ${GeneratedDir})
file(MAKE_DIRECTORY ${GeneratedTestDataDir})

# project configuration (keep this order)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})
include(pbbam-ccache)
include(pbbam-compilerflags)
include(pbbam-libtype)
include(pbbam-dependencies)

# project components (keep this order)
add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(docs)
add_subdirectory(tests)
