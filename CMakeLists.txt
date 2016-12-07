########################################################################
# CMake build script for PacBioBAM library.
########################################################################

cmake_policy(SET CMP0048 NEW)  # lets us set version in project()
project(PacBioBAM VERSION 0.7.4 LANGUAGES CXX C)
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

if (PacBioBAM_wrap_csharp OR PacBioBAM_wrap_r OR PacBioBAM_wrap_python)
    set(wrapping_swig TRUE)
else()
    set(wrapping_swig FALSE)
endif()

if(PacBioBAM_build_tests)
    enable_testing()
endif()

# project paths
set(PacBioBAM_RootDir       ${CMAKE_CURRENT_LIST_DIR})
set(PacBioBAM_DocsDir       ${PacBioBAM_RootDir}/docs)
set(PacBioBAM_IncludeDir    ${PacBioBAM_RootDir}/include)
set(PacBioBAM_SourceDir     ${PacBioBAM_RootDir}/src)
set(PacBioBAM_SwigSourceDir ${PacBioBAM_RootDir}/src/swig)
set(PacBioBAM_TestsDir      ${PacBioBAM_RootDir}/tests)
set(PacBioBAM_ThirdPartyDir ${PacBioBAM_RootDir}/third-party)
set(PacBioBAM_ToolsDir      ${PacBioBAM_RootDir}/tools)

if(NOT PacBioBAM_OutputDir)
    set(PacBioBAM_OutputDir ${CMAKE_CURRENT_BINARY_DIR})
else()
    if(${wrapping_swig})
        message(FATAL_ERROR "SWIG bindings not currently supported in modular build.")
    endif()
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
