
include(CheckCXXCompilerFlag)

# C++11 check & enabling
if (CMAKE_VERSION VERSION_LESS "3.1")
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")    # clang
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")        # gcc
    endif()
else() # 3.1+
    set(CMAKE_CXX_STANDARD          14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

# shared CXX flags for src & tests
if (MSVC)
    set(PacBioBAM_CXX_FLAGS "/Wall")
else()
    set(PacBioBAM_CXX_FLAGS "-Wall")
endif()

# NOTE: -Wno-unused-local-typedefs used to quash clang warnings w/ Boost
check_cxx_compiler_flag("-Wno-unused-local-typedefs" HAS_NO_UNUSED_LOCAL_TYPEDEFS)
if(HAS_NO_UNUSED_LOCAL_TYPEDEFS)
    set(PacBioBAM_CXX_FLAGS "${PacBioBAM_CXX_FLAGS} -Wno-unused-local-typedefs")
endif()

check_cxx_compiler_flag("-Wno-sign-compare" HAS_NO_SIGN_COMPARE)
if(HAS_NO_SIGN_COMPARE)
    set(PacBioBAM_CXX_FLAGS "${PacBioBAM_CXX_FLAGS} -Wno-sign-compare")
endif()

# Turn on windows-style filepath resolution.
# We need to add this #define early (not just in the C# SWIG wrapper)
if(WIN32)
    add_definitions(-DPBBAM_WIN_FILEPATHS)
endif()

# For now, keep @rpath out of install names on OS X, as it causes SWIG
# tests to fail.
if(APPLE)
    set(CMAKE_MACOSX_RPATH OFF)
endif()
