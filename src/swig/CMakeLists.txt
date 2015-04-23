# --------------------------------------------- @
# SWIG
# --------------------------------------------- @

# general SWIG
if(${wrapping_swig})
    find_package(SWIG 3.0.5 REQUIRED)
    include(${SWIG_USE_FILE})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    set(CMAKE_SWIG_FLAGS "")
    set_source_files_properties(PacBioBam.i PROPERTIES CPLUSPLUS ON)
endif()

# Python
if(PacBioBAM_wrap_python)
    include(WrapPython.cmake)
endif()

# R
if(PacBioBAM_wrap_r)
    include(WrapR.cmake)
endif()

# CSharp
if(PacBioBAM_wrap_csharp)
    include(WrapCSharp.cmake)
endif()
