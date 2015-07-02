# --------------------------------------------- @
# SWIG
# --------------------------------------------- @

# general SWIG
if(${wrapping_swig})

    find_package(SWIG 3.0.5 REQUIRED)

    include(${SWIG_USE_FILE})
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})


    check_cxx_compiler_flag("-Wno-unused-local-typedefs" HAS_NO_UNUSED_BUT_SET_VARIABLE)
    if(HAS_NO_UNUSED_BUT_SET_VARIABLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-but-set-variable")
    endif()

    check_cxx_compiler_flag("-Wno-dynamic-class-memaccess" HAS_NO_DYNAMIC_CLASS_MEMACCESS)
    if (HAS_NO_DYNAMIC_CLASS_MEMACCESS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-dynamic-class-memaccess")
    endif()

    set_source_files_properties(PacBioBam.i PROPERTIES CPLUSPLUS ON)

    if (APPLE)
    else()
        if(CMAKE_SIZEOF_VOID_P EQUAL 8)
            set_source_files_properties(PacBioBam.i PROPERTIES SWIG_FLAGS "-DSWIGWORDSIZE64")
        endif()
    endif()
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
