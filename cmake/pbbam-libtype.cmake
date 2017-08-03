
# determine if we need a shared lib
if(PacBioBAM_build_shared)
    set(BUILD_SHARED_LIBS ON)
    set(htslib_build_shared ON CACHE BOOL "force htslibConfig to export proper library name")
    set(PB_LIB_MODE SHARED)
    set(PB_LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
    set(BUILD_SHARED_LIBS OFF)
    set(PB_LIB_MODE STATIC)
    set(PB_LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

if(WIN32)
    # Limit the number of DLLs we will have to bundle
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libstdc++")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libstdc++")
endif()



