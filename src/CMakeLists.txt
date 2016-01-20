find_package(Threads)

if(WIN32)
  # Need winsock on windows
  set(SOCKET_LIBRARIES "ws2_32")

  # Limit the number of DLLs we will have to bundle
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libstdc++")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libstdc++")
endif()

# define PacBioBAM library dependencies
set(PacBioBAM_DependencyIncludes
    ${Boost_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
)
set(PacBioBAM_DependencyLibraries
    ${HTSLIB_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${SOCKET_LIBRARIES}
)

# set up library include dirs
include_directories(SYSTEM
    ${PacBioBAM_DependencyIncludes}
)

include_directories(
    ${PacBioBAM_IncludeDir}
)

# grab library source files
include(files.cmake)
set(SOURCES
    ${PacBioBAM_H}
    ${PacBioBAM_CPP}
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioBAM_CXX_FLAGS}")

# define actual library
add_definitions(-DPBBAM_LIBRARY)
if(PacBioBAM_wrap_r)
    # SWIG R does not support PBBAM_SHARED_PTR, but it does support boost::shared_ptr
    # So force boost if we're wrapping for R.
    add_definitions(-DPBBAM_USE_BOOST_SHARED_PTR)
endif()

add_library(pbbam ${SOURCES})
set_target_properties(pbbam PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
    RUNTIME_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
    LIBRARY_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
)
        
# link dependency libs to pbbam 
target_link_libraries(pbbam ${PacBioBAM_DependencyLibraries})
target_include_directories(pbbam 
    PUBLIC
    ${PacBioBAM_IncludeDir} 
    ${PacBioBAM_DependencyIncludes}
)

# define symbols for projects that use PacBioBAM
set(PacBioBAM_INCLUDE_DIRS
    ${PacBioBAM_IncludeDir} ${PacBioBAM_DependencyIncludes}
    CACHE INTERNAL
    "${PROJECT_NAME}: Include Directories"
    FORCE
)

set(PacBioBAM_LIBRARIES
    ${PacBioBAM_LibDir}/libpbbam${PB_LIB_SUFFIX} ${PacBioBAM_DependencyLibraries}
    CACHE INTERNAL
    "${PROJECT_NAME}: Libraries"
    FORCE
)

if(${wrapping_swig})
    # add SWIG directory
    add_subdirectory(swig)
endif()  # swig
