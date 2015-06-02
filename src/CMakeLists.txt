if(PacBioBAM_build_chef)
    set(boost_RootDir "${PacBioBAM_RootDir}/_output/deplinks/boost")
    set(boost_IncludeDir "${boost_RootDir}")
    set(zlib_RootDir "${PacBioBAM_RootDir}/_output/deplinks/zlib")
    set(zlib_IncludeDir "${zlib_RootDir}/include")
else()
    find_package(ZLIB REQUIRED)
    find_package(Boost REQUIRED)
endif()
find_package(Threads)

if(WIN32)
  # Need winsock on windows
  set(SOCKET_LIBRARIES "ws2_32")

  # Limit the number of DLLs we will have to bundle
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libstdc++")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libstdc++")
endif()

	
# define PacBioBAM library dependencies
if(PacBioBAM_build_chef)
    set(PacBioBAM_DependencyIncludes 
        ${htslib_IncludeDir}
        ${boost_IncludeDir}
        ${zlib_IncludeDir}
    )
    set(PacBioBAM_DependencyLibraries 
        ${htslib_LibDir}/libhts.a
        ${CMAKE_THREAD_LIBS_INIT}
   )
else()
    if(${building_shared})
        set(PacBioBAM_DependencyIncludes
            ${Htslib_IncludeDir}
            ${ZLIB_INCLUDE_DIRS}
            ${Boost_INCLUDE_DIRS}
        )
        set(PacBioBAM_DependencyLibraries
            ${Htslib_LibDir}/libhts${CMAKE_SHARED_LIBRARY_SUFFIX}
            ${CMAKE_THREAD_LIBS_INIT}
            ${ZLIB_LIBRARIES}
            ${SOCKET_LIBRARIES}
        )
    else()
        set(PacBioBAM_DependencyIncludes
            ${Htslib_IncludeDir}
            ${ZLIB_INCLUDE_DIRS}
            ${Boost_INCLUDE_DIRS}
        )
        set(PacBioBAM_DependencyLibraries
            ${Htslib_LibDir}/libhts.a
            ${CMAKE_THREAD_LIBS_INIT}
            ${ZLIB_LIBRARIES}
            ${SOCKET_LIBRARIES}
        )
    endif() # shared
endif() # chef


# set up library include dirs
include_directories(
	${PacBioBAM_DependencyIncludes}
	${PacBioBAM_IncludeDir}
)

# grab library source files
include(files.cmake)
set(SOURCES
    ${PacBioBAM_H}
    ${PacBioBAM_CPP}
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioBAM_CXX_FLAGS}")

if(PacBioBAM_build_chef)
else()
    # ensure we build self-contained htslib
    add_custom_target(
        htslib
        "make" ZLIB_INC=${ZLIB_INCLUDE_DIR} ZLIB_DIR=${ZLIB_LIBRARIES}
        COMMENT "building htslib ZLIB_INC=${ZLIB_INCLUDE_DIR} ZLIB_DIR=${ZLIB_LIBRARIES}"
        WORKING_DIRECTORY ${Htslib_RootDir}
        VERBATIM
    )

endif()

# define actual library
add_definitions(-DPBBAM_LIBRARY)
if(PacBioBAM_wrap_r)
    # SWIG R does not support PBBAM_SHARED_PTR, but it does support boost::shared_ptr
    # So force boost if we're wrapping for R.
    add_definitions(-DPBBAM_USE_BOOST_SHARED_PTR)
endif()

if(${building_shared})
    add_library(pbbam SHARED ${SOURCES})
    set_target_properties(pbbam PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
        LIBRARY_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
    )
else()
    add_library(pbbam STATIC ${SOURCES})
    set_target_properties(pbbam PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
        LIBRARY_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
    )
endif()

if(PacBioBAM_build_chef)
else()
    add_dependencies(pbbam htslib)
endif()
		
# link dependency libs to pbbam	
target_link_libraries(pbbam ${PacBioBAM_DependencyLibraries})

# define symbols for projects that use PacBioBAM
set(PacBioBAM_INCLUDE_DIRS
    ${PacBioBAM_IncludeDir} ${PacBioBAM_DependencyIncludes}
    CACHE INTERNAL
    "${PROJECT_NAME}: Include Directories"
    FORCE
)

if(${building_shared})
    set(PacBioBAM_LIBRARIES
        ${PacBioBAM_LibDir}/libpbbam${CMAKE_SHARED_LIBRARY_SUFFIX} ${PacBioBAM_DependencyLibraries}
        CACHE INTERNAL
        "${PROJECT_NAME}: Libraries"
        FORCE
    )
else()
    set(PacBioBAM_LIBRARIES
        ${PacBioBAM_LibDir}/libpbbam.a ${PacBioBAM_DependencyLibraries}
        CACHE INTERNAL
        "${PROJECT_NAME}: Libraries"
        FORCE
    )
endif()

# add SWIG directory
add_subdirectory(swig)
