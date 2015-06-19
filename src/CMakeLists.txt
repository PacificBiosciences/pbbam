if(PacBioBAM_use_modbuild)
    # boost
    set(boost_RootDir     "${PacBioBAM_OutputDir}/deplinks/boost")
    set(boost_IncludeDir  "${boost_RootDir}")
    # htslib
    set(Htslib_RootDir    "${PacBioBAM_OutputDir}/deplinks/htslib")
    set(Htslib_IncludeDir "${Htslib_RootDir}/include")
    set(Htslib_LibDir     "${Htslib_RootDir}/lib")
    # zlib
    set(zlib_RootDir      "${PacBioBAM_OutputDir}/deplinks/zlib")
    set(zlib_IncludeDir   "${zlib_RootDir}/include")
    set(zlib_LibDir       "${zlib_RootDir}/lib")

    set(Boost_INCLUDE_DIRS "${boost_IncludeDir}")
    set(ZLIB_INCLUDE_DIRS  "${zlib_IncludeDir}")

    if(${building_shared})
        set(Htslib_Libraries "${Htslib_LibDir}/libhts${CMAKE_SHARED_LIBRARY_SUFFIX}")
        set(ZLIB_LIBRARIES   "${zlib_LibDir}/libz${CMAKE_SHARED_LIBRARY_SUFFIX}")
    else()
        set(Htslib_Libraries "${Htslib_LibDir}/libhts.a")
        set(ZLIB_LIBRARIES   "${zlib_LibDir}/libhts.a")
    endif() # shared
else()
    find_package(ZLIB REQUIRED)
    find_package(Boost REQUIRED)

    if(${building_shared})
        set(Htslib_Libraries "${Htslib_LibDir}/libhts${CMAKE_SHARED_LIBRARY_SUFFIX}")
    else()
        set(Htslib_Libraries "${Htslib_LibDir}/libhts.a")
    endif() # shared
endif()  # use_modbuild
find_package(Threads)



if(WIN32)
  # Need winsock on windows
  set(SOCKET_LIBRARIES "ws2_32")

  # Limit the number of DLLs we will have to bundle
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libstdc++")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libstdc++")
endif()

if(PacBioBAM_use_modbuild)
  # Turn off setting of rpath in libs and exes, so we can use LD_LIBRARY_PATH
  # without getting overriden by the rpath setting if the rpath library 
  # happens to exist.
  SET(CMAKE_SKIP_BUILD_RPATH  TRUE)
endif()


# define PacBioBAM library dependencies
set(PacBioBAM_DependencyIncludes
    ${Boost_INCLUDE_DIRS}
    ${Htslib_IncludeDir}
    ${ZLIB_INCLUDE_DIRS}
)
set(PacBioBAM_DependencyLibraries
    ${Htslib_Libraries}
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${SOCKET_LIBRARIES}
)

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

if(PacBioBAM_use_modbuild)
else()
    # ensure we build self-contained htslib
    string(REPLACE "/libz.a" "" ZLIB_LIBRARY_REPLACED ${ZLIB_LIBRARY})
    add_custom_target(
        htslib
        "make" ZLIB_INC=${ZLIB_INCLUDE_DIR} ZLIB_DIR=${ZLIB_LIBRARY_REPLACED}
        COMMENT "building htslib ZLIB_INC=${ZLIB_INCLUDE_DIR} ZLIB_DIR=${ZLIB_LIBRARY_REPLACED}"
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

if(PacBioBAM_use_modbuild)
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

if(${wrapping_swig})
    # add SWIG directory
    add_subdirectory(swig)
endif()  # swig
