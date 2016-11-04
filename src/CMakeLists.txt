
# grab library source files
include(files.cmake)
set(SOURCES
    ${PacBioBAM_H}
    ${PacBioBAM_CPP}
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioBAM_CXX_FLAGS}")

# define actual library
add_library(pbbam ${SOURCES})

# library properties
target_compile_definitions(pbbam
    PRIVATE "-DPBBAM_LIBRARY"
)
set_target_properties(pbbam PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
    RUNTIME_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
    LIBRARY_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
)

if(PacBioBAM_wrap_r)
    # SWIG R does not support std::shared_ptr, but it does support boost::shared_ptr
    # So force boost if we're wrapping for R.
    target_compile_definitions(pbbam
        PUBLIC -DPBBAM_USE_BOOST_SHARED_PTR
    )
endif()

if(PacBioBAM_auto_validate)
    target_compile_definitions(pbbam
        PUBLIC "-DPBBAM_AUTOVALIDATE=1"
    )
endif()

# pbbam includes
target_include_directories(pbbam
    PUBLIC
    ${PacBioBAM_IncludeDir}
    ${HTSLIB_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
)

# set link dependencies
#   if htslib provided externally
if(HTSLIB_LIBRARIES)
    set(pbbam_all_dependency_libs
        ${HTSLIB_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${SOCKET_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
    )
# otherwise, use the "in-project" htslib target
else()
    set(pbbam_all_dependency_libs
        $<TARGET_FILE:hts>
        ${ZLIB_LIBRARIES}
        ${SOCKET_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
    )
endif()

target_link_libraries(pbbam
    PUBLIC
    ${pbbam_all_dependency_libs}
)

# define include paths for projects that use pbbam
set(PacBioBAM_INCLUDE_DIRS
    ${PacBioBAM_IncludeDir}
    ${HTSLIB_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
    CACHE INTERNAL
    "${PROJECT_NAME}: Include Directories"
    FORCE
)
set(PacBioBAM_LIBRARIES
    $<TARGET_FILE:pbbam>
    ${pbbam_all_dependency_libs}
    CACHE INTERNAL
    "${PROJECT_NAME}: Libraries"
    FORCE
)

# add SWIG directory
add_subdirectory(swig)
