
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
set_target_properties(pbbam PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
    RUNTIME_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
    LIBRARY_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
)

if(PacBioBAM_auto_validate)
    target_compile_definitions(pbbam
        PUBLIC "-DPBBAM_AUTOVALIDATE=1"
    )
endif()

# pbbam includes
target_include_directories(pbbam
    PUBLIC
    ${PacBioBAM_IncludeDir}
    ${Boost_INCLUDE_DIRS}
)

target_compile_options(pbbam
    PUBLIC
    ${HTSlib_CFLAGS}
    ${ZLIB_CFLAGS}
)

# set link dependencies
set(pbbam_all_dependency_libs
    ${HTSlib_LDFLAGS}
    ${ZLIB_LDFLAGS}
    ${SOCKET_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)

target_link_libraries(pbbam
    PUBLIC
    ${pbbam_all_dependency_libs}
)

# define include paths for projects that use pbbam
set(PacBioBAM_INCLUDE_DIRS
    ${PacBioBAM_IncludeDir}
    ${Boost_INCLUDE_DIRS}
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
