find_package(Threads REQUIRED)

# ensure output directories exists
file(MAKE_DIRECTORY ${PacBioBAM_TestsDir}/bin)

# generate paths/values used by for unit tests
configure_file(
    ${PacBioBAM_TestsDir}/src/TestData.h.in
    ${GeneratedDir}/TestData.h
)
configure_file(
    ${PacBioBAM_TestsDir}/data/group/group.fofn.in
    ${GeneratedDir}/group.fofn
)

# add PacBioBAM & GoogleTest includes
include_directories(
    ${GeneratedDir}
    ${PacBioBAM_INCLUDE_DIRS}
    ${gtest_SOURCE_DIR}/include 
    ${gtest_SOURCE_DIR}
)

# grab PacBioBAM unit test source files
include(files.cmake)
set(SOURCES
    ${PacBioBAMTest_H}
    ${PacBioBAMTest_CPP}
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioBAM_CXX_FLAGS}")

# define unit test executable
add_definitions(-DPBBAM_TESTING)
if(MSVC)
    # VS2012+ pooh-pooh's Derek's "#define private public" trick
    add_definitions(-D_ALLOW_KEYWORD_MACROS)
endif()

if(PacBioBAM_wrap_r)
    # SWIG R does not support std::shared_ptr, but it does support boost::shared_ptr
    # So force boost if we're wrapping for R.
    add_definitions(-DPBBAM_USE_BOOST_SHARED_PTR)
endif()
add_executable(test_pbbam ${SOURCES})
set_target_properties(test_pbbam PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PacBioBAM_TestsDir}/bin
)

# generate test data
add_custom_target(
    generate_test_data
    WORKING_DIRECTORY ${PacBioBAM_TestsDir}/scripts
    COMMAND "python" generate_data.py
        ${PacBioBAM_TestsDir}/data/
        ${GeneratedTestDataDir}
        baz
)

# add unit tests to test framework
add_test(
    NAME UnitTests
    WORKING_DIRECTORY ${PacBioBAM_TestsDir}/bin
    COMMAND test_pbbam
)
add_dependencies(test_pbbam pbbam generate_test_data)
target_link_libraries(test_pbbam
    ${PacBioBAM_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT} # quirky pthreads 
    gtest
    gtest_main
)
