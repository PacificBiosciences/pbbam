find_package(Threads REQUIRED)

# use self-contained GoogleTest (going to compile into text exe)
set(GTest_RootDir    ${PacBioBAM_ThirdPartyDir}/gtest-1.7.0)
set(GTest_IncludeDir ${GTest_RootDir}/include)
set(GTest_SourceDir  ${GTest_RootDir}/src)

# find samtools binary for comparing against 'gold standard' files
#set(Samtools_Dir "/Users/derek/development/samtools")
set(Samtools_Dir "${PacBioBAM_RootDir}/../../../../../prebuilt.out/samtools/samtools-0.1.19/ubuntu-1404/bin")
set(Samtools_Bin "${Samtools_Dir}/samtools")

# ensure tests directory exists
file(MAKE_DIRECTORY ${PacBioBAM_TestsDir}/bin)

# generate paths/values used by for test
configure_file(
    ${PacBioBAM_TestsDir}/src/TestData.h.in
    ${PacBioBAM_TestsDir}/src/TestData.h
)
configure_file(
    ${PacBioBAM_TestsDir}/data/test_group_query/group.fofn.in
    ${PacBioBAM_TestsDir}/data/test_group_query/group.fofn
)

# add PacBioBAM & GoogleTest includes
include_directories(
    ${PacBioBAM_INCLUDE_DIRS}
    ${GTest_RootDir}
    ${GTest_IncludeDir}
)

# grab PacBioBAM test source files
include(files.cmake)
set(SOURCES
    ${PacBioBAMTest_H}
    ${PacBioBAMTest_CPP}

    # GoogleTest compiled in to test exe - mostly to ensure same compile flags
    ${GTest_H}
    ${GTest_CPP}
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioBAM_CXX_FLAGS}")

# define test executable
add_definitions(-DPBBAM_TESTING)
if(MSVC)
    # VS2012+ pooh-pooh's Derek's "#define private public" trick
    add_definitions(-D_ALLOW_KEYWORD_MACROS)
endif()

if(PacBioBAM_wrap_r)
    # SWIG R does not support PBBAM_SHARED_PTR, but it does support boost::shared_ptr
    # So force boost if we're wrapping for R.
    add_definitions(-DPBBAM_USE_BOOST_SHARED_PTR)
endif()
add_executable(test_pbbam ${SOURCES})
set_target_properties(test_pbbam PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PacBioBAM_TestsDir}/bin
)

# set up unit test to run our executable
add_test(
    NAME UnitTests
    WORKING_DIRECTORY ${PacBioBAM_TestsDir}/bin
    COMMAND test_pbbam
)
add_dependencies(test_pbbam pbbam)
target_link_libraries(test_pbbam
    ${PacBioBAM_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT} # quirky pthreads 
)
