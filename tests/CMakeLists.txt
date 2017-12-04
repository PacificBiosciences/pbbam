
if(PacBioBAM_build_tests)

    # setup GoogleTest
    if (NOT GTEST_SRC_DIR)
        set(PREBUILT_GTEST_SRC ${PacBioBAM_RootDir}/../../../../prebuilt.tmpout/gtest/gtest_1.7.0/)
        if(EXISTS ${PREBUILT_GTEST_SRC})
            set(GTEST_SRC_DIR ${PREBUILT_GTEST_SRC})
        else()
            set(GTEST_SRC_DIR ${PacBioBAM_RootDir}/../gtest) # keep old fallback behavior for external builds, for now at least
        endif()
    endif()
    add_subdirectory(${GTEST_SRC_DIR} external/gtest/build)

    # generate paths/values used by for unit tests
    configure_file(
        ${PacBioBAM_TestsDir}/src/PbbamTestData.h.in
        ${CMAKE_BINARY_DIR}/generated/PbbamTestData.h
    )
    configure_file(
        ${PacBioBAM_TestsDir}/data/group/group.fofn.in
        ${CMAKE_BINARY_DIR}/generated/group.fofn
    )
    configure_file(
        ${PacBioBAM_TestsDir}/data/pbbamify/synthetic_movie_all.subreadset.xml.in
        ${CMAKE_BINARY_DIR}/generated/synthetic_movie_all.subreadset.xml
    )

    # grab PacBioBAM unit test source files
    include(files.cmake)
    set(SOURCES
        ${PacBioBAMTest_H}
        ${PacBioBAMTest_CPP}
    )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioBAM_CXX_FLAGS}")

    # define unit test executable
    if(MSVC)
        # VS2012+ pooh-pooh's Derek's "#define private public" trick
        add_definitions(-D_ALLOW_KEYWORD_MACROS)
    endif()

    add_executable(test_pbbam ${SOURCES})
    set_target_properties(test_pbbam PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PacBioBAM_BinDir}
    )
    target_include_directories(test_pbbam
        PUBLIC
        ${CMAKE_BINARY_DIR}/generated
        ${PacBioBAM_INCLUDE_DIRS}
        ${gtest_SOURCE_DIR}/include
        ${gtest_SOURCE_DIR}
    )

    # generate test data
    add_custom_target(
        generate_test_data
        WORKING_DIRECTORY ${PacBioBAM_TestsDir}/scripts
        COMMAND "python" generate_data.py
            ${PacBioBAM_TestsDir}/data/
            ${GeneratedTestDataDir}
    )

    # add unit tests to test framework
    add_test(
        NAME UnitTests
        WORKING_DIRECTORY ${PacBioBAM_BinDir}
        COMMAND test_pbbam
    )
    add_dependencies(test_pbbam generate_test_data)
    target_link_libraries(test_pbbam
        pbbam
        ${CMAKE_THREAD_LIBS_INIT} # quirky pthreads
        gtest
        gtest_main
    )
endif() # PacBioBAM_build_tests
