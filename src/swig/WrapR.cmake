# setup
set(R_INCLUDE_DIR_HINT /mnt/software/r/R/3.1.1/usr/share/R/include) # TODO: hard-coded hint for now, clean up later
find_package(R REQUIRED)
include_directories(${R_INCLUDE_DIR})
set(PacBioBAM_RLibDir ${PacBioBAM_LibDir}/R)
set(RTestRootDir ${PacBioBAM_TestsDir}/src/R)

# SWIG R does not support PBBAM_SHARED_PTR, but it does support boost::shared_ptr
# So force boost if we're wrapping for R.
add_definitions(-DPBBAM_USE_BOOST_SHARED_PTR)

# create wrapper & library
file(MAKE_DIRECTORY ${PacBioBAM_RLibDir})
set(CMAKE_SWIG_OUTDIR ${PacBioBAM_RLibDir})      # put PacBioBam.R wrapper in lib/R
swig_add_module(PacBioBam r PacBioBam.i)
swig_link_libraries(PacBioBam ${PacBioBAM_LIBRARIES})
if(R_LIBRARIES)
    swig_link_libraries(PacBioBam ${R_LIBRARIES})
endif()

# symlink htslib
add_custom_target(
    htslib_symlink
    ALL
    "${CMAKE_COMMAND}" -E create_symlink ${Htslib_Libraries} ${PacBioBAM_RLibDir}/libhts.1${CMAKE_SHARED_LIBRARY_SUFFIX}
    COMMENT "Symlinking htslib for R"
)

# make sure the library is named "PacBioBam.so" explicitly
# no "lib" prefix... that gets in the way of the name lookups between SWIG/R
# and make sure library ends up in lib/R
set_target_properties(
    ${SWIG_MODULE_PacBioBam_REAL_NAME}
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PacBioBAM_RLibDir}
    RUNTIME_OUTPUT_DIRECTORY ${PacBioBAM_RLibDir}
    SONAME PacBioBam.so
    PREFIX ""
)
add_dependencies(${SWIG_MODULE_PacBioBam_REAL_NAME} pbbam-shared)
add_dependencies(${SWIG_MODULE_PacBioBam_REAL_NAME} htslib_symlink)

# simple "wrapper worked" check
configure_file(
    ${RTestRootDir}/check_swig.R.in
    ${RTestRootDir}/check_swig.R
)

add_custom_target(
    check_swig_R
    ALL
    "R" --slave --no-save < ${RTestRootDir}/check_swig.R
    COMMENT "Checking R wrapper"
    WORKING_DIRECTORY ${PacBioBAM_RLibDir}
)
add_dependencies(check_swig_R ${SWIG_MODULE_PacBioBam_REAL_NAME})

# unit tests
if(PacBioBAM_build_tests)

    # configure script
    configure_file(
        ${RTestRootDir}/test_pbbam.sh.in
        ${RTestRootDir}/test_pbbam.sh
    )

    # test runner
    add_test(
        NAME RUnitTests
        COMMAND "sh" ${RTestRootDir}/test_pbbam.sh
        WORKING_DIRECTORY ${PacBioBAM_RLibDir}
    )
endif()
