
# setup
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})
set(PacBioBAM_PythonLibDir ${PacBioBAM_LibDir}/python)
set(PythonTestRootDir ${PacBioBAM_TestsDir}/src/python)

# create wrapper
file(MAKE_DIRECTORY ${PacBioBAM_PythonLibDir})
set(CMAKE_SWIG_OUTDIR ${PacBioBAM_PythonLibDir})  # put PacBioBam.py in lib/python

swig_add_module(PacBioBam python PacBioBam.i)
swig_link_libraries(PacBioBam ${PacBioBAM_LIBRARIES} ${PYTHON_LIBRARIES})
set_target_properties(
    ${SWIG_MODULE_PacBioBam_REAL_NAME}            # put _PacBioBam.so in lib/python
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PacBioBAM_PythonLibDir}
)
add_dependencies(${SWIG_MODULE_PacBioBam_REAL_NAME} pbbam)

# simple "wrapper worked" check
# this is run every build, to check importing from Python, but does NOT run full Python-side unit tests
add_custom_target(
    check_swig_python
    ALL
    "PYTHONPATH=${PacBioBAM_PythonLibDir}" python check_swig.py
    COMMENT "Checking Python wrapper"
    WORKING_DIRECTORY ${PythonTestRootDir}
)
add_dependencies(check_swig_python ${SWIG_MODULE_PacBioBam_REAL_NAME})

# unit tests
if(PacBioBAM_build_tests)

    # configure data directory info
    configure_file(
        ${PythonTestRootDir}/test/config.py.in
        ${PythonTestRootDir}/test/config.py
    )

    # test runner
    add_test(
        NAME PythonUnitTests
        WORKING_DIRECTORY ${PythonTestRootDir}
        COMMAND "python" test_pbbam.py
    )
    set_tests_properties(
        PythonUnitTests
        PROPERTIES
        ENVIRONMENT "PYTHONPATH=${PacBioBAM_PythonLibDir}"
    )

endif() # unit tests

