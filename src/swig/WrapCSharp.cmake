
find_package(CSharp REQUIRED)
include (${CSHARP_USE_FILE})

set(PacBioBAM_CSharpLibDir ${PacBioBAM_LibDir}/csharp)
set(PacBioBAM_CSharpDLL    ${PacBioBAM_LibDir}/csharp/bin/Debug/PacBio.BAM.dll)
set(CSharpTestRootDir      ${PacBioBAM_TestsDir}/src/csharp)
set(NativeLibraryPaths     ${PacBioBAM_CSharpLibDir}:${PacBioBAM_LibDir}:${Htslib_LibDir})

# create wrapper
file(MAKE_DIRECTORY ${PacBioBAM_CSharpLibDir})
set(CMAKE_SWIG_OUTDIR ${PacBioBAM_CSharpLibDir})   # ensure any swig files in lib/csharp
set_source_files_properties(
  PacBioBam.i PROPERTIES
  CPLUSPLUS ON
  SWIG_FLAGS "-namespace;PacBio.BAM")
swig_add_module(PacBioBam csharp PacBioBam.i)
swig_link_libraries(PacBioBam ${PacBioBAM_LIBRARIES}) # add any C# libs you need from <Find|Use>CSharp.cmake
set_target_properties(
    ${SWIG_MODULE_PacBioBam_REAL_NAME}            # ensure wrapper lib in lib/csharp
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PacBioBAM_CSharpLibDir}
)
add_dependencies(${SWIG_MODULE_PacBioBam_REAL_NAME} pbbam)

# create .csproj for IDEs.  Not strictly necessary to do it this way
# but prevents the problem of people trying to open the csproj before
# the .cs files have been generated.
configure_file(
  ${PacBioBAM_SwigSourceDir}/PacBio.BAM.csproj.in
  ${PacBioBAM_CSharpLibDir}/PacBio.BAM.csproj
)

#
# "wrapper worked" check
#
configure_file(
  ${CSharpTestRootDir}/check_swig.sh.in
  ${CSharpTestRootDir}/check_swig.sh
)
add_custom_target(
    check_swig_csharp
    ALL
    ./check_swig.sh
    COMMENT "Checking C# wrapper"
    WORKING_DIRECTORY ${CSharpTestRootDir}
)
add_dependencies(check_swig_csharp ${SWIG_MODULE_PacBioBam_REAL_NAME})

#
# Unit tests
#
if (PacBioBAM_build_tests)

  configure_file(
    ${CSharpTestRootDir}/TestPbbam.cs.in
    ${CSharpTestRootDir}/TestPbbam.cs)

  # Not clear how to steer these files where we want them.
  csharp_add_executable(
    TestPbbam
    "${PacBioBAM_CSharpDLL}"
    "${CSharpTestRootDir}/TestPbbam.cs;")
  add_dependencies(TestPbbam check_swig_csharp)

  add_test(
    NAME CSharpUnitTests
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/src/swig
    COMMAND ${CSHARP_INTERPRETER} TestPbbam.exe)

  set_tests_properties(
    CSharpUnitTests
    PROPERTIES
    ENVIRONMENT "MONO_PATH=${PacBioBAM_CSharpLibDir}/bin/Debug;LD_LIBRARY_PATH=${NativeLibraryPaths};DYLD_LIBRARY_PATH=${NativeLibraryPaths}")

endif() # unit tests
