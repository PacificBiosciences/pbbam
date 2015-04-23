
find_package(CSharp REQUIRED)
include (${CSHARP_USE_FILE})

set(PacBioBAM_CSharpLibDir ${PacBioBAM_LibDir}/csharp)
set(CSharpTestRootDir ${PacBioBAM_TestsDir}/src/csharp)

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


# "wrapper worked" check
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




# unit tests
