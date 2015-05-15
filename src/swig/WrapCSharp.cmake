
find_package(CSharp REQUIRED)
include (${CSHARP_USE_FILE})

set(PacBioBAM_CSharpLibDir  ${PacBioBAM_LibDir}/csharp/PacBio.BAM)
set(PacBioBAM_CSharpDLL     ${PacBioBAM_CSharpLibDir}/bin/Debug/PacBio.BAM.dll)
set(CSharpTestRootDir       ${PacBioBAM_TestsDir}/src/csharp)
set(NativeLibraryPaths      ${PacBioBAM_CSharpLibDir}:${PacBioBAM_LibDir}:${Htslib_LibDir})

#
# Create SWIG wrapper
#
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

#
# Write a csproj, then shell out to build and check the assembly---
# can't get it working nicely in CMake yet
#
configure_file(
   ${PacBioBAM_SwigSourceDir}/PacBio.BAM.csproj.in
   ${PacBioBAM_CSharpLibDir}/PacBio.BAM.csproj)
configure_file(
  ${CSharpTestRootDir}/TestPbbam.cs.in
  ${CSharpTestRootDir}/TestPbbam.cs)
configure_file(
   ${CSharpTestRootDir}/buildAssembly.sh.in
   buildAssembly.sh)
add_custom_command(
  OUTPUT ${PacBioBAM_CSharpDLL}
  DEPENDS ${SWIG_MODULE_PacBioBam_REAL_NAME}
  COMMAND ./buildAssembly.sh
)
add_custom_target(CSharpAssembly ALL DEPENDS ${PacBioBAM_CSharpDLL})
