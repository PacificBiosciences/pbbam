
set(ToolsCommonDir ${PacBioBAM_ToolsDir}/common)
set(PacBioBAM_CramTestsDir ${PacBioBAM_TestsDir}/src/cram)

# quash warning with OptionParser
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-unused-private-field" HAS_NO_UNUSED_PRIVATE_FIELD)
if(HAS_NO_UNUSED_PRIVATE_FIELD)
    set(PacBioBAM_CXX_FLAGS "${PacBioBAM_CXX_FLAGS} -Wno-unused-private-field")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioBAM_CXX_FLAGS}")

# tools
add_subdirectory(bam2sam)
add_subdirectory(pbindex)
add_subdirectory(pbindexdump)
add_subdirectory(pbmerge)
