
set(Bam2SamSrcDir ${PacBioBAM_ToolsDir}/bam2sam/src)

# create version header
set(Bam2Sam_VERSION ${PacBioBAM_VERSION})
configure_file(
    ${Bam2SamSrcDir}/Bam2SamVersion.h.in ${GeneratedDir}/Bam2SamVersion.h @ONLY
)

# list source files
set(BAM2SAM_SOURCES
    ${ToolsCommonDir}/OptionParser.cpp
    ${Bam2SamSrcDir}/main.cpp
    ${Bam2SamSrcDir}/Bam2Sam.cpp
)

# build bam2sam executable
include(PbbamTool)
create_pbbam_tool(
    TARGET  bam2sam
    SOURCES ${BAM2SAM_SOURCES}
)

# cram tests
if (PacBioBAM_build_tests)

    configure_file(
        ${PacBioBAM_CramTestsDir}/bam2sam.t.in
        ${GeneratedDir}/bam2sam.t
    )

    add_test(
        NAME bam2sam_CramTests
        WORKING_DIRECTORY ${PacBioBAM_TestsDir}/scripts
        COMMAND "python" cram.py
            ${GeneratedDir}/bam2sam.t
    )

endif()
