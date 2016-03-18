
set(PbmergeSrcDir ${PacBioBAM_ToolsDir}/pbmerge/src)

# create version header
set(PbMerge_VERSION ${PacBioBAM_VERSION})
configure_file(
    ${PbmergeSrcDir}/PbMergeVersion.h.in PbMergeVersion.h @ONLY
)

# list source files
set(PBMERGE_SOURCES
    ${ToolsCommonDir}/BamFileMerger.h
    ${ToolsCommonDir}/OptionParser.cpp
    ${PbmergeSrcDir}/main.cpp
)

# build pbmerge executable
include(PbbamTool)
create_pbbam_tool(
    TARGET  pbmerge
    SOURCES ${PBMERGE_SOURCES}
)

# cram tests
if (PacBioBAM_build_tests)
    if(PacBioBAM_auto_validate)
        # skip for now til we clean up merge tests under autovalidate, too
    else()
        add_test(
        NAME pbmerge_CramTests
        WORKING_DIRECTORY ${PacBioBAM_TestsDir}/scripts
        COMMAND "python" cram.py
            ${PacBioBAM_CramTestsDir}/pbmerge_pacbio_ordering.t
            ${PacBioBAM_CramTestsDir}/pbmerge_aligned_ordering.t
            ${PacBioBAM_CramTestsDir}/pbmerge_mixed_ordering.t
            ${PacBioBAM_CramTestsDir}/pbmerge_dataset.t
            ${PacBioBAM_CramTestsDir}/pbmerge_fofn.t
    )
    endif()
endif()
