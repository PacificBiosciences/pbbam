
set(PbmergeSrcDir ${PacBioBAM_ToolsDir}/pbmerge/src)

# create version header
set(PbMerge_VERSION ${PacBioBAM_VERSION})
configure_file(
    ${PbmergeSrcDir}/PbMergeVersion.h.in ${GeneratedDir}/PbMergeVersion.h @ONLY
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

        configure_file(
            ${PacBioBAM_CramTestsDir}/pbmerge_pacbio_ordering.t.in
            ${GeneratedDir}/pbmerge_pacbio_ordering.t
            @ONLY
        )
 
        configure_file(
            ${PacBioBAM_CramTestsDir}/pbmerge_aligned_ordering.t.in
            ${GeneratedDir}/pbmerge_aligned_ordering.t
            @ONLY
        )

        configure_file(
            ${PacBioBAM_CramTestsDir}/pbmerge_mixed_ordering.t.in
            ${GeneratedDir}/pbmerge_mixed_ordering.t
            @ONLY
        )

        configure_file(
            ${PacBioBAM_CramTestsDir}/pbmerge_dataset.t.in
            ${GeneratedDir}/pbmerge_dataset.t
            @ONLY
        )

        configure_file(
            ${PacBioBAM_CramTestsDir}/pbmerge_fofn.t.in
            ${GeneratedDir}/pbmerge_fofn.t
            @ONLY
        )

        add_test(
            NAME pbmerge_CramTests
            WORKING_DIRECTORY ${PacBioBAM_TestsDir}/scripts
            COMMAND "python" cram.py
                ${GeneratedDir}/pbmerge_pacbio_ordering.t
                ${GeneratedDir}/pbmerge_aligned_ordering.t
                ${GeneratedDir}/pbmerge_mixed_ordering.t
                ${GeneratedDir}/pbmerge_dataset.t
                ${GeneratedDir}/pbmerge_fofn.t
        )

    endif()
endif()
