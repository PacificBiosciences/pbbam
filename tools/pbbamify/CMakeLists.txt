
set(PbbamifySrcDir ${PacBioBAM_ToolsDir}/pbbamify/src)

# create version header
set(PbBamify_VERSION ${PacBioBAM_VERSION})
configure_file(
    ${PbbamifySrcDir}/PbBamifyVersion.h.in ${GeneratedDir}/PbBamifyVersion.h @ONLY
)

# list source files
set(PBBAMIFY_SOURCES
    ${ToolsCommonDir}/OptionParser.cpp
    ${PbbamifySrcDir}/main.cpp
    ${PbbamifySrcDir}/PbBamify.cpp
    ${PbbamifySrcDir}/QueryLookup.cpp
)

# build pbbamify executable
include(PbbamTool)
create_pbbam_tool(
    TARGET  pbbamify
    SOURCES ${PBBAMIFY_SOURCES}
)

# cram tests
if (PacBioBAM_build_tests AND PacBioBAM_permissive_cigar)

    configure_file(
        ${PacBioBAM_CramTestsDir}/pbbamify.t.in
        ${GeneratedDir}/pbbamify.t
    )

    add_test(
        NAME pbbamify_CramTests
        WORKING_DIRECTORY ${PacBioBAM_TestsDir}/scripts
        COMMAND "python" cram.py
            ${GeneratedDir}/pbbamify.t
    )

endif()
