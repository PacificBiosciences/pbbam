
set(PbindexdumpSrcDir ${PacBioBAM_ToolsDir}/pbindexdump/src)

# create version header
set(PbIndexDump_VERSION ${PacBioBAM_VERSION})
configure_file(
    ${PbindexdumpSrcDir}/PbIndexDumpVersion.h.in ${GeneratedDir}/PbIndexDumpVersion.h @ONLY
)

# list source files
set(PBINDEXDUMP_SOURCES
    ${ToolsCommonDir}/OptionParser.cpp
    ${PbindexdumpSrcDir}/CppFormatter.cpp
    ${PbindexdumpSrcDir}/JsonFormatter.cpp
    ${PbindexdumpSrcDir}/PbIndexDump.cpp
    ${PbindexdumpSrcDir}/main.cpp
)

# build pbindexdump executable
include(PbbamTool)
create_pbbam_tool(
    TARGET  pbindexdump
    SOURCES ${PBINDEXDUMP_SOURCES}
)

# cram tests
if (PacBioBAM_build_tests)

    configure_file(
        ${PacBioBAM_CramTestsDir}/pbindexdump_json.t.in
        ${GeneratedDir}/pbindexdump_json.t
    )

    configure_file(
        ${PacBioBAM_CramTestsDir}/pbindexdump_cpp.t.in
        ${GeneratedDir}/pbindexdump_cpp.t
    )

    add_test(
        NAME pbindexdump_CramTests
        WORKING_DIRECTORY ${PacBioBAM_TestsDir}/scripts
        COMMAND "python" cram.py
            ${GeneratedDir}/pbindexdump_json.t
            ${GeneratedDir}/pbindexdump_cpp.t
    )

endif()
