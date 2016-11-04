include(CMakeParseArguments)

function(create_pbbam_tool)

    # parse args
    set(oneValueArgs TARGET)
    set(multiValueArgs SOURCES)
    cmake_parse_arguments(create_pbbam_tool "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # create executable
    include_directories(
        ${ToolsCommonDir}           # shared tool code
        ${GeneratedDir}             # generated version headers
        ${PacBioBAM_INCLUDE_DIRS}   # pbbam/htslib includes
    )
    add_executable(${create_pbbam_tool_TARGET} ${create_pbbam_tool_SOURCES})
    set_target_properties(
        ${create_pbbam_tool_TARGET} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${PacBioBAM_BinDir}
    )
    target_link_libraries(${create_pbbam_tool_TARGET} pbbam)

endfunction(create_pbbam_tool)
