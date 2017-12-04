

if(DEFINED PacBioBAM_build_pbindex)

    # Deprecating the "PacBioBAM_build_pbindex" command line option in favor of more
    # general "PacBioBAM_build_tools", as we're starting to add new utilities.
    #
    # That said, I don't want to break current auto tests/builds, so I'm providing a
    # warning message so devs are aware.
    #
    # construct warning message
    set(pbindex_warning "\nDeprecated:\n-DPacBioBAM_build_pbindex\n")
    if (PacBioBAM_build_pbindex)
        set(pbindex_warning "${pbindex_warning} Building as requested,")
    else()
        set(pbindex_warning "${pbindex_warning} Skipping as requested,")
    endif()
    set(pbindex_warning "${pbindex_warning} but support for this option will be removed at some point in the future.\n")
    message(AUTHOR_WARNING "${pbindex_warning} ** Use -DPacBioBAM_build_tools instead. **\n")

    # force PacBioBAM_build_tools option
    set(PacBioBAM_build_tools
        ${PacBioBAM_build_pbindex} CACHE BOOL
        "Build PacBioBAM with add'l utilities (e.g. pbindex, pbindexdump)." FORCE)
endif()

if (PacBioBAM_build_tools)

    # tools directory
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
    add_subdirectory(pbbamify)

endif()
