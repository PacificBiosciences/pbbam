find_package(Doxygen)

if(DOXYGEN_FOUND)
    configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${PacBioBAM_DocsDir}/Doxyfile @ONLY )
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${PacBioBAM_DocsDir}/Doxyfile
        WORKING_DIRECTORY ${PacBioBAM_DocsDir}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
endif()
