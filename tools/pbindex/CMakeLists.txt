
# create version header
set(PbIndex_VERSION ${PacBioBAM_VERSION})
configure_file(
    ${PacBioBAM_RootDir}/tools/pbindex/src/PbIndexVersion.h.in
    ${PacBioBAM_RootDir}/tools/pbindex/src/PbIndexVersion.h
)

#pbindex sources
set(PBINDEX_SOURCES
    ${PacBioBAM_RootDir}/tools/pbindex/src/main.cpp
    ${PacBioBAM_RootDir}/tools/pbindex/src/OptionParser.h
    ${PacBioBAM_RootDir}/tools/pbindex/src/OptionParser.cpp
    ${PacBioBAM_RootDir}/tools/pbindex/src/PbIndex.h
    ${PacBioBAM_RootDir}/tools/pbindex/src/PbIndex.cpp
    ${PacBioBAM_RootDir}/tools/pbindex/src/PbIndexVersion.h
)

# build pbindex executable
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PacBioBAM_CXX_FLAGS}")
include_directories(${PacBioBAM_RootDir}/tools/pbindex/src ${PacBioBAM_INCLUDE_DIRS})
add_executable(pbindex ${PBINDEX_SOURCES})
set_target_properties(pbindex PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PacBioBAM_BinDir}
)
target_link_libraries(pbindex pbbam)

