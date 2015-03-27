# test case headers
set( PacBioBAMTest_H

)

# test case sources
set( PacBioBAMTest_CPP

    ${PacBioBAM_TestsDir}/src/test_Accuracy.cpp
    ${PacBioBAM_TestsDir}/src/test_BamHeader.cpp
    ${PacBioBAM_TestsDir}/src/test_BamReader.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecord.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecordBuilder.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecordClipping.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecordImplCore.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecordImplTags.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecordImplVariableData.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecordMapping.cpp
    ${PacBioBAM_TestsDir}/src/test_BamWriter.cpp
    ${PacBioBAM_TestsDir}/src/test_Cigar.cpp
    ${PacBioBAM_TestsDir}/src/test_EndToEnd.cpp
    ${PacBioBAM_TestsDir}/src/test_EntireFileQuery.cpp
    ${PacBioBAM_TestsDir}/src/test_GenomicIntervalQuery.cpp
    ${PacBioBAM_TestsDir}/src/test_IndexedFastaReader.cpp
    ${PacBioBAM_TestsDir}/src/test_Intervals.cpp
    ${PacBioBAM_TestsDir}/src/test_QualityValues.cpp
    ${PacBioBAM_TestsDir}/src/test_SequenceUtils.cpp
#    ${PacBioBAM_TestsDir}/src/test_SamHeader.cpp
    ${PacBioBAM_TestsDir}/src/test_Tags.cpp
    ${PacBioBAM_TestsDir}/src/test_UnmappedReadsQuery.cpp
)

# GoogleTest headers
set( GTest_H

    ${GTest_IncludeDir}/gtest/gtest-death-test.h
    ${GTest_IncludeDir}/gtest/gtest-message.h
    ${GTest_IncludeDir}/gtest/gtest-param-test.h
    ${GTest_IncludeDir}/gtest/gtest-printers.h
    ${GTest_IncludeDir}/gtest/gtest-spi.h
    ${GTest_IncludeDir}/gtest/gtest-test-part.h
    ${GTest_IncludeDir}/gtest/gtest-typed-test.h
    ${GTest_IncludeDir}/gtest/gtest.h
    ${GTest_IncludeDir}/gtest/gtest_pred_impl.h
    ${GTest_IncludeDir}/gtest/gtest_prod.h
    ${GTest_IncludeDir}/gtest/internal/gtest-death-test-internal.h
    ${GTest_IncludeDir}/gtest/internal/gtest-filepath.h
    ${GTest_IncludeDir}/gtest/internal/gtest-internal.h
    ${GTest_IncludeDir}/gtest/internal/gtest-linked_ptr.h
    ${GTest_IncludeDir}/gtest/internal/gtest-param-util-generated.h
    ${GTest_IncludeDir}/gtest/internal/gtest-param-util.h
    ${GTest_IncludeDir}/gtest/internal/gtest-port.h
    ${GTest_IncludeDir}/gtest/internal/gtest-string.h
    ${GTest_IncludeDir}/gtest/internal/gtest-tuple.h
    ${GTest_IncludeDir}/gtest/internal/gtest-type-util.h

    ${GTest_SourceDir}/gtest-internal-inl.h
)

# GoogleTest sources
set( GTest_CPP
    ${GTest_SourceDir}/gtest-all.cc
    ${GTest_SourceDir}/gtest_main.cc
)
