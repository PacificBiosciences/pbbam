pbbam_test_cpp_sources = files([
  'src/test_Accuracy.cpp',
  'src/test_AlignmentPrinter.cpp',
  'src/test_BamFile.cpp',
  'src/test_BamHeader.cpp',
  'src/test_BamRecord.cpp',
  'src/test_BamRecordBuilder.cpp',
  'src/test_BamRecordClipping.cpp',
  'src/test_BamRecordImplCore.cpp',
  'src/test_BamRecordImplTags.cpp',
  'src/test_BamRecordImplVariableData.cpp',
  'src/test_BamRecordMapping.cpp',
  'src/test_BamWriter.cpp',
  'src/test_BarcodeQuery.cpp',
  'src/test_Cigar.cpp',
  'src/test_Compare.cpp',
  'src/test_DataSetCore.cpp',
  'src/test_DataSetIO.cpp',
  'src/test_DataSetQuery.cpp',
  'src/test_DataSetXsd.cpp',
  'src/test_EndToEnd.cpp',
  'src/test_EntireFileQuery.cpp',
  'src/test_Fasta.cpp',
  'src/test_Fastq.cpp',
  'src/test_FileUtils.cpp',
  'src/test_Frames.cpp',
  'src/test_GenomicIntervalQuery.cpp',
  'src/test_IndexedFastaReader.cpp',
  'src/test_Intervals.cpp',
  'src/test_PacBioIndex.cpp',
  'src/test_PbiFilter.cpp',
  'src/test_PbiFilterQuery.cpp',
  'src/test_QNameQuery.cpp',
  'src/test_QualityValues.cpp',
  'src/test_Pulse2BaseCache.cpp',
  'src/test_ReadAccuracyQuery.cpp',
  'src/test_ReadGroupInfo.cpp',
  'src/test_SamWriter.cpp',
  'src/test_SequenceUtils.cpp',
  'src/test_StringUtils.cpp',
  'src/test_SubreadLengthQuery.cpp',
  'src/test_Tags.cpp',
  'src/test_TimeUtils.cpp',
#  'src/test_UnmappedReadsQuery.cpp',
  'src/test_Validator.cpp',
  'src/test_Version.cpp',
  'src/test_WhitelistedZmwReadStitcher.cpp',
  'src/test_ZmwReadStitcher.cpp',
  'src/test_ZmwQuery.cpp'])

if get_option('enable-pbbam-tests')
  pbbam_gtest_dep = dependency('gtest_main', fallback : ['gtest', 'gtest_with_main_dep'])

  pbbam_PbbamTestData_h_config = configuration_data()
  pbbam_PbbamTestData_h_config.set('PacBioBAM_BinDir', join_paths([meson.build_root(), 'tools']))
  pbbam_PbbamTestData_h_config.set('PacBioBAM_TestsDir', meson.current_source_dir())
  pbbam_PbbamTestData_h_config.set('CMAKE_CURRENT_BINARY_DIR', meson.current_build_dir())
  pbbam_PbbamTestData_h_config.set('GeneratedDir', meson.current_build_dir())
  pbbam_PbbamTestData_h_config.set('GeneratedTestDataDir', meson.current_build_dir())

  configure_file(input : 'src/PbbamTestData.h.in', output : 'PbbamTestData.h', configuration : pbbam_PbbamTestData_h_config)
  configure_file(input : 'data/group/group.fofn.in', output : 'group.fofn', configuration : pbbam_PbbamTestData_h_config)

  pbbam_test = executable(
    'pbbam_test',
    pbbam_test_cpp_sources,
    dependencies : [pbbam_gtest_dep, pbbam_boost_dep, pbbam_htslib_dep],
    include_directories : pbbam_include_directories,
    link_with : pbbam_lib,
    cpp_args : [pbbam_warning_flags, '-DPBBAM_TESTING'],
    install : false)

  custom_target('pbbam_generate_data',
    output : 'input.fa',
    command : [
      pbbam_python,
      join_paths([meson.current_source_dir(), 'scripts/generate_data.py']),
      join_paths([meson.current_source_dir(), 'data']),
      meson.current_build_dir()],
    build_by_default : true,
    install : false)

  test('pbbam test', pbbam_test)
endif
