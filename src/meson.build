###########
# sources #
###########

pbbam_cpp_sources = files([
  'Accuracy.cpp',
  'AlignmentPrinter.cpp',
  'BaiIndexedBamReader.cpp',
  'BamFile.cpp',
  'BamHeader.cpp',
  'BamReader.cpp',
  'BamRecord.cpp',
  'BamRecordBuilder.cpp',
  'BamRecordImpl.cpp',
  'BamRecordTags.cpp',
  'BamTagCodec.cpp',
  'BamWriter.cpp',
  'BarcodeQuery.cpp',
  'ChemistryTable.cpp',
  'Cigar.cpp',
  'CigarOperation.cpp',
  'Compare.cpp',
  'Config.cpp',
  'DataSet.cpp',
  'DataSetBaseTypes.cpp',
  'DataSetElement.cpp',
  'DataSetIO.cpp',
  'DataSetTypes.cpp',
  'DataSetXsd.cpp',
  'EntireFileQuery.cpp',
  'FastaReader.cpp',
  'FastaSequenceQuery.cpp',
  'FastqReader.cpp',
  'FileProducer.cpp',
  'FileUtils.cpp',
  'FofnReader.cpp',
  'Frames.cpp',
  'GenomicInterval.cpp',
  'GenomicIntervalQuery.cpp',
  'IndexedBamWriter.cpp',
  'IndexedFastaReader.cpp',
  'MD5.cpp',
  'MemoryUtils.cpp',
  'PbiBuilder.cpp',
  'PbiFile.cpp',
  'PbiFilter.cpp',
  'PbiFilterQuery.cpp',
  'PbiFilterTypes.cpp',
  'PbiIndexedBamReader.cpp',
  'PbiIndexIO.cpp',
  'PbiRawData.cpp',
  'ProgramInfo.cpp',
  'QNameQuery.cpp',
  'QualityValue.cpp',
  'ReadAccuracyQuery.cpp',
  'ReadGroupInfo.cpp',
  'SamTagCodec.cpp',
  'SamWriter.cpp',
  'SequenceInfo.cpp',
  'SubreadLengthQuery.cpp',
  'Validator.cpp',
  'ValidationErrors.cpp',
  'ValidationException.cpp',
  'Version.cpp',
  'VirtualZmwBamRecord.cpp',
  'VirtualZmwCompositeReader.cpp',
  'VirtualZmwReader.cpp',
  'VirtualRegionTypeMap.cpp',
  'XmlReader.cpp',
  'XmlWriter.cpp',
  'WhitelistedZmwReadStitcher.cpp',
  'ZmwGroupQuery.cpp',
  'ZmwReadStitcher.cpp',
  'ZmwQuery.cpp',
  'ZmwTypeMap.cpp',

  # vcf
  'vcf/VcfVariant.cpp',
  'vcf/VcfFormat.cpp',
  'vcf/VcfHeader.cpp',
  'vcf/VcfReader.cpp',
  'vcf/VcfSort.cpp',
  'vcf/VcfQuery.cpp',
  'vcf/VcfWriter.cpp',

  # XML I/O
  'pugixml/pugixml.cpp'
])

pbbam_extra_flags = []
if get_option('auto-validate')
  pbbam_extra_flags += '-DPBBAM_AUTOVALIDATE=1'
endif

# install library if
# - either running as a proper project
# - or using shared libraries
pbbam_lib_install = (not meson.is_subproject()) or (get_option('default_library') == 'shared')

pbbam_lib = library(
  'pbbam',
  pbbam_cpp_sources,
  # use boost SONAME practice:
  #   cause ld startup issues before
  #   you even have the chance of running
  #   into ABI issues.
  soversion : meson.project_version(),
  version : meson.project_version(),
  install : pbbam_lib_install,
  dependencies : [pbbam_thread_dep, pbbam_boost_dep, pbbam_zlib_dep, pbbam_htslib_dep],
  include_directories : pbbam_include_directories,
  cpp_args : [pbbam_extra_flags, pbbam_warning_flags, pbbam_macros])
