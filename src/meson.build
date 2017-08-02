###########
# sources #
###########

pbbam_cpp_sources = [
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
  'IndexedFastaReader.cpp',
  'IRecordWriter.cpp',
  'MD5.cpp',
  'MemoryUtils.cpp',
  'PbiBuilder.cpp',
  'PbiFile.cpp',
  'PbiFilter.cpp',
  'PbiFilterQuery.cpp',
  'PbiFilterTypes.cpp',
  'PbiIndex.cpp',
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
  'Tag.cpp',
  'TagCollection.cpp',
#  'UnmappedReadsQuery.cpp',
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

  # XML I/O
  'pugixml/pugixml.cpp'
]

pbbam_extra_flags = ['-DPBBAM_LIBRARY']

if get_option('enable-pbbam-auto-validate')
  pbbam_extra_flags += '-DPBBAM_AUTOVALIDATE=1'
endif

pbbam_lib = library(
  'pbbam',
  pbbam_cpp_sources,
  soversion : 0,
  install : true,
  dependencies : [pbbam_thread_dep, pbbam_boost_dep, pbbam_zlib_dep, pbbam_htslib_dep],
  include_directories : pbbam_include_directories,
  cpp_args : [pbbam_extra_flags, pbbam_warning_flags])

###################
# dependency info #
###################

import('pkgconfig').generate(
  libraries : pbbam_lib,
  version : meson.project_version(),
  name : 'pbbam',
  requires : 'htslib',
  filebase : 'pbbam',
  description : 'Library for accessing PacBio-compatible BAM files')

declare_dependency(
  include_directories : pbbam_include_directories,
  link_with : pbbam_lib,
  version : meson.project_version(),
  sources : pbbam_cpp_sources)
