###########
# bam2sam #
###########

pbbam_Bam2SamVersion_h_config = configuration_data()
pbbam_Bam2SamVersion_h_config.set('Bam2Sam_VERSION', meson.project_version())
pbbam_Bam2SamVersion_h_config.set('PacBioBAM_BinDir', meson.current_build_dir())
pbbam_Bam2SamVersion_h_config.set('PacBioBAM_TestsDir', join_paths([meson.current_source_dir(), '../tests']))

pbbam_Bam2SamVersion_h = configure_file(
  input : files('bam2sam/src/Bam2SamVersion.h.in'),
  output : 'Bam2SamVersion.h',
  configuration : pbbam_Bam2SamVersion_h_config)

pbbam_bam2sam_cpp_sources = [pbbam_Bam2SamVersion_h]
pbbam_bam2sam_cpp_sources += files([
  'bam2sam/src/main.cpp',
  'bam2sam/src/Bam2SamSettings.cpp',
  'bam2sam/src/Bam2SamWorkflow.cpp'])

pbbam_bam2sam = executable(
  'bam2sam',
  pbbam_bam2sam_cpp_sources,
  dependencies : [pbbam_htslib_dep, pbbam_thread_dep, pbbam_zlib_dep, pbbam_boost_dep, pbbam_pbcopper_dep],
  include_directories : [pbbam_include_directories, include_directories('bam2sam')],
  install : not meson.is_subproject(),
  cpp_args : pbbam_warning_flags)

# tests
if pbbam_run_tests
  bam2sam_t = configure_file(
    input : pbbam_cram_bam2sam_t_in,
    output : 'bam2sam.t',
    configuration : pbbam_Bam2SamVersion_h_config)

  test(
    'bam2sam_CramTests',
    pbbam_cram_script,
    args : [
      '--xunit-file=' + join_paths(meson.project_build_root(), 'pbbam-cram-bam2sam.xml'),
      '--verbose'] + [
        bam2sam_t],
    timeout : 1800)
endif

###########
# pbindex #
###########

pbbam_PbIndexVersion_h_config = configuration_data()
pbbam_PbIndexVersion_h_config.set('PbIndex_VERSION', meson.project_version())
pbbam_PbIndexVersion_h = configure_file(
  input : files('pbindex/src/PbIndexVersion.h.in'),
  output : 'PbIndexVersion.h',
  configuration : pbbam_PbIndexVersion_h_config)

pbbam_pbindex_cpp_sources = [pbbam_PbIndexVersion_h]
pbbam_pbindex_cpp_sources += files([
  'pbindex/src/main.cpp',
  'pbindex/src/PbIndexSettings.cpp',
  'pbindex/src/PbIndexWorkflow.cpp'])

pbbam_pbindex = executable(
  'pbindex',
  pbbam_pbindex_cpp_sources,
  dependencies : [pbbam_htslib_dep, pbbam_zlib_dep, pbbam_boost_dep, pbbam_pbcopper_dep],
  include_directories : [pbbam_include_directories, include_directories('pbindex')],
  link_with : pbbam_lib,
  install : not meson.is_subproject(),
  cpp_args : pbbam_warning_flags)

###############
# pbindexdump #
###############

pbbam_PbIndexDumpVersion_h_config = configuration_data()
pbbam_PbIndexDumpVersion_h_config.set('PbIndexDump_VERSION', meson.project_version())
pbbam_PbIndexDumpVersion_h_config.set('PacBioBAM_BinDir', meson.current_build_dir())
pbbam_PbIndexDumpVersion_h_config.set('PacBioBAM_TestsDir', join_paths([meson.current_source_dir(), '../tests']))
pbbam_PbIndexDumpVersion_h = configure_file(
  input : files('pbindexdump/src/PbIndexDumpVersion.h.in'),
  output : 'PbIndexDumpVersion.h',
  configuration : pbbam_PbIndexDumpVersion_h_config)

pbbam_pbindexdump_cpp_sources = [pbbam_PbIndexDumpVersion_h]
pbbam_pbindexdump_cpp_sources += files([
  'pbindexdump/src/CppFormatter.cpp',
  'pbindexdump/src/JsonFormatter.cpp',
  'pbindexdump/src/PbIndexDumpSettings.cpp',
  'pbindexdump/src/PbIndexDumpWorkflow.cpp',
  'pbindexdump/src/main.cpp'])

pbbam_pbindexdump = executable(
  'pbindexdump',
  pbbam_pbindexdump_cpp_sources,
  dependencies : [pbbam_htslib_dep, pbbam_zlib_dep, pbbam_boost_dep, pbbam_pbcopper_dep],
  include_directories : [pbbam_include_directories, include_directories('pbindexdump')],
  link_with : pbbam_lib,
  install : not meson.is_subproject(),
  cpp_args : pbbam_warning_flags)

# tests
if pbbam_run_tests
  pbindexdump_json_t = configure_file(
    input : pbbam_cram_pbindexdump_json_t_in,
    output : 'pbindexdump_json.t',
    configuration : pbbam_PbIndexDumpVersion_h_config)
  pbindexdump_cpp_t = configure_file(
    input : pbbam_cram_pbindexdump_cpp_t_in,
    output : 'pbindexdump_cpp.t',
    configuration : pbbam_PbIndexDumpVersion_h_config)

  test(
    'pbindexdump_CramTests',
    pbbam_cram_script,
    args : [
      '--xunit-file=' + join_paths(meson.project_build_root(), 'pbbam-cram-pbindexdump.xml'),
      '--verbose'] + [
        pbindexdump_json_t,
        pbindexdump_cpp_t],
    timeout : 1800)
endif

###########
# pbmerge #
###########

pbbam_PbMergeVersion_h_config = configuration_data()
pbbam_PbMergeVersion_h_config.set('PbMerge_VERSION', meson.project_version())
pbbam_PbMergeVersion_h_config.set('PacBioBAM_VERSION', meson.project_version())
pbbam_PbMergeVersion_h_config.set('PacBioBAM_BinDir', meson.current_build_dir())
pbbam_PbMergeVersion_h_config.set('PacBioBAM_TestsDir', join_paths([meson.current_source_dir(), '../tests']))
pbbam_PbMergeVersion_h_config.set('GeneratedTestDataDir', join_paths(meson.current_build_dir()))
pbbam_PbMergeVersion_h = configure_file(
  input : files('pbmerge/src/PbMergeVersion.h.in'),
  output : 'PbMergeVersion.h',
  configuration : pbbam_PbMergeVersion_h_config)

pbbam_pbmerge_cpp_sources = [pbbam_PbMergeVersion_h]
pbbam_pbmerge_cpp_sources += files([
  'pbmerge/src/main.cpp',
  'pbmerge/src/PbMergeSettings.cpp',
  'pbmerge/src/PbMergeWorkflow.cpp'])

pbbam_pbmerge = executable(
  'pbmerge',
  pbbam_pbmerge_cpp_sources,
  dependencies : [pbbam_boost_dep, pbbam_htslib_dep, pbbam_zlib_dep, pbbam_pbcopper_dep],
  include_directories : [pbbam_include_directories, include_directories('pbmerge')],
  link_with : pbbam_lib,
  install : not meson.is_subproject(),
  cpp_args : pbbam_warning_flags)

# tests
if pbbam_run_tests and (not get_option('auto-validate'))
  pbmerge_pacbio_ordering_t = configure_file(
    input : pbbam_cram_pbmerge_pacbio_ordering_t_in,
    output : 'pbmerge_pacbio_ordering.t', configuration : pbbam_PbMergeVersion_h_config)
  pbmerge_aligned_ordering_t = configure_file(
    input : pbbam_cram_pbmerge_aligned_ordering_t_in,
    output : 'pbmerge_aligned_ordering.t', configuration : pbbam_PbMergeVersion_h_config)
  pbmerge_mixed_ordering_t = configure_file(
    input : pbbam_cram_pbmerge_mixed_ordering_t_in,
    output : 'pbmerge_mixed_ordering.t', configuration : pbbam_PbMergeVersion_h_config)
  pbmerge_dataset_t = configure_file(
    input : pbbam_cram_pbmerge_dataset_t_in,
    output : 'pbmerge_dataset.t', configuration : pbbam_PbMergeVersion_h_config)
  pbmerge_fofn_t = configure_file(
    input : pbbam_cram_pbmerge_fofn_t_in,
    output : 'pbmerge_fofn.t', configuration : pbbam_PbMergeVersion_h_config)

  test(
    'pbmerge_CramTests',
    pbbam_cram_script,
    args : [
      '--xunit-file=' + join_paths(meson.project_build_root(), 'pbbam-cram-pbmerge.xml'),
      '--verbose'] + [
        pbmerge_pacbio_ordering_t,
        pbmerge_aligned_ordering_t,
        pbmerge_mixed_ordering_t,
        pbmerge_dataset_t,
        pbmerge_fofn_t],
    timeout : 1800)
endif

############
# pbbamify #
############

pbbam_PbBamifyVersion_h_config = configuration_data()
pbbam_PbBamifyVersion_h_config.set('PbBamify_VERSION', meson.project_version())
pbbam_PbBamifyVersion_h_config.set('PacBioBAM_VERSION', meson.project_version())
pbbam_PbBamifyVersion_h_config.set('PacBioBAM_BinDir', meson.current_build_dir())
pbbam_PbBamifyVersion_h_config.set('PacBioBAM_TestsDir', join_paths([meson.current_source_dir(), '../tests']))
pbbam_PbBamifyVersion_h_config.set('GeneratedTestDataDir', join_paths(meson.current_build_dir()))
pbbam_PbBamifyVersion_h_config.set('GeneratedDir', join_paths(meson.current_build_dir(), '../tests'))
pbbam_PbBamifyVersion_h = configure_file(
  input : files('pbbamify/src/PbBamifyVersion.h.in'),
  output : 'PbBamifyVersion.h',
  configuration : pbbam_PbBamifyVersion_h_config)

pbbam_pbbamify_cpp_sources = [pbbam_PbBamifyVersion_h]
pbbam_pbbamify_cpp_sources += files([
  'pbbamify/src/main.cpp',
  'pbbamify/src/PbBamifySettings.cpp',
  'pbbamify/src/PbBamifyWorkflow.cpp',
  'pbbamify/src/QueryLookup.cpp'
  ])

pbbam_pbbamify = executable(
  'pbbamify',
  pbbam_pbbamify_cpp_sources,
  dependencies : [pbbam_boost_dep, pbbam_htslib_dep, pbbam_zlib_dep, pbbam_pbcopper_dep],
  include_directories : [pbbam_include_directories, include_directories('pbbamify')],
  link_with : pbbam_lib,
  install : not meson.is_subproject(),
  cpp_args : pbbam_warning_flags)


# tests
if pbbam_run_tests
  pbbam_test_samtools = find_program('samtools')

  pbbamify_t = configure_file(
    input : pbbam_cram_pbbamify_t_in,
    output : 'pbbamify.t',
    configuration : pbbam_PbBamifyVersion_h_config)

  test(
    'pbbamify_CramTests',
    pbbam_cram_script,
    args : [
      '--xunit-file=' + join_paths(meson.project_build_root(), 'pbbam-cram-pbbamify.xml'),
      '--verbose'] + [
        pbbamify_t],
    timeout : 1800)
endif

############################
# ccs-kinetics-bystrandify #
############################

pbbam_CcsKineticsBystrandifyVersion_h_config = configuration_data()
pbbam_CcsKineticsBystrandifyVersion_h_config.set('CcsKineticsBystrandify_VERSION', meson.project_version())
pbbam_CcsKineticsBystrandifyVersion_h_config.set('PacBioBAM_VERSION', meson.project_version())
pbbam_CcsKineticsBystrandifyVersion_h_config.set('PacBioBAM_BinDir', meson.current_build_dir())
pbbam_CcsKineticsBystrandifyVersion_h_config.set('PacBioBAM_TestsDir', join_paths([meson.current_source_dir(), '../tests']))
pbbam_CcsKineticsBystrandifyVersion_h_config.set('GeneratedTestDataDir', join_paths(meson.current_build_dir()))
pbbam_CcsKineticsBystrandifyVersion_h_config.set('GeneratedDir', join_paths(meson.current_build_dir(), '../tests'))
pbbam_CcsKineticsBystrandifyVersion_h = configure_file(
  input : files('ccs-kinetics-bystrandify/src/CcsKineticsBystrandifyVersion.h.in'),
  output : 'CcsKineticsBystrandifyVersion.h',
  configuration : pbbam_CcsKineticsBystrandifyVersion_h_config)

pbbam_ccs_kinetics_bystrandify_cpp_sources = [pbbam_CcsKineticsBystrandifyVersion_h]
pbbam_ccs_kinetics_bystrandify_cpp_sources += files([
  'ccs-kinetics-bystrandify/src/main.cpp',
  'ccs-kinetics-bystrandify/src/CcsKineticsBystrandifySettings.cpp',
  'ccs-kinetics-bystrandify/src/CcsKineticsBystrandifyWorkflow.cpp',
  ])

pbbam_ccs_kinetics_bystrandify = executable(
  'ccs-kinetics-bystrandify',
  pbbam_ccs_kinetics_bystrandify_cpp_sources,
  dependencies : [pbbam_boost_dep, pbbam_htslib_dep, pbbam_zlib_dep, pbbam_pbcopper_dep],
  include_directories : [pbbam_include_directories, include_directories('ccs-kinetics-bystrandify')],
  link_with : pbbam_lib,
  install : not meson.is_subproject(),
  cpp_args : pbbam_warning_flags)

# tests
if pbbam_run_tests
  test(
    'ccs-kinetics-bystrandify_CramTests',
    pbbam_cram_script,
    args : [
      '--xunit-file=' + join_paths(meson.project_build_root(), 'pbbam-cram-ccs-kinetics-bystrandify-check.xml'),
      '--verbose'] + files('../tests/src/cram/ccs-kinetics-bystrandify.t'),
    env : [
      'CCS_KINETICS_BYSTRANDIFY=' + pbbam_ccs_kinetics_bystrandify.full_path(),
    ])
endif
