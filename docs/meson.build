#################
# documentation #
#################

doxygen = find_program('doxygen', required : true)

pbbam_doxygen_config = configuration_data()
pbbam_doxygen_config.set('PacBioBAM_NAME', meson.project_name())
pbbam_doxygen_config.set('PacBioBAM_VERSION', meson.project_version())
pbbam_doxygen_config.set('PacBioBAM_DocsDir', '.')
pbbam_doxygen_config.set('PacBioBAM_IncludeDir', join_paths([meson.current_source_dir(), '../include']))

doxyfile = configure_file(
  input : 'Doxyfile.in',
  output : 'Doxyfile',
  configuration : pbbam_doxygen_config,
  install : false)

custom_target('docs',
  input : doxyfile,
  output : 'docs',
  command : [doxygen, doxyfile],
  build_by_default : true,
  install : false)
