project(
  'PacBioBAM',
  'cpp',
  version : '0.9.0',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++14'],
  license : 'BSD-3',
  meson_version : '>= 0.41.1')

############
# CXXFLAGS #
############

pbbam_warning_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wno-unused-result',
  '-Wno-unused-parameter',
  '-Wno-sign-compare',
  '-Wno-non-virtual-dtor']
    if cpp.has_argument(cflag)
      pbbam_warning_flags += cflag
    endif
endforeach

if cpp.get_id() == 'clang'
  foreach cflag: [
    '-Wno-unused-private-field',
    '-Wno-unused-command-line-argument',
    '-Wno-keyword-macro']
      if cpp.has_argument(cflag)
        pbbam_warning_flags += cflag
      endif
  endforeach
endif

################
# dependencies #
################

# threads
pbbam_thread_dep = dependency('threads', required : true)

# boost
pbbam_boost_dep = dependency('boost', required : true)

# TODO(dseifert): Add test for Winsock on Windows
# Winsock for htslib on Windows
#if(WIN32)
#    set(SOCKET_LIBRARIES "ws2_32")
#endif()

# zlib
pbbam_zlib_dep = dependency('zlib', required : true)

# htslib
pbbam_htslib_dep = dependency('htslib', required : true)

###########
# headers #
###########

subdir('include')

#####################
# sources + library #
#####################

subdir('src')

#########
# tools #
#########

pbbam_cram_py = files('tests/scripts/cram.py')
if get_option('enable-build-tools') or get_option('enable-tests')
  pbbam_python = find_program('python')
  subdir('tools')
endif

#########
# tests #
#########

subdir('tests')

#################
# documentation #
#################

if get_option('enable-build-docs')
  subdir('docs')
endif

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

pbbam_dep = declare_dependency(
  include_directories : pbbam_include_directories,
  link_with : pbbam_lib,
  dependencies : pbbam_htslib_dep,
  version : meson.project_version(),
  sources : pbbam_cpp_sources,
  compile_args : pbbam_warning_flags)
