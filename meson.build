project(
  'PacBioBAM',
  'cpp',
  version : '0.18.0',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++11',
    'b_ndebug=if-release'],
  license : 'BSD-3',
  meson_version : '>= 0.46.0')

############
# CXXFLAGS #
############

pbbam_warning_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wduplicated-cond',
  '-Wduplicated-branches',
  '-Wlogical-op',
  '-Wrestrict',
  '-Wnull-dereference',
  '-Wuseless-cast',
  '-Wdouble-promotion',
  '-Wshadow',
  '-Wformat=1']
    if cpp.has_argument(cflag)
      pbbam_warning_flags += cflag
    endif
endforeach

pbbam_macros = []
if get_option('permissive-cigar')
  pbbam_macros += ['-DPBBAM_PERMISSIVE_CIGAR']
  warning('**********************************************')
  warning('* You have enabled allowing "M" in BAM files *')
  warning('*    This is an unsupported combination!     *')
  warning('**********************************************')
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
pbbam_htslib_dep = dependency('htslib', required : true, version : '>=1.4')

###########
# headers #
###########

subdir('include')

#####################
# sources + library #
#####################

subdir('src')

#########
# tests #
#########

if not meson.is_subproject()
  if get_option('build-tools') or get_option('tests')
    pbbam_python = find_program('python')

    if get_option('tests')
      pbbam_clang_formatter = find_program('tools/check-formatting')
      subdir('tests')
    endif
  endif
endif

#########
# tools #
#########

if not meson.is_subproject()
  if get_option('build-tools') or get_option('tests')
    subdir('tools')
  endif
endif

#################
# documentation #
#################

if get_option('build-docs')
  subdir('docs')
endif

###################
# dependency info #
###################

if not meson.is_subproject()
  import('pkgconfig').generate(
    libraries : pbbam_lib,
    version : meson.project_version(),
    name : 'pbbam',
    requires : 'htslib',
    filebase : 'pbbam',
    description : 'Library for accessing PacBio-compatible BAM files')
endif

pbbam_dep = declare_dependency(
  include_directories : pbbam_include_directories,
  link_with : pbbam_lib,
  dependencies : pbbam_htslib_dep,
  version : meson.project_version(),
  compile_args : pbbam_macros)
