project(
  'PacBioBAM',
  'cpp',
  version : '0.7.4',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++14'],
  license : 'BSD-3',
  meson_version : '>= 0.41.1')

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

############
# CXXFLAGS #
############

pbbam_warning_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
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
if get_option('enable-pbbam-build-tools') or get_option('enable-pbbam-tests')
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

if get_option('enable-pbbam-build-docs')
  subdir('docs')
endif
