#!/bin/bash -l

set -e

################
# DEPENDENCIES #
################

# buildsystem dependencies
module load meson/0.41.1
module load ninja/1.7.2

module load gcc/6.4.0
module load ccache/3.3.4
module load git/2.8.3

module load zlib/1.2.11
module load htslib/1.5
module load gtest/1.8.0_p20170819
module load boost/1.60
# remove trailing "/include", because CMake is brain-damaged
BOOST_ROOT=${BOOST_ROOT%/include}
# unset these variables to have meson discover all
# boost-dependent variables from BOOST_ROOT alone
unset BOOST_INCLUDEDIR
unset BOOST_LIBRARYDIR

#########
# FLAGS #
#########

export LDFLAGS="-static-libstdc++ -static-libgcc"

##########
# BUILDS #
##########

for i in static shared; do
  for j in on off; do
    CURRENT_BUILD_DIR="build_libs=${i}_unity=${j}"
    CURRENT_CONFIG="with libs=${i^^} and unity=${j^^}"

    echo "======================"
    echo "Current configuration:"
    echo "  Libs:  ${i^^}"
    echo "  Unity: ${j^^}"
    echo "======================"

    # 1. configure
    # '--wrap-mode nofallback' prevents meson from downloading
    # stuff from the internet or using subprojects.
    meson \
      --wrap-mode nofallback \
      --default-library "${i}" \
      --unity "${j}" \
      -Denable-build-tools=true \
      -Denable-tests=true \
      "${CURRENT_BUILD_DIR}" .

    # 2. build
    ninja -v -C "${CURRENT_BUILD_DIR}"

    # 3. tests
    ninja -v -C "${CURRENT_BUILD_DIR}" test
  done
done
