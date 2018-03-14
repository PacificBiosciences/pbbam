#!/bin/bash

set -e
type module >& /dev/null || . /mnt/software/Modules/current/init/bash
################
# DEPENDENCIES #
################

# buildsystem dependencies
module load meson
module load ninja

module load gcc
module load ccache
if [[ $USER == "bamboo" ]]; then
  export CCACHE_DIR=/mnt/secondary/Share/tmp/bamboo.${bamboo_shortPlanKey}.ccachedir
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi
module load git

module load zlib
module load htslib
module load gtest
module load boost
module load cram
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

export CCACHE_BASEDIR=$PWD
for i in static shared; do
  for j in on off; do
    CURRENT_BUILD_DIR="build_libs=${i}_unity=${j}"
    CURRENT_CONFIG="with libs=${i^^} and unity=${j^^}"
    mkdir -p ${CURRENT_BUILD_DIR}/test-reports

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
      -Dbuild-tools=true \
      -Dtests=true \
      "${CURRENT_BUILD_DIR}" .

    # 2. build
    ninja -v -C "${CURRENT_BUILD_DIR}"

    # 3. tests
    GTEST_OUTPUT="xml:${CURRENT_BUILD_DIR}/test-reports/pbbam_results.xml" ARGS=-V VERBOSE=1 \
    ninja -v -C "${CURRENT_BUILD_DIR}" test
    cram --xunit-file=${CURRENT_BUILD_DIR}/test-reports/pbbam_cramunit.xml ${CURRENT_BUILD_DIR}/tools
  done
done
