#!/usr/bin/env bash
set -vex
# This script contains pacbio-specific build details that we do not want to push to github.
# In a purely internal project, many of these functions would be performed by cmake/make.

echo "################"
echo "# DEPENDENCIES #"
echo "################"

echo "## Load modules"
type module >& /dev/null || . /mnt/software/Modules/current/init/bash

set +vx
#module load pacbio-devtools
module purge
module load gcc
module load ccache

module load meson
module load ninja

module load zlib
module load htslib
module load samtools

module load boost

module load gtest
module load cram
set -vx

BOOST_ROOT="${BOOST_ROOT%/include}"
# unset these variables to have meson discover all
# boost-dependent variables from BOOST_ROOT alone
unset BOOST_INCLUDEDIR
unset BOOST_LIBRARYDIR

export CC="ccache gcc"
export CXX="ccache g++"
export CCACHE_BASEDIR="${PWD}"

if [[ $USER == bamboo ]]; then
  export CCACHE_DIR=/mnt/secondary/Share/tmp/bamboo.${bamboo_shortPlanKey}.ccachedir
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi


echo "#########"
echo "# BUILD #"
echo "#########"

case "${bamboo_planRepository_branchName}" in
  develop|master)
    PREFIX_ARG="/mnt/software/p/pbbam/${bamboo_planRepository_branchName}"
    BUILD_NUMBER="${bamboo_globalBuildNumber:-0}"
    ;;
  *)
    BUILD_NUMBER="0"
    ;;
esac

# in order to make shared libraries consumable
# by conda and other package managers
export LDFLAGS="-static-libstdc++ -static-libgcc"

# i : unity build
for i in "on" "off"; do
  CURRENT_BUILD_DIR="build_unity=${i^^}"
  mkdir -p "${CURRENT_BUILD_DIR}"/test-reports

  echo "=============================="
  echo "Current configuration:"
  echo "  Unity:             ${i^^}"
  echo "=============================="

  # 1. configure
  # '--wrap-mode nofallback' prevents meson from downloading
  # stuff from the internet or using subprojects.
  echo "## Configuring source (${CURRENT_BUILD_DIR})"
  meson \
    --wrap-mode nofallback \
    --backend ninja \
    --buildtype release \
    -Db_ndebug=true \
    --strip \
    --default-library shared \
    --warnlevel 3 \
    --libdir lib \
    --unity "${i}" \
    --prefix "${PREFIX_ARG:-/usr/local}" \
    -Dbuild-tools=true \
    -Dtests=true \
    -Dpermissive-cigar=false \
    "${CURRENT_BUILD_DIR}" .

  # 2. build
  echo "## Building source (${CURRENT_BUILD_DIR})"
  ninja -C "${CURRENT_BUILD_DIR}" -v

  # 3. tests
  echo "## Tests (${CURRENT_BUILD_DIR})"
  GTEST_OUTPUT="xml:${CURRENT_BUILD_DIR}/test-reports/pbbam_results.xml" ARGS=-V VERBOSE=1 \
  ninja -C "${CURRENT_BUILD_DIR}" -v test
  cram --xunit-file=${CURRENT_BUILD_DIR}/test-reports/pbbam_cramunit.xml ${CURRENT_BUILD_DIR}/tools
done

if [[ -z ${PREFIX_ARG+x} ]]; then
  echo "Not installing anything (branch: ${bamboo_planRepository_branchName}), exiting."
  exit 0
fi

echo "###########"
echo "# INSTALL #"
echo "###########"

echo "## Cleaning out old installation from /mnt/software"
rm -rf "${PREFIX_ARG}"/*

echo "## Installing to /mnt/software"
ninja -C "${CURRENT_BUILD_DIR}" -v install

if [[ ${BUILD_NUMBER} == 0 ]]; then
  echo "Build number is 0, hence not creating artifact"
  exit 0
fi

echo "## Creating artifact"
# install into staging dir with --prefix /usr/local
# in order to sidestep all the artifact policy
rm -rf staging
meson configure -Dprefix=/usr/local "${CURRENT_BUILD_DIR}"
DESTDIR="${PWD}/staging" ninja -C "${CURRENT_BUILD_DIR}" -v install

if [[ ${BUILD_NUMBER} = 0 ]]; then
  exit 0
elif [[ $bamboo_planRepository_branchName == master ]]; then
  VERSION="$(${CURRENT_BUILD_DIR}/tools/bam2sam --version)".${BUILD_NUMBER}
  NEXUS_REPO=maven-releases
elif [[ $bamboo_planRepository_branchName == develop ]]; then
  VERSION="$(${CURRENT_BUILD_DIR}/tools/bam2sam --version)".SNAPSHOT${BUILD_NUMBER}
  NEXUS_REPO=maven-snapshots
else
  exit 0
fi
NEXUS_VERSION="$(${CURRENT_BUILD_DIR}/tools/bam2sam --version)".${BUILD_NUMBER}

( cd staging && tar zcf ../pbbam-${VERSION}-x86_64.tgz . )
md5sum  pbbam-${VERSION}-x86_64.tgz | awk -e '{print $1}' >| pbbam-${VERSION}-x86_64.tgz.md5
sha1sum pbbam-${VERSION}-x86_64.tgz | awk -e '{print $1}' >| pbbam-${VERSION}-x86_64.tgz.sha1

NEXUS_URL=http://ossnexus.pacificbiosciences.com/repository/${NEXUS_REPO}/pacbio/sat/pbbam/pbbam/${NEXUS_VERSION}
curl -vn --upload-file pbbam-${VERSION}-x86_64.tgz      ${NEXUS_URL}/gcc-6.4.0/pbbam-${VERSION}-x86_64.tgz
curl -vn --upload-file pbbam-${VERSION}-x86_64.tgz.md5  ${NEXUS_URL}/gcc-6.4.0/pbbam-${VERSION}-x86_64.tgz.md5
curl -vn --upload-file pbbam-${VERSION}-x86_64.tgz.sha1 ${NEXUS_URL}/gcc-6.4.0/pbbam-${VERSION}-x86_64.tgz.sha1
