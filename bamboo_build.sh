#!/usr/bin/env bash
set -vex

################
# DEPENDENCIES #
################

## Load modules
type module >& /dev/null || . /mnt/software/Modules/current/init/bash

module purge

module load meson
module load ninja

module load zlib
module load htslib
module load samtools

module load boost

module load cram


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

case "${bamboo_planRepository_branchName}" in
  develop|master)
    export PREFIX_ARG="/mnt/software/p/pbbam/${bamboo_planRepository_branchName}"
    export BUILD_NUMBER="${bamboo_globalBuildNumber:-0}"
    ;;
  *)
    export BUILD_NUMBER="0"
    ;;
esac

# in order to make shared libraries consumable
# by conda and other package managers
export LDFLAGS="-static-libstdc++ -static-libgcc"

# i : unity build
for i in "on" "off"; do
  for j in "system-gcc" "gcc/8.1.0" "gcc"; do
    # 1. load either current MOBS GCC or RHEL-default GCC
    if [[ ${j} == system-gcc ]]; then
      module load gtest/gcc48
    else
      module load ${j} gtest
    fi
    module load ccache

    export CURRENT_BUILD_DIR="build_unity=${i^^}_gcc=${j/\//_}"
    export ENABLED_TESTS="true"
    export ENABLED_UNITY_BUILD="${i}"

    bash scripts/ci/build.sh
    bash scripts/ci/test.sh

    module unload ccache gtest
    [[ ${j} != system-gcc ]] && module unload gcc
  done
done

# create symlink so Bamboo can find the xunit output
ln -s "${CURRENT_BUILD_DIR}" build

if [[ -z ${PREFIX_ARG+x} ]]; then
  echo "Not installing anything (branch: ${bamboo_planRepository_branchName}), exiting."
  exit 0
fi

bash scripts/ci/install.sh

if [[ ${BUILD_NUMBER} == 0 ]]; then
  echo "Build number is 0, hence not creating artifact"
  exit 0
fi

echo "## Creating artifact"
# install into staging dir with --prefix /usr/local
# in order to sidestep all the artifact policy
rm -rf staging
meson configure -Dprefix=/usr/local -Dtests=false "${CURRENT_BUILD_DIR}"
DESTDIR="${PWD}/staging" ninja -C "${CURRENT_BUILD_DIR}" -v install

if [[ ${BUILD_NUMBER} = 0 ]]; then
  exit 0
elif [[ $bamboo_planRepository_branchName == master ]]; then
  VERSION="$(${CURRENT_BUILD_DIR}/tools/bam2sam --version)".${BUILD_NUMBER}
  NEXUS_REPO=maven-releases
elif [[ $bamboo_planRepository_branchName == develop ]]; then
  VERSION="$(${CURRENT_BUILD_DIR}/tools/bam2sam --version)".SNAPSHOT${BUILD_NUMBER}
  NEXUS_REPO=maven-snapshots
  rm -rf /mnt/secondary/builds/unsupported/pbbam.previous
  if [[ -e /mnt/secondary/builds/unsupported/pbbam ]]; then
    mv /mnt/secondary/builds/unsupported/pbbam \
       /mnt/secondary/builds/unsupported/pbbam.previous
  fi
  DESTDIR="/mnt/secondary/builds/unsupported/pbbam/" ninja -C "${CURRENT_BUILD_DIR}" -v install
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
