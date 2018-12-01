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

# in order to make shared libraries consumable
# by conda and other package managers
export LDFLAGS="-static-libstdc++ -static-libgcc"

source scripts/ci/setup.sh
source scripts/ci/build.sh
source scripts/ci/test.sh

if [[ ${BUILD_NUMBER} == 0 ]]; then
  echo "Not installing anything (branch: ${bamboo_planRepository_branchName}), exiting."
  exit 0
fi

source scripts/ci/install.sh

## Creating artifact
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
