#!/usr/bin/env bash
set -vex

############
# ARTIFACT #
############

if [[ ${_create_artifact} != true ]]; then
  echo "Not creating artifact (branch: ${bamboo_planRepository_branchName}), returning."
  return 0
fi

# *never* create artifacts with ASAN enabled
meson configure -Dprefix=/ -Db_sanitize=none "${CURRENT_BUILD_DIR:-build}"

NEXUS_VERSION="$(${CURRENT_BUILD_DIR:-build}/tools/bam2sam --version)".${BUILD_NUMBER}
case "${bamboo_planRepository_branchName}" in
  develop)
    VERSION="$(${CURRENT_BUILD_DIR:-build}/tools/bam2sam --version)".SNAPSHOT${BUILD_NUMBER}
    NEXUS_REPO=maven-snapshots
    ;;
  master)
    VERSION="${NEXUS_VERSION}"
    NEXUS_REPO=maven-releases
    ;;
  *)
    echo "You can only create artifacts from 'develop' or 'master' branches"
    exit 1
    ;;
esac

DESTDIR="${PWD}/staging" ninja -C "${CURRENT_BUILD_DIR:-build}" -v install
if [[ ${_artifact_versionprepend:-false} == true ]]; then
  ( cd staging && tar --transform "s,^,${VERSION}," zcf ../pbbam-${VERSION}-x86_64.tgz . )
else
  ( cd staging && tar zcf ../pbbam-${VERSION}-x86_64.tgz . )
fi

md5sum  pbbam-${VERSION}-x86_64.tgz | awk -e '{print $1}' >| pbbam-${VERSION}-x86_64.tgz.md5
sha1sum pbbam-${VERSION}-x86_64.tgz | awk -e '{print $1}' >| pbbam-${VERSION}-x86_64.tgz.sha1

NEXUS_URL=http://ossnexus.pacificbiosciences.com/repository/${NEXUS_REPO}/${NEXUS_PROJECT:-pacbio/sat/pbbam/pbbam}/${NEXUS_VERSION:-gcc-6.4.0}${NEXUS_TC}
curl -vn --upload-file pbbam-${VERSION}-x86_64.tgz      ${NEXUS_URL}/pbbam-${VERSION}-x86_64.tgz
curl -vn --upload-file pbbam-${VERSION}-x86_64.tgz.md5  ${NEXUS_URL}/pbbam-${VERSION}-x86_64.tgz.md5
curl -vn --upload-file pbbam-${VERSION}-x86_64.tgz.sha1 ${NEXUS_URL}/pbbam-${VERSION}-x86_64.tgz.sha1
