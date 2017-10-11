#!/bin/bash -vex

# This script contains pacbio-specific build details that we do not want to push to github.
# In a purely internal project, many of these functions would be performed by cmake/make.


BUILD_NUMBER=0
if [ -n "$bamboo_planRepository_branchName" ]
  then
  if [ "$bamboo_planRepository_branchName" = "master" -o "$bamboo_planRepository_branchName" = "develop" ]
	then
	BUILD_NUMBER=${bamboo_globalBuildNumber:-0}
  fi
fi

# tool deps
CMAKE_VERSION=3

mkdir -p build
cd build
if [ $(basename `pwd`) != build ]; then
  echo $0 must be run from the build directory.  Current directory is `pwd`
  exit 1
fi

module () { eval `modulecmd bash $*`; }
cmake --version | grep "cmake version ${CMAKE_VERSION}" || {
  echo cmake ${CMAKE_VERSION} not on the path
  module use /pbi/dept/primary/modulefiles
  module load pacbio-devtools 
  module load gcc/6.4.0
  module load zlib/1.2.11
  cmake --version # print version found
  cmake --version | grep "cmake version ${CMAKE_VERSION}"
}

function resolve_dep {
  DEP_GROUP=$1 # use / not .
  DEP=$2
  DEP_VERSION=$3
  DEP_CLASSIFIER=${4:-x86_64}
  rm -rf ${DEP}-${DEP_VERSION}
  # TODO use https once it is available
  curl http://nexus.pacificbiosciences.com/repository/maven-public/${DEP_GROUP}/${DEP}/${DEP_VERSION}/${DEP}-${DEP_VERSION}-${DEP_CLASSIFIER}.tgz | tar xz
}
  
# external lib deps
resolve_dep boost boost 1.58
resolve_dep google gtest 1.7.0 src
#resolve_dep pacbio/sat/htslib htslib 4

# mhsieh invasion
CFLAGS="-fPIC"
CXXFLAGS="-fPIC"
LDFLAGS="-static-libstdc++ -static-libgcc"
export CFLAGS CXXFLAGS LDFLAGS

# TODO use standard names in cmake for this project, so these -D options can be added by resolve_dep rather than by hand
#cmake -DBOOST_ROOT=`pwd`/boost-1.58 -DGTEST_SRC_DIR=`pwd`/gtest-1.7.0 -Dhtslib_DIR=`pwd`/htslib-4 ..
cmake -DBOOST_ROOT=`pwd`/boost-1.58 -DGTEST_SRC_DIR=`pwd`/gtest-1.7.0 -DHTSLIB_INCLUDE_DIRS=/pbi/dept/secondary/builds/develop/current_thirdpartyall-release_installdir/htslib/htslib_1.3.2/include -DHTSLIB_LIBRARIES=/pbi/dept/secondary/builds/develop/current_thirdpartyall-release_installdir/htslib/htslib_1.3.2/lib/libhts.so ..
make -j
DIR=${bamboo_build_working_directory}
if [ -z "${DIR}" ]
  then
  DIR=`pwd`
fi  
mkdir -p ${DIR}/test-reports
GTEST_OUTPUT="xml:${DIR}/test-reports/pbbam_results.xml" ARGS=-V VERBOSE=1 make test
module load /mnt/software/modulefiles/cram/0.7
cram --xunit-file=${DIR}/test-reports/pbbam_cramunit.xml generated
rm -f /tmp/pbbam_*

if [ "$bamboo_planRepository_branchName" = "master" ]
  then
  VERSION=`bin/bam2sam --version`.${BUILD_NUMBER}
else
  VERSION=`bin/bam2sam --version`.SNAPSHOT${BUILD_NUMBER}
fi
NEXUS_VERSION=`bin/bam2sam --version`.${BUILD_NUMBER}

TAR_ROOT=pbbam-${VERSION}
mkdir -p ${TAR_ROOT}
cp -r lib ${TAR_ROOT}/
cp -r bin ${TAR_ROOT}/
#cp external/htslib/*.a ${TAR_ROOT}/lib/
cp -r ../include ${TAR_ROOT}/
#cp -r ../third-party/htslib/htslib/{cram,htslib} ${TAR_ROOT}/include/

f=pbbam-${VERSION}-x86_64.tgz
tar cfz $f ${TAR_ROOT}

if [ ${BUILD_NUMBER} != 0  ]; then
  if [ "$bamboo_planRepository_branchName" = "master" ]
    then
    NEXUS_REPO=maven-releases
  else
    NEXUS_REPO=maven-snapshots
  fi
  # gradle or maven could do this...
  md5sum $f | awk -e '{print $1}' >| ${f}.md5
  sha1sum $f | awk -e '{print $1}' >| ${f}.sha1

  for ext in "" .md5 .sha1 ; do
    # TODO use https once available
    NEXUS_URL=http://ossnexus.pacificbiosciences.com/repository/${NEXUS_REPO}/pacbio/sat/pbbam/pbbam/${NEXUS_VERSION}/$f${ext}
    curl -fv -n --upload-file $f${ext} ${NEXUS_URL}
  done
fi

