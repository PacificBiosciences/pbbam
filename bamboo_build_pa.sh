#!/bin/bash -vex
type module >& /dev/null || . /mnt/software/Modules/current/init/bash

module use /mnt/software/modulefiles
module use /pbi/dept/primary/modulefiles

module load cmake
module load ccache
export CCACHE_DIR="/mnt/secondary/Share/tmp/bamboo.${bamboo_shortPlanKey}.ccache"
# module load composer_xe/2017.4.196

HTSLIB_VERSION=$(/bin/ls -d src/htslib-*|sed -e 's/.*htslib-//'|sort -V|tail -1)
PBBAM_VERSION=$(grep 'PacBioBAM VERSION ' src/pbbam/CMakeLists.txt|sed -e 's/.*VERSION //'|awk '{print $1}')
# project(PacBioBAM VERSION 0.14.0 LANGUAGES CXX C)
BUILD_NUMBER=0
if [ -n "$bamboo_planRepository_branchName" ]; then
  BUILD_NUMBER=${bamboo_globalBuildNumber:-0}
fi

rm -f *.tgz || true
rm -rf pbbam* || true
# rm -rf prefix && mkdir -p prefix
# cd src/htslib-${HTSLIB_VERSION}
# export CCACHE_BASEDIR=$PWD
# CC=icc CXX=icpc \
# CFLAGS='-fPIC -Os' bash ./configure --prefix=$PWD/../../prefix
# VERBOSE=1 make CC='ccache icc' install
# rm -rf $PWD/../../prefix/lib/pkgconfig
# 
# cd -
# cd src/pbbam
# export CCACHE_BASEDIR=$PWD
# rm -rf build && mkdir -p build
# cd build
# cmake \
#   -DCMAKE_CXX_COMPILER="ccache" \
#   -DCMAKE_CXX_COMPILER_ARG1="icpc -fPIC" \
#   -DCMAKE_C_COMPILER="ccache" \
#   -DCMAKE_C_COMPILER_ARG1="icc -fPIC" \
#   -DPacBioBAM_build_shared=OFF \
#   -DPacBioBAM_build_docs=OFF \
#   -DPacBioBAM_build_tests=OFF \
#   -DHTSLIB_INCLUDE_DIRS=$PWD/../../../prefix/include \
#   -DHTSLIB_LIBRARIES=$PWD/../../../prefix/lib/libhts.a \
#   -DBoost_INCLUDE_DIRS=/mnt/software/b/boost/1.60/include \
#   -DCMAKE_BUILD_TYPE=RelWithDebInfo \
#   -DCMAKE_SKIP_BUILD_RPATH=FALSE ..
# VERBOSE=1 make -j
# tar c bin lib | tar xv -C ../../../prefix/
# cd ..
# tar c include/pbbam | tar xv -C ../../prefix/
# cd ../..
# 
# if [ ! -n "$bamboo_planRepository_branchName" ]; then
#   SNAPSHOT="_branch_"
# fi
# if [ "$bamboo_planRepository_branchName" = "develop" ]; then
#   SNAPSHOT="SNAPSHOT"
# elif [ "$bamboo_planRepository_branchName" = "master" ]; then
#   SNAPSHOT=""
# else
#   SNAPSHOT="_branch_"
# fi
# 
# rsync -ax --delete prefix/ pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}/
# 
# tar zcf pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}
# sha1sum pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz | awk -e '{print $1}' >| pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz.sha1
# md5sum  pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz | awk -e '{print $1}' >| pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz.md5
# if [ "$bamboo_planRepository_branchName" = "develop" -o "$bamboo_planRepository_branchName" = "master" ]; then
#   NEXUS_URL=http://ossnexus.pacificbiosciences.com/repository/maven-snapshots/pacbio/seq/pa/pbbam/${PBBAM_VERSION}.${BUILD_NUMBER}
#   curl -L -fvn --upload-file pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz.md5  $NEXUS_URL/pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz.md5
#   curl -L -fvn --upload-file pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz.sha1 $NEXUS_URL/pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz.sha1
#   curl -L -fvn --upload-file pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz      $NEXUS_URL/pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz
# else
#   echo "[INFO] pbbam-${PBBAM_VERSION}.SNAPSHOT${BUILD_NUMBER}-x86_64.tgz if the branch is develop"
# fi

# GCC BUILD

# module unload composer_xe/2017.4.196
# module load gcc/4.9.2
rm -rf prefix && mkdir -p prefix
cd src/htslib-${HTSLIB_VERSION}
export CCACHE_BASEDIR=$PWD
make distclean
CC=gcc CFLAGS='-fPIC -O' bash ./configure --prefix=$PWD/../../prefix --disable-bz2 --disable-lzma --disable-libcurl
VERBOSE=1 make install
rm -rf $PWD/../../prefix/lib/pkgconfig

cd -
cd src/pbbam
export CCACHE_BASEDIR=$PWD
rm -rf build && mkdir -p build
cd build
curl -sL -O http://nexus/repository/maven-thirdparty/libquadmath/4.8.5-11/libquadmath-devel-4.8.5-11.el7.x86_64.rpm
rpm2cpio libquadmath-devel-4.8.5-11.el7.x86_64.rpm | cpio -vid
CXXFLAGS="-fPIC -I$PWD/usr/lib/gcc/x86_64-redhat-linux/4.8.5/include" \
CFLAGS="-fPIC" \
cmake \
  -DCMAKE_CXX_COMPILER="g++" \
  -DCMAKE_C_COMPILER="gcc" \
  -DPacBioBAM_build_shared=OFF \
  -DPacBioBAM_build_docs=OFF \
  -DPacBioBAM_build_tests=OFF \
  -DHTSLIB_INCLUDE_DIRS=$PWD/../../../prefix/include \
  -DHTSLIB_LIBRARIES=$PWD/../../../prefix/lib/libhts.a \
  -DBoost_INCLUDE_DIRS=/mnt/software/b/boost/1.60/include \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_SKIP_BUILD_RPATH=FALSE ..
VERBOSE=1 make -j
tar c bin lib | tar xv -C ../../../prefix/
cd ..
tar c include/pbbam | tar xv -C ../../prefix/
cd ../..

if [ ! -n "$bamboo_planRepository_branchName" ]; then
  SNAPSHOT="_branch_"
fi
if [ "$bamboo_planRepository_branchName" = "develop" ]; then
  SNAPSHOT="SNAPSHOT"
elif [ "$bamboo_planRepository_branchName" = "master" ]; then
  SNAPSHOT=""
else
  SNAPSHOT="_branch_"
fi

rsync -ax --delete prefix/ pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}/

tar zcf pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}
sha1sum pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz | awk -e '{print $1}' >| pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz.sha1
md5sum  pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz | awk -e '{print $1}' >| pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz.md5
if [ "$bamboo_planRepository_branchName" = "develop" ]; then
  NEXUS_URL=http://ossnexus.pacificbiosciences.com/repository/maven-snapshots/pacbio/seq/pa/pbbam/${PBBAM_VERSION}.${BUILD_NUMBER}
elif [ "$bamboo_planRepository_branchName" = "master" ]; then
  NEXUS_URL=http://ossnexus.pacificbiosciences.com/repository/maven-releases/pacbio/seq/pa/pbbam/${PBBAM_VERSION}.${BUILD_NUMBER}
else
  echo "[INFO] pbbam-${PBBAM_VERSION}.SNAPSHOT${BUILD_NUMBER}-x86_64.tgz if the branch was develop"
  echo "[INFO] pbbam-${PBBAM_VERSION}.${BUILD_NUMBER}-x86_64.tgz if the branch was master"
  exit
fi
curl -L -fvn --upload-file pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz.md5  $NEXUS_URL/pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz.md5
curl -L -fvn --upload-file pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz.sha1 $NEXUS_URL/pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz.sha1
curl -L -fvn --upload-file pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz      $NEXUS_URL/pbbam-${PBBAM_VERSION}.${SNAPSHOT}${BUILD_NUMBER}-x86_64.tgz
