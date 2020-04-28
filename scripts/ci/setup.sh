#!/usr/bin/env bash
set -vex

export ENABLED_TESTS="true"

case "${bamboo_planRepository_branchName}" in
  master)
    _pbcopper_module="pbcopper/master"
    ;;
  *)
    _pbcopper_module="pbcopper/develop"
    ;;
esac

case "${GCC_VERSION}" in
  next)
    module load gcc/8.1.0
    module load gtest
    module load ${_pbcopper_module}
    ;;

  PA)
    # have to build htslib for PA
    module unload htslib
    module load zlib
    module load gtest/gcc48

    # load SCL GCC
    source /opt/rh/devtoolset-6/enable

    export NEXUS_PROJECT=pacbio/seq/pa/pbbam
    export NEXUS_TC=""
    export _artifact_versionprepend="true"
    ;;

  ICC)
    module load devtoolset/6
    module load composer_xe/2017.4.196
    module load gtest/gcc48

    CC="icc"
    CXX="icpc"
    ;;

  clang)
    module load gtest/gcc48

    source /opt/rh/llvm-toolset-6.0/enable
    CC="clang"
    CXX="clang++"
    ;;

  *)
    module load gcc

    if [[ ${ENABLED_WRAP_MODE:-nofallback} == nofallback ]]; then
      module load gtest
      module load ${_pbcopper_module}
    else
      # don't want to rely on internet connectivity for CI
      mkdir -p subprojects/packagecache
      cp /mnt/software/m/meson/packagecache/gtest* subprojects/packagecache/
    fi
    ;;
esac

module load ccache

export CC="ccache ${CC:-gcc}"
export CXX="ccache ${CXX:-g++}"
export CCACHE_BASEDIR="${PWD}"

if [[ -z ${bamboo_planRepository_branchName+x} ]]; then
  : #pass
elif [[ ! -d /pbi/flash/bamboo/ccachedir ]]; then
  echo "[WARNING] /pbi/flash/bamboo/ccachedir is missing"
elif [[ $bamboo_planRepository_branchName == develop ]]; then
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}.develop
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
elif [[ $bamboo_planRepository_branchName == master ]]; then
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}.master
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
elif [[ $USER == bamboo ]]; then
  _shortPlanKey=$(echo ${bamboo_shortPlanKey}|sed -e 's/[0-9]*$//')
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}
  if [[ -d /pbi/flash/bamboo/ccachedir/${_shortPlanKey}.${bamboo_shortJobKey}.develop ]]; then
    ( cd /pbi/flash/bamboo/ccachedir/
      cp -a ${_shortPlanKey}.${bamboo_shortJobKey}.develop $CCACHE_DIR
    )
  fi
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi
