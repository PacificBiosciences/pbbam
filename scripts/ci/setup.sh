#!/usr/bin/env bash
set -vex

export BUILD_NUMBER="0"
export ENABLED_TESTS="true"

case "${GCC_VERSION}" in
  4.8)
    module load gtest/gcc48
    ;;

  next)
    module load gcc/8.1.0
    module load gtest
    ;;

  *)
    case "${bamboo_planRepository_branchName}-${BUILDTYPE:-release}-${ENABLED_UNITY_BUILD:-off}-${ENABLED_COVERAGE:-false}" in
      develop-release-off-false|master-release-off-false)
        export PREFIX_ARG="/mnt/software/p/pbbam/${bamboo_planRepository_branchName}"
        export BUILD_NUMBER="${bamboo_globalBuildNumber:-0}"
        ;;
    esac

    module load gcc
    module load gtest
    ;;
esac

module load ccache

export CC="ccache gcc"
export CXX="ccache g++"
export CCACHE_BASEDIR="${PWD}"

if [[ $USER == bamboo ]]; then
  export CCACHE_DIR=/mnt/secondary/Share/tmp/bamboo.${bamboo_shortPlanKey}.${bamboo_shortJobKey}.ccachedir
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi
