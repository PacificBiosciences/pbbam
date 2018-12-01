#!/usr/bin/env bash
set -vex

########
# TEST #
########

type module >& /dev/null || . /mnt/software/Modules/current/init/bash

ninja -C "${CURRENT_BUILD_DIR:-build}" -v test

############
# COVERAGE #
############

if [[ ${ENABLED_COVERAGE} == true ]]; then
  pushd "${CURRENT_BUILD_DIR:-build}"
  find . -type f -iname '*.o' | xargs gcov -acbrfu {} \; >/dev/null && \
    mkdir coverage && pushd coverage && mv ../*.gcov . && \
    sed -i -e 's@Source:@Source:../@' *.gcov && \
    sed -i -e 's@Graph:@Graph:../@' *.gcov && \
    sed -i -e 's@Data:@Data:../@' *.gcov && \
    rm pugixml* && popd
  popd
fi
