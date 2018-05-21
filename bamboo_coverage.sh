#!/usr/bin/env bash
set -xve

source /mnt/software/Modules/current/init/bash
module load gcc meson ccache ninja zlib htslib samtools cram boost gtest gcov

echo "#####################"
echo "# BUILD & RUN TESTS #"
echo "#####################"

rm -rf build
mkdir build 
cd build

meson \
  --werror \
  --backend ninja \
  --buildtype debug \
  --default-library shared \
  --libdir lib \
  --wrap-mode nofallback \
  --prefix "${PREFIX_ARG:-/usr/local}" \
  -Db_coverage=true \
  ..

ninja test

echo "################"
echo "# COVERAGE     #"
echo "################"

find . -type f -iname '*.o' | xargs gcov -acbrfu {} \; >/dev/null && \
mkdir coverage && pushd coverage && mv ../*.gcov . && \
sed -i -e 's@Source:@Source:../@' *.gcov && \
sed -i -e 's@Graph:@Graph:../@' *.gcov && \
sed -i -e 's@Data:@Data:../@' *.gcov && \
rm pugixml*

