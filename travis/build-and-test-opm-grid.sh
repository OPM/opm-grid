#!/usr/bin/env bash
set -ex

pushd . > /dev/null
opm-grid/travis/build-opm-grid.sh
cd opm-grid/build
ctest --output-on-failure
popd > /dev/null
