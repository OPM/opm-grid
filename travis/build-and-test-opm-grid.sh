#!/usr/bin/env bash
set -ex

pushd . > /dev/null
dune-cornerpoint/travis/build-dune-cornerpoint.sh
cd dune-cornerpoint/build
ctest --output-on-failure
popd > /dev/null
