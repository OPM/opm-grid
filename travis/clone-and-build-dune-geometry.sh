#!/usr/bin/env bash
set -ex

pushd . > /dev/null
git clone https://github.com/dune-project/dune-geometry.git
cd dune-geometry
git checkout tags/v2.3.1
mkdir build
cd build
cmake ../
make
popd > /dev/null
