#!/bin/bash

source `dirname $0`/build-opm-grid.sh

declare -a upstreams
upstreams=(opm-parser
           opm-output
           opm-material
           opm-core)

declare -A upstreamRev
upstreamRev[opm-parser]=master
upstreamRev[opm-output]=master
upstreamRev[opm-material]=master
upstreamRev[opm-core]=master

ERT_REVISION=master
OPM_COMMON_REVISION=master

build_opm_grid
test $? -eq 0 || exit 1

cp serial/build-opm-grid/testoutput.xml .
