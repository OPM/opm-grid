#!/bin/bash

source `dirname $0`/build-opm-grid.sh

ERT_REVISION=master
OPM_COMMON_REVISION=master
OPM_PARSER_REVISION=master
OPM_MATERIAL_REVISION=master
OPM_CORE_REVISION=master

build_opm_grid
test $? -eq 0 || exit 1

cp serial/build-opm-grid/testoutput.xml .
