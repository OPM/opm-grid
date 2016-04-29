#!/bin/bash

source `dirname $0`/build-opm-grid.sh

declare -a upstreams
upstreams=(opm-parser
           opm-material
           opm-core)

declare -A upstreamRev
upstreamRev[opm-parser]=master
upstreamRev[opm-material]=master
upstreamRev[opm-core]=master

ERT_REVISION=master
OPM_COMMON_REVISION=master

if grep -q "ert=" <<< $ghprbCommentBody
then
  ERT_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*ert=([0-9]+).*/\1/g'`/merge
fi

if grep -q "opm-common=" <<< $ghprbCommentBody
then
  OPM_COMMON_REVISION=pull/`echo $ghprbCommentBody | sed -r 's/.*opm-common=([0-9]+).*/\1/g'`/merge
fi

for upstream in $upstreams
do
  if grep -q "$upstream=" <<< $ghprbCommentBody
  then
    upstreamRev[$upstream]=pull/`echo $ghprbCommentBody | sed -r "s/.*$upstream=([0-9]+).*/\1/g"`/merge
  fi
done

echo "Building with ert=$ERT_REVISION opm-common=$OPM_COMMON_REVISION opm-parser=${upstreamRev[opm-parser]} opm-material=${upstreamRev[opm-material]} opm-core=${upstreamRev[opm-core]} opm-grid=$sha1"

build_opm_grid
test $? -eq 0 || exit 1

# If no downstream builds we are done
if ! grep -q "with downstreams" <<< $ghprbCommentBody
then
  cp serial/build-opm-grid/testoutput.xml .
  exit 0
fi

# Downstreams and default revisions
declare -a downstreams
downstreams=(opm-output
             opm-simulators
             opm-upscaling
             ewoms)

declare -A downstreamRev
downstreamRev[opm-output]=master
downstreamRev[opm-simulators]=master
downstreamRev[opm-upscaling]=master
downstreamRev[ewoms]=master

source $WORKSPACE/deps/opm-common/jenkins/setup-opm-data.sh
source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh

build_downstreams opm-grid

test $? -eq 0 || exit 1
