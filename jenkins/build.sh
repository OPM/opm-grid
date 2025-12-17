#!/bin/bash

declare -a upstreams
upstreams=(opm-common)

declare -A upstreamRev
upstreamRev[opm-common]=master

if grep -q "opm-common=" <<< $ghprbCommentBody
then
  upstreamRev[opm-common]=pull/`echo $ghprbCommentBody | sed -r 's/.*opm-common=([0-9]+).*/\1/g'`/merge
fi

# Downstreams and revisions
declare -a downstreams
downstreams=(opm-simulators
             opm-upscaling)

declare -A downstreamRev
downstreamRev[opm-simulators]=master
downstreamRev[opm-upscaling]=master

# Clone opm-common
mkdir -p $WORKSPACE/deps/opm-common
pushd $WORKSPACE/deps/opm-common
git init .
git remote add origin https://github.com/OPM/opm-common
git fetch --depth 1 origin ${upstreamRev[opm-common]}:branch_to_build
test $? -eq 0 || exit 1
git checkout branch_to_build
popd

source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh

parseRevisions
printHeader opm-grid

clone_repositories opm-grid

# Setup opm-data
if grep -q "with downstreams" <<< $ghprbCommentBody
then
  source $WORKSPACE/deps/opm-common/jenkins/setup-opm-tests.sh
fi

build_module_full opm-grid
