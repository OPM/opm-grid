# opm-grid jenkins build scripts:

**build-opm-grid.sh**:
This is a helper script which contains a function for building,
testing and cloning opm-grid and its dependencies.

**build.sh**:
This script will build dependencies, then build opm-grid and execute its tests.
It is intended for post-merge builds of the master branch.

**build-pr.sh**:
This script will build dependencies, then build opm-grid and execute its tests.
It inspects the $ghbPrBuildComment environmental variable to obtain a pull request
to use for ert, opm-common, opm-parser, opm-material and
opm-core (defaults to master) and then builds $sha1 of opm-grid.

It is intended for pre-merge builds of pull requests.

You can optionally specify a given pull request to use for opm-common, opm-parser, opm-material and opm-core through the trigger.
The trigger line needs to contain ert=&lt;pull request number&gt; and/or
opm-common=&lt;pull request number&gt; and/or opm-parser=&lt;pull request number&gt;
and/or opm-material=&lt;pull request number&gt; and/or opm-core=&lt;pull request number&gt;.
