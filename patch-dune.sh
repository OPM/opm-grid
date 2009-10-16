#!/bin/bash --norc

d=${BASH_SOURCE%/*}

for src in $(find ${d}/changed_from_dune/ -type f | grep -v .svn); do
	f="${src##*dune/}"
	dest="dune-${f%%/*}/dune/${f}"
	cmd="cp ${src} ${dest}"
	echo "${cmd}"
	${cmd}
done
