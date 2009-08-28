#!/bin/bash --norc
for src in $(find dune-cornerpoint/changed_from_dune/ -type f | grep -v .svn); do
	f="${src##*dune/}"
	dest="dune-${f%%/*}/${f}"
	cmd="cp ${src} ${dest}"
	echo "${cmd}"
	${cmd}
done
