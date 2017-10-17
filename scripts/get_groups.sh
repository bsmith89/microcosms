#!/usr/bin/env bash
#
# See http://www.mothur.org/wiki/Make.contigs#file and
# http://www.mothur.org/wiki/Group_file

for file in "$@"; do
    groupname=$(basename $file | sed 's:\([^.]\+\)\..*:\1:')
    grep '^>' $file \
        | sed 's:^>\(.*\)\s*:\1:' \
        | awk -v group=$groupname '{print $1 "\t" group}' \
    || exit 1
done
