#!/bin/bash

set -euo pipefile

script=$(basename $0)

function usage() {
cat <<EOF
USAGE: bash $script <query> <target>

<query>    :        File of reads to use as the query
<target>   :        A FastA database to use as the target

EOF
}

function print_error() {
cat <<ERR

ERROR: Command line not parsed correctly. Check input.

ERR
}

if [ $# -gt 2 ]; then
    print_error
    usage
    exit 1
fi

blast_cmd=/usr/local/wublast/latest/blastn

query=$1
target=$2
q_base=${query%.*}
t_base=${target%.*}
db=${t_base}_xddb
out=${q_base}_${t_base}_wublast.bln

## set manually, for now
cpus=4
format=2
##

xdformat -n -o $db $target

$blast_cmd $db $query M=1 N=-3 -Q 3 -R 1 -notes -errors -warnings -cpus $cpus -mformat $format -o $out 2>&1 > /dev/null
