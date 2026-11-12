#!/bin/bash

# arguments: pairs of an input file and its output folder, one argument each
n=$(( $# / 2 ))
i=0
while [ $# -ge 2 ]; do
    i=$(( i + 1 ))
    echo "$1"
    echo "$2"
    if [ ! -d "$2" ]; then
        mkdir -p "$2"
    fi

    rsync -aP "$1" "$2/"
    echo $i / $n
    shift 2
done
