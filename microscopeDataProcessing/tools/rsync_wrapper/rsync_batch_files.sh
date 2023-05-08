#!/bin/bash

# input and output are separated by comma
IN="$1"
OUT="$2"

# inputFullpaths=($(echo $IN | tr "," "\n"))
# outputPaths=($(echo $OUT | tr "," "\n"))

IFS=',' read -r -a inputFullpaths <<< "$IN"
IFS=',' read -r -a outputPaths <<< "$OUT"

for i in "${!inputFullpaths[@]}"; do 
    # echo $i
    echo ${inputFullpaths[$i]}
    echo ${outputPaths[$i]}
    if [ ! -d "${outputPaths[$i]}" ]; then
        mkdir -p "${outputPaths[$i]}"
    fi

    rsync -aP "${inputFullpaths[$i]}" "${outputPaths[$i]}/"
    echo $(( $i + 1 )) / $(( ${#inputFullpaths[@]}))
done