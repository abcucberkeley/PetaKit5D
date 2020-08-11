#!/bin/bash
export MCR_CACHE_ROOT=/scratch/shaol/mcr_cache_root.$JOB_ID

umask 0002

resizeFactor=1
scalingThresh=1

# -t: resize (down or up-sampling) factor (default 1)
# -l: rescaling threshold (ask Wes; only has effect when -t is used and not 1) (default 1)

while getopts "t:l:" opt; do
  case $opt in
    t)
      resizeFactor=$OPTARG
      ;;
    l)
      scalingThresh=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))

$HOME/matlab/rescale_resample $1 $resizeFactor $scalingThresh

rm -rf /scratch/shaol/mcr_cache_root.$JOB_ID
