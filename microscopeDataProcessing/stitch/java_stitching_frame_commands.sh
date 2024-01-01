#!/bin/bash
#
# Currently the script only support single channel stitching

# parse options

AxisOrder=''
ImagePath=''
channel=''
Flatfield=''
ImageListFileName=''
Resolution=''

print_usage() {
  printf "Usage: Stitching for a single frame and a single channel"
}

while getopts 'a:b:c:fi:r:' flag; do
  case "${flag}" in
    a) AxisOrder="${OPTARG}" ;;
    b) ImagePath="${OPTARG}" ;;
    c) channel="${OPTARG}" ;;
    f) Flatfield='true' ;;
    i) ImageListFileName="${OPTARG}" ;;
    r) Resolution="${OPTARG}" ;;
    *) print_usage
       exit 1 ;;
  esac
done

Ch=$channel

# save intermediate results in the folder of the image list csv file
ResultDir=$(dirname "${ImageListFileName}")

# move to the result dir
cd $ResultDir

# load module
# module load stitching-spark/01102020

STITCHING_DIR=/home/xruan/Projects/tmp/stitching-spark

# generate JSON file
python $STITCHING_DIR/startup-scripts/spark-local/parse-imagelist-metadata.py -b $ImagePath -a $AxisOrder -i $ImageListFileName  -r $Resolution


# convert tiff to n5 
python $STITCHING_DIR/startup-scripts/spark-local/convert-tiff-tiles-n5.py -i ${Ch}nm-n5.json -o $ResultDir


# flatfield correction
if [ $Flatfield ]
then
    python $STITCHING_DIR/startup-scripts/spark-local/flatfield.py -i ${Ch}nm-n5.json
fi


# stitching
# first remove old optimation results if there is
for d in iter* ; do
    if [ -d "$ResultDir/$d" ]
    then
    rm -r $ResultDir/$d
    fi
done

export _JAVA_OPTIONS="-Xms51200m -Xmx102400m"
python $STITCHING_DIR/startup-scripts/spark-local/stitch.py -i ${Ch}nm-n5.json -m full


# export to n5
python $STITCHING_DIR/startup-scripts/spark-local/export.py -i ${Ch}nm-n5-final.json


# convert n5 to tiff
python $STITCHING_DIR/startup-scripts/spark-local/n5-slice-tiff.py -i export.n5


