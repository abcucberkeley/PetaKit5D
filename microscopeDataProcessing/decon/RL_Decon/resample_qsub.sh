
#how many cores to use; default to 1
ncores=1

#parse option flag "-c ncores" to see if more or less than 4 cores are requested
while getopts ":c:" opt; do
  case $opt in
    c)
      ncores=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))

# remember the data folder
datafolder=$1

umask 0002
# make a subfolder matlab_decon in $datafolder
mkdir "$datafolder/downsampled_data" 2>/dev/null

# $pat stores the identifying pattern in the name of the files to be processed
pat=$2

shift 2

find "$datafolder" -maxdepth 1 -name "*${pat}*.tif" -exec qsub -j y -b y -cwd -V -pe batch $ncores $HOME/matlab/resample.sh $@ \"{}\" \;
