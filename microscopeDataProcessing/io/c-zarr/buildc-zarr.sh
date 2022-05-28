#!/bin/bash
read -p "This build script will restart the current shell and set PATH VARS in your bashrc. Continue (y/N)?" answer
if [ "${answer,,}" == "y" ]
then
	echo "Building cBlosc, cBlosc2, and cJSON";
else
	exit 0;
fi
BASEDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

CBLOSC="c-blosc-1.21.0"
CBLOSC2="c-blosc2-2.0.4"
CJSON="cJSON-1.7.15"

CBLOSCPREFIX="${BASEDIR}/cBlosc"
CBLOSC2PREFIX="${BASEDIR}/cBlosc2"
CJSONPREFIX="${BASEDIR}/cJSON"

while getopts CBLOSC2PREFIX:CJSONPREFIX flag; do
	case "${flag}" in
		CBLOSCPREFIX) CBLOSCPREFIX="${OPTARG}" ;;
		CBLOSC2PREFIX) CBLOSC2PREFIX="${OPTARG}" ;;
		CJSONPREFIX) CJSONPREFIX="${OPTARG}" ;;
	esac
done

wget https://github.com/Blosc/c-blosc/archive/refs/tags/v1.21.0.tar.gz
wget https://github.com/DaveGamble/cJSON/archive/refs/tags/v1.7.15.tar.gz
wget https://github.com/Blosc/c-blosc2/archive/refs/tags/v2.0.4.tar.gz

tar -xvf v1.21.0.tar.gz
tar -xvf v1.7.15.tar.gz
tar -xvf v2.0.4.tar.gz

rm v1.21.0.tar.gz
rm v1.7.15.tar.gz
rm v2.0.4.tar.gz

mkdir $BASEDIR/$CBLOSC/build
mkdir $BASEDIR/$CBLOSC2/build
mkdir $BASEDIR/$CJSON/build

cd $BASEDIR/$CBLOSC/build
cmake -DCMAKE_INSTALL_PREFIX=$CBLOSCPREFIX ..
make -j
make -j install

cd ../../$CBLOSC2/build
cmake -DCMAKE_INSTALL_PREFIX=$CBLOSC2PREFIX ..
make -j
make -j install

cd ../../$CJSON/build
cmake -DCMAKE_INSTALL_PREFIX=$CJSONPREFIX ..
make -j
make -j install

cd ../..

rm -rf $CBLOSC
rm -rf $CBLOSC2
rm -rf $CJSON

export PATH=":${CBLOSCPREFIX}/lib:${CBLOSC2PREFIX}/lib64:${CJSONPREFIX}/lib64:${CBLOSC2PREFIX}/lib:${CJSONPREFIX}/lib:${CBLOSCPREFIX}:${CBLOSC2PREFIX}:${CJSONPREFIX}:${PATH}"
export LD_LIBRARY_PATH=":${CBLOSCPREFIX}/lib:${CBLOSC2PREFIX}/lib64:${CJSONPREFIX}/lib64:${CBLOSC2PREFIX}/lib:${CJSONPREFIX}/lib:${CBLOSCPREFIX}:${CBLOSC2PREFIX}:${CJSONPREFIX}:${LD_LIBRARY_PATH}"

printf "\n" >> ~/.bashrc
printf "# Add paths for c-zarr\n" >> ~/.bashrc
printf "export PATH=\":${CBLOSCPREFIX}/lib:${CBLOSC2PREFIX}/lib64:${CJSONPREFIX}/lib64:${CBLOSC2PREFIX}/lib:${CJSONPREFIX}/lib:${CBLOSCPREFIX}:${CBLOSC2PREFIX}:${CJSONPREFIX}:\${PATH}\"\n" >> ~/.bashrc
printf "export LD_LIBRARY_PATH=\":${CBLOSCPREFIX}/lib:${CBLOSC2PREFIX}/lib64:${CJSONPREFIX}/lib64:${CBLOSC2PREFIX}/lib:${CJSONPREFIX}/lib:${CBLOSCPREFIX}:${CBLOSC2PREFIX}:${CJSONPREFIX}:\${LD_LIBRARY_PATH}\"\n" >> ~/.bashrc
exec "$SHELL"
