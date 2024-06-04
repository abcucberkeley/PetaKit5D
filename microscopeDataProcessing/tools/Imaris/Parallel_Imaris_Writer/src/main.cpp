/***************************************************************************
 *   Copyright (c) 2020-present Bitplane AG Zuerich                        *
 *                                                                         *
 *   Licensed under the Apache License, Version 2.0 (the "License");       *
 *   you may not use this file except in compliance with the License.      *
 *   You may obtain a copy of the License at                               *
 *                                                                         *
 *       http://www.apache.org/licenses/LICENSE-2.0                        *
 *                                                                         *
 *   Unless required by applicable law or agreed to in writing, software   *
 *   distributed under the License is distributed on an "AS IS" BASIS,     *
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or imp   *
 *   See the License for the specific language governing permissions and   *
 *   limitations under the License.                                        *
 ***************************************************************************/

/*
gcc -I. -L. -o bpImarisWriter96TestProgram.exe bpImarisWriter96TestProgram.c -lbpImarisWriter96
./bpImarisWriter96TestProgram.exe
*/
#include "bpImageConverterInterfaceC.h"
#include "helperfunctions.h"
#include "parallelreadtiff.h"
#include "zarr.h"
#include "parallelreadzarr.h"
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>
#include <stdint.h>
#include <sys/stat.h>
#include <omp.h>

#ifndef _WIN32
#include <wordexp.h>
// Expand the tilde to the home directory on Mac and Linux
const char* expandTilde(const char* path) {
    if(strchr(path,'~')){
        wordexp_t expPath;
        wordexp(path, &expPath, 0);
        return expPath.we_wordv[0];
    }
    else return path;
}
#endif

typedef struct
{
    unsigned int mImageIndex;
    int mProgress;
} bpCallbackData;

static int compare(const void* str1, const void* str2) {
   if(strcmp(*(const char**)str1, *(const char**)str2) >= 0)
      return 1;
   else return 0;
}

void listFiles(char *basePath, char* pattern, uint8_t build, uint64_t cPattN, char** files, uint8_t count, int64_t* timepoints)
{
    //char path[10000];
    struct dirent *dp;
    DIR *dir = opendir(basePath);
    uint64_t baseSize = strlen(basePath)+1;
    uint64_t cFile = *timepoints*cPattN;

    // Unable to open directory stream
    if (!dir)
        return;

    while ((dp = readdir(dir)) != NULL)
    {
        if (strcmp(dp->d_name, ".") != 0 && strcmp(dp->d_name, "..") != 0 && strrchr(dp->d_name, '.'))
        {
            if(((!strcmp(strrchr(dp->d_name, '.'),".tif")) || (!strcmp(strrchr(dp->d_name, '.'),".zarr"))) && strstr(dp->d_name,pattern) != NULL){
                if(count) *timepoints += 1;
                if(build){
                    files[cFile] = (char*)malloc(baseSize+strlen(dp->d_name)+1);
                    files[cFile][0] = '\0';
                    strcat(files[cFile],basePath);
                    #ifdef _WIN32
                    strcat(files[cFile],"\\");
                    #else
                    strcat(files[cFile],"/");
                    #endif
                    strcat(files[cFile],dp->d_name);
                    cFile++;
                }
            }
        }
    }

    closedir(dir);
}

void ProgressCallback(bpConverterTypesC_Float aProgress, bpConverterTypesC_UInt64 aTotalBytesWritten, void* aUserData)
{
    bpCallbackData* vCallbackData = (bpCallbackData*)aUserData;

    int vProgress = (int)(aProgress * 100);
    if (vProgress - vCallbackData->mProgress < 5) {
        return;
    }

    unsigned int vImageIndex = vCallbackData->mImageIndex;
    if (aTotalBytesWritten < 10 * 1024 * 1024) {
        printf("Progress image %u: %d%% [%llu KB]\n", vImageIndex, vProgress, aTotalBytesWritten / 1024);
    }
    else {
        printf("Progress image %u: %d%% [%llu MB]\n", vImageIndex, vProgress, aTotalBytesWritten / (1024 * 1024));
    }
    vCallbackData->mProgress = vProgress;
}

unsigned int NumBlocks(unsigned int aSize, unsigned int aBlockSize)
{
    return (aSize + aBlockSize - 1) / aBlockSize;
}

void CheckErrors(bpImageConverterCPtr aConverter)
{
    const char* vException = bpImageConverterC_GetLastException(aConverter);
    if (vException) {
        printf("%s", vException);
        exit(1);
    }
}

void convertToImaris(int argc, char **argv)
{
    int64_t timepoints = -1;
    int64_t channels = -1;
    int option_index = 0;
    char* vVoxelSizes = NULL;
    float vVoxelSizeX = -1.0;
    float vVoxelSizeY = -1.0;
    float vVoxelSizeZ = -1.0;
    char* fileNames = NULL;
    char* folderNames = NULL;
    char* patterns = NULL;
    DIR* outputDir = NULL;
    char* outputDirString = NULL;
    char* outName = NULL;
	char* reader = NULL;
	char* blockSizes = NULL;
	
	uint8_t crop = 0;
	char* boundingBoxString = NULL;
	uint64_t boundingBox[6] = {0,0,0,0,0,0};
    while (( option_index = getopt(argc, argv, ":c:P:t:f:F:o:r:v:n:b:B:")) != -1){
        switch (option_index) {
        case 'c':
            channels = atoi(optarg);
            break;
        case 'P':
            patterns = strdup(optarg);
            break;
        case 't':
            timepoints = atoi(optarg);
            break;
        case 'f':
            fileNames = strdup(optarg);
            break;
        case 'F':
            folderNames = strdup(optarg);
            break;
        case 'o':
            outputDirString = strdup(optarg);
            break;
        case 'r':
            reader = strdup(optarg);
            break;
        case 'v':
            vVoxelSizes = strdup(optarg);
            break;
		case 'n':
			outName = strdup(optarg);
			break;
		case 'b':
			blockSizes = strdup(optarg);
			break;
		case 'B':
			crop = 1;
			boundingBoxString = strdup(optarg);
			break;
        default:
            printf("Option incorrect\n");
        }
    }
    const char* delim = ",";

    // Expand the tilde to the home directory on Mac and Linux
    #ifndef _WIN32
    std::stringstream folderNamesStrean(folderNames);
    std::string token;
    std::string folderNamesString;

    while(std::getline(folderNamesStrean, token, delim[0]))
    {
        folderNamesString.append(expandTilde(token.c_str()));
        folderNamesString.append(delim);
    }
    folderNamesString.pop_back();
    free(folderNames);
    folderNames = strdup(folderNamesString.c_str());
    #endif

    if(!reader){
        printf("Reader type not specified. Defaulting to tiff.\n");
        reader = "tiff";
    }
    if(strcmp(reader,"tiff") && strcmp(reader,"zarr")){
        printf("ERROR: Reader type \"%s\" is not a supported reader type.\n",reader);
        return;
    }

    if(!outputDirString){
        if(!folderNames){
            printf("ERROR: Output Directory is Required and there are no Folder names provided to use as a backup output folder. Use argument -o path/to/dir for a custom output folder.\n");
            return;
        }
        char* backupFolders = strdup(folderNames);
        char *saveptr0 = NULL;
        char* backupFolder = strtok_r(backupFolders,delim,&saveptr0);
        outputDirString = strdup(backupFolder);
        printf("Output Directory not provided. Using %s as the output directory.\n",outputDirString);
    }
    #ifdef WIN32
    int madeDir = mkdir(outputDirString);
    #else
    // Expand the tilde to the home directory on Mac and Linux
    const char* expandedOutputDirString = expandTilde(outputDirString);
    free(outputDirString);
    outputDirString = strdup(expandedOutputDirString);
    int madeDir = mkdir(outputDirString, 0775);
    #endif
    outputDir = opendir(outputDirString);
    if(!outputDir){
        printf("ERROR: Given Output Directory \"%s\" does not exist and was unable to be created.\n",outputDirString);
        return;
    }
	// Check if we can write to the output dir
	if(access(outputDirString, W_OK)){
		printf("ERROR: Given Output Directory \"%s\" exists but cannot be written to.\n",outputDirString);
		return;
	}
	if (!outName) outName = "output";
    char* timepointsFile = (char*)malloc(strlen(outputDirString)+16);
    char* outputFile = (char*)malloc(strlen(outputDirString)+strlen(outName)+6);
    timepointsFile[0] = '\0';
    outputFile[0] = '\0';

    strcat(outputFile,outputDirString);
    #ifdef _WIN32
    strcat(outputFile,"\\");
    #else
    strcat(outputFile,"/");
    #endif

	strcat(outputFile,outName);
	strcat(outputFile,".ims");

    strcat(timepointsFile,outputDirString);
    #ifdef _WIN32
    strcat(timepointsFile,"\\timepoints.txt");
    #else
    strcat(timepointsFile,"/timepoints.txt");
    #endif


    char* cPatt = NULL;
    char* cFolder = NULL;
    uint64_t nChannels = 0;
    uint64_t nTimepoints = 0;
	uint64_t lnTimePoints = 0;
    char* pDup = strdup(patterns);
    char* fDup = NULL;
    char *saveptr1 = NULL, *saveptr2 = NULL;
    cPatt = strtok_r(pDup,delim,&saveptr1);
    while(cPatt){
        free(fDup);
        char* fDup = strdup(folderNames);
        cFolder = strtok_r(fDup,delim,&saveptr2);
        lnTimePoints = nTimepoints;
		while(cFolder){
            listFiles(cFolder,cPatt,0,0,NULL,1,(int64_t*)&nTimepoints);
            cFolder = strtok_r(NULL,delim,&saveptr2);
        }
		if(lnTimePoints == nTimepoints){
            printf("ERROR: Given pattern \"%s\" was not found in any of the provided folders.\n",cPatt);
            return;
        }
        nChannels++;
        cPatt = strtok_r(NULL,delim,&saveptr1);
    }
    if(channels == -1) channels = nChannels;
    if(channels != nChannels){
        printf("Number of channels given \"%lld\" is incorrect",channels);
        return;
    }

    if(timepoints == -1) timepoints = nTimepoints/nChannels;
    if(nTimepoints/nChannels != timepoints){
        printf("Number of timepoints given \"%lld\" is incorrect",timepoints);
        return;
    }

    // Save Files and Patterns
    char** files = (char**)malloc(sizeof(char*)*timepoints*channels);
    char** patternArr = (char**)malloc(sizeof(char*)*channels);

    free(pDup);
    free(fDup);

    pDup = strdup(patterns);
    fDup = NULL;
    char* saveptr3 = NULL, *saveptr4 = NULL;
    cPatt = strtok_r(pDup,delim,&saveptr3);
    uint64_t cPattN = 0;
    while(cPatt){
        patternArr[cPattN] = strdup(cPatt);
        free(fDup);
        char* fDup = strdup(folderNames);
        cFolder = strtok_r(fDup,delim,&saveptr4);
        while(cFolder){
            listFiles(cFolder,cPatt,1,cPattN,files,0,&timepoints);
            cFolder = strtok_r(NULL,delim,&saveptr4);
        }
        cPatt = strtok_r(NULL,delim,&saveptr3);
        cPattN++;
    }

    // Sort Channel Patterns
    for(int i = 0; i < channels; i++){
        qsort(&files[i*timepoints], timepoints, sizeof(char*), compare);
    }

    FILE* fp;
    fp = fopen(timepointsFile,"a");
    fprintf(fp,"Folder Names: %s\n",folderNames);
    for(int i = 0; i < timepoints*channels; i++){
        if(!(i%timepoints)) fprintf(fp,"Channel Pattern: %s\n",patternArr[i/timepoints]);
        fprintf(fp,"%s\n",files[i]);
    }
    fclose(fp);

    /*
    for(int i = 0; i < timepoints*channels; i++){
        printf("%s\n", files[i]);
    }
    return;
    */

    uint64_t shapeX = 0;
    uint64_t shapeY = 0;
    uint64_t shapeZ = 0;
    uint64_t chunkXSize = 0;
    uint64_t chunkYSize = 0;
    uint64_t chunkZSize = 0;
    char dtype[4];
    char order;
    uint64_t bits = 1;
    std::string cname;
	zarr Zarr;

    if(!strcmp(reader,"tiff")){
        uint64_t *size = getImageSize(files[0]);
        shapeX = size[0];
        shapeY = size[1];
        shapeZ = size[2];
        chunkXSize = 64;
        chunkYSize = 64;
        chunkZSize = 78;
        bits = getDataType(files[0]);
    }
    else{
		Zarr = zarr(files[0]);
		chunkXSize = Zarr.get_chunks(0);
		chunkYSize = Zarr.get_chunks(1);
		chunkZSize = Zarr.get_chunks(2);
		bits = Zarr.dtypeBytes()*8;
		order = Zarr.get_order()[0];
		shapeX = Zarr.get_shape(0);
		shapeY = Zarr.get_shape(1);
		shapeZ = Zarr.get_shape(2);
		cname = Zarr.get_cname();
        //setValuesFromJSON(files[0],&chunkXSize,&chunkYSize,&chunkZSize,dtype,&order,&shapeX,&shapeY,&shapeZ,&cname);
        //bits = dTypeToBits(dtype);
    }
	// Use custom chunk sizes if provided
    if(blockSizes){
        char* saveptrB = NULL;
        char* cB = NULL;
        cB = strtok_r(blockSizes,delim,&saveptrB);
        chunkXSize = strtof(cB,NULL);
        cB = strtok_r(NULL,delim,&saveptrB);
        chunkYSize = strtof(cB,NULL);
        cB = strtok_r(NULL,delim,&saveptrB);
        chunkZSize = strtof(cB,NULL);
    }
	printf("Block Sizes: %d,%d,%d\n",chunkXSize,chunkYSize,chunkZSize);
    // TODO: 64 BIT NOT YET SUPPORTED
    if(bits == 64){
        printf("64 bit data is not yet supported by the Imaris Converter");
        return;
    }

    bpConverterTypesC_DataType aDataType;
    switch(bits){
        case 8:
            aDataType = bpConverterTypesC_UInt8Type;
            break;
        case 16:
            aDataType = bpConverterTypesC_UInt16Type;
            break;
        case 32:
            aDataType = bpConverterTypesC_FloatType;
            break;
        case 64:
            //aDataType = bpConverterTypesC_DoubleType;
            break;
    }

    bpConverterTypesC_Size5D aImageSize;

	// Bounding Box
	if(crop){
		char* saveptrV = NULL;
        char* cV = NULL;
        cV = strtok_r(boundingBoxString,delim,&saveptrV);
        boundingBox[0] = strtof(cV,NULL);
		for(uint8_t i = 1; i < 6; i++){
			cV = strtok_r(NULL,delim,&saveptrV);
        	boundingBox[i] = strtof(cV,NULL);
		}
		shapeX = boundingBox[3]-boundingBox[0];
		shapeY = boundingBox[4]-boundingBox[1];
		shapeZ = boundingBox[5]-boundingBox[2];
		printf("Bounding Box: %d %d %d %d %d %d\n",boundingBox[0],boundingBox[1],boundingBox[2],boundingBox[3],boundingBox[4],boundingBox[5]);
		printf("New Shape: %d %d %d\n",shapeX,shapeY,shapeZ);
	}
	// If there is no bounding box then read the entire image
	else{
		boundingBox[3] = shapeX;
		boundingBox[4] = shapeY;
		boundingBox[5] = shapeZ;
	}

    if (shapeX % chunkXSize != 0){
        aImageSize.mValueX = shapeX+(shapeX % chunkXSize);
    }
    else aImageSize.mValueX = shapeX;
    if (shapeY % chunkYSize != 0){
        aImageSize.mValueY = shapeY+(shapeY % chunkYSize);
    }
    else aImageSize.mValueY = shapeY;
    if (shapeZ % chunkZSize != 0){
        aImageSize.mValueZ = shapeZ+(shapeZ % chunkZSize);
    }
    else aImageSize.mValueZ = shapeZ;

    aImageSize.mValueT = timepoints;
    aImageSize.mValueC = channels;


    vVoxelSizeX = 1.08f;
    vVoxelSizeY = 1.08f;
    vVoxelSizeZ = 1.08f;
    if(vVoxelSizes){
        char* saveptrV = NULL;
        char* cV = NULL;
        cV = strtok_r(vVoxelSizes,delim,&saveptrV);
        vVoxelSizeX = strtof(cV,NULL);
        cV = strtok_r(NULL,delim,&saveptrV);
        vVoxelSizeY = strtof(cV,NULL);
        cV = strtok_r(NULL,delim,&saveptrV);
        vVoxelSizeZ = strtof(cV,NULL);
    }

    printf("Voxel Sizes: %f,%f,%f\n",vVoxelSizeX,vVoxelSizeY,vVoxelSizeZ);


    bpConverterTypesC_ImageExtent aImageExtent = {
        0, 0, 0,
        aImageSize.mValueX * vVoxelSizeX,
        aImageSize.mValueY * vVoxelSizeY,
        aImageSize.mValueZ * vVoxelSizeZ
    };
    bpConverterTypesC_Size5D aSample = { 1, 1, 1, 1, 1 };

    bpConverterTypesC_DimensionSequence5D aDimensionSequence = {
        bpConverterTypesC_DimensionX,
        bpConverterTypesC_DimensionY,
        bpConverterTypesC_DimensionZ,
        bpConverterTypesC_DimensionC,
        bpConverterTypesC_DimensionT
    };

    bpConverterTypesC_Size5D aBlockSize = {
        chunkXSize, chunkYSize, chunkZSize, 1, 1
    };

    bpConverterTypesC_Options aOptions;
    aOptions.mThumbnailSizeXY = 256;
    aOptions.mFlipDimensionX = false;
    aOptions.mFlipDimensionY = false;
    aOptions.mFlipDimensionZ = false;
    aOptions.mForceFileBlockSizeZ1 = false;
    aOptions.mEnableLogProgress = true;
    aOptions.mNumberOfThreads = omp_get_max_threads();
    aOptions.mCompressionAlgorithmType = eCompressionAlgorithmLShuffleLZ4;

    bpConverterTypesC_String aApplicationName = "TestC";
    bpConverterTypesC_String aApplicationVersion = "1.0.0";
    bpConverterTypesC_ProgressCallback aProgressCallback = ProgressCallback;

    bpCallbackData aCallbackUserData;
    aCallbackUserData.mImageIndex = 0;
    aCallbackUserData.mProgress = -5;

    bpImageConverterCPtr vConverter = bpImageConverterC_Create(
                aDataType, &aImageSize, &aSample,
                &aDimensionSequence, &aBlockSize,
                outputFile, &aOptions,
                aApplicationName, aApplicationVersion,
                aProgressCallback, &aCallbackUserData
                );
    CheckErrors(vConverter);

    unsigned long long vBlockSize =
            (unsigned long long)aBlockSize.mValueX *
            aBlockSize.mValueY * aBlockSize.mValueZ *
            aBlockSize.mValueC * aBlockSize.mValueT;

    unsigned int vNBlocksX = NumBlocks(aImageSize.mValueX, aBlockSize.mValueX);
    unsigned int vNBlocksY = NumBlocks(aImageSize.mValueY, aBlockSize.mValueY);
    unsigned int vNBlocksZ = NumBlocks(aImageSize.mValueZ, aBlockSize.mValueZ);
    unsigned int vNBlocksC = NumBlocks(aImageSize.mValueC, aBlockSize.mValueC);
    unsigned int vNBlocksT = NumBlocks(aImageSize.mValueT, aBlockSize.mValueT);

    bpConverterTypesC_Index5D aBlockIndex = {
        0, 0, 0, 0, 0
    };
    //uint16_t* buffer = (uint16_t*)malloc(vBlockSize*sizeof(uint16_t));
    void* buffer = malloc(vBlockSize*(bits/8));

    //uint64_t skipN = 0;
    for (uint64_t vC = 0; vC < vNBlocksC; ++vC) {
        aBlockIndex.mValueC = vC;
        //skipN = 0;
        for (uint64_t vT = 0; vT < vNBlocksT; ++vT) {
            aBlockIndex.mValueT = vT;
            char fileName[10000];
            //sprintf(fileName,"D:\\imarisData\\ch%llu_z0_t%04llu.tif",vC,vT+1);
            //sprintf(fileName,files,vC,vT+1);
            //sprintf(fileName,"X:\\Data\\20220115_Korra_LLCPK_LFOV_0p1PSAmpKan\\DSR_RotateXY_Cropped_filtered_combined\\ch%llu_z0_t%04llu.tif",vC,vT+1+skipN);
            //skipN += 7;
            sprintf(fileName,files[vT+(vC*vNBlocksT)]);
            printf("Current File: %s\n",fileName);
            fflush(stdout);
            //if(strcmp("i", "i"))
            //uint16_t* vData = (uint16_t*)readTiffParallelWrapper(fileName);
            void* vData;
            if(!strcmp(reader,"tiff")) vData = readTiffParallelWrapper(fileName);
            else vData = readZarrParallelHelper(fileName,boundingBox[0],boundingBox[1],boundingBox[2],boundingBox[3],boundingBox[4],boundingBox[5],0);

            for (uint64_t vZ = 0; vZ < vNBlocksZ; ++vZ) {
                aBlockIndex.mValueZ = vZ;
                for (uint64_t vY = 0; vY < vNBlocksY; ++vY) {
                    aBlockIndex.mValueY = vY;
                    for (uint64_t vX = 0; vX < vNBlocksX; ++vX) {
                        aBlockIndex.mValueX = vX;
                        #pragma omp parallel for collapse(3)
                        for(uint64_t z = 0; z < aBlockSize.mValueZ; z++){
                            for(uint64_t y = 0; y < aBlockSize.mValueY; y++){
                                for(uint64_t x = 0; x < aBlockSize.mValueX; x++){
                                    if(x+(vX*aBlockSize.mValueX) >= shapeX || y+(vY*aBlockSize.mValueY) >= shapeY || z+(vZ*aBlockSize.mValueZ) >= shapeZ){
                                        switch(bits){
                                            case 8:
                                                ((uint8_t*)buffer)[x+(y*aBlockSize.mValueX)+(z*aBlockSize.mValueX*aBlockSize.mValueY)] = 0;
                                                break;
                                            case 16:
                                                ((uint16_t*)buffer)[x+(y*aBlockSize.mValueX)+(z*aBlockSize.mValueX*aBlockSize.mValueY)] = 0;
                                                break;
                                            case 32:
                                                ((float*)buffer)[x+(y*aBlockSize.mValueX)+(z*aBlockSize.mValueX*aBlockSize.mValueY)] = 0;
                                                break;
                                            case 64:
                                                ((double*)buffer)[x+(y*aBlockSize.mValueX)+(z*aBlockSize.mValueX*aBlockSize.mValueY)] = 0;
                                                break;
                                        }
                                        //buffer[x+(y*aBlockSize.mValueX)+(z*aBlockSize.mValueX*aBlockSize.mValueY)] = 0;
                                    }
                                    else{
                                        uint64_t nX = (x+(vX*aBlockSize.mValueX));
                                        uint64_t nY = ((y+(vY*aBlockSize.mValueY))*shapeX);
                                        uint64_t nZ = ((z+(vZ*aBlockSize.mValueZ))*shapeX*shapeY);
                                        switch(bits){
                                            case 8:
                                                ((uint8_t*)buffer)[x+(y*aBlockSize.mValueX)+(z*aBlockSize.mValueX*aBlockSize.mValueY)] = ((uint8_t*)vData)[nX+nY+nZ];
                                                break;
                                            case 16:
                                                ((uint16_t*)buffer)[x+(y*aBlockSize.mValueX)+(z*aBlockSize.mValueX*aBlockSize.mValueY)] = ((uint16_t*)vData)[nX+nY+nZ];
                                                break;
                                            case 32:
                                                ((float*)buffer)[x+(y*aBlockSize.mValueX)+(z*aBlockSize.mValueX*aBlockSize.mValueY)] = ((float*)vData)[nX+nY+nZ];
                                                break;
                                            case 64:
                                                ((double*)buffer)[x+(y*aBlockSize.mValueX)+(z*aBlockSize.mValueX*aBlockSize.mValueY)] = ((double*)vData)[nX+nY+nZ];
                                                break;
                                        }
                                        //buffer[x+(y*aBlockSize.mValueX)+(z*aBlockSize.mValueX*aBlockSize.mValueY)] = vData[nX+nY+nZ];
                                    }
                                }
                            }
                        }
                        switch(bits){
                            case 8:
                                bpImageConverterC_CopyBlockUInt8(vConverter, (uint8_t*)buffer, &aBlockIndex);
                                break;
                            case 16:
                                bpImageConverterC_CopyBlockUInt16(vConverter, (uint16_t*)buffer, &aBlockIndex);
                                break;
                            case 32:
                                bpImageConverterC_CopyBlockFloat(vConverter, (float*)buffer, &aBlockIndex);
                                break;
                            case 64:
                                //bpImageConverterC_CopyBlockDouble(vConverter, (double*)buffer, &aBlockIndex);
                                break;
                        }
                        //bpImageConverterC_CopyBlockUInt16(vConverter, buffer, &aBlockIndex);
                        CheckErrors(vConverter);
                    }
                }
            }
            free(vData);
        }
    }
    free(buffer);

    unsigned int vNumberOfOtherSections = 1; // Image
    unsigned int vNumberOfSections = vNumberOfOtherSections + aImageSize.mValueC;
    bpConverterTypesC_ParameterSection* vParameterSections = (bpConverterTypesC_ParameterSection*)malloc(vNumberOfSections * sizeof(bpConverterTypesC_ParameterSection));

    bpConverterTypesC_Parameter vUnitParameter = { "Unit", "um" };
    bpConverterTypesC_ParameterSection* vImageSection = &vParameterSections[0];
    vImageSection->mName = "Image";
    vImageSection->mValuesCount = 1;
    vImageSection->mValues = &vUnitParameter;

    char vChannelNamesBuffer[1024]; // will this be enough?
    char* vChannelNameBuffer = vChannelNamesBuffer;

    unsigned int vNumberOfParametersPerChannel = 3;
    bpConverterTypesC_Parameter* vChannelParameters = (bpConverterTypesC_Parameter*)malloc(aImageSize.mValueC * vNumberOfParametersPerChannel * sizeof(bpConverterTypesC_Parameter));
    for (unsigned int vC = 0; vC < aImageSize.mValueC; ++vC) {
        bpConverterTypesC_Parameter* vThisChannelParameters = &vChannelParameters[vNumberOfParametersPerChannel * vC];
        vThisChannelParameters[0].mName = "Name";
        vThisChannelParameters[0].mValue = vC == 0 ? "First channel" : vC == 1 ? "Second channel" : vC == 2 ? "Third channel" : "Other channel";
        vThisChannelParameters[1].mName = "LSMEmissionWavelength";
        vThisChannelParameters[1].mValue = "700";
        vThisChannelParameters[2].mName = "OtherChannelParameter";
        vThisChannelParameters[2].mValue = "OtherChannelValue";
        bpConverterTypesC_ParameterSection* vChannelSection = &vParameterSections[vNumberOfOtherSections + vC];
        int vChannelNameLength = sprintf(vChannelNameBuffer, "Channel %i", vC);
        vChannelSection->mName = vChannelNameBuffer;
        vChannelNameBuffer += vChannelNameLength + 1;
        vChannelSection->mValues = vThisChannelParameters;
        vChannelSection->mValuesCount = vNumberOfParametersPerChannel;
    }

    bpConverterTypesC_Parameters aParameters;
    aParameters.mValuesCount = vNumberOfSections;
    aParameters.mValues = vParameterSections;

    bpConverterTypesC_TimeInfo* vTimeInfos = (bpConverterTypesC_TimeInfo*)malloc(aImageSize.mValueT * sizeof(bpConverterTypesC_TimeInfo));
    for (unsigned int vT = 0; vT < aImageSize.mValueT; ++vT) {
        vTimeInfos[vT].mJulianDay = 2458885; // 5 feb 2020
        unsigned long long vSeconds = vT + 4 + 60 * (27 + 60 * 15); // 3:27.04 PM + 1 sec per time point
        vTimeInfos[vT].mNanosecondsOfDay = vSeconds * 1000000000;
    }
    bpConverterTypesC_TimeInfos aTimeInfoPerTimePoint;
    aTimeInfoPerTimePoint.mValuesCount = aImageSize.mValueT;
    aTimeInfoPerTimePoint.mValues = vTimeInfos;

    bpConverterTypesC_ColorInfo* vColorInfos = (bpConverterTypesC_ColorInfo*)malloc(aImageSize.mValueC * sizeof(bpConverterTypesC_ColorInfo));
    for (unsigned int vC = 0; vC < aImageSize.mValueC; ++vC) {
        bpConverterTypesC_ColorInfo* vColor = &vColorInfos[vC];
        vColor->mIsBaseColorMode = true;
        vColor->mBaseColor.mRed = (vC % 3) == 0 ? 1 : 0;
        vColor->mBaseColor.mGreen = (vC % 3) == 1 ? 1 : 0;
        vColor->mBaseColor.mBlue = (vC % 3) == 2 ? 1 : 0;
        vColor->mBaseColor.mAlpha = 1;
        vColor->mColorTableSize = 0;
        vColor->mOpacity = 0;
        vColor->mRangeMin = 0;
        vColor->mRangeMax = 255;
        vColor->mGammaCorrection = 1;
    }
    bpConverterTypesC_ColorInfos aColorInfoPerChannel;
    aColorInfoPerChannel.mValuesCount = aImageSize.mValueC;
    aColorInfoPerChannel.mValues = vColorInfos;

    bool aAutoAdjustColorRange = true;

    bpImageConverterC_Finish(vConverter,
                             &aImageExtent, &aParameters, &aTimeInfoPerTimePoint,
                             &aColorInfoPerChannel, aAutoAdjustColorRange);
    CheckErrors(vConverter);

    free(vTimeInfos);
    free(vColorInfos);

    free(vParameterSections);
    free(vChannelParameters);

    bpImageConverterC_Destroy(vConverter);
}

int main(int argc, char **argv)
{
    //listFiles("C:/Users/Matt/Desktop/test","CamA");
    //return 0;

    convertToImaris(argc,argv);
    return 0;
}
