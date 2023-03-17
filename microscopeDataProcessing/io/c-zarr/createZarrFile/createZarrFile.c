#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <errno.h>
#include <cjson/cJSON.h>
#include <omp.h>
#ifdef _WIN32
#include <stdarg.h> 
#include <sys/time.h>
#else
#include <uuid/uuid.h>
#endif
#include <sys/stat.h>
#include "mex.h"

//compile
//mex -v COPTIMFLAGS="-O3 -DNDEBUG" LDOPTIMFLAGS="-Wl',-rpath='''$ORIGIN'''' -O3 -DNDEBUG" CFLAGS='$CFLAGS -O3 -fopenmp' LDFLAGS='$LDFLAGS -O3 -fopenmp' '-I/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/include/' -L'/global/home/groups/software/sl-7.x86_64/modules/cJSON/1.7.15/lib64' -lcjson -luuid createZarrFile.c

// Handle the tilde character in filenames on Linux/Mac
#ifndef _WIN32
#include <wordexp.h>
char* expandTilde(char* path) {
    wordexp_t expPath;
    wordexp(path, &expPath, 0);
    return expPath.we_wordv[0];
}
#endif

const char fileSep =
#ifdef _WIN32
    '\\';
#else
    '/';
#endif
    
#ifdef _WIN32
char* strndup (const char *s, size_t n)
{
  size_t len = strnlen (s, n);
  char *new = (char *) malloc (len + 1);
  if (new == NULL)
    return NULL;
  new[len] = '\0';
  return (char *) memcpy (new, s, len);
}
int _vscprintf_so(const char * format, va_list pargs) {
    int retval;
    va_list argcopy;
    va_copy(argcopy, pargs);
    retval = vsnprintf(NULL, 0, format, argcopy);
    va_end(argcopy);
    return retval;
}

int vasprintf(char **strp, const char *fmt, va_list ap) {
    int len = _vscprintf_so(fmt, ap);
    if (len == -1) return -1;
    char *str = malloc((size_t) len + 1);
    if (!str) return -1;
    int r = vsnprintf(str, len + 1, fmt, ap); /* "secure" version of vsprintf */
    if (r == -1) return free(str), -1;
    *strp = str;
    return r;
}

int asprintf(char *strp[], const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    int r = vasprintf(strp, fmt, ap);
    va_end(ap);
    return r;
}
#endif

static void mkdirRecursive(const char *dir) {
    char tmp[8192];
    char *p = NULL;
    size_t len;
    #ifdef _WIN32
    char fileSep = '\\';
    #else
    char fileSep = '/';
    #endif
    int status;
    snprintf(tmp, sizeof(tmp),"%s",dir);
    len = strlen(tmp);
    if (tmp[len - 1] == fileSep)
        tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++){
        if (*p == fileSep) {
            *p = 0;

            #ifdef _WIN32
            mkdir(tmp);
            #else
            mkdir(tmp, 0775);
            #endif

            chmod(tmp, 0775);
            *p = fileSep;
        }
    }
    #ifdef _WIN32
    mkdir(tmp);
    #else
    mkdir(tmp, 0775);
    #endif
    chmod(tmp, 0775);
}

void setJSONValues(char* fileName,uint64_t *chunkXSize,uint64_t *chunkYSize,uint64_t *chunkZSize,char* dtype,char* order,uint64_t *shapeX,uint64_t *shapeY,uint64_t *shapeZ, char* cname, uint64_t level, uint64_t *subfolderSizeX, uint64_t *subfolderSizeY, uint64_t *subfolderSizeZ){
    
    // Overflows for ints greater than 32 bit for chunkSizes and shape
    cJSON* zArray = cJSON_CreateObject();
    const int chunkSizes[3] = {*chunkXSize,*chunkYSize,*chunkZSize};
    cJSON* chunks = cJSON_CreateIntArray(chunkSizes, 3);
    cJSON_AddItemToObject(zArray, "chunks", chunks);

    cJSON* compressor = cJSON_CreateObject();
    cJSON_AddItemToObject(zArray, "compressor", compressor);
    if(!strcmp(cname,"lz4") || !strcmp(cname,"blosclz") || !strcmp(cname,"lz4hc") || !strcmp(cname,"zlib") || !strcmp(cname,"zstd")){
        cJSON_AddNumberToObject(compressor, "blocksize", 0);
        cJSON_AddNumberToObject(compressor, "clevel", level);
        cJSON_AddStringToObject(compressor, "cname", cname);
        cJSON_AddStringToObject(compressor, "id", "blosc");
        cJSON_AddNumberToObject(compressor, "shuffle", 1);
    }
    else if(!strcmp(cname,"gzip")){
        cJSON_AddStringToObject(compressor, "id", cname);
        cJSON_AddNumberToObject(compressor, "level", level);
    }
    else mexErrMsgIdAndTxt("zarr:zarrayError","Compressor: \"%s\" is not currently supported\n",cname);

    cJSON_AddStringToObject(zArray, "dtype", dtype);
    cJSON_AddNumberToObject(zArray, "fill_value", 0);
    cJSON_AddNullToObject(zArray, "filters");
    cJSON_AddStringToObject(zArray, "order", order);

    const int shapeSizes[3] = {*shapeX,*shapeY,*shapeZ};
    cJSON* shape = cJSON_CreateIntArray(shapeSizes, 3);
    cJSON_AddItemToObject(zArray, "shape", shape);
    cJSON_AddNumberToObject(zArray, "zarr_format", 2);

    const int subfolderSizeSizes[3] = {*subfolderSizeX,*subfolderSizeY,*subfolderSizeZ};
    cJSON* subfolderSize = cJSON_CreateIntArray(subfolderSizeSizes, 3);
    cJSON_AddItemToObject(zArray, "subfolders", subfolderSize);
    
    uint64_t uuidLen;
    #ifdef _WIN32
    uuidLen = 5;
    char *uuid = malloc(uuidLen+1);
    char *seedArr = malloc(1000);
    struct timeval cSeed;
    gettimeofday(&cSeed,NULL);
    int nChars = sprintf(seedArr,"%d%d",cSeed.tv_sec,cSeed.tv_usec);
    int aSeed = 0;
    char* ptr;
    if(nChars > 9)
        aSeed = strtol(seedArr+nChars-10, &ptr, 9);
    else aSeed = strtol(seedArr, &ptr, 9);
    srand(aSeed);
    sprintf(uuid,"%.5d",rand() % 99999);
    #else
    uuidLen = 36;
    uuid_t binuuid;
    uuid_generate_random(binuuid);
    char *uuid = malloc(uuidLen+1);
    uuid_unparse(binuuid, uuid);
    #endif

    char* zArrayS = ".zarray";
    char* fnFull = (char*)malloc(strlen(fileName)+8+uuidLen+1);
    fnFull[0] = '\0';
    char fileSepS[2];
    fileSepS[0] = fileSep;
    fileSepS[1] = '\0';

    strcat(fnFull,fileName);
    strcat(fnFull,fileSepS);
    strcat(fnFull,zArrayS);
    char* fileNameFinal = strndup(fnFull,strlen(fileName)+1+8);
    strcat(fnFull,uuid);
    
    char* string = cJSON_Print(zArray);

    FILE *fileptr = fopen(fnFull, "w+");
    if(!fileptr) mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s\n",fnFull);
    fprintf(fileptr,"%s",string);
    fclose(fileptr);
    
    rename(fnFull,fileNameFinal);
    //file = fopen(fileNameFinal, "r");
    //if(!file) rename(fileName,fileNameFinal);
    
    cJSON_Delete(zArray);
    free(fnFull);
    free(uuid);
    free(fileNameFinal);
}

uint64_t fastCeilDiv(uint64_t num, uint64_t denom){
    return 1 + ((num - 1) / denom);
}

void createSubfolders(char* folderName, uint64_t shapeX, uint64_t shapeY, uint64_t shapeZ, uint64_t chunkXSize, uint64_t chunkYSize, uint64_t chunkZSize, uint64_t subfolderSizeX, uint64_t subfolderSizeY, uint64_t subfolderSizeZ){
    if(subfolderSizeX == 0 && subfolderSizeY == 0 && subfolderSizeZ == 0) return;
    uint64_t nChunksX = fastCeilDiv(shapeX,chunkXSize);
    uint64_t nChunksY = fastCeilDiv(shapeY,chunkYSize);
    uint64_t nChunksZ = fastCeilDiv(shapeZ,chunkZSize);

    uint64_t nSubfoldersX = 1;
    uint64_t nSubfoldersY = 1;
    uint64_t nSubfoldersZ = 1;

    if(subfolderSizeX > 0) nSubfoldersX = fastCeilDiv(nChunksX,subfolderSizeX);
    if(subfolderSizeY > 0) nSubfoldersY = fastCeilDiv(nChunksY,subfolderSizeY);
    if(subfolderSizeZ > 0) nSubfoldersZ = fastCeilDiv(nChunksZ,subfolderSizeZ);
    
    #pragma omp parallel for collapse(3)
    for(uint64_t x = 0; x < nSubfoldersX; x++){
        for(uint64_t y = 0; y < nSubfoldersY; y++){
            for(uint64_t z = 0; z < nSubfoldersZ; z++){
                char* currName = NULL;
                asprintf(&currName,"%s/%llu_%llu_%llu",folderName,x,y,z);
                mkdirRecursive(currName);
                free(currName);
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    // Needed variables for .zarray file
    uint64_t chunkXSize = 256;
    uint64_t chunkYSize = 256;
    uint64_t chunkZSize = 256;
    char* cname = "lz4";
    char* dtype = "<u2";
    char* order = "F";
    uint64_t shapeX = 0;
    uint64_t shapeY = 0;
    uint64_t shapeZ = 0;
    uint64_t level = 5;
    uint64_t subfolderSizeX = 0;
    uint64_t subfolderSizeY = 0;
    uint64_t subfolderSizeZ = 0;
    if(nrhs < 1) mexErrMsgIdAndTxt("zarr:inputError","This functions requires at least 1 argument\n");
    if(nrhs == 2) mexErrMsgIdAndTxt("zarr:inputError","This functions does not accept only 2 arguments\n");
    if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("zarr:inputError","The first argument must be a string\n");
    char* folderName = mxArrayToString(prhs[0]);

    // Handle the tilde character in filenames on Linux/Mac
    #ifndef _WIN32
    if(strchr(folderName,'~')) folderName = expandTilde(folderName);
    #endif

    for(int i = 1; i < nrhs; i+=2){
        if(i+1 == nrhs) mexErrMsgIdAndTxt("zarr:inputError","Mismatched argument pair for input number %d\n");
        if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("zarr:inputError","The argument in input location %d is not a string\n",i+1);
        char* currInput = mxArrayToString(prhs[i]);

        if(!strcmp(currInput,"chunks")){
            if(mxGetN(prhs[i+1]) != 3) mexErrMsgIdAndTxt("zarr:inputError","chunks must be an array of 3 numbers\n");
            chunkXSize = (uint64_t)*(mxGetPr(prhs[i+1]));
            chunkYSize = (uint64_t)*((mxGetPr(prhs[i+1])+1));
            chunkZSize = (uint64_t)*((mxGetPr(prhs[i+1])+2));
        }
        else if(!strcmp(currInput,"cname")){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","cname must be a string\n");
            cname = mxArrayToString(prhs[i+1]);
            if(!strcmp(currInput,cname) || !strcmp(currInput,cname)) mexErrMsgIdAndTxt("zarr:inputError","cname supported compressors are lz4 and zlib\n");
        }
        else if(!strcmp(currInput,"dtype")){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","dtype must be a string\n");
            dtype = mxArrayToString(prhs[i+1]);
            int dtypeLen = strlen(dtype);

            if(dtypeLen == 2){
                char* dtypeCopy = strdup(dtype);
                dtype = (char*)malloc(4);
                dtype[0] = '<';
                dtype[1] = dtypeCopy[0];
                dtype[2] = dtypeCopy[1];
                dtype[3] = '\0';
                free(dtypeCopy);
            }
            else if(dtypeLen != 3) mexErrMsgIdAndTxt("zarr:inputError","dtype must be of length 2 or 3\n");
            
            if(dtype[0] != '<' && dtype[0] != '>'){
                mexErrMsgIdAndTxt("zarr:inputError","The first character of dtype must be \"<\" or \">\"\n");
            }

            if(dtype[1] == 'u' || dtype[1] == 'f'){
                if(dtype[1] == 'u'){
                    if(dtype[2] != '1' && dtype[2] != '2'){
                        mexErrMsgIdAndTxt("zarr:inputError","For dtype u, the third character of dtype must be \"1\" or \"2\"\n");
                    }
                }
                else{
                    if(dtype[2] != '4' && dtype[2] != '8'){
                        mexErrMsgIdAndTxt("zarr:inputError","For dtype f, the third character of dtype must be \"4\" or \"8\"\n");
                    }
                }
            }
            else mexErrMsgIdAndTxt("zarr:inputError","The second character of dtype must be \"u\" or \"f\"\n");
        }
        else if(!strcmp(currInput,"order")){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","order must be a string\n");
            order = mxArrayToString(prhs[i+1]);
            if(strlen(order) > 1) mexErrMsgIdAndTxt("zarr:inputError","order must be of length 1\n");
            if(order[0] != 'F' || order[0] != 'C' || order[0] != 'f' || order[0] != 'c'){
                if(order[0] == 'f') order = "F";
                else if(order [0] == 'c') order = "C";
            }
            else mexErrMsgIdAndTxt("zarr:inputError","order must be \"F\" or \"C\"\n");
        }
        else if(!strcmp(currInput,"shape")){
            if(mxGetN(prhs[i+1]) != 3) mexErrMsgIdAndTxt("zarr:inputError","shape must be an array of 3 numbers\n");
            shapeX = (uint64_t)*(mxGetPr(prhs[i+1]));
            shapeY = (uint64_t)*((mxGetPr(prhs[i+1])+1));
            shapeZ = (uint64_t)*((mxGetPr(prhs[i+1])+2));
        }
        else if(!strcmp(currInput,"level")){
            if(!mxIsNumeric(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","level must be a numerical value\n");
            level = (uint64_t)*(mxGetPr(prhs[i+1]));
        }
        else if(!strcmp(currInput,"subfolders")){
            if(mxGetN(prhs[i+1]) != 3) mexErrMsgIdAndTxt("zarr:inputError","subfolders must be an array of 3 numbers\n");
            subfolderSizeX = (uint64_t)*(mxGetPr(prhs[i+1]));
            subfolderSizeY = (uint64_t)*((mxGetPr(prhs[i+1])+1));
            subfolderSizeZ = (uint64_t)*((mxGetPr(prhs[i+1])+2));
        }
        else{
            mexErrMsgIdAndTxt("zarr:inputError","The argument \"%s\" does not match the name of any supported input name.\n \
                              Currently Supported Names: chunks, cname, dtype, order, shape\n",currInput);
        }
    }
    



    char* zArray = ".zarray";
    char* fnFull = (char*)malloc(strlen(folderName)+9);
    fnFull[0] = '\0';
    char fileSepS[2];
    fileSepS[0] = '/';
    fileSepS[1] = '\0';
    
    strcat(fnFull,folderName);
    strcat(fnFull,fileSepS);
    strcat(fnFull,zArray);
    
    FILE* f = fopen(fnFull,"r");
    if(f) fclose(f);
    else{
        mkdirRecursive(folderName);
        chmod(folderName, 0775);
    }
    createSubfolders(folderName,shapeX,shapeY,shapeZ,chunkXSize,chunkYSize,chunkZSize,subfolderSizeX,subfolderSizeY,subfolderSizeZ);
    setJSONValues(folderName,&chunkXSize,&chunkYSize,&chunkZSize,dtype,order,&shapeX,&shapeY,&shapeZ,cname,level,&subfolderSizeX,&subfolderSizeY,&subfolderSizeZ);


    /////////
}