#include <cstdint>
#include <cstring>
#include <sys/stat.h>
#include <string>
// For std::replace
#ifdef _WIN32
#include <algorithm>
#endif
#include "tiffio.h"
#include "helperfunctions.h"
// Handle the tilde character in filenames on Linux/Mac
#ifndef _WIN32
#include <wordexp.h>
char* expandTilde(char* path) {
    wordexp_t expPath;
    wordexp(path, &expPath, 0);
    return expPath.we_wordv[0];
}
#endif

bool folderExists(const std::string &folderName){
    struct stat info;

    if(stat(folderName.c_str(), &info) != 0)
        return false;
    else if(info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}

// Recursively make a directory
void mkdirRecursive(const char *dir) {
    if(folderExists(dir)) return;

    std::string dirPath(dir);
    if(!dirPath.size()) return;

    // Convert all \\ to / if on Windows
    #ifdef _WIN32
    std::replace(dirPath.begin(), dirPath.end(), '\\', '/');
    #endif

    // If there is a slash at the end, remove it
    if(dirPath.back() == '/') dirPath.pop_back();
    
    for(size_t i = 0; i < dirPath.size(); i++){
        if(dirPath[i] == '/'){
            dirPath[i] = '\0';

            #ifdef _WIN32
            mkdir(dirPath.c_str());
            #else
            mkdir(dirPath.c_str(), 0777);
            #endif

            dirPath[i] = '/';
        }
    }

    #ifdef _WIN32
    mkdir(dirPath.c_str());
    #else
    mkdir(dirPath.c_str(), 0777);
    #endif
}

void DummyHandler(const char* module, const char* fmt, va_list ap)
{
    // ignore errors and warnings
}

// The check for if it is an ImageJ image needs to be improved in the future
uint8_t isImageJIm(const char* fileName){
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) return 0;
    char* tiffDesc = NULL;
    if(TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &tiffDesc)){
        if(strstr(tiffDesc, "ImageJ")){
            uint64_t* size = getImageSize(fileName);
            if(size[2] > 1){
                if(TIFFSetDirectory(tif,1)){
                    free(size);
                    return 0;
                }
            }
            free(size);
            uint16_t compressed = 1;
            TIFFGetField(tif, TIFFTAG_COMPRESSION, &compressed);
            TIFFClose(tif);
            if(compressed != 1) return 0;
            else return 1;
        }
    }
    TIFFClose(tif);
    return 0;
}

uint64_t imageJImGetZ(const char* fileName){
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) return 0;
    char* tiffDesc = NULL;
    if(TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &tiffDesc)){
        if(strstr(tiffDesc, "ImageJ")){
            char* nZ = strstr(tiffDesc,"images=");
            if(nZ){
                TIFFClose(tif);
                nZ+=7;
                char* temp;
                return strtol(nZ,&temp,10);
            }
        }
    }
    TIFFClose(tif);
    return 0;
}

uint64_t* getImageSize(const char* fileName){
    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) printf("File \"%s\" cannot be opened",fileName);

    uint64_t x = 1,y = 1,z = 1;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &x);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &y);
    uint16_t s = 0, m = 0, t = 1;
    while(TIFFSetDirectory(tif,t)){
        s = t;
        t *= 8;
        if(s > t){
            t = 65535;
            printf("Number of slices > 32768\n");
            break;
        }
    }
    while(s != t){
        m = (s+t+1)/2;
        if(TIFFSetDirectory(tif,m)){
            s = m;
        }
        else{
            if(m > 0) t = m-1;
            else t = m;
        }
    }
    z = s+1;

    TIFFClose(tif);
    uint64_t* dims = (uint64_t*)malloc(3*sizeof(uint64_t));
    dims[0] = y;
    dims[1] = x;
    dims[2] = z;
    return dims;
}

// Returns number of bits the tiff file is.
uint64_t getDataType(const char* fileName){
    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) printf("File \"%s\" cannot be opened",fileName);

    uint64_t bits = 1;
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits);
    TIFFClose(tif);

    return bits;
}

