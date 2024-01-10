#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
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

void mkdirRecursive(const char *dir) {
    char tmp[8192];
    char *p = NULL;
    size_t len;
    #ifdef _WIN32
    char fileSep = '\\';
    #else
    char fileSep = '/';
    #endif
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
            // Check if the image was written by tifffile
            char* tiffSoftware = NULL;
            if(TIFFGetField(tif, TIFFTAG_SOFTWARE, &tiffSoftware)){
                if(!strcmp(tiffSoftware,"tifffile.py")) return 0;
            }
            uint16_t compressed = 1;
            TIFFGetField(tif, TIFFTAG_COMPRESSION, &compressed);
            if(compressed != 1) return 0;
            else return 1;
        }
    }
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
                nZ+=7;
                char* temp;
                return strtol(nZ,&temp,10);
            }
        }
    }
    return 0;
}


// The following are unused maybe?

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

