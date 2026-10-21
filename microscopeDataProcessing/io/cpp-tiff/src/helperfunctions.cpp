#include <cstdint>
#include <cstring>
#include <sys/stat.h>
#include <string>
// For std::replace
#ifdef _WIN32
#include <algorithm>
#endif
#include <omp.h>
#include "tiffio.h"
#include "helperfunctions.h"
#ifdef __linux__
#include <fcntl.h>
#endif
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

// Open a file for walking its directory chain. Counting directories touches a
// few bytes per IFD (entry count + next-IFD offset), and the IFDs of a stack
// are spread through the whole file, so the walk is one dependent disk read per
// IFD. With libtiff's default memory-mapped open each of those reads is a page
// fault that also pulls in the ~128 KB read-ahead window around the IFD; on a
// cold cache that made a 10,000-directory file take ~8 s to count (1.3 GB read
// for ~0.1 MB of IFDs). On Linux, opening without mmap ("rm") plus a random-
// access hint reads only the touched pages: ~1.6 s, ~40 MB, warm walks
// unchanged. Windows and macOS measured slower this way (the extra read() per
// hop costs more than the mapped fault path there), so they keep the default
// mmap open. Only the count uses this; the readers keep their mmap handles.
static TIFF* openForDirectoryWalk(const char* fileName){
#ifdef __linux__
    TIFF* tif = TIFFOpen(fileName, "rm");
    if(!tif) return NULL;
  #ifdef POSIX_FADV_RANDOM
    posix_fadvise(TIFFFileno(tif), 0, 0, POSIX_FADV_RANDOM);
  #endif
    return tif;
#else
    return TIFFOpen(fileName, "r");
#endif
}

// The check for if it is an ImageJ image needs to be improved in the future
uint8_t isImageJIm(const char* fileName){
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif) return 0;
    char* tiffDesc = NULL;
    if(TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &tiffDesc)){
        if(strstr(tiffDesc, "ImageJ")){
            // The ImageJ fast layout is a single uncompressed IFD followed by
            // the raw slices. Read the compression from directory 0 first, then
            // probe for a second directory (one hop) instead of counting the
            // whole chain, which walked every IFD of a normal multi-page file.
            uint16_t compressed = 1;
            TIFFGetField(tif, TIFFTAG_COMPRESSION, &compressed);
            uint8_t multiDir = TIFFSetDirectory(tif,1) ? 1 : 0;
            TIFFClose(tif);
            if(multiDir || compressed != 1) return 0;
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
                // Parse before TIFFClose: tiffDesc points into libtiff's
                // directory storage, which TIFFClose frees. Reading it after
                // the close only worked by allocator luck on Linux/Windows and
                // returned garbage on macOS.
                char* temp;
                uint64_t z = strtol(nZ+7,&temp,10);
                TIFFClose(tif);
                return z;
            }
        }
    }
    TIFFClose(tif);
    return 0;
}

// The number of directories, found by stepping ahead 8x at a time and then
// bisecting; the handle starts on directory 0.
static uint32_t countDirectories(TIFF* tif){
    uint32_t s = 0, m = 0, t = 1;
    while(TIFFSetDirectory(tif,t)){
        s = t;
        t *= 8;
        if(s > t){
            t = 4294967295;
            printf("Number of slices > 1073741824\n");
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
    return s+1;
}

uint32_t getImageSizeZ(const char* fileName){
    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = openForDirectoryWalk(fileName);
    if(!tif){ printf("File \"%s\" cannot be opened\n",fileName); return 0; }
    const uint32_t z = countDirectories(tif);
    TIFFClose(tif);
    return z;
}

uint8_t getTiffInfo(const char* fileName, TiffInfo &info){
    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = openForDirectoryWalk(fileName);
    if(!tif) return 1;
    info = TiffInfo();
    uint32_t x = 1, y = 1, rowsPerStrip = 0;
    uint16_t bits = 1, spp = 1, planarConfig = PLANARCONFIG_CONTIG, sampleFormat = SAMPLEFORMAT_UINT, compression = COMPRESSION_NONE;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &x);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &y);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits);
    TIFFGetFieldDefaulted(tif, TIFFTAG_ROWSPERSTRIP, &rowsPerStrip);
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFGetFieldDefaulted(tif, TIFFTAG_PLANARCONFIG, &planarConfig);
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &sampleFormat);
    TIFFGetFieldDefaulted(tif, TIFFTAG_COMPRESSION, &compression);
    info.x = x; info.y = y; info.bits = bits; info.samples = spp;
    info.sampleFormat = sampleFormat; info.stripSize = rowsPerStrip; info.compression = compression;
    info.planar = spp > 1 && planarConfig == PLANARCONFIG_SEPARATE;
    info.tiled = TIFFIsTiled(tif);
    // The description belongs to directory 0, which the walk below frees, so
    // it is parsed first.
    bool imageJDescription = false;
    char* desc = NULL;
    if(TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &desc) && desc && strstr(desc, "ImageJ")){
        imageJDescription = true;
        const char* images = strstr(desc, "images=");
        if(images){
            char* end;
            info.imageJImages = strtol(images+7, &end, 10);
        }
    }
    info.directories = countDirectories(tif);
    TIFFClose(tif);
    info.z = info.directories;
    info.imageJ = imageJDescription && info.directories == 1 && compression == COMPRESSION_NONE && !info.tiled && !info.planar;
    if(info.imageJ && info.imageJImages) info.z = info.imageJImages;
    return 0;
}

int callThreads(){
    return omp_get_max_threads();
}

void setCallThreads(int n){
    if(n > 0) omp_set_num_threads(n);
}

uint64_t* getImageSize(const char* fileName){
    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif){ printf("File \"%s\" cannot be opened\n",fileName); return NULL; }

    uint64_t x = 1,y = 1,z = 1;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &x);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &y);
    z = getImageSizeZ(fileName);
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
    if(!tif){ printf("File \"%s\" cannot be opened\n",fileName); return 0; }

    uint64_t bits = 1;
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits);
    TIFFClose(tif);

    return bits;
}

// Returns the samples per pixel (1 = grayscale, 3 = RGB, 4 = RGBA).
// Uses the defaulted getter so files that omit the tag report 1.
uint64_t getSamplesPerPixel(const char* fileName){
    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif){ printf("File \"%s\" cannot be opened",fileName); return 1; }

    uint16_t spp = 1;
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFClose(tif);

    return spp;
}

// Returns the sample format (1 = unsigned int, 2 = signed int, 3 = IEEE float).
// Uses the defaulted getter so files that omit the tag report unsigned int.
uint64_t getSampleFormat(const char* fileName){
    TIFFSetWarningHandler(DummyHandler);
    TIFF* tif = TIFFOpen(fileName, "r");
    if(!tif){ printf("File \"%s\" cannot be opened",fileName); return 1; }

    uint16_t sf = SAMPLEFORMAT_UINT;
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &sf);
    TIFFClose(tif);

    return sf;
}

