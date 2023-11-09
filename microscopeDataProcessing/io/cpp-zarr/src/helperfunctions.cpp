#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <omp.h>
#ifdef _WIN32
#include <stdarg.h> 
#include <sys/time.h>
#else
#include <uuid/uuid.h>
#endif
#include <sys/stat.h>
#include <filesystem>
#include "helperfunctions.h"
#include "zarr.h"

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
    
#ifdef _WIN32
char* strndup (const char *s, size_t n)
{
  size_t len = strnlen (s, n);
  char *newS = (char *) malloc (len + 1);
  if (newS == NULL)
    return NULL;
  newS[len] = '\0';
  return (char *) memcpy (newS, s, len);
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
    char *str = (char*)malloc((size_t) len + 1);
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

std::string generateUUID(){
    #ifdef _WIN32
    // uuid of length 5 for windows
    char uuidC[6];
    char *seedArr = (char*)malloc(1000);
    struct timeval cSeed;
    gettimeofday(&cSeed,NULL);
    int nChars = sprintf(seedArr,"%d%d",cSeed.tv_sec,cSeed.tv_usec);
    int aSeed = 0;
    char* ptr;
    if(nChars > 9)
        aSeed = strtol(seedArr+nChars-10, &ptr, 9);
    else aSeed = strtol(seedArr, &ptr, 9);
    srand(aSeed);
    sprintf(uuidC,"%.5d",rand() % 99999);
    free(seedArr);
    #else
    uuid_t binuuid;
    uuid_generate_random(binuuid);
    char uuidC[37];
    uuid_unparse(binuuid, uuidC);
    #endif
    return std::string(uuidC);
}

bool folderExists(const std::string &folderName){
    struct stat info;

    if(stat(folderName.c_str(), &info) != 0)
        return false;
    else if(info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}

// Recursively make a directory and set its permissions to 775
void mkdirRecursive(const char *dir) {
    if(folderExists(dir)) return;
    char tmp[8192];
    char *p = NULL;
    size_t len;
    snprintf(tmp, sizeof(tmp),"%s",dir);
    len = strlen(tmp);
    if (tmp[len - 1] == '/')
        tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++){
        if (*p == '/') {
            *p = 0;

            #ifdef _WIN32
            mkdir(tmp);
            #else
            mkdir(tmp, 0775);
            #endif

            chmod(tmp, 0775);
            *p = '/';
        }
    }
    #ifdef _WIN32
    mkdir(tmp);
    #else
    mkdir(tmp, 0775);
    #endif
    chmod(tmp, 0775);
}

bool fileExists(const std::string &fileName){
    if (FILE *file = fopen(fileName.c_str(), "r")){
        fclose(file);
        return true;
    }
    else return false;
}

void makeDimensionFolders(const std::string &fileName){
    size_t lastSlash = fileName.find_last_of("/");
    mkdirRecursive(fileName.substr(0,lastSlash).c_str());
}