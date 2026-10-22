#include <cstdint>
#ifdef _WIN32
#include <stdarg.h>
#include <sys/time.h>
#endif
#include <sys/stat.h>
#include "mex.h"
#include "../src/zarr.h"
#include "../src/helperfunctions.h"
#include "zarrmexhelpers.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
#if defined(__APPLE__)
    // keep the mex file loaded until MATLAB exits; unloading it crashes MATLAB on exit
    if (!mexIsLocked()) {
        mexLock();
    }
#endif
    if(nrhs < 1) mexErrMsgIdAndTxt("zarr:inputError","This functions requires at least 1 argument\n");
    if(nrhs == 2) mexErrMsgIdAndTxt("zarr:inputError","This functions does not accept only 2 arguments\n");
    if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("zarr:inputError","The first argument must be a string\n");
    zarr Zarr;
    Zarr.set_fileName(mxArrayToString(prhs[0]));
    // Per-axis options are fit to the number of dimensions after all options are read
    std::vector<uint64_t> chunksVals, subfoldersVals, chunkShapeVals;
    bool hasChunks = false, hasSubfolders = false, hasChunkShape = false;
    // Applied after the Zarr format, which sets their Zarr v3 defaults
    uint64_t zarrFormat = 0;
    std::string order, dimension_separator;

    for(int i = 1; i < nrhs; i+=2){
        if(i+1 == nrhs) mexErrMsgIdAndTxt("zarr:inputError","Mismatched argument pair for input number %d\n",i+1);
        if(!mxIsChar(prhs[i])) mexErrMsgIdAndTxt("zarr:inputError","The argument in input location %d is not a string\n",i+1);
        std::string currInput = mxArrayToString(prhs[i]);

        if(currInput == "chunks"){
            chunksVals = mexVector(prhs[i+1]);
            hasChunks = true;
        }
        else if(currInput == "cname"){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","cname must be a string\n");
            Zarr.set_cname(mxArrayToString(prhs[i+1]));
        }
        else if(currInput == "dtype"){
            Zarr.set_dtype(mexDtype(prhs[i+1]));
        }
        else if(currInput == "order"){
            order = mexOrder(prhs[i+1]);
        }
        else if(currInput == "dimension_separator"){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","dimension_separator must be a string\n");
            dimension_separator = mxArrayToString(prhs[i+1]);
            if(dimension_separator != "." && dimension_separator != "/") mexErrMsgIdAndTxt("zarr:inputError","dimension_separator must be a . or /\n");
        }
        else if(currInput == "zarr_format"){
            zarrFormat = mexZarrFormat(prhs[i+1]);
        }
        else if(currInput == "shape"){
            const std::vector<uint64_t> shape = mexVector(prhs[i+1]);
            // An empty shape is a 0-dimensional array (a single element)
            Zarr.set_shape(shape);
        }
        else if(currInput == "clevel"){
            Zarr.set_clevel(mexClevel(prhs[i+1]));
        }
        else if(currInput == "subfolders"){
            subfoldersVals = mexVector(prhs[i+1]);
            hasSubfolders = true;
        }
        else if(currInput == "chunk_shape"){
            chunkShapeVals = mexVector(prhs[i+1]);
            hasChunkShape = true;
        }
        else{
            mexErrMsgIdAndTxt("zarr:inputError","The argument \"%s\" does not match the name of any supported input name.\n \
            Currently Supported Names: chunks, cname, dtype, order, dimension_separator, shape, clevel, subfolders, chunk_shape, zarr_format\n",currInput.c_str());
        }
    }

    // The Zarr format: as given, or that of the array already at this path (2 for
    // a new one). Zarr v3 defaults to C order and "/" (chunk files c/0/0/0).
    // (set_fileName expanded a leading ~)
    const std::string &folderName = Zarr.get_fileName();
    if(!zarrFormat && !fileExists(folderName+"/.zarray") && fileExists(folderName+"/zarr.json")) zarrFormat = 3;
    if(zarrFormat) Zarr.set_zarr_format(zarrFormat);
    if(!order.empty()) Zarr.set_order(order);
    if(!dimension_separator.empty()) Zarr.set_dimension_separator(dimension_separator);

    const uint64_t nDims = Zarr.get_ndims();
    if(hasSubfolders) Zarr.set_subfolders(mexFitAxes(subfoldersVals, nDims, "subfolders"));
    if(hasChunks) Zarr.set_chunks(mexFitAxes(chunksVals, nDims, "chunks"));
    if(hasChunkShape){
        Zarr.set_shard(true);
        Zarr.set_chunk_shape(mexFitAxes(chunkShapeVals, nDims, "chunk_shape"));
    }

    try{
        Zarr.write_zarray();
    }
    catch(const std::string &e){
        if(e == "unsupportedCompressor"){
            mexErrMsgIdAndTxt("zarr:zarrayError","Compressor: \"%s\" is not currently supported\n",Zarr.get_cname().c_str());
        }
        else if(e.find("cannotOpenZarray") != std::string::npos){
            mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s for writing. Try checking permissions and path.\n",e.substr(e.find(':')+1).c_str());
        }
        else mexZarrError(e);
    }
}