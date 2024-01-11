#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <errno.h>
#include <omp.h>
#ifdef _WIN32
#include <stdarg.h>
#include <sys/time.h>
#else
#include <uuid/uuid.h>
#endif
#include <sys/stat.h>
#include "mex.h"
#include "../src/zarr.h"
#include "../src/helperfunctions.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if(nrhs < 1) mexErrMsgIdAndTxt("zarr:inputError","This functions requires at least 1 argument\n");
    if(nrhs == 2) mexErrMsgIdAndTxt("zarr:inputError","This functions does not accept only 2 arguments\n");
    if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("zarr:inputError","The first argument must be a string\n");
    zarr Zarr;
    Zarr.set_fileName(mxArrayToString(prhs[0]));

    for(int i = 1; i < nrhs; i+=2){
        if(i+1 == nrhs) mexErrMsgIdAndTxt("zarr:inputError","Mismatched argument pair for input number %d\n",i+1);
        if(!mxIsChar(prhs[i])) mexErrMsgIdAndTxt("zarr:inputError","The argument in input location %d is not a string\n",i+1);
        std::string currInput = mxArrayToString(prhs[i]);

        if(currInput == "chunks"){
            if(mxGetN(prhs[i+1]) != 3) mexErrMsgIdAndTxt("zarr:inputError","chunks must be an array of 3 numbers\n");
            Zarr.set_chunks({(uint64_t)*(mxGetPr(prhs[i+1])),
                             (uint64_t)*((mxGetPr(prhs[i+1])+1)),
                             (uint64_t)*((mxGetPr(prhs[i+1])+2))});
        }
        else if(currInput == "cname"){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","cname must be a string\n");
            Zarr.set_cname(mxArrayToString(prhs[i+1]));
        }
        else if(currInput == "dtype"){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","dtype must be a string\n");
            std::string dtype = mxArrayToString(prhs[i+1]);
            if(dtype.size() == 2){
                dtype.insert(0,1,'<');
            }
            else if(dtype.size() != 3) mexErrMsgIdAndTxt("zarr:inputError","dtype must be of length 2 or 3\n");

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
            Zarr.set_dtype(dtype);
        }
        else if(currInput == "order"){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","order must be a string\n");
            std::string order = mxArrayToString(prhs[i+1]);
            if(order.size() > 1) mexErrMsgIdAndTxt("zarr:inputError","order must be of length 1\n");
            if(order[0] != 'F' || order[0] != 'C' || order[0] != 'f' || order[0] != 'c'){
                if(order[0] == 'f') order = "F";
                else if(order [0] == 'c') order = "C";
            }
            else mexErrMsgIdAndTxt("zarr:inputError","order must be \"F\" or \"C\"\n");
            Zarr.set_order(order);
        }
        else if(currInput == "dimension_separator"){
            if(!mxIsChar(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","dimension_separator must be a string\n");
            const std::string dimension_separator(mxArrayToString(prhs[i+1]));
            if(dimension_separator != "." && dimension_separator != "/") mexErrMsgIdAndTxt("zarr:inputError","dimension_separator must be a . or /\n");
            Zarr.set_dimension_separator(dimension_separator);
        }
        else if(currInput == "shape"){
            if(mxGetN(prhs[i+1]) != 3) mexErrMsgIdAndTxt("zarr:inputError","shape must be an array of 3 numbers\n");
            Zarr.set_shape({(uint64_t)*(mxGetPr(prhs[i+1])),
                           (uint64_t)*((mxGetPr(prhs[i+1])+1)),
                           (uint64_t)*((mxGetPr(prhs[i+1])+2))});
        }
        else if(currInput == "clevel"){
            if(!mxIsNumeric(prhs[i+1])) mexErrMsgIdAndTxt("zarr:inputError","clevel must be a numerical value\n");
            Zarr.set_clevel((uint64_t)*(mxGetPr(prhs[i+1])));
        }
        else if(currInput == "subfolders"){
            if(mxGetN(prhs[i+1]) != 3) mexErrMsgIdAndTxt("zarr:inputError","subfolders must be an array of 3 numbers\n");
            Zarr.set_subfolders({(uint64_t)*(mxGetPr(prhs[i+1])),
                               (uint64_t)*((mxGetPr(prhs[i+1])+1)),
                               (uint64_t)*((mxGetPr(prhs[i+1])+2))});
        }
        else if(currInput == "chunk_shape"){
            if(mxGetN(prhs[i+1]) != 3) mexErrMsgIdAndTxt("zarr:inputError","chunk_shape must be an array of 3 numbers\n");
            Zarr.set_shard(true);
            Zarr.set_chunk_shape({(uint64_t)*(mxGetPr(prhs[i+1])),
                               (uint64_t)*((mxGetPr(prhs[i+1])+1)),
                               (uint64_t)*((mxGetPr(prhs[i+1])+2))});
        }
        else{
            mexErrMsgIdAndTxt("zarr:inputError","The argument \"%s\" does not match the name of any supported input name.\n \
            Currently Supported Names: chunks, cname, dtype, order, shape, clevel, subfolders, chunk_shape\n",currInput.c_str());
        }
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
        else mexErrMsgIdAndTxt("zarr:zarrayError","Unknown error occurred\n");
    }
}