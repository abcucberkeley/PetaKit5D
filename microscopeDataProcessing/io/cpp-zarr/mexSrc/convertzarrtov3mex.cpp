#include <string>
#include "mex.h"
#include "../src/helperfunctions.h"
#include "../src/zarr.h"
#include "zarrmexhelpers.h"

// convertZarrToV3(folderName): turn a Zarr v2 array into a Zarr v3 array in place.
// A zarr.json replaces its .zarray (and takes in its .zattrs attributes), and its
// chunk files are kept as they are, so no data is rewritten.
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if(nrhs != 1) mexErrMsgIdAndTxt("zarr:inputError","This function requires 1 argument, the path of the array\n");
    if(!mxIsChar(prhs[0])) mexErrMsgIdAndTxt("zarr:inputError","The first argument must be a string\n");
    std::string folderName(mxArrayToString(prhs[0]));
    // Handle the tilde character in filenames on Linux/Mac
    #ifndef _WIN32
    folderName = expandTilde(folderName.c_str());
    #endif
    if(!fileExists(folderName+"/.zarray")){
        mexErrMsgIdAndTxt("zarr:zarrayError","%s is not a Zarr v2 array (no .zarray metadata file was found)\n",folderName.c_str());
    }

    zarr Zarr;
    try{
        Zarr = zarr(folderName);
    }
    catch(const std::string &e){
        if(e.find("metadataFileMissing") != std::string::npos){
            mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s for reading. Try checking permissions or the file path.\n",e.substr(e.find(':')+1).c_str());
        }
        else if(e == "metadataIncomplete"){
            mexErrMsgIdAndTxt("zarr:zarrayError","Metadata is incomplete. Check the .zarray file");
        }
        else mexZarrError(e);
    }

    try{
        Zarr.convert_to_v3();
    }
    catch(const std::string &e){
        if(e.rfind("v3Unsupported:", 0) == 0){
            mexErrMsgIdAndTxt("zarr:zarrayError","Cannot convert %s to Zarr v3: %s\n",folderName.c_str(),e.substr(e.find(':')+1).c_str());
        }
        else if(e.find("cannotOpenZarray") != std::string::npos){
            mexErrMsgIdAndTxt("zarr:zarrayError","Cannot open %s for writing. Try checking permissions and path.\n",e.substr(e.find(':')+1).c_str());
        }
        else mexZarrError(e);
    }
}
