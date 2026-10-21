#ifndef CPPZARR_ZARRMEXHELPERS_H
#define CPPZARR_ZARRMEXHELPERS_H
// Argument handling shared by the MEX functions for arrays of any number of
// dimensions. Existing 3D calls behave exactly as before; lower-dimensional
// arrays also accept the old 3D-style forms (extra trailing singleton axes).
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>
#include "mex.h"
#include "../src/helperfunctions.h"

// Numeric vector argument (read as double, like mxGetPr)
static std::vector<uint64_t> mexVector(const mxArray* a){
    const size_t n = mxGetNumberOfElements(a);
    std::vector<uint64_t> v(n);
    const double* p = mxGetPr(a);
    for(size_t i = 0; i < n; i++) v[i] = (uint64_t)p[i];
    return v;
}

// Fit a per-axis option (chunks, subfolders, chunk_shape) to nDims axes. Fewer
// values are fine (the library fills in its defaults for the remaining axes); 3
// values are accepted for 1D/2D arrays (the old 3D form); otherwise more values
// than axes is an error.
static std::vector<uint64_t> mexFitAxes(const std::vector<uint64_t> &v, const uint64_t nDims, const char* name){
    if(v.size() <= nDims) return v;
    if(v.size() == 3 && nDims < 3) return std::vector<uint64_t>(v.begin(), v.begin()+nDims);
    mexErrMsgIdAndTxt("zarr:inputError","%s has %d values but the array has %d dimensions\n",
                      name, (int)v.size(), (int)nDims);
    return v;
}

// Parse a region given as [starts ends] (1-based starts, inclusive ends) for an
// array with nDims axes: 2*nDims values; or 2*M values with M > nDims whose extra
// axes are singleton (start 1, end 1); or, when the array's shape is given, 2*M
// values with M < nDims if the array's remaining axes have size 1 (e.g. a 2D
// region on a [m n 1] array). Returns false if the length does not fit.
static bool mexParseRegion(const std::vector<uint64_t> &b, const uint64_t nDims,
                           std::vector<uint64_t> &start, std::vector<uint64_t> &end,
                           const std::vector<uint64_t> &shape = std::vector<uint64_t>()){
    if(b.empty() || b.size()%2) return false;
    const uint64_t m = b.size()/2;
    for(uint64_t d = nDims; d < m; d++){
        if(b[d] != 1 || b[m+d] != 1) return false;
    }
    for(uint64_t d = m; d < nDims; d++){
        if(d >= shape.size() || shape[d] != 1) return false;
    }
    start.assign(nDims,0);
    end.assign(nDims,1);
    for(uint64_t d = 0; d < nDims && d < m; d++){
        if(b[d] < 1) mexErrMsgIdAndTxt("zarr:inputError","Lower bounds must be at least 1");
        start[d] = b[d]-1;
        end[d] = b[m+d];
    }
    return true;
}
// The order option: "F" or "C" (either case)
static std::string mexOrder(const mxArray* a){
    if(!mxIsChar(a)) mexErrMsgIdAndTxt("zarr:inputError","order must be a string\n");
    const std::string order(mxArrayToString(a));
    if(order == "F" || order == "f") return "F";
    if(order == "C" || order == "c") return "C";
    mexErrMsgIdAndTxt("zarr:inputError","order must be \"F\" or \"C\"\n");
    return order;
}

// The zarr_format option: 2 or 3
static uint64_t mexZarrFormat(const mxArray* a){
    if(!mxIsNumeric(a) || mxGetNumberOfElements(a) != 1) mexErrMsgIdAndTxt("zarr:inputError","zarr_format must be 2 or 3\n");
    const double v = mxGetScalar(a);
    if(v != 2 && v != 3) mexErrMsgIdAndTxt("zarr:inputError","zarr_format must be 2 or 3\n");
    return (uint64_t)v;
}

// The clevel option: a whole number from 0 to 9 (the levels of blosc's compressors
// and of gzip)
static uint64_t mexClevel(const mxArray* a){
    if(!mxIsNumeric(a) || mxGetNumberOfElements(a) != 1) mexErrMsgIdAndTxt("zarr:inputError","clevel must be a whole number from 0 to 9\n");
    const double v = mxGetScalar(a);
    if(!(v >= 0 && v <= 9) || v != (double)(uint64_t)v) mexErrMsgIdAndTxt("zarr:inputError","clevel must be a whole number from 0 to 9\n");
    return (uint64_t)v;
}

// The dtype option: a NumPy typestr ("<u2", "f4", "|b1", ...; little-endian when
// it has no byte order mark), a type name ("uint16", "single", "float16", "half",
// "complex64", "logical", ...), or "bfloat16"
static std::string mexDtype(const mxArray* a){
    if(!mxIsChar(a)) mexErrMsgIdAndTxt("zarr:inputError","dtype must be a string\n");
    std::string s = mxArrayToString(a);
    static const char* const names[][2] = {
        {"bool", "|b1"}, {"logical", "|b1"}, {"b1", "|b1"},
        {"int8", "<i1"}, {"uint8", "<u1"}, {"int16", "<i2"}, {"uint16", "<u2"},
        {"int32", "<i4"}, {"uint32", "<u4"}, {"int64", "<i8"}, {"uint64", "<u8"},
        {"float16", "<f2"}, {"half", "<f2"}, {"float32", "<f4"}, {"single", "<f4"},
        {"float64", "<f8"}, {"double", "<f8"}, {"complex64", "<c8"}, {"complex128", "<c16"}};
    for(const auto &n : names){
        if(s == n[0]){
            s = n[1];
            break;
        }
    }
    if(!s.empty() && s[0] != '<' && s[0] != '>' && s[0] != '|' && s != "bfloat16") s.insert(0, 1, '<');
    if(!dtypeInfo(s).kind){
        mexErrMsgIdAndTxt("zarr:inputError","dtype \"%s\" is not supported. Supported: bool (b1), int8 to int64 (i1 to i8), "
                          "uint8 to uint64 (u1 to u8), float16 (f2), single (f4), double (f8), bfloat16, complex64 (c8) "
                          "and complex128 (c16)\n",mxArrayToString(a));
    }
    return s;
}

// The data type of a MATLAB array's own class: logical is bool, complex single and
// double are complex64 and complex128. Empty for a class Zarr has no type for.
static std::string mexClassDtype(const mxArray* a){
    const bool cplx = mxIsComplex(a);
    switch(mxGetClassID(a)){
        case mxLOGICAL_CLASS: return "|b1";
        case mxUINT8_CLASS:  return cplx ? "" : "<u1";
        case mxINT8_CLASS:   return cplx ? "" : "<i1";
        case mxUINT16_CLASS: return cplx ? "" : "<u2";
        case mxINT16_CLASS:  return cplx ? "" : "<i2";
        case mxUINT32_CLASS: return cplx ? "" : "<u4";
        case mxINT32_CLASS:  return cplx ? "" : "<i4";
        case mxUINT64_CLASS: return cplx ? "" : "<u8";
        case mxINT64_CLASS:  return cplx ? "" : "<i8";
        case mxSINGLE_CLASS: return cplx ? "<c8" : "<f4";
        case mxDOUBLE_CLASS: return cplx ? "<c16" : "<f8";
        default: return "";
    }
}

// MATLAB's conversions to an integer class: rounded to the nearest integer (halves
// away from zero) and saturated at the class's limits, NaN to 0
template <typename D> static inline D mexRoundToInt(const double v){
    if(std::isnan(v)) return 0;
    const double r = std::round(v);
    if(r <= (double)std::numeric_limits<D>::min()) return std::numeric_limits<D>::min();
    if(r >= (double)std::numeric_limits<D>::max()) return std::numeric_limits<D>::max();
    return (D)r;
}
template <typename D> static inline D mexSignedToInt(const int64_t v){
    if(std::numeric_limits<D>::is_signed){
        if(v < (int64_t)std::numeric_limits<D>::min()) return std::numeric_limits<D>::min();
        if(v > (int64_t)std::numeric_limits<D>::max()) return std::numeric_limits<D>::max();
        return (D)v;
    }
    if(v < 0) return 0;
    return (uint64_t)v > (uint64_t)std::numeric_limits<D>::max() ? std::numeric_limits<D>::max() : (D)v;
}
template <typename D> static inline D mexUnsignedToInt(const uint64_t v){
    return v > (uint64_t)std::numeric_limits<D>::max() ? std::numeric_limits<D>::max() : (D)v;
}
template <typename S, typename D> static inline D mexToInt(const S v){
    if(!std::numeric_limits<S>::is_integer) return mexRoundToInt<D>((double)v);
    if(std::numeric_limits<S>::is_signed) return mexSignedToInt<D>((int64_t)v);
    return mexUnsignedToInt<D>((uint64_t)v);
}

template <typename S, typename D>
static void mexConvertInts(const S* re, const uint64_t stride, const int64_t n, void* out){
    D* o = (D*)out;
    #pragma omp parallel for
    for(int64_t i = 0; i < n; i++) o[i] = mexToInt<S, D>(re[i*stride]);
}

// n elements (real parts re, imaginary parts im or NULL, stride elements apart)
// converted to the data type t, as MATLAB converts between classes; any nonzero
// value is true
template <typename S>
static void mexConvert(const S* re, const S* im, const uint64_t stride, const uint64_t n, const DtypeInfo &t, void* out){
    const int64_t N = (int64_t)n;
    if(t.kind == 'b'){
        uint8_t* o = (uint8_t*)out;
        #pragma omp parallel for
        for(int64_t i = 0; i < N; i++) o[i] = re[i*stride] != 0 || (im && im[i*stride] != 0);
    }
    else if(t.kind == 'u'){
        if(t.bytes == 1) mexConvertInts<S, uint8_t>(re, stride, N, out);
        else if(t.bytes == 2) mexConvertInts<S, uint16_t>(re, stride, N, out);
        else if(t.bytes == 4) mexConvertInts<S, uint32_t>(re, stride, N, out);
        else mexConvertInts<S, uint64_t>(re, stride, N, out);
    }
    else if(t.kind == 'i'){
        if(t.bytes == 1) mexConvertInts<S, int8_t>(re, stride, N, out);
        else if(t.bytes == 2) mexConvertInts<S, int16_t>(re, stride, N, out);
        else if(t.bytes == 4) mexConvertInts<S, int32_t>(re, stride, N, out);
        else mexConvertInts<S, int64_t>(re, stride, N, out);
    }
    else if(t.kind == 'f' && t.bytes == 2){
        uint16_t* o = (uint16_t*)out;
        const bool bf = t.bfloat16;
        #pragma omp parallel for
        for(int64_t i = 0; i < N; i++) o[i] = bf ? floatToBfloat16((float)re[i*stride]) : doubleToHalf((double)re[i*stride]);
    }
    else if(t.kind == 'f' && t.bytes == 4){
        float* o = (float*)out;
        #pragma omp parallel for
        for(int64_t i = 0; i < N; i++) o[i] = (float)re[i*stride];
    }
    else if(t.kind == 'f'){
        double* o = (double*)out;
        #pragma omp parallel for
        for(int64_t i = 0; i < N; i++) o[i] = (double)re[i*stride];
    }
    else if(t.kind == 'c' && t.bytes == 8){
        float* o = (float*)out;
        #pragma omp parallel for
        for(int64_t i = 0; i < N; i++){
            o[2*i] = (float)re[i*stride];
            o[2*i+1] = im ? (float)im[i*stride] : 0.0f;
        }
    }
    else if(t.kind == 'c'){
        double* o = (double*)out;
        #pragma omp parallel for
        for(int64_t i = 0; i < N; i++){
            o[2*i] = (double)re[i*stride];
            o[2*i+1] = im ? (double)im[i*stride] : 0.0;
        }
    }
}

template <typename S>
static void mexConvertArray(const mxArray* a, const uint64_t n, const DtypeInfo &t, void* out){
    const S* re = (const S*)mxGetData(a);
    const S* im = NULL;
    uint64_t stride = 1;
    if(mxIsComplex(a)){
#if (MX_HAS_INTERLEAVED_COMPLEX)
        // (real and imaginary parts alternate)
        im = re+1;
        stride = 2;
#else
        im = (const S*)mxGetImagData(a);
#endif
    }
    mexConvert<S>(re, im, stride, n, t, out);
}

// A MATLAB array's n elements as the data type (a malloc'd buffer the caller
// frees), converted as MATLAB converts between classes. NULL for classes Zarr has
// no type for.
static void* mexConvertToDtype(const mxArray* a, const std::string &dtype, const uint64_t n){
    const DtypeInfo t = dtypeInfo(dtype);
    if(!t.kind) return NULL;
    void* out = malloc(n*t.bytes ? n*t.bytes : 1);
    if(!out) return NULL;
    switch(mxGetClassID(a)){
        case mxLOGICAL_CLASS: mexConvertArray<mxLogical>(a, n, t, out); break;
        case mxUINT8_CLASS:   mexConvertArray<uint8_t>(a, n, t, out); break;
        case mxINT8_CLASS:    mexConvertArray<int8_t>(a, n, t, out); break;
        case mxUINT16_CLASS:  mexConvertArray<uint16_t>(a, n, t, out); break;
        case mxINT16_CLASS:   mexConvertArray<int16_t>(a, n, t, out); break;
        case mxUINT32_CLASS:  mexConvertArray<uint32_t>(a, n, t, out); break;
        case mxINT32_CLASS:   mexConvertArray<int32_t>(a, n, t, out); break;
        case mxUINT64_CLASS:  mexConvertArray<uint64_t>(a, n, t, out); break;
        case mxINT64_CLASS:   mexConvertArray<int64_t>(a, n, t, out); break;
        case mxSINGLE_CLASS:  mexConvertArray<float>(a, n, t, out); break;
        case mxDOUBLE_CLASS:  mexConvertArray<double>(a, n, t, out); break;
        default: free(out); return NULL;
    }
    return out;
}

// Report an error the zarr class threw (as a string) that the caller does not
// handle itself
static void mexZarrError(const std::string &e){
    const std::string detail = e.substr(e.find(':')+1);
    if(e.rfind("metadataUnsupported:", 0) == 0){
        mexErrMsgIdAndTxt("zarr:zarrayError","This zarr array cannot be read: %s\n",detail.c_str());
    }
    if(e.rfind("zarrV3NotWritable:", 0) == 0){
        mexErrMsgIdAndTxt("zarr:zarrayError","%s is a Zarr v3 array. Write it with 'zarr_format', 3, or delete it first\n",detail.c_str());
    }
    if(e.rfind("v3Unsupported:", 0) == 0){
        mexErrMsgIdAndTxt("zarr:zarrayError","This array cannot be written as Zarr v3: %s\n",detail.c_str());
    }
    mexErrMsgIdAndTxt("zarr:zarrayError","Unknown error occurred (%s)\n",e.c_str());
}
#endif
