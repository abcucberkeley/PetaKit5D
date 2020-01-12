/* [] = binterp(signal, coords);
 *
 * Francois Aguet, May 2012 (last modified May 11, 2012)
 *
 * Compilation:
 * Mac/Linux: mex -I/usr/local/include -I../../mex/include binterp.cpp
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT"  -I"..\..\mex\include" binterp.cpp
 */

#include <iostream>
#include <cstring> // memcpy
#include <algorithm> // transform
#include <math.h>
#include "mex.h"
#include "interpolator.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //===========================================================
    // Input checks 
    //===========================================================
    /* The required inputs are:
        input: input signal (1D or 2D)
        xi: x interpolation coordinates
        yi: y interpolation coordinates
       Other relevant variables:
        N: number of points in input    
        ni: number of interpolation points
        vi: interpolated values
    */
    
    if (nrhs<2) {
        mexErrMsgTxt("Insufficient inputs.");
    }
    
    // check image/signal
    if (!mxIsDouble(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 2) // Matlab 1D signals have 2 dimensions
        mexErrMsgTxt("Input must be a 1D or 2D signal (double array).");
    int nx = (int)mxGetN(prhs[0]); // #cols
    int ny = (int)mxGetM(prhs[0]);
    int N = nx*ny;
    double* input = mxGetPr(prhs[0]);
    
    int dims;
    if (nx==1 || ny==1) {
        dims = 1;
        if (nx==1) { // use x coordinate as index
            nx = ny;
            ny = 1;
        }
    } else {
        dims = 2;
    }
    
    // check input number and type
    if (dims==1 && nrhs!=2 && !(nrhs==3 && mxIsChar(prhs[2]))) {
        mexErrMsgTxt("Incompatible input. For 1D signals, inputs should be:\nbinterp(signal, coordinates, {'mirror'|'periodic'})");
    }
    if (dims==2 && nrhs!=3 && !(nrhs==4 && mxIsChar(prhs[3]))) {
        mexErrMsgTxt("Incompatible input. For 2D signals, inputs should be:\nbinterp(signal, x-coordinates, y-coordinates, {'mirror'|'periodic'})");
    }
    
    // check for NaNs in input, as these will result in a crash
    for (int i=0;i<N;++i) {
        if (mxIsNaN(input[i])) {
            mexErrMsgTxt("Input contains NaNs.");
            break;
        }
    }
    
    // get x interpolation coordinates
    if (!mxIsDouble(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("Interpolation coordinates must match dimensionality of input.");
    int nx_X = (int)mxGetN(prhs[1]);
    int ny_X = (int)mxGetM(prhs[1]);
    int ni = nx_X*ny_X;
    
    const double* xin = mxGetPr(prhs[1]);
    double* xi = new double[ni];
    
    // verify that coordinates are within correct range ([0 ... nx+1])
    for (int i=0;i<ni;++i) {
        if (mxIsNaN(xin[i]) || xin[i]<0.0 || xin[i]>nx+1.0) {
            mexErrMsgTxt("x-interpolation coordinates must be real-valued in interval [0..nx+1] (data range: [1..nx]).");
            break;
        } else {
            xi[i] = xin[i] - 1.0; // internal index starts at 0, Matlab at 1
        }
    }
    
    // get y interpolation coordinates
    double* yi = NULL;
    if (dims==2) {
        if (!mxIsDouble(prhs[2]) || mxGetNumberOfDimensions(prhs[2]) != 2)
        mexErrMsgTxt("Interpolation coordinates must match dimensionality of input.");
        int nx_Y = (int)mxGetN(prhs[2]);
        int ny_Y = (int)mxGetM(prhs[2]);
        if (nx_X!=nx_Y || ny_X!=ny_Y) {
            mexErrMsgTxt("Sizes of interpolation coordinate arrays must match.");
        }
        const double* yin = mxGetPr(prhs[2]);
        yi = new double[ni];
        // verify that coordinates are within signal range
        for (int i=0;i<ni;++i) {
            if (mxIsNaN(yin[i]) || yin[i]<0.0 || yin[i]>ny+1.0) {
                mexErrMsgTxt("y-interpolation coordinates must be real-valued in interval [0..ny+1] (data range: [1..ny]).");
                break;
            } else {
                yi[i] = yin[i] - 1.0;
            }
        }
    }
    
    
    // check border condition input
    int bcIdx = 0;
    if (dims==1 && nrhs==3) {
        bcIdx = 2;
    } else if (dims==2 && nrhs==4) {
        bcIdx = 3;
    }
    
    // default border conditions
    int borderCondition = Interpolator::MIRROR;
    if (bcIdx!=0) {
        if (!mxIsChar(prhs[bcIdx])) { // border condition
            mexErrMsgTxt("Border condition must be identifier string: 'mirror' or 'periodic'.");
        }        
        // parse border conditions
        size_t nchar = mxGetNumberOfElements(prhs[bcIdx])+1;
        char *ch = new char[nchar];
        int f = mxGetString(prhs[bcIdx], ch,  nchar);
        if (f!=0) {
            mexErrMsgTxt("Error parsing border condition.");
        }
        string str = ch;
        delete ch;
        int (*pf)(int) = tolower;
        transform(str.begin(), str.end(), str.begin(), pf);
        if (str.compare("mirror")==0 || str.compare("symmetric")==0) {
            borderCondition = Interpolator::MIRROR;
        } else if (str.compare("periodic")==0) {
            borderCondition = Interpolator::PERIODIC;
        } else {
            mexErrMsgTxt("Unsupported border conditions.");
        }
    }
    
    //===========================================================
    // Perform interpolation
    //===========================================================
    // Transpose for row-major index (Matlab uses column-major)
    double* pixels = new double[N];
    div_t divRes;
    for (int i=0;i<N;++i) {
        divRes = div(i, ny);
        pixels[divRes.quot+divRes.rem*nx] = input[i];
    }
    
    // The 'xi' and 'yi' matrices do not need transposition, they contain independent indexes.
    // As a result, the interpolated array 'vi' is already in the correct order for Matlab.
    Interpolator ip = Interpolator(pixels, nx, ny, borderCondition);
    double* vi = new double[ni]; // interpolated values
    if (dims==1) {
        ip.interp(xi, ni, vi);
    } else if (dims==2) {
        ip.interp(xi, yi, ni, vi);
    }

    
    //===========================================================
    // Return results to Matlab
    //===========================================================
    if (nlhs > 0) {
        plhs[0] = mxCreateDoubleMatrix(ny_X, nx_X, mxREAL);
        memcpy(mxGetPr(plhs[0]), vi, ni*sizeof(double));
    }
    
    if (nlhs>1 && dims==1) {

        double* xip = new double[ni];
        double* xim = new double[ni];
        double* vip = new double[ni];
        double* vim = new double[ni];
        for (int i=0;i<ni;++i) {
            xip[i] = xi[i]+0.5; // some values will fall outside of signal range, 'wrap()' in interp takes care of that
            xim[i] = xi[i]-0.5;
        }
        ip.interp(xip, ni, vip, 2); // compute using degree 2 spline (derivative of 3 -> 2)
        ip.interp(xim, ni, vim, 2);
        
        //B1(x+0.5) - B1(x-0.5);
        for (int i=0;i<ni;++i) {
            vip[i] -= vim[i];
        }
        plhs[1] = mxCreateDoubleMatrix(ny_X, nx_X, mxREAL);
        memcpy(mxGetPr(plhs[1]), vip, ni*sizeof(double));
        
        delete vim;
        delete vip;
        delete xim;
        delete xip;
    }
    
    if (nlhs>2 && dims==1) {
                
        double* xip = new double[ni];
        double* xim = new double[ni];
        double* vip = new double[ni];
        double* vi0 = new double[ni];
        double* vim = new double[ni];
        
        for (int i=0;i<ni;++i) {
            xip[i] = xi[i]+1.0; // some values will fall outside of signal range, 'wrap()' in interp takes care of that
            xim[i] = xi[i]-1.0;
        }
        ip.interp(xip, ni, vip, 1); // compute using degree 1 spline (2nd derivative of 3 -> 1)
        ip.interp(xim, ni, vim, 1);
        ip.interp(xi, ni, vi0, 1);
        
        //B1(x-1) - 2*B1(x) + B1(x+1);
        for (int i=0;i<ni;++i) {
            vip[i] += vim[i] - 2.0*vi0[i];
        }
        
        plhs[2] = mxCreateDoubleMatrix(ny_X, nx_X, mxREAL);
        memcpy(mxGetPr(plhs[2]), vip, ni*sizeof(double));
        
        delete vim;
        delete vi0;
        delete vip;
        delete xim;
        delete xip;
    }
    
    if (dims==2) {    
        if (nlhs>1) {

            double* xip = new double[ni];
            double* xim = new double[ni];
            double* vip = new double[ni];
            double* vim = new double[ni];
            
            // dx
            for (int i=0;i<ni;++i) {
                xip[i] = xi[i]+0.5; // some values will fall outside of signal range, 'wrap()' in interp takes care of that
                xim[i] = xi[i]-0.5;
            }
            ip.interp(xip, yi, ni, vip, 2); // compute using degree 2 spline (derivative of 3 -> 2)
            ip.interp(xim, yi, ni, vim, 2);
            for (int i=0;i<ni;++i) {
                vip[i] -= vim[i];
            }
            plhs[1] = mxCreateDoubleMatrix(ny_X, nx_X, mxREAL);
            memcpy(mxGetPr(plhs[1]), vip, ni*sizeof(double));
            
            if (nlhs>2) {
                // dy
                for (int i=0;i<ni;++i) {
                    xip[i] = yi[i]+0.5;
                    xim[i] = yi[i]-0.5;
                }
                ip.interp(xi, xip, ni, vip, 2);
                ip.interp(xi, xim, ni, vim, 2);
                for (int i=0;i<ni;++i) {
                    vip[i] -= vim[i];
                }
                plhs[2] = mxCreateDoubleMatrix(ny_X, nx_X, mxREAL);
                memcpy(mxGetPr(plhs[2]), vip, ni*sizeof(double));
            }
            
            delete vim;
            delete vip;
            delete xim;
            delete xip;
        }
        if (nlhs>3) {
            double* xip = new double[ni];
            double* xim = new double[ni];
            double* vip = new double[ni];
            double* vi0 = new double[ni];
            double* vim = new double[ni];
            
            // d2x
            for (int i=0;i<ni;++i) {
                xip[i] = xi[i]+1.0;
                xim[i] = xi[i]-1.0;
            }
            ip.interp(xip, yi, ni, vip, 1);
            ip.interp(xim, yi, ni, vim, 1);
            ip.interp(xi, yi, ni, vi0, 1);
            for (int i=0;i<ni;++i) {
                vip[i] += vim[i] - 2.0*vi0[i];
            }
            plhs[3] = mxCreateDoubleMatrix(ny_X, nx_X, mxREAL);
            memcpy(mxGetPr(plhs[3]), vip, ni*sizeof(double));
            
            if (nlhs>4) {
                // d2y
                for (int i=0;i<ni;++i) {
                    xip[i] = yi[i]+1.0;
                    xim[i] = yi[i]-1.0;
                }
                ip.interp(xi, xip, ni, vip, 1);
                ip.interp(xi, xim, ni, vim, 1);
                ip.interp(xi, yi, ni, vi0, 1);
                for (int i=0;i<ni;++i) {
                    vip[i] += vim[i] - 2.0*vi0[i];
                }
                plhs[4] = mxCreateDoubleMatrix(ny_X, nx_X, mxREAL);
                memcpy(mxGetPr(plhs[4]), vip, ni*sizeof(double));
            }
            
            delete vim;
            delete vi0;
            delete vip;
            delete xim;
            delete xip;
            
        }
    }
    
    //===========================================================
    // De-allocate
    //===========================================================
    delete xi;
    delete vi;
    
    if (dims>1) {
        delete yi;
    }
    
}

// compiled with:
// export DYLD_LIBRARY_PATH=/Applications/MATLAB_R2012a.app/bin/maci64 && g++ -Wall -g -DARRAY_ACCESS_INLINING -I. -L/Applications/MATLAB_R2012a.app/bin/maci64 -I../../mex/include/ -I/Applications/MATLAB_R2012a.app/extern/include binterp.cpp -lmx -lmex
// tested with:
// valgrind --tool=memcheck --leak-check=full --show-reachable=yes ./a.out 2>&1 | grep binterp

/*int main(void) {
    double f[] = {1.0, 2.0, 3.0};
    double xi[] = {1.0, 1.5, 2.0, 2.5, 3.0};
    
    Interpolator ip = Interpolator(f, 3, 1);
    int nxi = 5;
    double* vi = new double[nxi];

    ip.interp(xi, nxi, vi);
    //ip.interp(xi, yi, ni, vi);
}*/
