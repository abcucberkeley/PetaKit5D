#include <matrix.h>
#include <mex.h>
#include <math.h>

int detectAndRemoveHotPixels(unsigned short *buf, int nx, int ny, float background, int harshness, int border_width=1)
{
  if (border_width<1)
    throw -1;

  float factor = 4 + harshness;

  float numPixels = (2*border_width+1) * (2*border_width+1) - 1;
  for (int i=border_width; i<ny-border_width; i++)
    for (int j=border_width; j<nx-border_width; j++) {
      int ii, jj;
      float stddev;

      stddev = sqrt(buf[i*nx+j] - background); // use the square root of current pixel as standard deviation
      // decide if all neighbouring pixels are different from (i, j) significantly
      int how_many_different = 0;
      for (ii=-1; ii<=1; ii++)
        for (jj=-1; jj<=1; jj++)
          if (!(ii==0 && jj==0)) {
            if (buf[i*nx+j] - buf[(i+ii)*nx + (j+jj)] > factor * stddev)
              how_many_different ++;
          }
      // If all neighbouring pixel are significantly different, replace it with the mean
      if (how_many_different == 8) {
        // Calculate mean and stddev of the 9 neighboring pixels
        float mean=0;
        for (ii=-border_width; ii<=border_width; ii++)
          for (jj=-border_width; jj<=border_width; jj++)
            if (!(ii==0 && jj==0))
              mean += buf[(i+ii)*nx + (j+jj)];

        // mexPrintf("mean=%.3f at (%d, %d)\n", mean, j, i);

        buf[i*nx + j] = (unsigned short) rint(mean/numPixels);
      }
    }
  return 0;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mwSize *dims;
  int nx, ny, nz, nxy;

  dims = mxGetDimensions(prhs[0]);

  nx = dims[0];
  ny = dims[1];
  nz = dims[2];
  // mexPrintf("nx=%d, ny=%d, nz=%d\n", nx, ny, nz);

  nxy = nx*ny;

  float background =  (float) mxGetScalar(prhs[1]);
  int harshness =  (int) mxGetScalar(prhs[2]);
  int border_width =  (int) mxGetScalar(prhs[3]);

  plhs[0] = mxDuplicateArray(prhs[0]);
  unsigned short * ptr_out = (unsigned short*) mxGetPr(plhs[0]);

  for (int z=0; z<nz; z++) {
    detectAndRemoveHotPixels(ptr_out, nx, ny, background, harshness);
    ptr_out += nxy;
  }
}
