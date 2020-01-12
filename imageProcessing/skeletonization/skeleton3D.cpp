/*
 *
 * Compilation:
 * Mac/Linux: mex -I/usr/local/include  skeleton3D.cpp
 * Windows: mex COMPFLAGS="$COMPFLAGS /TP /MT" -output skeleton3D skeleton3D.cpp
 */

#include "mex.h"
#include "matrix.h"

#include "skeleton3D.h"
#include <stdlib.h>

int volNeighbors[27] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0
};



inline void CopyNeighborhoodInBuffer(unsigned char *vol, int L, int M, int N, 
				     int idx, bool nb[3][3][3])
{
  int nidx;
  char i, j, k, ii;
  
  ii = 0;
  for(k=0; k < 3; k++) {
    for(j=0; j < 3; j++) {
      for(i=0; i < 3; i++) {
	nidx = idx + volNeighbors[ii];
	
	if(vol[nidx] != 0) {
	  nb[i][j][k] = true;
	}
	else {
	  nb[i][j][k] = false;
	}
	
	ii++;
      }
    }
  }
  
return;
}

// Mex entrypoint for matlab use
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    // Check input and output arguments
    
    if (nrhs != 1)
        mexErrMsgTxt("This function accepts exactly 1 input argument - a binary 3D matrix!");
    
    if (nlhs != 1)
        mexErrMsgTxt("This function has exactly 1 output argument - a binary 3D matrix!");    
    
    if (mxGetNumberOfDimensions(prhs[0]) != 3)
        mexErrMsgTxt("The input matrix must be 3D!");
    
    if (!mxIsLogical(prhs[0]))
        mexErrMsgTxt("The input matrix must be logical!");
    
    //Get size of input matrix
    int const *matSize = mxGetDimensions(prhs[0]);
    
    //Initialize output array and in/out pointers
    plhs[0] = mxCreateLogicalArray(3,matSize);    
    bool* matIn = (bool*) mxGetPr(prhs[0]);                        
    bool* skel = (bool*) mxGetPr(plhs[0]);      
    
    //Get matrix sizes for readability
    int L = matSize[0];
    int M = matSize[1];
    int N = matSize[2];
    int i, j, k;        
    int sz, slsz, idx;
    sz = L * M * N;
    slsz = L * M;
        
    unsigned char *volmat = (unsigned char*) malloc(sz * sizeof(unsigned char)); //Volume for labelling borders, deleted voxels etc. Stored as char to minimize memory use
    unsigned char *vol = volmat;//Pointer to volume matrix 
    char dir; //Direction of deletion
    //int nrPasses; //Number of deletion passes           
    int nrDel; //number of deleted voxels
    
    bool nb[3][3][3];  //Neighborhood of current point
    bool USn[3][3][3]; //Transformed neighborhoods
        
    
    // initialize global neighbors array
    // lower plane
    volNeighbors[0] = (-slsz -L -1);
    volNeighbors[1] = (-slsz -L +0);
    volNeighbors[2] = (-slsz -L +1);
    volNeighbors[3] = (-slsz +0 -1);
    volNeighbors[4] = (-slsz +0 +0);
    volNeighbors[5] = (-slsz +0 +1);
    volNeighbors[6] = (-slsz +L -1);
    volNeighbors[7] = (-slsz +L +0);
    volNeighbors[8] = (-slsz +L +1);
    // same plane
    volNeighbors[9]  = (+0 -L -1);
    volNeighbors[10] = (+0 -L +0);
    volNeighbors[11] = (+0 -L +1);
    volNeighbors[12] = (+0 +0 -1);
    volNeighbors[13] = (+0 +0 +0);
    volNeighbors[14] = (+0 +0 +1);
    volNeighbors[15] = (+0 +L -1);
    volNeighbors[16] = (+0 +L +0);
    volNeighbors[17] = (+0 +L +1);
    // upper plane
    volNeighbors[18] = (+slsz -L -1);
    volNeighbors[19] = (+slsz -L +0);
    volNeighbors[20] = (+slsz -L +1);
    volNeighbors[21] = (+slsz +0 -1);
    volNeighbors[22] = (+slsz +0 +0);
    volNeighbors[23] = (+slsz +0 +1);
    volNeighbors[24] = (+slsz +L -1);
    volNeighbors[25] = (+slsz +L +0);
    volNeighbors[26] = (+slsz +L +1);       

    //Label the starting object voxels in the volume using the input matrix
    int nnz = 0;
    for(idx=0; idx < sz; idx++) {
        if(matIn[idx] != 0){
            vol[idx] = OBJECT;
            //nnz++;
        }
        else{
            vol[idx] = 0;
        }
    }
    
    //std::cout << "Number of object voxels = " << nnz << "\n";

    nrDel = 1;
    //nrPasses = 1;    

    //Loop through thinning until no more points can be deleted
    while(nrDel > 0) {
        
        //std::cout << "Pass" << nrPasses << "\n";
        nrDel = 0;
        
        //Go through each direction and mark boundary points in that direction
        for(dir = 0; dir < 12; dir++) {
            
            switch(dir) {
                
                case UP_SOUTH:
                    // UP
                    markBoundaryInDirection(vol, L, M, N, UP);
                    // SOUTH
                    markBoundaryInDirection(vol, L, M, N, SOUTH);
                    break;
                case NORTH_EAST:
                    // NORTH
                    markBoundaryInDirection(vol, L, M, N, NORTH);
                    // EAST
                    markBoundaryInDirection(vol, L, M, N, EAST);
                    break;
                case WEST_DOWN:
                    // WEST
                    markBoundaryInDirection(vol, L, M, N, WEST);
                    // DOWN
                    markBoundaryInDirection(vol, L, M, N, DOWN);
                    break;
                case EAST_SOUTH:
                    // EAST
                    markBoundaryInDirection(vol, L, M, N, EAST);
                    // SOUTH
                    markBoundaryInDirection(vol, L, M, N, SOUTH);
                    break;
                case UP_WEST:
                    // UP
                    markBoundaryInDirection(vol, L, M, N, UP);
                    // WEST
                    markBoundaryInDirection(vol, L, M, N, WEST);
                    break;
                case NORTH_DOWN:
                    // NORTH
                    markBoundaryInDirection(vol, L, M, N, NORTH);
                    // DOWN
                    markBoundaryInDirection(vol, L, M, N, DOWN);
                    break;
                case SOUTH_WEST:
                    // SOUTH
                    markBoundaryInDirection(vol, L, M, N, SOUTH);
                    // WEST
                    markBoundaryInDirection(vol, L, M, N, WEST);
                    break;
                case UP_NORTH:
                    // UP
                    markBoundaryInDirection(vol, L, M, N, UP);
                    // NORTH
                    markBoundaryInDirection(vol, L, M, N, NORTH);
                    break;
                case EAST_DOWN:
                    // EAST
                    markBoundaryInDirection(vol, L, M, N, EAST);
                    // DOWN
                    markBoundaryInDirection(vol, L, M, N, DOWN);	
                    break;
                case NORTH_WEST:
                    // NORTH
                    markBoundaryInDirection(vol, L, M, N, NORTH);
                    // WEST
                    markBoundaryInDirection(vol, L, M, N, WEST);
                    break;
                case UP_EAST:
                    // UP
                    markBoundaryInDirection(vol, L, M, N, UP);
                    // EAST
                    markBoundaryInDirection(vol, L, M, N, EAST);
                    break;
                case SOUTH_DOWN:
                    // SOUTH
                    markBoundaryInDirection(vol, L, M, N, SOUTH);
                    // DOWN
                    markBoundaryInDirection(vol, L, M, N, DOWN);	
                    break;
            }                        
            
            //check each boundary point and remove it if it matches a template
            for(k=1; k < (N-1); k++) {
                for(j=1; j < (M-1); j++) {
                    for(i=1; i < (L-1); i++) {
                        
                        //Convert 3D index to 1D index
                        idx = k * slsz + j * L + i;
	    
                        //If it is a border point
                        if(vol[idx] == D_BORDER) {                                                        
                            
                            // copy this point's neighborhood into buffer                              
                            CopyNeighborhoodInBuffer(vol, L, M, N, idx, nb);
                            
                            //Get appropriate transformation of neighborhood for this direction
                            TransformNeighborhood(nb, dir, USn);	      
                            if(MatchesATemplate(USn)) {		  
                                
                                //If it matches a template, it is simple, so mark it as such
                                vol[idx] = SIMPLE;
                                nrDel++;
                            }
                        }
                    }
                }
            }
            
            //std::cout << "# border voxels = " << nN << "\n";
            
            //Delete all simple points, and re-label all the remaining points as object
            for(idx=0; idx < sz; idx++) {
                if(vol[idx] == SIMPLE) 
                    vol[idx] = 0;
                if(vol[idx] != 0) 
                    vol[idx] = OBJECT;
            }
            
            //printf("done.\n");
        }//End direction loop
        
    //std::cout << "Number of deleted voxels : " <<nrDel << "\n";
    //nrPasses++;
        
    }//End thinning loop
    
//The remaining points are all non-simple and are therefore define object's skeleton
//Return these points as output
for(idx=0; idx < sz; idx++) {        
    if(vol[idx] != 0) 
        skel[idx] = 1;
}
    free(volmat);
}//End mexFunction

