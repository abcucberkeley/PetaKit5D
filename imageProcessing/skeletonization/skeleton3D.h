// Skeletonization of 3D binary images via thinning, as described in:
// [1] Palagyi and Kuba in "A Parallel 3D 12-Subiteration Thinning Algorithm", 1999
//
// Adapted to mex from C++ code from:
// [2] N. Cornea et al, Curve Skeleton Properties, Applications and Algorithms, 
// IEEE Trans. on Vis. and Comp. Graphics, Vol 13, No 3 (2007)
//
// Adapted by:
// Hunter Elliott
// 6/2010
//

#include <mex.h>

#include <memory.h>
#include <stdio.h>

// #define ROTATIONS_ONLY 1

//Label for object voxels
#define OBJECT 1 
//Label for border voxels
#define D_BORDER 100 
//Label for simple points
#define SIMPLE 200   

//List and order of deletion directions
enum Direction {
  UP_SOUTH = 0, 
  NORTH_EAST,
  WEST_DOWN, 
  
  EAST_SOUTH,
  UP_WEST, 
  NORTH_DOWN, 
  
  SOUTH_WEST,
  UP_NORTH, 
  EAST_DOWN, 
  
  NORTH_WEST,
  UP_EAST, 
  SOUTH_DOWN,
  
  UP,       //12
  DOWN, 
  EAST, 
  WEST, 
  NORTH, 
  SOUTH
};

//Subrouting for checking if a point's neighborhood matches any of the 
//templates from [2] Figure 3

bool MatchesATemplate(bool n[3][3][3]) {
  // T1
  if((n[1][1][1] && n[1][0][1]) 
     &&
     (!n[0][2][0] && !n[1][2][0] && !n[2][2][0] &&
      !n[0][2][1] && !n[1][2][1] && !n[2][2][1] &&
      !n[0][2][2] && !n[1][2][2] && !n[2][2][2])
     &&
     (n[0][0][0] || n[1][0][0] || n[2][0][0] ||
      n[0][0][1]               || n[2][0][1] ||
      n[0][0][2] || n[1][0][2] || n[2][0][2] ||
      n[0][1][0] || n[1][1][0] || n[2][1][0] ||
      n[0][1][1]               || n[2][1][1] ||
      n[0][1][2] || n[1][1][2] || n[2][1][2])
     )
  {
    return true;
  }
  
  // T2
  if((n[1][1][1] && n[1][1][2]) 
     &&
     (!n[0][0][0] && !n[1][0][0]  && !n[2][0][0] &&
      !n[0][1][0] && !n[1][1][0]  && !n[2][1][0] &&
      !n[0][2][0] && !n[1][2][0]  && !n[2][2][0])
     &&
     (n[0][0][1] || n[1][0][1] || n[2][0][1] ||
      n[0][1][1]               || n[2][1][1] ||
      n[0][2][1] || n[1][2][1] || n[2][2][1] ||
      n[0][0][2] || n[1][0][2] || n[2][0][2] ||
      n[0][1][2]               || n[2][1][2] ||
      n[0][2][2] || n[1][2][2] || n[2][2][2])
     )
  {
    return true;
  }

  // T3
  if((n[1][1][1] && n[1][0][2]) 
     &&
     (!n[0][0][0] && !n[1][0][0] && !n[2][0][0] &&
      !n[0][1][0] && !n[1][1][0] && !n[2][1][0] &&
      !n[0][2][0] && !n[1][2][0] && !n[2][2][0] &&
      !n[0][2][1] && !n[1][2][1] && !n[2][2][1] &&
      !n[0][2][2] && !n[1][2][2] && !n[2][2][2])
     &&
     (n[0][0][1] || n[2][0][1] || n[0][0][2] || n[2][0][2] ||
      n[0][1][1] || n[2][1][1] || n[0][1][2] || n[2][1][2])
     )
  {
    return true;
  }

  // T4
  if((n[1][1][1] && n[1][0][1] && n[1][1][2])
     &&
     (!n[1][1][0] && !n[0][2][0] && !n[1][2][0] && !n[2][2][0] && !n[1][2][1])
     &&
     !(n[0][1][0] && n[0][2][1])
     &&
     !(n[2][1][0] && n[2][2][1])
     )
  {
    return true;
  }

  // T5
  if((n[1][1][1] && n[1][0][1] && n[1][1][2] && n[2][2][0])
     &&
     (!n[1][1][0] && !n[0][2][0] && !n[1][2][0] && !n[1][2][1])
     &&
     !(n[0][1][0] && n[0][2][1]) 
     &&
     ((!n[2][1][0] && n[2][2][1]) || (n[2][1][0] && !n[2][2][1]))
     )
  {
    return true;
  }

  // T6
  if((n[1][1][1] && n[1][0][1] && n[1][1][2] && n[0][2][0]) 
     &&
     (!n[1][1][0] && !n[1][2][0] && !n[2][2][0] && !n[1][2][1])
     && 
     ((!n[0][1][0] && n[0][2][1]) || (n[0][1][0] && !n[0][2][1]))
     &&
     !(n[2][1][0] && n[2][2][1])
     )
  {
    return true;
  }

  // T7
  if((n[1][1][1] && n[1][0][1] && n[2][1][1] && n[1][1][2])
     &&
     (!n[1][1][0] && !n[0][2][0] && !n[1][2][0] && !n[1][2][1])
     &&
     !(n[0][1][0] && n[0][2][1])
     )
  {
    return true;
  }
  
  // T8
  if((n[1][1][1] && n[1][0][1] && n[0][1][1] && n[1][1][2])
     &&
     (!n[1][1][0] && !n[1][2][0] && !n[2][2][0] && !n[1][2][1])
     &&
     !(n[2][1][0] && n[2][2][1])
     )
  {
    return true;
  }

  // T9
  if((n[1][1][1] && n[1][0][1] && n[2][1][1] && n[1][1][2] && n[0][2][0])
     &&
     (!n[1][1][0] && !n[1][2][0] && !n[1][2][1])
     &&
     ((n[0][1][0] && !n[0][2][1]) || (!n[0][1][0] && n[0][2][1]))
     )
  {
    return true;
  }

  // T10
  if((n[1][1][1] && n[1][0][1] && n[0][1][1] && n[1][1][2] && n[2][2][0])
     &&
     (!n[1][1][0] && !n[1][2][0] && !n[1][2][1])
     &&
     ((n[2][1][0] && !n[2][2][1]) || (!n[2][1][0] && n[2][2][1]))
     )
  {
    return true;
  }

  // T11
  if((n[1][1][1] && n[2][0][1] && n[1][0][2])
     &&
     (!n[0][0][0] && !n[1][0][0] && 
      !n[0][1][0] && !n[1][1][0] && 
      !n[0][2][0] && !n[1][2][0] && !n[2][2][0] &&
      !n[0][2][1] && !n[1][2][1] && !n[2][2][1] &&
      !n[0][2][2] && !n[1][2][2] && !n[2][2][2])
     )
  {
    return true;
  }

  // T12
  if((n[1][1][1] && n[0][0][1] && n[1][0][2])
     &&
     (!n[1][0][0] && !n[2][0][0] && 
      !n[1][1][0] && !n[2][1][0] &&
      !n[0][2][0] && !n[1][2][0] && !n[2][2][0] &&
      !n[0][2][1] && !n[1][2][1] && !n[2][2][1] &&
      !n[0][2][2] && !n[1][2][2] && !n[2][2][2])
     )
  {
    return true;
  }

  // T13
  if((n[1][1][1] && n[1][0][2] && n[2][1][2])
     &&
     (!n[0][0][0] && !n[1][0][0] && !n[2][0][0] &&
      !n[0][1][0] && !n[1][1][0] && !n[2][1][0] &&
      !n[0][2][0] && !n[1][2][0] && !n[2][2][0] &&
      !n[0][2][1] && !n[1][2][1] && 
      !n[0][2][2] && !n[1][2][2])
     )
  {
    return true;
  }

  // T14
  if((n[1][1][1] && n[1][0][2] && n[0][1][2])
     &&
     (!n[0][0][0] && !n[1][0][0] && !n[2][0][0] &&
      !n[0][1][0] && !n[1][1][0] && !n[2][1][0] &&
      !n[0][2][0] && !n[1][2][0] && !n[2][2][0] &&
      !n[1][2][1] && !n[2][2][1] &&
      !n[1][2][2] && !n[2][2][2])
     )
  {
    return true;
  }

  return false;
}



// transform neighborhood from a different direction
bool TransformNeighborhood(bool n[3][3][3], char direction, 
			   bool USn[3][3][3]) 
{
  char i, j, k;
  bool tmp[3][3][3];
  
  switch(direction) {
  case 0:  //UP_SOUTH = 0, 
    // just copy
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[i][j][k] = n[i][j][k];
	}
      }
    }
    break;

    // The following cases are clearly rotations only
  case 3: //   EAST_SOUTH,
    
    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[2-j][i][k] = n[i][j][k];
	}
      }
    }
    break;

  case 6: //  SOUTH_WEST,
    
    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[j][2-i][k] = n[i][j][k];
	}
      }
    }
    break;

  case 10: //   UP_EAST, 

    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[k][j][2-i] = n[i][j][k];
	}
      }
    }
    break;
  case 4: //  UP_WEST, 
    
    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[2-k][j][i] = n[i][j][k];
	}
      }
    }
    break;  


    // the next 3 are opposite sides - we can do either 2 or 3 rotations 
    // or a single reflection
#ifdef ROTATIONS_ONLY
  case 11: //    SOUTH_DOWN  
    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  tmp[2-j][i][k] = n[i][j][k];
	}
      }
    }
    // 2
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[2-j][i][k] = tmp[i][j][k];
	}
      }
    }
    break;

  case 7: //  UP_NORTH, 
    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  tmp[k][j][2-i] = n[i][j][k];
	}
      }
    }
    // 2
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[k][j][2-i] = tmp[i][j][k];
	}
      }
    }
    break;


  case 5: //   NORTH_DOWN, 
    
    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  tmp[k][j][2-i] = n[i][j][k];
	}
      }
    }
    // 2
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  n[i][2-k][j] = tmp[i][j][k];
	}
      }
    }
    // 3
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[2-j][i][k] = n[i][j][k];
	}
      }
    }
    break;

#else // not ROTATIONS_ONLY
    // use reflection
  case 11: //    SOUTH_DOWN  
    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[i][2-j][k] = n[i][j][k];
	}
      }
    }
    break;
    
  case 7: //  UP_NORTH, 

    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[i][j][2-k] = n[i][j][k];
	}
      }
    }
    break;

  case 5: //   NORTH_DOWN, 
    /*
    // 1 reflection
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[i][k][j] = n[i][j][k];
	}
      }
    }
    break;
    */
    
    // OR - two reflections
    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  tmp[i][j][2-k] = n[i][j][k];
	}
      }
    }
    // 2
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[i][2-j][k] = tmp[i][j][k];
	}
      }
    }
    break;
    
#endif // ROTATIONS_ONLY




    // The following 4 cases need
    // either 2 rotations or one reflection and one rotation
#ifdef ROTATIONS_ONLY
  case 8: //   EAST_DOWN, 
    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  tmp[i][2-k][j] = n[i][j][k];
	}
      }
    }
    // 2
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[2-j][i][k] = tmp[i][j][k];
	}
      }
    }
    break;

  case 2: // WEST_DOWN, 
    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  tmp[i][2-k][j] = n[i][j][k];
	}
      }
    }
    // 2
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[j][2-i][k] = tmp[i][j][k];
	}
      }
    }
    break;

  case 1:  //NORTH_EAST,

    // 1 - rotation
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  tmp[i][k][2-j] = n[i][j][k];
	}
      }
    }
    // 2 - rotation
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[k][j][2-i] = tmp[i][j][k];
	}
      }
    }
    break;

  case 9: //   NORTH_WEST,

    // 1 - rotation
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  tmp[i][k][2-j] = n[i][j][k];
	}
      }
    }
    // 2 - rotation
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[2-k][j][i] = tmp[i][j][k];
	}
      }
    }
    break;

#else // not ROTATIONS_ONLY
  case 8: //   EAST_DOWN, 
  
     // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  tmp[k][j][2-i] = n[i][j][k];
	}
      }
    }
    // 2
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[i][2-j][k] = tmp[i][j][k];
	}
      }
    }
    break;

  case 2: // WEST_DOWN, 
    // 1
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  tmp[2-k][j][i] = n[i][j][k];
	}
      }
    }
    // 2
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[i][2-j][k] = tmp[i][j][k];
	}
      }
    }
    break;

  case 1:  //NORTH_EAST,

    // 1 - reflection
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  tmp[i][j][2-k] = n[i][j][k];
	}
      }
    }
    // 2 - rotation
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[2-j][i][k] = tmp[i][j][k];
	}
      }
    }
    break;

  case 9: //   NORTH_WEST,

    // 1 - reflection
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  tmp[i][j][2-k] = n[i][j][k];
	}
      }
    }
    // 2 - rotation
    for(k=0; k < 3; k++) {
      for(j=0; j < 3; j++) {
	for(i=0; i < 3; i++) {
	  USn[j][2-i][k] = tmp[i][j][k];
	}
      }
    }
    break;
#endif // ROTATIONS_ONLY

  }
  
  return true; 
}


bool markBoundaryInDirection(unsigned char *vol, int L, int M, int N,
			      char direction)
{
  
  int slsz = L*M;
  int idx;
  int i, j, k;

  
  // neighbor index in 18 directions (only last 6 used)
  int nb[18] = {
    +L - slsz,  // UP_SOUTH,   0
    +slsz + 1,// NORTH_EAST     1
    -1 - L,     // WEST_DOWN,  2 

    +1 -slsz, // EAST_SOUTH    3
    +L - 1,     // UP_WEST,    4
    +slsz - L,  // NORTH_DOWN, 5 

    -slsz - 1,// SOUTH_WEST    6
    +L + slsz,  // UP_NORTH,   7
    +1 - L,     // EAST_DOWN,  8

    +slsz - 1,// NORTH_WEST     9
    +L + 1,     // UP_EAST,   10
    -slsz - L,  // SOUTH_DOWN 11

    +L, // UP                 12
    -L, // DOWN               13
    +1,  // EAST,             14
    -1,  // WEST,             15
    +slsz,  // NORTH,         16
    -slsz  // SOUTH,          17
  };
  
  //printf("Marking boundary in direction: %d (offset: %d)\n", direction, nb[direction]);

  for(k=1; k < (N-1); k++) {
    for(j=1; j < (M-1); j++) {
      for(i=1; i < (L-1); i++) {
	idx = k*slsz + j*L + i;
	
	if((vol[idx] == OBJECT) && (vol[idx + nb[direction]] == 0)) {
	  vol[idx] = D_BORDER;
	}
      }
    }
  }
  
  return true;
}
