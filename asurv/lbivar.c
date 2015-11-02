#include <stdio.h>
 
int correlation(argc, argv) int argc; void *argv[];{ 
  
  extern void mulvar_l_();  /* function prototype for the fortran funct */
  
  // Declare variables
  
  double      *X,*Y; // FORTRAN REAL*8
  short       *IND,*NTOT, *ND, *NYC;  // FORTRAN INTEGER
  double      *RETURNED_PROB,*RETURNED_ALPHA,*RETURNED_ALPHA_SIGMA; // FORTRAN REAL*8
  
  //Insure that the correct number of arguments were passed in (argc = 9)
 
  if(argc != 9) {    //Print an error message and return
    fprintf(stderr,"correlation : Incorrect number of arguments\n");
    return(0);  /* Signal an error */
  }
  
  // Cast the pointer in argv to the pointer variables
  
  X                    = (double *) argv[0];
  Y                    = (double *) argv[1];
  IND                  = (short *) argv[2];
  NTOT                 = (short *) argv[3];
  ND                   = (short *) argv[4];
  NYC                  = (short *) argv[5];
  RETURNED_PROB        = (double *) argv[6];
  RETURNED_ALPHA       = (double *) argv[7];
  RETURNED_ALPHA_SIGMA = (double *) argv[8];

  mulvar_l_(X,Y,IND,NTOT,ND,NYC,RETURNED_PROB,RETURNED_ALPHA,RETURNED_ALPHA_SIGMA);

  return(1);
}
