#include <stdio.h>
 
int two_sample_test(argc, argv) int argc; void *argv[];{ 
  
  extern void twost_l_();  /* function prototype for the fortran funct */
  
  // Declare variables
  
  int i;
  float      *Z; // FORTRAN REAL*4
  short      *IND,*ISTA,*IS,*NG1,*NG2,*NTOT,*MVAR,*NDAT; // FORTRAN INTEGER*2
  float      *RETURNED_PROB; // FORTRAN REAL*4
  
  //Insure that the correct number of arguments were passed in (argc = 10)
 
  if(argc != 10) {    //Print an error message and return
    fprintf(stderr,"two_sample_test: Incorrect number of arguments\n");
    return(0);  /* Signal an error */
  }
  
  // Cast the pointer in argv to the pointer variables
  
  Z             = (float *) argv[0];
  IND           = (short *) argv[1];
  ISTA          = (short *) argv[2];
  IS            = (short *) argv[3];
  NG1           = (short *) argv[4];
  NG2           = (short *) argv[5];
  NTOT          = (short *) argv[6];
  MVAR          = (short *) argv[7];
  NDAT          = (short *) argv[8];
  RETURNED_PROB = (float *) argv[9];

  twost_l_(Z,IND,ISTA,IS,NG1,NG2,NTOT,MVAR,NDAT,RETURNED_PROB);

  return(1);
}
 
 
