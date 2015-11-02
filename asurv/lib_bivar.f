C      **********************************************************************
C      ********************** SUBROUTINE BIVAR ******************************
C      **********************************************************************
C
C       SUBROUTINE BIVAR(IBACK)
          
C
C      *   CORRELATION AND REGRESSION                                       *
C      *   PARAMETERS                                                       *
C      *      MVAR     :  THE MAXIMUM NUMBER OF VARIABLES (COLUMNS)         *
C      *                   ALLOWED IN A DATA SET.                           *
C      *      NDAT     :  THE MAXIMUM NUMBER OF DATA POINTS (ROWS)          *
C      *                   ALLOWED IN A DATA SET.                           *
C      *      IBIN     :  THE DIMENSION SIZE FOR BINS - USED IN SCHMITT'S   *
C      *                   PROCEDURE COMPUTATIONS.                          *
C      *      LENG     :  (MVAR+1)+NTOT - USED IN EM ALGORITHM COMPUTATIONS *
C      *      LEGWRK   :  (MVAR+1)*NTOT - USED IN EM ALGORITHM COMPUTATIONS *
C      *   INPUT                                                            *
C      *      FILE     :  NAME OF DATA FILE (9 LETTERS)                     *
C      *      TITLE    :  TITLE OF THE PROBLEM (80 LETTERS)                 *
C      *      NVAR     :  THE ACTUAL NUMBER OF VARIABLES IN THE DATA SET    *
C      *      NTOT     :  THE ACTUAL NUMBER OF DATA POINTS IN THE DATA SET  *
C      *      ICOL     :  INDICATOR OF VARIABLE (<=NVAR)                    *
C      *                  IF A MULTIVARIATE PROBLEM IS NEEDED, SET ICOL=0   *
C      *      COLM     :  NAME OF THE INDEPENDENT VARIABLE                  *
C      *      YNAME    :  NAME OF THE DEPENDENT VARIABLE                    *
C      *      COMMAND  :  NAME OF THE "COMMAND" FILE                        *
C      *      OUTPUT   :  NAME OF THE OUTPUT FILE                           *
C      *      IND(1,I) :  INDICATOR OF CENSORING                            *
C      *                    IF =0,  DETECTED                                *
C      *                       =1,  Y LOWER LIMIT                           *
C      *                       =2,  X LOWER LIMIT                           *
C      *                       =3,  DOUBLE LOWER LIMIT                      *
C      *                       =4,  X UPPER LIMIT AND Y LOWER LIMIT         *
C      *                  FOR THE UPPER LIMITS, CHANGE THE SIGN             *
C      *                  2, 3, AND 4 CAN BE USED ONLY IN GEN. KENDALL'S    *
C      *                  TAU, GEN. SPEARMAN'S RHO, AND SCHMITT'S METHOD    *
C      *      X(J,I)   :  INDEPENDENT VARIABLES                             *
C      *      Y(I)     :  DEPENDENT VARIABLE                                *
C      *     IPROG(I)  :  INDICATOR OF METHODS                              *
C      *     NOTEST    :  NUMBERS OF TEST                                   *
C      *  INPUT FOR EM ALGORITHM                                            *
C      *      TOL      :  TOLERANCE (DEFAULT 1.0E-5)                        *
C      *      MAX      :  MAXIMUM ITERATION (DEFAULT 20)                    *
C      *      IBET     :  IF 0, NO DEPENDENT VARIABLE IS CONFINED BETWEEN   *
C      *                        TWO VALUES                                  *
C      *                     1, THERE ARE SOME DEPENDENT VARIABLE WHICH     *
C      *                        ARE CONFINED BETWEEN TWO VALUES             *
C      *    ALPHA(K)   :  INITIAL ESTIMATE OF REGRESSION COEFFICIENTS       *
C      *  INPUTS FOR BUCKLEY-JAMES METHOD                                   *
C      *      TOL1     :  TOLERANCE (DEFAULT 1.0E-5)                        *
C      *      MAX1     :  MAXIMUM ITERATION (DEFAULT 20)                    *
C      *  INPUTS FOR SCHMITT'S BINNING METHOD                               *
C      *      MX       :  BIN NUMBER OF X AXES                              *
C      *      MY       :  BIN NUMBER OF Y AXES                              *
C      *      TOL3     :  TOLERANCE LEVEL (DEFAULT 1.0E-5)                  *
C      *      MAX3     :  MAXIMUM ITERATION (DEFAULT 20)                    *
C      *      XBIN     :  BIN SIZE FOR X AXES                               *
C      *      YBIN     :  BIN SIZE FOR Y AXES                               *
C      *      XORG     :  ORIGIN OF X AXES                                  *
C      *      YORG     :  ORIGIN OF Y AXES                                  *
C      *      ISKIP    :  IF 0, THE PROGRAM WILL PROVIDE XBIN, YBIN, XORG,  *
C      *                        AND YORG.                                   *
C      *                    >0, THESE VALUES MUST BE PROVIDED BY THE USER   *
C      *      IPIRNT   :  IF 0, NO TWO DIMENSIONAL K-M ESTIMATOR WILL BE    *
C      *                        PRINTED                                     *
C      *                    >0, TWO DIMENSIONAL K-M ESTIMATOR WILL BE       *
C      *                        PRINTED                                     *
C      *                                                                    *
C      *    WORKING VARIABLES AND ARRAYS:                                   *
C      *      NTOT     :  NUMBER OF DATA POINTS                             *
C      *      ND       :  NUMBER OF DETECTED POINTS                         *
C      *      NC1      :  NUMBER OF Y LOWER LIMITS                          *
C      *      NC2      :  NUMBER OF X LOWER LIMITS                          * 
C      *      NC3      :  NUMBER OF DOUBLE LOWER LIMITS                     * 
C      *      NC4      :  NUMBER OF Y LOWER AND X UPPER LIMITS              *
C      *      NC5      :  NUMBER OF Y UPPER LIMITS                          *
C      *      NC6      :  NUMBER OF X UPPER LIMITS                          *
C      *      NC7      :  NUMBER OF DOUBLE UPPER LIMITS                     *
C      *      NC8      :  NUMBER OF Y UPPER AND X LOWER LIMITS              *
C      *      ICENS    :  IF 0, CENSORING IS MIXED                          *
C      *                     1, CENSORING IS Y LOWER LIMITS ONLY            *
C      *                    -1, CENSORING IS Y UPPER LIMITS ONLY            *
C      *      NYC      :  NC1+NC2                                           *
C      *      NXC      :  NC3+NC4                                           *
C      *      NBC      :  NC5+NC6+NC7+NC8                                   *
C      *      IDO      :  NXC+NBC                                           *
C      *      IMUL     :  INDICATOR OF MULTIVARIATE PROBLEM                 *
C      *      XX(J,I)  :  =X(ICOL,I), EXCEPT FOR MULTI INDEPENDENT VARIABLE *
C      *                  CASE (J=1,NVAR).                                  *
C      *      IND2(I)  :  =IND(1,I)                                         *
C      *                                                                    *
C      *  OUTPUT                                                            *
C      *     COXREG                                                         *
C      *      CHI      : GLOBAL CHI-SQUARE                                  *
C      *      PROB     : PROBABILITY FOR NULL                               *
C      *     BHK (GNERALIZED KENDALL'S TAU)                                 *
C      *       Z       : DEVIATION                                          *
C      *      PROB     : PROBABILITY FOR NULL                               *
C      *     EM ALGORITHM                                                   *
C      *      ALPHA(K) : LINEAR REGRESSION COEFFICIENTS  (K=1,NVAR+1)       *
C      *     ALPHA(K+2): STANDARD DEVIATION                                 *
C      *     SIGMAA(K) : ERROR                                              *
C      *     ITE       : NUMBER OF ITERATION                                *
C      *     BUCKLEY-JAMES                                                  *
C      *      ALPHA(K) : LINEAR REGRESSION COEFFICIENTS (K=1,NVAR+1)        *
C      *     ALPHA(K+2): STANDARD DEVIATION                                 *
C      *     SIGMAA(K) : ERROR                                              *
C      *     SCHMITT                                                        *
C      *     ALPHA     : INTERCEPT COEFFICIENT                              *
C      *     BETA      : SLOPE COEFFICIENT                                  *
C      *   *****  ALL OUTPUTS ARE INSIDE OF EACH SUBROUTINE                 *
C      *                                                                    *
C      *   SUBROUTINES                                                      *
C      *     DATA1, DATREG, DATA2, MULVAR                                   *
C      *                                                                    *

      IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)

C      *   THIS PARAMETER STATEMENT AND THE ONE IN UNIVAR.F ARE THE ONLY    *
C      *   STATEMENTS THAT NEED TO BE ADJUSTED IF THE USER WISHES TO        *
C      *   ANALYZE DATA SETS OF MORE THAN 500 OBSERVATIONS OR MORE THAN     *
C      *   VARIABLES.                                                       *

C  **************************************************************************
       PARAMETER(MVAR=4, NDAT=2000, IBIN=50)
C  **************************************************************************

       CHARACTER*1 CHECK,CHAR(4,10)
       CHARACTER*7 BB(10),YY
       CHARACTER*9 FILE,COMMAND,OUTPUT,COLM
       CHARACTER*9 YNAME,DUMMY1
       CHARACTER*80 TITLE 
       
       DIMENSION IND(MVAR,NDAT),X(MVAR,NDAT),Y(NDAT)
       DIMENSION IPROG(6),IIND(NDAT)

C       DIMENSION DWRK1(MVAR,NDAT),DWRK2(MVAR,NDAT)
C       DIMENSION DWRK3(MVAR,NDAT),DWRK4(MVAR,NDAT)
C       DIMENSION DWRK5(MVAR,NDAT),DWRK6(MVAR,NDAT)
C       DIMENSION DWRK8(MVAR,NDAT)
       
C       DIMENSION EWRK1(4,4),RWRK1(NDAT,MVAR)

C       DIMENSION AWRK(5,IBIN)
C       DIMENSION WWRK1((MVAR+1)+NDAT)
C       DIMENSION WWRK2((MVAR+1)+NDAT)
C       DIMENSION VWRK1((MVAR+1)*NDAT)
C       DIMENSION VWRK2((MVAR+1)*NDAT)

C       DIMENSION WRK1(NDAT),WRK2(NDAT),WRK3(NDAT),WRK4(NDAT)
C       DIMENSION WRK5(NDAT),WRK6(NDAT),WRK7(NDAT),WRK8(NDAT)
C       DIMENSION WRK9(NDAT),WRK10(NDAT),WRK11(NDAT)
C       DIMENSION WRK12(NDAT)

C       DIMENSION SWRK1(MVAR),SWRK2(MVAR),SWRK3(MVAR)
C       DIMENSION SWRK4(MVAR),SWRK5(MVAR),SWRK6(MVAR)
C       DIMENSION SWRK7(MVAR),SWRK8(MVAR),SWRK9(MVAR)
C       DIMENSION SWRK10(MVAR),SWRK11(MVAR),SWRK17(MVAR)

C       DIMENSION LWRK1(MVAR,NDAT), LWRK2(MVAR,NDAT)
C       DIMENSION LWRK3(MVAR,NDAT)
C       DIMENSION IWRK1(NDAT),IWRK2(NDAT),IWRK3(NDAT)
C       DIMENSION IWRK4(NDAT),IWRK5(NDAT),IWRK6(NDAT)
C       DIMENSION IWRK7(NDAT),IWRK8(NDAT)

C       DIMENSION IWRK9(NDAT),CWRK1(IBIN),CWRK2(IBIN)

C       DIMENSION IBWRK1(IBIN,IBIN),IBWRK2(IBIN,IBIN)
C       DIMENSION IBWRK3(IBIN,IBIN),IBWRK4(IBIN,IBIN)
C       DIMENSION IBWRK5(IBIN,IBIN),IBWRK6(IBIN,IBIN)
C       DIMENSION IBWRK7(IBIN,IBIN),IBWRK8(IBIN,IBIN)
C       DIMENSION IBWRK9(IBIN,IBIN)
C       DIMENSION BWRK1(IBIN,IBIN),BWRK2(IBIN,IBIN)

       DIMENSION RETURNED_PROB(3), RETURNED_ALPHA(2,3)
       DIMENSION RETURNED_ALPHA_SIGMA(2,3)
       
C       LENG = (MVAR+1)+NDAT
C       LEGWRK = (MVAR+1)*NDAT

       DO 5000 K=1,10
       BB(K)='     X '
 5000  CONTINUE
       YY='     Y '
C
C 6000  PRINT *
C       PRINT *
C       PRINT *,'      CORRELATION AND REGRESSION CALCULATIONS'
C       PRINT *
C       PRINT *,' CORRELATION OPTIONS        LINEAR REGRESSION OPTIONS'
C       PRINT *,' 1. COX HAZARD MODEL        4. EM ALGORITHM WITH '
C       PRINT *,'                               NORMAL DISTRIBUTION'
C       PRINT *,' 2. GEN. KENDALL`S TAU      5. BUCKLEY-JAMES METHOD'
C       PRINT *,' 3. GEN. SPEARMAN`S RHO     6. SCHMITT`S BINNING METHOD'
C       PRINT *
C       PRINT *,' DATA SETS WITH CEDNSORING IN ONLY ONE DIRECTION OF THE'
C       PRINT *,' DEPENDENT VARIABLE CAN USE ALL METHODS.'
C       PRINT *
C       PRINT *,' DATA SETS WITH SEVERAL INDEPENDENT AND ONE DEPENDENT'
C       PRINT *,' VARIABLE CAN USE ONLY THE COX PROPORTIONAL HAZARD'
C       PRINT *,' MODEL,EM ALGORITHM, OR BUCKLEY-JAMES METHOD.  ONLY'
C       PRINT *,' ONE TYPE OF CENSORING IN THE DEPENDENT VARIABLE IS'
C       PRINT *,' ALLOWED.'
C       PRINT *
C       PRINT *,' IF YOUR DATA SET HAS CENSORED POINTS IN THE '
C       PRINT *,' INDEPENDENT VARIABLE AND/OR DUAL CENSORED POINTS,'
C       PRINT *,' YOU CAN USE ONLY THE GEN. KENDALL`S TAU OR GEN.'
C       PRINT *,' SPEARMAN`S RHO CORRELATION COEFFICIENT, OR'
C       PRINT *,' SCHMITT`S BINNED LINEAR REGRESSION.'
C       PRINT *
C 6010  PRINT *
CC
C      *  CHECK WHETHER THE USER WANTS TO USE COMMAND FILE INPUTS. IF SO,   *
C      *  GO TO 6660                                                        *
C
   50  FORMAT(A1)
 1380  FORMAT(A9)
C
       OUTPUT='         '
       ICOMM=0
       ICOL=1
C       PRINT *,'DO YOU WANT TO READ ALL INFORMATION'
C       WRITE(6,6020)
C 6020  FORMAT('      FROM A COMMAND FILE (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6660
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6025
C       GOTO 6010
C
C      *          READ FROM THE TERMINAL                                    *
C
C 6025  PRINT *
C       PRINT *,'          OK, LET US READ FROM THE TERMINAL '
C
C      *           READ TITLE                                               *
C
C 6030  PRINT *
C       WRITE(6,6040)
C 6040  FORMAT('WHAT IS THE TITLE OF THE PROBLEM ? ')
C       READ(5,6050) TITLE
C 6050  FORMAT(A80)
C
C      *           READ DATA FILE NAME                                      *
C
C 6051  PRINT *
C       WRITE(6,6052)
C 6052  FORMAT('WHAT IS THE DATA FILE NAME ? ')
C       READ(5,1380) FILE
C
C      *           READ NUMBER OF INDEPENDENT VARIABLES                     *
C
C 6060  PRINT *
C       WRITE(6,6070)
C 6070  FORMAT('HOW MANY INDEPENDENT VARIABLES DO YOU HAVE ? ')
C       CALL DATA1(NVAR)
C       IF((NVAR.GE.1).AND.(NVAR.LE.MVAR-2)) GOTO 6080
C       PRINT *
C       PRINT *,'    YOUR CHOICE IS NOT ACCEPTABLE. PLEASE TYPE IT AGAIN'
C       GOTO 6060
C
C      *    CALL SUBROUTINE "DATREG" TO READ DATA                           *
C
C 6080  CALL DATREG(NVAR,IND,X,Y,FILE,NTOT,ND,NC1,NC2,NC3,NC4,NC5,
C     +             NC6,NC7,NC8,ICENS,NYC,NXC,NBC,MVAR,NDAT)
C       DO 6090 I = 1, NTOT
C            IIND(I) = IND(1,I)    
C 6090  CONTINUE
CC
CC      *        CHECK WHICH METHODS THE USER CAN USE                        *
CC
C       IDC=NXC+NBC
C       IF((NVAR.EQ.1).AND.(IDC.EQ.0)) GOTO 6530
C       IF((NVAR.NE.1).AND.(IDC.NE.0)) GOTO 6340
C       IF((NVAR.EQ.1).AND.(IDC.NE.0)) GOTO 6400
C 6170  PRINT *
C       WRITE(6,6180)
C 6180  FORMAT('IS THIS A MULTIVARIATE PROBLEM  (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6220
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6340
C       GOTO 6170
CC
CC      *      DATA SET WITH MORE THAN ONE INDEPENDENT VARIABLES             *
CC
C 6220  PRINT *
C       PRINT *,'           YOU CAN USE THE NEXT METHODS   '
C       PRINT *
C       PRINT *,'         1. COX HAZARD METHOD'
C       PRINT *,'         4. EM ALGORITHM WITH NORMAL DISTRIBUTION'
C       PRINT *,'         5. BUCKLEY-JAMES METHOD'
CC
C       ICOL=0
C       J=1
C 6230  PRINT *
C       WRITE(6,6240)
C 6240  FORMAT('WHICH METHOD DO YOU WANT TO USE ? ')
C       CALL DATA1(IPROG(J))
C       IF((IPROG(J).EQ.1).OR.(IPROG(J).EQ.4).OR.(IPROG(J).EQ.5))
C     +                                                    GOTO 6245
C       GOTO 6230
C 6245  IF(J.EQ.1) GOTO 6260
C       J1=J-1
C       DO 6250 K=1,J1
C       IF(IPROG(K).NE.IPROG(J)) GOTO 6250
C       PRINT *
C       PRINT *,' YOU ALREADY CHOSE THAT METHOD.'
C       PRINT *,'  PLEASE CHOOSE ANOTHER ONE'
C       GOTO 6230
C 6250  CONTINUE
C 6260  IF(J.GE.3) GOTO 6280
C 6265  PRINT *
C       WRITE(6,6270)
C 6270  FORMAT('DO YOU WANT TO USE ANY OTHER METHOD (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') J=J+1
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6230
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6280
C       GOTO 6265
CC
C 6280  NOTEST=J
C       GOTO 6604
CC
CC      *    FOR THE CASE  THAT THE DATA SET CONTAINS MIXED CENSORING        *
CC      *    (THAT IS, UPPER AND LOWER LIMITS SIMULTANEOUSLY AND/OR          *
CC      *    CENSORING IN BOTH VARIABLES).                                   *
CC
C 6340  PRINT *
C       WRITE(6,6350)
C 6350  FORMAT('WHICH INDEPENDENT VARIABLE DO YOU WANT TO USE ? ')
C       CALL DATA1(ICOL)
C       IF(ICOL.GT.NVAR) GOTO 6340
C       IF(ICOL.LE.0) GOTO 6340
CC
C 6400  IF(NBC.EQ.0) GOTO 6530
C       J=1
C       PRINT *
C       PRINT *,'          YOU CAN USE THE FOLLOWING METHODS'
C       PRINT *,'          2. GEN. KENDALL`S TAU METHOD'
C       PRINT *,'          3. GEN. SPEARMAN`S RHO METHOD'
C       PRINT *,'          6. SCHMITT`S BINNING METHOD'
C 6410  PRINT *
C       WRITE(6,6420)
C 6420  FORMAT('WHICH METHOD DO YOU WANT TO USE ? ')
C       CALL DATA1(IPROG(J))
C       IF((IPROG(J).EQ.2).OR.(IPROG(J).EQ.3).OR.(IPROG(J).EQ.6)) 
C     +                                                     GOTO 6425
C       GOTO 6410
C 6425  IF(J.EQ.1) GOTO 6440
C       J1=J-1
C       DO 6430 K=1,J1
C       IF(IPROG(K).NE.IPROG(J)) GOTO 6430
C       PRINT *
C       PRINT *,'   YOU ALREADY CHOSE THAT METHOD.'
C       PRINT *,'    PLEASE CHOOSE THE OTHER ONE'
C       GOTO 6410
C 6430  CONTINUE
C 6440  IF(J.EQ.3) GOTO 6600
C 6450  PRINT *
C       WRITE(6,6460)
C 6460  FORMAT('DO YOU WANT TO USE THE OTHER METHOD, TOO (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') J=J+1
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6410
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6600
C       GOTO 6450
CC
CC      *  FOR THE CASE THAT THE DATA SET CONTAINS ONE INDEPENDENT AND ONE   *
CC      *  DEPENDENT VARIABLES AND ONE KIND OF CENSORING IN THE DEPENDENT    *
CC      *  VARIABLE.                                                         *
CC
C 6530  PRINT *
C       PRINT *,'     YOU CAN USE THE FOLLOWING METHODS'
C       PRINT *
C       PRINT *,'     1. COX HAZARD MODEL    4. EM ALGORITHM WITH'
C       PRINT *,'                               NORMAL DISTRIBUTION'
C       PRINT *,'     2. KENDALL`S TAU       5. BUCKLEY-JAMES REGRESSION'
C       PRINT *,'     3. SPEARMAN`S RHO',
C     +         '      6. SCHMITT`S BINNED REGRESSION'
C       PRINT *
C       J=1
C 6540  PRINT *
C       WRITE(6,6550)
C 6550  FORMAT('WHICH METHOD DO YOU WANT TO USE ? ')
C       CALL DATA1(IPROG(J))
C       IF(IPROG(J).LT.1) GOTO 6540
C       IF(IPROG(J).GT.6) GOTO 6570
C       IF(J.EQ.1) GOTO 6580
C       J1=J-1
C       DO 6560 K=1,J1
C       IF(IPROG(K).NE.IPROG(J)) GOTO 6560
C       PRINT *
C       PRINT *,'   YOU ALREADY CHOSE THAT METHOD.'
C       PRINT *,'    PLEASE CHOOSE ANOTHER ONE.'
C       GOTO 6540
C 6560  CONTINUE
C 6570  IF(J.GE.6) GOTO 6600
C 6580  PRINT *
C       WRITE(6,6590)
C 6590  FORMAT('DO YOU WANT TO USE ANOTHER METHOD (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') J=J+1
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6540
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 6600
C       GOTO 6580
CC
C 6600  NOTEST=J
CC
CC      *        NAME THE VARIABLES                                          *
CC
C 6601  PRINT *
C       PRINT *,'          PLEASE NAME THE VARIABLES : '
C       PRINT *
C       WRITE(6,6602)
C 6602  FORMAT('WHAT IS THE NAME OF THE INDEPENDENT VARIABLE ? ')
C       READ(5,1380) COLM
CC
C       PRINT *
C       WRITE(6,6603)
C 6603  FORMAT('WHAT IS THE NAME OF THE DEPENDENT VARIABLE ? ')
C       READ(5,1380) YNAME
CC
C 6604  PRINT *
C       WRITE(6,6605)
C 6605  FORMAT('DO YOU WANT TO PRINT THE ORIGINAL DATA  (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IDATA=1
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') IDATA=0
C       IF((CHECK.EQ.'Y').OR.(CHECK.EQ.'N')) GOTO 6609
C       IF((CHECK.EQ.'y').OR.(CHECK.EQ.'n')) GOTO 6609
C       GOTO 6604
CC
CC      *   CHECK WHETHER THE USER WANT TO SAVE THE RESULT IN AN OUTPUT FILE *
CC
C 6609  PRINT *
C       PRINT *,'DO YOU WANT TO SAVE THE RESULT '
C       WRITE(6,6610)
C 6610  FORMAT('     IN AN OUTPUT FILE (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 6620
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 7116
C       GOTO 6609
C 6620  PRINT *
C       WRITE(6,6630)
C 6630  FORMAT('WHAT IS THE NAME OF THE OUTPUT FILE ? ')
C       READ(5,1380) OUTPUT
C       GOTO 7116
CC
C
C
C      *    USE "COMMAND" FILE FOR INPUTS                                   *
C
C 6660  PRINT *
C       WRITE(6,6670)
C 6670  FORMAT('WHAT IS THE NAME OF THE COMMAND FILE ? ')
C       READ(5,1380) COMMAND
C
       COMMAND="test5.dat"

 6700  OPEN(UNIT=50, FILE=COMMAND, STATUS='OLD', FORM='FORMATTED')
       ICOMM=1
C
C      *   READ TITLE OF THE PROBLEM ; NAME OF THE DATA FILE                *
C
       READ(50,6710) TITLE
 6710  FORMAT(A80)
       READ(50,1380) FILE
C
C      * READ NUMBER OF VARIABLES; WHICH VARIABLE WILL BE USED; AND HOW     *
C      * MANY METHODS THE USER WANTS TO USE.                                *
C
       READ(50,6720) ((CHAR(I,J),I=1,4),J=1,3)
 6720  FORMAT(20A1)
       CALL DATA2(CHAR,1,3,NVAR,LIND)
       IF(LIND.EQ.0) GOTO 6750
 6730  PRINT *
       PRINT *,'   TOTAL NUMBER OF INDEPENDENT VARIABLES IS NOT CLEAR.'
       STOP
 6750  IF(NVAR.LT.1) GOTO 6730
       IF(NVAR.NE.1) GOTO 6760
       ICOL=1
       GOTO 6915
C 
 6760  CALL DATA2(CHAR,50,3,ICOL,LIND)
       IF(LIND.EQ.0) GOTO 6900
 6860  PRINT *
       PRINT *,'     THE CHOICE OF THE VARIABLE IS NOT CLEAR'
       STOP
 6900  IF(ICOL.LE.0) GOTO 6860
       IF(ICOL.GT.NVAR) GOTO 6860
C
C      *         CHOICE OF THE METHODS                                      *
C
 6915  CALL DATA2(CHAR,3,3,NOTEST,LIND)
       IF(LIND.EQ.0) GOTO 6950
 6930  PRINT *
       PRINT *,'    IT IS NOT CLEAR HOW MANY METHODS YOU WANT TO USE '
       STOP
 6950  IF(NOTEST.LE.0) GOTO 6930
       IF(NOTEST.GT.6) GOTO 6930
C
       READ(50,6960) ((CHAR(I,J),I=1,4),J=1,NOTEST)
 6960  FORMAT(30A1)
       DO 7020 I=1,NOTEST
       CALL DATA2(CHAR,I,NOTEST,IPROG(I),LIND)
       IF(LIND.EQ.0) GOTO 7010
 6970  PRINT *
       IF(I.EQ.1) PRINT *,'     FIRST PROGRAM NUMBER IS NOT CLEAR'
       IF(I.EQ.2) PRINT *,'     SECOND PROGRAM NUMBER IS NOT CLEAR'
       IF(I.GE.3) WRITE(6,6780) I
 6780  FORMAT(5X,I4,'-TH PROGRAM NUMBER IS NOT CLEAR')
       STOP
 7010  IF(IPROG(I).LE.0) GOTO 6970
       IF(IPROG(I).GT.6) GOTO 6970
 7020  CONTINUE
C
C      *  READ NAMES OF THE INDEPENDENT AND DEPENDENT VARIABLES. IF THE     *
C      *  PROBLEM HAS MULTI-INDEPENDENT VARIABLES, THESE NAMES WIL BE       *
C      *  IGNORED.                                                          *
C
       READ(50,7022) COLM,YNAME
 7022  FORMAT(2A9)
C
       CLOSE(UNIT=50,STATUS='KEEP')
C
C      *      CALL SUBROUTINE "DATREG" TO READ DATA                         *
C
       CALL DATREG(NVAR,IND,X,Y,FILE,NTOT,ND,NC1,NC2,NC3,NC4,NC5,
     +             NC6,NC7,NC8,ICENS,NYC,NXC,NBC,MVAR,NDAT)
C
       OPEN(UNIT=50,FILE=COMMAND,STATUS='OLD',FORM='FORMATTED')
C
C      * THE NEXT SEVERAL LINES READ IN DUMMY VALUES TO PREVENT READING     *
C      * THE COMMANDS A SECOND TIME.                                        *
C
       READ(50,1380) DUMMY1
       READ(50,1380) DUMMY1
       READ(50,7029) IDUMMY
       READ(50,7029) IDUMMY
       READ(50,1380) DUMMY1
 7029  FORMAT(I4)
C
C      *  CHECK WHETHER THE ASSIGNED METHODS CAN BE USED FOR THE DATA       *
C
       IF(NVAR.GE.2) GOTO 7070
       IF((NXC.EQ.0).AND.(NBC.EQ.0)) GOTO 7110
C
C      *   THE CASE WITH MIXED CENSORING IN DATA                            *
C
       I=1
 7030  IF(IPROG(I).NE.1) GOTO 7040
       PRINT *
       PRINT *,'      YOU CANNOT USE COX HAZARD MODEL FOR THIS DATA SET'
       PRINT *,'      THIS METHOD WILL BE IGNORED'
       IPROG(I)=-9
 7040  IF(IPROG(I).NE.4) GOTO 7050
       PRINT *
       PRINT *,'       YOU CANNOT USE EM ALGORITHM FOR THIS DATA SET'
       PRINT *,'       THIS METHOD WILL BE IGNORED'
       IPROG(I)=-9
 7050  IF(IPROG(I).NE.5) GOTO 7060
       PRINT *
       PRINT *,'       YOU CANNOT USE BUCKLEY-JAMES METHOD FOR THIS'
       PRINT *,'        DATA SET. THIS METHOD WILL BE IGNORED'
       IPROG(I)=-9
 7060  IF(I.GE.NOTEST) GOTO 7110
       I=I+1
       GOTO 7030
C
C      *     THE CASE WITH MORE THAN ONE INDEPENDENT VARIABLES              *
C
 7070  I=1
 7080  IF(IPROG(I).NE.2) GOTO 7085
       PRINT *
       PRINT *,'         YOU CANNOT USE THE KENDALL`S TAU METHOD FOR'
       PRINT *,'         THIS DATA SET'
       PRINT *,'         THIS METHOD WILL BE IGNORED'
       IPROG(I)=-9
 7085  IF(IPROG(I).NE.3) GOTO 7090
       PRINT *
       PRINT *,'       YOU CANNOT USE SPEARMAN`S RHO FOR THIS DATA SET'
       PRINT *,'       THIS METHOD WILL BE IGNORED'
       IPROG(I)=-9
 7090  IF(IPROG(I).NE.6) GOTO 7100
       PRINT *
       PRINT *,'         YOU CANNOT USE SCHMITT`S BINNED REGRESSION'
       PRINT *,'         THIS METHOD WILL BE IGNORED'
       IPROG(I)=-9
 7100  IF(I.EQ.NOTEST) GOTO 7110
       I=I+1
       GOTO 7080
C
C      *        READ PRINT OUT INDICATOR FOR THE DATA                       *
C
 7110  READ(50,6960) (CHAR(I,1),I=1,4)
       CALL DATA2(CHAR,1,1,IDATA,LIND)
       IF(LIND.EQ.0) GOTO 7114
 7112  PRINT *
       PRINT *,'     THE VALUE FOR "IDATA" IS NOT ACCEPTABLE'
       STOP
 7114  IF((IDATA.EQ.0).OR.(IDATA.EQ.1)) GOTO 7115
       GOTO 7112
C
C      *       READ OUTPUT FILE NAME                                        *
C
 7115  READ(50,1380) OUTPUT
C
C      * CALL SUBROUTINE "MULVAR" TO COMPUTE CORRELATION/REGRESSION PROBLEMS*
C
C     7116 IF(OUTPUT .NE. '  ') OPEN(UNIT=60,FILE=OUTPUT, + STATUS='NEW'
C     THE FOLLOWING LINE SHOULD BE COMMENTED OUT ON ALL MACHINES OTHER
C     THAN VAX/VMS MACHINES.  + ,CARRIAGECONTROL='LIST' + ) * ADJUST THE
C     INPUT DATA FORMAT *
C     
       CALL  MULVAR_L(X,Y,IND,NTOT, ND, NYC, RETURNED_PROB,
     $      RETURNED_ALPHA,RETURNED_ALPHA_SIGMA
C     +              LENG,LEGWRK)
C     +              NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8,
C     +              FILE,YNAME,TITLE
C     +              COLM,IPROG,NOTEST,ICOMM, DWRK1,IWRK9,SWRK17,DWRK2,
C     +              DWRK3,DWRK4,DWRK5,DWRK6,DWRK8,RWRK1,
C     +              EWRK1,AWRK,WWRK1,WWRK2,
C     +              VWRK1,VWRK2,WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,
C     +              WRK7,WRK8,WRK9,WRK10,WRK11,WRK12,
C     +              SWRK1,SWRK2,SWRK3,SWRK4,SWRK5,SWRK6,SWRK7,
C     +              SWRK8,SWRK9,SWRK10,SWRK11,LWRK1,LWRK2,LWRK3,
C     +              IWRK1,IWRK2,IWRK3,IWRK4,IWRK5,IWRK6,IWRK7,
C     +              IWRK8,CWRK1,CWRK2,IBWRK1,IBWRK2,IBWRK3,
C     +              IBWRK4,IBWRK5,IBWRK6,IBWRK7,IBWRK8,IBWRK9,
C     +              BWRK1,BWRK2
     +              )


       PRINT *,"Proba of non correlation ", RETURNED_PROB
       DO 20 J=1,2
          DO 10 I=1,3
             PRINT *,"ALPHA(",J,",",I,") : ", RETURNED_ALPHA(J,I), " +-"
     $            ,RETURNED_ALPHA_SIGMA(J,I)
 10       CONTINUE
 20    CONTINUE

C
C       IF(IDATA.EQ.0) GOTO 7219
C       IF(OUTPUT.NE.'         ') WRITE(60,7140)
C       IF(OUTPUT.NE.'         ') WRITE(60,7117) FILE
C       IF(OUTPUT.EQ.'         ') PRINT 7140
C       IF(OUTPUT.EQ.'         ') PRINT 7117, FILE
C 7117  FORMAT(5X,'INPUT DATA FILE : ',A9)
C       IF(ICOL.NE.0) GOTO 7130
C       IF(OUTPUT.NE.'         ') WRITE(60,7118) (BB(K),K,K=1,NVAR),YY
C       IF(OUTPUT.EQ.'         ') PRINT 7118,(BB(K),K,K=1,NVAR),YY
C 7118  FORMAT(4X,'CENSORSHIP',12(A7,I2,1X))
C       DO 7119 I=1,NTOT
C       IF(OUTPUT.NE.'         ') WRITE(60,7120) IIND(I),
C     +                                      (X(J,I),J=1,NVAR),Y(I)
C       IF(OUTPUT.EQ.'         ') PRINT 7120,IIND(I),
C     +                                      (X(J,I),J=1,NVAR),Y(I)
C 7119  CONTINUE
C 7120  FORMAT(7X,I4,3X,10F10.3)
C       GOTO 7219
C 7130  IF(OUTPUT.NE.'         ') WRITE(60,7133)
C       IF(OUTPUT.EQ.'         ') PRINT 7133
C 7133  FORMAT(5X,' CENSORSHIP        X        Y')
C       DO 7134 I=1,NTOT
C       IF(OUTPUT.NE.'         ') WRITE(60,7135) IIND(I),X(ICOL,I),Y(I)
C       IF(OUTPUT.EQ.'         ') PRINT 7135,IIND(I),X(ICOL,I),Y(I)
C 7134  CONTINUE
C 7135  FORMAT(8X,I4,5X,2F10.3)
C 7140  FORMAT('     ')
CC
C 7219  IF(OUTPUT .NE. '         ') CLOSE(UNIT=60)
C       PRINT *
C       PRINT *
C       PRINT *,'    COMPUTATIONS FOR CORRELATION/REGRESSION'
C       PRINT *,'                          PROBLEMS ARE FINISHED'
C 7220  PRINT *
C       WRITE(6,7230)
C 7230  FORMAT('DO YOU WANT TO DO ANY OTHER ANALYSIS (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IBACK=1 
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') RETURN
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') STOP
C       GOTO 7220

       END

