C
C      **********************************************************************
C      ********************** SUBROUTINE UNIVAR *****************************
C      **********************************************************************
C
C       SUBROUTINE UNIVAR(IBACK)
C
C      *      NDAT     : DIMENSION DECLARATION                              *
C      *                                                                    *
C      *   UNIVARIATE PROBLEMS                                              *
C      *   PARAMETERS                                                       *
C      *       MVAR    :  THE MAXIMUM NUMBER OF VARIABLES ALLOWED IN A DATA *
C      *                    SET.                                            *
C      *       NDAT    :  THE MAXIMUM NUMBER OF DATA POINTS ALLOWED IN A    *
C      *                    DATA SET.                                       *
C      *       IBIN    :  THE DIMENSION SIZE FOR BINS - USED IN THE KAPLAN- *
C      *                    MEIER ESTIMATION PROCEDURE.                     *
C      *   COMMON FOR KAPLAN-MEIER AND TWO SAMPLE TESTS :                   *
C      *       FILE    :  NAME OF DATA FILE    (9 LETTERS)                  *
C      *       TITLE   :  TITLE OF THE PROBLEM (80 LETTERS)                 *
C      *       IUNI    :  INDICATOR OF PROBLEM                              *
C      *                    IF 1 : KAPLAN-MEIER ESTIMATOR                   *
C      *                    IF 2 : TWO-SAMPLE TESTS                         *
C      *       NTOT    :  THE ACTUAL NUMBER OF VARIABLES IN THE DATA SET    *
C      *       NVAR    :  THE ACTUAL NUMBER OF VARIABLES IN THE DATA SET    * 
C      *       ITEST   :  NUMBER OF VARIABLES TO BE TESTED (<=NVAR)         *
C      *       ICOL    :  INDICATOR OF SAMPLES                              *
C      *       COLM    :  NAME OF THE SAMPLE SETS                           *
C      *       COMMAND :  NAME OF COMMAND FILE                              *
C      *       OUTPUT  :  NAME OF OUTPUT FILE                               *
C      *       IND(I,J):  INDICATOR OF CENSORING (J-TH DATA POINTS OF I-TH  *
C      *                  VARIABLE)                                         *
C      *       X(I,J)  :  DATA POINTS                                       *
C      *                                                                    *
C      *    KAPLAN-MEIER ESTIMATOR                                          *
C      *       IKM     :  INDICATOR OF PRINTOUT                             *
C      *                    IF 0 : MEAN AND ERROR ONLY                      *
C      *                    IF 1 : MEAN, ERROR, SURVIVAL DISTRIBUTION       *
C      *                                                                    *
C      *    TWO-SAMPLE TESTS                                                *
C      *       IPROG(I):  INDICATOR OF TEST (I=1,5; OR 6 FOR EXIT)          *
C      *       NGROUP  :  NUMBER OF GROUPS                                  *
C      *      IGROUP(I):  INDICATOR OF GROUP                                *
C      *       IFIRST  :  INDICATOR OF FIRST GROUP                          *
C      *       ISECON  :  INDICATOR OF SECOND GROUP                         *
C      *       IFULL   :  INDICATOR OF PRINTOUT FOR THE I-TH COMBINATION    *
C      *       GROUP(I):  NAME OF THE GROUPS                                *
C      *        LKM    :  INDICATOR OF K-M ESTIMATOR                        *
C      *        IKM    :  INDICATOR OF PRINTOUT.                            *
C      *       NOTEST  :  NUMBER OF TESTS TO BE USED                        *
C      *       ISTA(I) :  INDICATOR OF GROUPS                               *
C      *                                                                    *
C      *    WRK VARIABLES AND ARRAYS                                        *
C      *       LCOMM   :  INDICATOR OF USE OF COMMAND FILE                  *
C      *                      IF 0, READ INFORMATION FROM THE TERMINAL      *
C      *                      IF 1, READ INFORMATION FROM THE COMMAND FILE  *
C      *       CHECK   :  READER OF Y/N QUESTIONS                           *
C      *       NTOT    :  NUMBER OF DATA POINTS                             *
C      *       NCHANGE :                                                    *
C      *       ICHANGE :                                                    *
C      *       FGROUP  :                                                    *
C      *       SGROUP  :                                                    *
C      *      CHAR(I,J):  READ-IN FORM OF SEVERAL INPUTS                    *
C      *     CTEST(I,1):  READ-IN FORM OF NOTEST                            *
C      *     CPROG(I,J):  READ-IN FORM OF IPROG(J)                          *
C      *      CCOL(I,J):  READ-IN FORM OF ICOL(J)                           *
C      *      IFIRST   :                                                    *
C      *      ISECON   :                                                    *
C      *      CF(I,1)  :  READ-IN FORM OF IFIRST                            *
C      *      CS(I,1)  :  READ-IN FORM OF ISECON                            *
C      *     CIKM1(I,J):  READ-IN FORM OF IKM1(J)                           *
C      *       N1      :  NUMBER OF DATA POINTS IN THE FIRST GROUP          *
C      *       N2      :  NUMBER OF DATA POINTS IN THE SECOND GROUP         *
C      *      JIND(I,J):  IND(I,J) IN FIRST OR SECOND GROUP                 *
C      *       Z(I,J)  :  X(I,J) IN FIRST OR SECOND GROUP                   *
C      *                                                                    *
C      *   NOTE:                                       *
C      *     ALL VARIABLES WRK_, IWRK_, DWRK1, AND SWRK1 ARE USED TO COLLECT*
C      *     ALL DIMENSION DECLARATIONS TO THE TOP LEVEL SUBROUTINE.        *
C      *     THEY DO NOT DIRECTLY AFFECT ANY OTHER SUBROUTINES. IF ONE NEEDS*
C      *     TO KNOW WHAT KIND OF WORK A VARIABLE DOES, GO TO THE LOWER     *
C      *     SUBROUTINES AND READ THE DESCRIPTION OF THE VARIABLE.          *
C      *                                                                    *
C      *   SUBROUTINES                                                      *
C      *     DATA1, DATA2, DATAIN, KMESTM, TWOST                            *
C      *                                                                    *
C
C
C      *                                                                    *
C      *  START UNIVARIATE PROBLEM: K-M ESTIMATOR OR TWO-SAMPLE TESTS       *
C
       IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
       
C      *   THIS PARAMETER STATEMENT AND THE ONE IN BIVAR.F ARE ALL THAT     *
C      *   NEEDS TO BE CHANGED IF THE USER WISHES TO WORK ON DATA SETS      *
C      *   OF MORE THAN 500 OBSERVATIONS OR WITH MORE THAN FOUR VARIABLES.  *

C       INTEGER*2 MVAR,NDAT,IBIN
C  **************************************************************************
       PARAMETER(MVAR=4, NDAT=2000, IBIN=50)
C  **************************************************************************

       INTEGER*2 IS,NG1,NG2,NTOT

       CHARACTER*1 CHECK,CHAR(4,10)
       CHARACTER*1 CF(4,1),CS(4,1),CIKM(4,1),CIS4(4,1)
       CHARACTER*9 FILE,COMMAND,OUTPUT,COLM,GROUP(10)
       CHARACTER*80 TITLE 
       DIMENSION RETURNED_PROB(5)
       DIMENSION IND(MVAR,NDAT),ISTA(NDAT)

       DIMENSION Z(MVAR,NDAT)
       DIMENSION X(MVAR,NDAT)
       DIMENSION JIND(MVAR,NDAT)
       DIMENSION JGROUP(MVAR),IGROUP(MVAR)

C     *     THE DIMENSIONS BELOW WERE ADDED TO COLLECT ALL DIMENSION        *
C     *     DECLARATIONS IN THE TOP SUBROUTNE                        *
C
       DIMENSION WRK1(NDAT),WRK2(NDAT),WRK3(NDAT),WRK4(NDAT)
       DIMENSION WRK5(NDAT),WRK6(NDAT),WRK7(NDAT),WRK8(NDAT)
       DIMENSION WRK9(NDAT),WRK10(NDAT),WRK11(NDAT)
       DIMENSION WRK12(NDAT),WRK13(NDAT),WRK14(NDAT)
       DIMENSION WRK15(NDAT),WRK16(NDAT),WRK17(NDAT)
       DIMENSION WRK18(NDAT),WRK19(NDAT),WRK20(NDAT)
       DIMENSION WRK21(NDAT),WRK22(NDAT),WRK23(NDAT)
       DIMENSION WRK24(NDAT),WRK25(NDAT),WRK26(NDAT)

       DIMENSION IWRK1(NDAT), IWRK2(NDAT), IWRK3(NDAT)
       DIMENSION BWRK1(IBIN),BWRK2(IBIN), BWRK3(IBIN)
       DIMENSION DWRK1(MVAR,NDAT), SWRK1(MVAR)
C
C


   50  FORMAT(A1)
C
C 1000  PRINT *
C       PRINT *,'    SELECT PROBLEM: ' 
C       PRINT *,'     1 KAPLAN-MEIER DISTRIBUTION '
C       PRINT *,'     2 TWO SAMPLE TESTS'
C       PRINT *,'     3 EXIT '
C       PRINT *
C       PRINT *,'    (IF YOU CHOOSE OPTION 2, YOU CAN STILL DO 1 LATER) '
C       PRINT *
C 1010  WRITE(6,1020)
C
C      *         SELECT PROBLEM                                             *
C
C 1020  FORMAT(' CHOICE? ')
C
C       CALL DATA1(IUNI)
       IUNI=2
C       IF((IUNI.EQ.1).OR.(IUNI.EQ.2).OR.(IUNI.EQ.3)) GOTO 1030
C       PRINT *
C       PRINT *,'      PLEASE TYPE ONCE MORE'
C       GOTO 1010
C
C 1030  IF(IUNI.EQ.3) STOP 
C
C 1120  IF(IUNI.EQ.2) GOTO 1140
C
C       PRINT *
C       PRINT *,'   ***  KAPLAN-MEIER ESTIMATOR  ***'
C       PRINT *
C       GOTO 1330
C
C
C      *     DISPLAY THE INFORMATION ABOUT TWO SAMPLE TESTS                 *
C
C 1140  PRINT *
C       PRINT *,'   ***     TWO-SAMPLE TESTS     ***'
C       PRINT *
C       PRINT *
C
C       LCOMM=1
C
C      *   CHECK WHETHER THE DATA NEEDS TO BE READ FROM A FILE              *
C
C 1330  PRINT *
C       PRINT *,'DO YOU WANT TO READ THE INPUTS'
C       WRITE(6,1340)
C 1340  FORMAT('     FROM A COMMAND FILE (Y/N)? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 2680
C       GOTO 2680
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 1350
C       GOTO 1330
C
C      *  READ INFORMATION COMMON TO K-M ESTIMATOR AND TWO SAMPLE TESTS     *
C
C 1350  PRINT *
C
C      *               READ NAME OF THE DATA FILE                           *
C
C 1360  PRINT *
C       WRITE(6,1370)
C 1370  FORMAT('WHAT IS THE DATA FILE NAME ? ')
C       READ(5,1380) FILE
 1380  FORMAT(A9)
C       PRINT *
C       WRITE(6,1400) FILE
C 1400  FORMAT(5X,'THE FILE NAME IS ',A9)
C       PRINT *
C
C      *              READ TITLE OF THE PROBLEM                             *
C
C 1410  PRINT *
C       WRITE(6,1420)
C 1420  FORMAT('WHAT IS THE PROBLEM TITLE? ')
C       READ(5,1430) TITLE 
 1430  FORMAT(A80)
C
C      *            READ THE NUMBER OF VARIABLES                            *
C
C       PRINT *
C       WRITE(6,1480)
C 1480  FORMAT('HOW MANY VARIABLES DO YOU HAVE? ')
C       CALL DATA1(NVAR)
C
C      *          CHECK WHICH VARIABLE SHOULD BE TESTED                     *
C
C       ICOL=1
C       IF(NVAR.EQ.1) GOTO 1840 
C 1550  PRINT *
C       PRINT *,'      WHICH VARIABLE DO YOU WANT TO TEST? '
C 1560  WRITE(6,1570)
C 1570  FORMAT(' VARIABLE NUMBER: ')
C       READ(5,1580) (CHAR(I,1),I=1,4)
C 1580  FORMAT(4A1)
C
C      *     CHECK IF THE CHOICE IS CORRECT                                 *
C
C       CALL DATA2(CHAR,1,1,ICOL,LIND)
C       IF(LIND.NE.0) PRINT *,
C     +               '    PLEASE TYPE IN THE VARIABLE NUMBER AGAIN'
C       IF(LIND.NE.0) GOTO 1550
C       IF(ICOL.LE.NVAR) GOTO 1840
C       PRINT *
C       PRINT *,'   THE NUMBER IS LARGER THAN THE NUMBER OF VARIABLES'
C       GOTO 1560
C
C      *                READ NAME OF THE VARIABLE                           *
C
C 1840  PRINT *
C       WRITE(6,1850) ICOL
C 1850  FORMAT('VARIABLE ',I4,' IS NAMED')
C       READ(5,1380) COLM
C
C      *   THE NEXT FEW LINES CONCERN ONLY 2-SAMPLE TESTS                   *
C      *   IF THE PROBLEM IS K-M ESTIMATION, GO TO LINE 2630                *
C
C      *          READ THE NUMBER OF GROUPS                                 *
C
C       IF(IUNI.EQ.1) GOTO 2630
C 2030  WRITE(6,2040)
C 2040  FORMAT(/'HOW MANY GROUPS DO YOU HAVE? ')
C       CALL DATA1(NGROUP)
C       IF(NGROUP.LT.2) THEN
C            PRINT *,'      NUMBER OF GROUPS MUST BE TWO OR MORE'
C            GOTO 2030
C       ENDIF
C
C       IF(NGROUP.EQ.2) GOTO 2180C
C
C
C      *  IF THE NUMBER OF GROUPS IS MORE THAN TWO, SPECIFY COMBINATIONS    *
C
C 2170  PRINT *
C       PRINT *,'     WHICH COMBINATION DO YOU WANT TO TEST? '
C 2180  PRINT *
C       WRITE(6,2190)
C 2190  FORMAT('FIRST GROUP INDICATOR  ')
C       CALL DATA1(IFIRST)
C       IGROUP(1) = IFIRST
C 2210  PRINT *
C       WRITE(6,2220)
C 2220  FORMAT('SECOND GROUP INDICATOR  ')
C       CALL DATA1(ISECON)
C 2240  IF(IFIRST.EQ.ISECON) THEN
C            PRINT *,' YOU CHOSE THE SAME GROUP.'
C            PRINT *,' PLEASE CHANGE THE SECOND GROUP.'
C            GOTO 2210
C       ENDIF
C       IGROUP(2) = ISECON
C
C      *         READ THE NAME OF THE GROUPS                                *
C
C 2250  PRINT *
C       WRITE(6,2255) IFIRST   
C 2255  FORMAT('WHAT IS THE NAME OF GROUP ',I4,' ? ')
C       READ(5,1380) GROUP(1)
C       PRINT *
C       WRITE(6,2258) ISECON   
C 2258  FORMAT('WHAT IS THE NAME OF GROUP ',I4,' ? ')
C       READ(5,1380) GROUP(2)
C       PRINT *
C
C      *      READ WHETHER TO PRINTOUT ONLY RESULTS OR TO GIVE FULL DETAILS *
C
C 2312  WRITE(6,2314)
C 2314  FORMAT('DO YOU WANT PRINTOUTS OF COMPUTATIONAL',
C     +       ' DETAILS (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.NE.'Y'.AND.CHECK.NE.'y') GOTO 2316
C       IFULL=1
C       GOTO 2400
C 2316  IF(CHECK.NE.'N'.AND.CHECK.NE.'n') GOTO 2312
C       IFULL=0
C
C      *     LKM IS SET TO ONE IN THE TWO-SAMPLE TEST ROUTINE, SO THAT      *
C      *     THE KAPLAN-MEIER PERCENTILES AND MEAN FOR EACH GROUP ARE       *
C      *     AUTOMATICALLY PROVIDED.                                        *
C
C
C 2400  LKM=1
C
C
C      *  CHECK WHETHER THE FULL K-M ESTIMATOR IS NEEDED.                   *
C      *  FROM THE NEXT LINE, THE INPUTS ARE COMMON FOR BOTH KAPLAN-MEIER   *
C      *  ESTIMATION AND TWO SAMPLE TESTS.                                  *
C
C 2630  PRINT *
C       WRITE(6,2640)
C       IKM=0
C 2640  FORMAT('DO YOU WANT TO PRINT OUT THE FULL K-M ',
C     +        'ESTIMATE (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') THEN
C          IKM=1
C
C       *    INFORMATION ABOUT THE DIFFERENTIAL KM ESTIMATE                 *
C
C 4644    PRINT *
C         PRINT *,'DO YOU WANT TO PRINT OUT '
C         WRITE(6,2644)
C 2644    FORMAT('           THE DIFFERENTIAL FORM (Y/N)?')
C         READ(5,50) CHECK
C         IF(CHECK .EQ. 'Y'.OR.CHECK.EQ.'y') THEN
C
C           KDIFF = 1
C 2646      WRITE(6,2647)
C 2647      FORMAT(' SPECIFY STARTING VALUE : ')
C           READ(5,2648) START
C 2648      FORMAT(F10.3)
C
C 4502      WRITE(6,4503) 
C 4503      FORMAT('HOW MANY BINS DO YOU WANT? : ')
C           READ(5,4554) LSTEP
C 4554      FORMAT(I4)
C           IF(LSTEP .LE. 0) GOTO 4502
C
C 4506      WRITE(6,4507)
C 4507      FORMAT('SPECIFY BIN SIZE :')
C           READ(5,2648) BINSIZ
C           IF(BINSIZ .LE. 0.0) GOTO 4506
C
C         ELSEIF(CHECK .EQ. 'N'.OR. CHECK .EQ. 'n') THEN
C           KDIFF = 0
C         ELSE
C           GOTO 4644
C         ENDIF
C       ELSEIF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') THEN
C         IKM=0
C       ELSE
C         GOTO 2630
C       ENDIF
C       
C 2650  PRINT *
C        WRITE(6,2655)
C 2655  FORMAT('DO YOU WANT TO PRINT THE ORIGINAL DATA (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IDATA=1
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') IDATA=0
C       IF((CHECK.EQ.'Y').OR.(CHECK.EQ.'N')) GOTO 3230
C       IF((CHECK.EQ.'y').OR.(CHECK.EQ.'n')) GOTO 3230
C       GOTO 2650
C
C      *   IF ALL INFORMATION IS TO BE READ FROM THE TERMINAL, GOTO 3230    *
C      *   FROM THE NEXT LINE, THE INPUTS ARE FROM "COMMAND" FILE           *
C
C
C      *   READ THE NAME OF "COMMAND" FILE                                  *
C
C 2680  PRINT *
C       WRITE(6,2690)
C 2690  FORMAT('WHAT IS THE NAME OF YOUR COMMAND FILE? ')
C       READ(5,1380) COMMAND
C       WRITE(6,2710) COMMAND
C 2710  FORMAT(5X,'YOUR COMMAND FILE IS CALLED ',A9)
       COMMAND="test2.dat"
C
       OPEN(UNIT=50,FILE=COMMAND,STATUS='OLD',FORM='FORMATTED')
C
C      *              READ THE DATA FILE NAME                               *
C
 2820  READ(50,1380) FILE
C
C      *              READ THE TITLE OF THE PROBLEM                         *
C
       READ(50,1430) TITLE 
C
C      * READ THE NUMBER OF VARIABLES.                                      *
C
       READ(50,2830) (CHAR(I,1),I=1,4)
 2830  FORMAT(12A1)
C
C      *      CHECK THE NUMBER OF VARIABLES                                 *
C
       CALL DATA2(CHAR,1,3,NVAR,LIND)
       IF(LIND.EQ.0) GOTO 2845 
 2835  PRINT *
       PRINT *,'     NUMBER OF VARIABLES IS NOT READABLE'
       CLOSE(UNIT=50)
       STOP
 2845  IF(NVAR.LE.0) GOTO 2835
C
C      *          READ WHICH VARIABLE NEEDS TO BE TESTED                    *
C
       READ(50,2910) (CHAR(I,1),I=1,4)
 2910  FORMAT(4A1)
       CALL DATA2(CHAR,1,1,ICOL,IN)
       IF(IN.EQ.0) GOTO 2935
 2915  PRINT *
       WRITE(6,2920) 
 2920  FORMAT(5X,'THE VARIABLE IS NOT READABLE')
       CLOSE(UNIT=50)
       STOP
 2935  IF((ICOL.LE.0).OR.(ICOL.GT.NVAR)) GOTO 2915
 2940  CONTINUE
C
C
C      *      READ THE NAME OF THE VARIABLE                                 *
C
 2962  READ(50,2963) COLM
 2963  FORMAT(10A9)
       IF(IUNI.EQ.1) GOTO 3180
C
C      *   FROM THE NEXT LINE, INPUTS ARE ONLY FOR TWO SAMPLE TESTS         *
C      *   IF IT IS A K-M ESTIMATOR PROBLEM, GO TO 3180                     *
C
C      *   READ THE NUMBER OF GROUPS                                        *
C
       READ(50,2970) (CHAR(I,1),I=1,4)
 2970  FORMAT(4A1)
       CALL DATA2(CHAR,1,1,NGROUP,LIND)
       IF(LIND.EQ.0) GOTO 3000
 2980  PRINT *,'     IT IS NOT CLEAR HOW MANY GROUPS YOU HAVE'
       CLOSE(UNIT=50)
       STOP
 3000  IF(NGROUP.GT.1) GOTO 3005
       GOTO 2980
C
C      *     READ THE INDICATOR OF THE GROUPS                               *
C
 3005  READ(50,3010) ((CHAR(I,J),I=1,4),J=1,NGROUP)
 3010  FORMAT(60A1)
 3020  DO 3050 J=1,NGROUP
       CALL DATA2(CHAR,J,NGROUP,JGROUP(J),LIND)
       IF(LIND.EQ.0) GOTO 3050
 3025  PRINT *
       WRITE(6,3030) J
 3030  FORMAT(5X,'THE INDICATOR OF ',I4,'TH GROUP IS NOT CLEAR')
       CLOSE(UNIT=50)
       STOP
 3050  CONTINUE
C
C      *  READ NUMBER OF FIRST GROUP,SECOND GROUP,                          *
C      *       WHETHER PRINT OUT ALL OR RESULTS ONLY                        *
C      *       WHETHER K-M ESTIMATOR IS NEEDED                              *
C      *       WHETHER PRINT OUT ALL OR RESULTS ONLY FOR K-M                *
C
       READ(50,3085) (CF(I1,1),I1=1,4),(CS(I2,1),I2=1,4),
     +              (CIS4(I3,1),I3=1,4),(CIKM(I4,1),I4=1,4)
 3085  FORMAT(16A1)
C
       CALL DATA2(CF,1,1,IFIRST,LIND)
       IF(LIND.EQ.0) GOTO 3087
 3086  PRINT *
       PRINT *,'   THE VALUE FOR "IFIRST" IS NOT ACCEPTABLE'
       CLOSE(UNIT=50)
       STOP
 3087  IF((IFIRST.LT.0).OR.(IFIRST.GT.NGROUP)) GOTO 3086
C
       CALL DATA2(CS,1,1,ISECON,LIND)
       IF(LIND.EQ.0) GOTO 3089
 3088  PRINT *
       PRINT *,'   THE VALUE FOR "ISECON" IS NOT ACCEPTABLE'
       CLOSE(UNIT=50)
       STOP
 3089  IF((ISECON.LT.0).OR.(ISECON.GT.NGROUP)) GOTO 3087
       IF(ISECON.EQ.IFIRST) GOTO 3087
       IGROUP(1) = IFIRST
       IGROUP(2) = ISECON
C
       CALL DATA2(CIS4,1,1,IFULL,LIND)
       IF(LIND.EQ.0) GOTO 3091
 3090  PRINT *
       PRINT *,'    THE VALUE FOR "IFULL" IS NOT ACCEPTABLE'
       CLOSE(UNIT=50)
       STOP
 3091  IF((IFULL.EQ.0).OR.(IFULL.EQ.1)) GOTO 3092
       GOTO 3090
C
 3092  LKM = 1
C
 3095  CALL DATA2(CIKM,1,1,IKM,LIND)
       IF(LIND.EQ.0) GOTO 3097
 3096  PRINT *
       PRINT *,'    THE VALUE FOR "IKM" IS NOT ACCEPTABLE'
       CLOSE(UNIT=50)
       STOP
 3097  IF((IKM.EQ.0).OR.(IKM.EQ.1)) GOTO 5190
       GOTO 3096
C
C      *    INFORMATION ABOUT THE DIFFERENTIAL KM ESTIMATOR                  *
C
 5190  READ(50,2970) (CHAR(I,1),I=1,4)
       CALL DATA2(CHAR,1,1,KDIFF,LIND)
       IF(LIND.EQ.0) GOTO 5201
 5200  PRINT *,'    IT IS NOT CLEAR WHETHER YOU WANT TO PRINT'
       PRINT *,'    DIFFERENTIAL KM ESTIMATOR'
       CLOSE(UNIT=50)
       STOP
 5201  IF(KDIFF.EQ.1) GOTO 5202
       IF(KDIFF.EQ.0) GOTO 3102
       GOTO 5200
C
 5202  READ(50,4203) START
 5203  FORMAT(F10.3)
       READ(50,4204) LSTEP
 5204  FORMAT(I4)
       READ(50,4203) BINSIZ
C
 3102  READ(50,2970) (CHAR(I,1),I=1,4)
       CALL DATA2(CHAR,1,1,IDATA,LIND)
       IF(LIND.EQ.0) GOTO 3110
 3100  PRINT *
       PRINT *,'    THE VALUE FOR "IDATA" IS NOT ACCEPTABLE'
       CLOSE(UNIT=50)
       STOP
 3110  IF((IDATA.EQ.0).OR.(IDATA.EQ.1)) GOTO 3168
       GOTO 3100
C
C      *            READ THE NAME OF THE GROUPS                             *
C
 3168  READ(50,1380) GROUP(1)
       READ(50,1380) GROUP(2)
C
C      *   READ NAME OF THE OUTPUT FILE. IF THE NAME IS BLANK, THE RESULTS  *
C      *   WILL BE SHOWN ON THE TERMINAL.                                   *
C
       READ(50,1380) OUTPUT
       IF(OUTPUT.NE.'         ') GOTO 3300
       GOTO 3230
C
C      *      FROM THE NEXT LINE, INPUTS ARE ONLY FOR THE K-M ESTIMATOR     *
C
 3180  READ(50,3200) (CHAR(I,1),I=1,4)
 3200  FORMAT(4A1)
       CALL DATA2(CHAR,1,1,IKM,LIND)
       IF(LIND.EQ.0) GOTO 3205
 3203  PRINT *,'     IT IS NOT CLEAR WHETHER YOU WANT TO PRINT OUT ALL'
       PRINT *,'     KM ESTIMATORS'
       CLOSE(UNIT=50)
       STOP
 3205  IF((IKM.EQ.0).OR.(IKM.EQ.1)) GOTO 4190
       GOTO 3203
C
C      *    INFORMATION ABOUT THE DIFFERENTIAL KM ESTIMATOR                  *
C
 4190  READ(50,2970) (CHAR(I,1),I=1,4)
       CALL DATA2(CHAR,1,1,KDIFF,LIND)
       IF(LIND.EQ.0) GOTO 4201
 4200  PRINT *,'    IT IS NOT CLEAR WHETHER YOU WANT TO PRINT'
       PRINT *,'    DIFFERENTIAL KM ESTIMATOR'
       CLOSE(UNIT=50)
       STOP
 4201  IF((KDIFF.EQ.0).OR.(KDIFF.EQ.1)) GOTO 4202
       GOTO 4200
C
 4202  READ(50,4203) START
 4203  FORMAT(F10.3)
       READ(50,4204) LSTEP
 4204  FORMAT(I4)
       READ(50,4203) BINSIZ
C
C      *    INFORMATION ABOUT PRINTOUT                                       *
C
 3210  READ(50,2970) (CHAR(I,1),I=1,4)
       CALL DATA2(CHAR,1,1,IDATA,LIND)
       IF(LIND.EQ.0) GOTO 3220
 3215  PRINT *
       PRINT *,'    THE VALUE FOR "IDATA" IS NOT ACCEPTABLE'
       CLOSE(UNIT=50)
       STOP
 3220  IF((IDATA.EQ.0).OR.(IDATA.EQ.1)) GOTO 3225
       GOTO 3215
C
C
 3225  READ(50,1380) OUTPUT
       CLOSE(UNIT=50)
       IF(OUTPUT.NE.'         ') GOTO 3300
C
C
C      *         LEAVE THE "COMMAND" FILE                                   *
C      *         CHECK OUTPUT FILE                                          *
C
 3230  OUTPUT='         '
C 3240  PRINT *
C       WRITE(6,3250)
C 3250  FORMAT('DO YOU WANT TO SAVE THE RESULTS IN A FILE (Y/N)? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') GOTO 3260
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') GOTO 3300
C       GOTO 3240
C 3260  PRINT *
C       WRITE(6,3270)
C 3270  FORMAT('WHAT IS THE NAME OF THE FILE?  ')
C       READ(5,1380) OUTPUT
C
C      *          READ IN DATA THOUGH THE SUBROUTINE "DATAIN"               *
C
 3300  CALL DATAIN(IUNI,FILE,NVAR,ISTA,IND,X,NTOT,NDAT,MVAR)
C
C       IF(IUNI.EQ.2) GOTO 3360
C
C      *              COMPUTE THE K-M ESTIMATOR                             *
C
C
C       IF(OUTPUT.NE.'         ') 
C     +                OPEN(UNIT=60,FILE=OUTPUT,STATUS='NEW'
C  THE FOLLOWING LINE SHOULD BE COMMENTED OUT ON ALL MACHINES OTHER THAN
C  VAX/VMS MACHINES.
C     +                ,CARRIAGECONTROL='LIST'
C     +                )
C
C       CALL KMESTM(IND,X,NTOT,ICOL,IKM,TITLE,COLM,OUTPUT,IBIN,0,
C     +             KDIFF,START,BINSIZ,LSTEP,FILE,
C     +             WRK1,WRK2,WRK3,WRK4,DWRK1,IWRK1,IWRK2,
C     +             WRK7,WRK8,BWRK1,BWRK2,BWRK3,IWRK3,SWRK1,
C     +             WRK9,MVAR)
C
C       IF(IDATA.EQ.0) GOTO 3335
C       IF(OUTPUT.NE.'         ') WRITE(60,3331) FILE
C       IF(OUTPUT.EQ.'         ') WRITE(6,3331)  FILE
C       IF(OUTPUT.NE.'         ') WRITE(60,3332) 
C       IF(OUTPUT.EQ.'         ') PRINT 3332
 3331  FORMAT(7X,' INPUT DATA FILE: ',A9)
C 3332  FORMAT(5X,'   CENSORSHIP     X ')
C       DO 3333 I=1,NTOT
C       IF(OUTPUT.NE.'         ') WRITE(60,3334) IND(ICOL,I),X(ICOL,I)
C       IF(OUTPUT.EQ.'         ') PRINT 3334,IND(ICOL,I),X(ICOL,I)
C 3333  CONTINUE
C 3334  FORMAT(12X,I4,F10.3)
C
C
C       IF(OUTPUT.NE.'         ') CLOSE(UNIT = 60)
C
C 3335  PRINT *
C       PRINT *,'    K-M ESTIMATOR COMPUTATIONS ARE FINISHED'
C
C      *      CHECK WHETHER THE USER WANTS TO USE OTHER METHODS            *
C
C 3340  WRITE(6,3350)
C 3350  FORMAT('DO YOU WANT ANY OTHER ANALYSIS (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IBACK=1
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') RETURN
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') STOP
C       GOTO 3340
C
C
C      *                COMPUTE TWO SAMPLE TESTS                            *
C
C 3360  IF(OUTPUT.EQ.'         ') GOTO 3370
C       OPEN(UNIT=60,FILE=OUTPUT,STATUS='NEW'
C  THE FOLLOWING LINE SHOULD BE COMMENTED OUT ON ALL MACHINES OTHER THAN
C  VAX/VMS MACHINES.
C     +    ,CARRIAGECONTROL='LIST'
C     +    )
C       WRITE(60,3380)
C       WRITE(60,3390)
C       WRITE(60,3400) TITLE  
C       WRITE(60,3390)
C       GOTO 3410
C 3370  WRITE(6,3380)
C       PRINT *
C       WRITE(6,3400) TITLE  
C       PRINT *
C 3380  FORMAT(8X,'   ******* TWO SAMPLE TEST ******')
C 3390  FORMAT('    ')
C 3400  FORMAT(8X,'TITLE : ',A80)    
C
C 3410  IF(OUTPUT.EQ.'         ') GOTO 3420
C       WRITE(60,3430) FILE
C       WRITE(60,3435) COLM
C       GOTO 3440
C 3420  WRITE(6,3430) FILE
C       WRITE(6,3435) COLM
C 3430  FORMAT(8X,'DATA SET : ',A9)
C 3435  FORMAT(8X,'VARIABLE : ',A9)
C
C 3440  IF(OUTPUT.EQ.'         ') GOTO 3450
C       WRITE(60,3460) GROUP(1),GROUP(2)
C       GOTO 3470
C 3450  WRITE(6,3460) GROUP(1),GROUP(2)
C 3460  FORMAT(8X,'TESTED GROUPS ARE ',A9,' AND ',A9)
C
C
C 3470   DO 3480 M=1,NOTEST
C

 3470  CALL TWOST_L(X,IND,ISTA,ICOL,IGROUP(1),IGROUP(2),NTOT,
     +            MVAR,NDAT,RETURNED_PROB)
C     FILE)
C     +            WRK1,WRK2,WRK3,WRK4,WRK5,WRK6,WRK7,WRK8,
C     +            WRK9,WRK10,WRK11,WRK12,IWRK1,IWRK2,
C     +            WRK13,WRK14,WRK15,WRK16,WRK17,WRK18,WRK19,
C     +            WRK20,WRK21,WRK22,WRK23,WRK24,WRK25,WRK26,IWRK3)
C
C 3480  CONTINUE


       PRINT *,RETURNED_PROB

C
C      *         IF K-M ESTIMATOR IS NOT REQUESTED, GOTO 3510               *
C
C       IF(LKM.EQ.0) GOTO 3510
C       N1=0
C       N2=0
C       DO 3490 I=1,NTOT
C          IF(ISTA(I).NE.IGROUP(1)) GOTO 3490
C          N1=N1+1
C          JIND(ICOL,N1)=IND(ICOL,I)
C          Z(ICOL,N1)=X(ICOL,I)
C 3490  CONTINUE
C
C       CALL KMESTM(JIND,Z,N1,ICOL,IKM,TITLE,GROUP(1),OUTPUT,
C     +             IBIN,1,KDIFF,START,BINSIZ,LSTEP,FILE,
C     +             WRK1,WRK2,WRK3,WRK4,DWRK1,IWRK1,IWRK2,
C     +             WRK7,WRK8,BWRK1,BWRK2,BWRK3,IWRK3,SWRK1,
C     +             WRK9,MVAR)
CC
C       DO 3500 I=1,NTOT
C          IF(ISTA(I).NE.IGROUP(2)) GOTO 3500
C          N2=N2+1
C          JIND(ICOL,N2)=IND(ICOL,I)
C          Z(ICOL,N2)=X(ICOL,I)
C 3500  CONTINUE
CC
C       CALL KMESTM(JIND,Z,N2,ICOL,IKM,TITLE,GROUP(2),OUTPUT,
C     +             IBIN,1,KDIFF,START,BINSIZ,LSTEP,FILE,
C     +             WRK1,WRK2,WRK3,WRK4,DWRK1,IWRK1,IWRK2,
C     +             WRK7,WRK8,BWRK1,BWRK2,BWRK3,IWRK3,SWRK1,
C     +             WRK9,MVAR)

 3510  IF(IDATA.EQ.0) GOTO 3525

       IF(OUTPUT.NE.'         ') WRITE(60,3331) FILE
       IF(OUTPUT.EQ.'         ') PRINT 3331, FILE
       IF(OUTPUT.NE.'         ') WRITE(60,3520) 
       IF(OUTPUT.EQ.'         ') PRINT 3520
 3520  FORMAT(5X,'  CENSORSHIP     GROUP      X ')
       DO 3521 I=1,NTOT
       IF(OUTPUT.NE.'         ') WRITE(60,3522) IND(ICOL,I),ISTA(I),
     +                                                        X(ICOL,I)
       IF(OUTPUT.EQ.'         ') PRINT 3522,IND(ICOL,I),ISTA(I),
     +                                                        X(ICOL,I)
 3521  CONTINUE
 3522  FORMAT(12X,I4,6X,I4,F10.3)

       IF(OUTPUT.NE.'         ') CLOSE(UNIT = 60)

 3525  PRINT *
C       PRINT *,'     THE TWO SAMPLE TESTS ARE FINISHED'
C 3530  PRINT *
C       WRITE(6,3540)
C 3540  FORMAT('DO YOU WANT ANY OTHER ANALYSIS (Y/N) ? ')
C       READ(5,50) CHECK
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') IBACK=1
C       IF(CHECK.EQ.'Y'.OR.CHECK.EQ.'y') RETURN
C       IF(CHECK.EQ.'N'.OR.CHECK.EQ.'n') STOP
C       GOTO 3530
       STOP
       END
