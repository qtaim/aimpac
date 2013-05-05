
      PROGRAM GRDVEC
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*80 LINE,TITLE
      CHARACTER*40 WVEC,WFN,WGVP
      CHARACTER*7 WORD
      CHARACTER*4 FVEC /'.vec'/, FWFN /'.wfn'/, FGVP /'.gvp'/
C
      DIMENSION A(3,3),C(3,60),CR(2),CT(3),CX(3),IAR(60),IR(60)
      DIMENSION NUP(60),NDN(60),ORG(3,60),SV(3),XS(3,60),XY(2)
      DIMENSION H(3,3),EUL(3),F(3)
C
      DATA IWFN /10/, IVEC /15/, IOUT /13/
C
      COMMON CO(30000),IC(30000),MODE,NCENT,NMO,NPRIMS
      COMMON /OFFSET/ ITYPE,ICENT,KATOM,IEORB,IE,IMO,ICHARG,IXC,IYC,IZC,
     +  IXX, IYY, IZZ,IXS,IYS,IZS,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,ID2,IGXX,
     +  IGXY,IGXZ,IGYY,IGYZ,IGZZ
      COMMON /TRANS/ A,CR,CX,CT,SCAL,IPROJ
      COMMON /UNITS/ ISRF,IVEC,IOUT,IWFN,IDBG
C
      CALL MAKNAME (1,WVEC,ILEN,FVEC)
      IF (ILEN .EQ. 0) STOP 'usage: grdvec vecfile wfnfile'
      CALL MAKNAME (1,WGVP,ILEN,FGVP)
      IF (ILEN .EQ. 0) STOP 'usage: grdvec vecfile wfnfile'
      CALL MAKNAME (2,WFN,ILEN,FWFN)
      IF (ILEN .EQ. 0) STOP 'usage: grdvec vecfile wfnfile'
C
      OPEN (IVEC,FILE=WVEC)
      OPEN (IOUT,FILE=WGVP)
      OPEN (IWFN,FILE=WFN)
C
      CALL RDWFN
C
C    TITLE
C
      READ (IVEC,1000) TITLE
C
C    PLOT:
C
      READ (IVEC,1010) LINE
      LPST = 8
      DO 90 I = 1,2
        IF (NUMBER(LINE,LPST,NUM,XY(I)) .GT. 0) GOTO 1990
90    CONTINUE
C 
C    DETERMINE CENTER OF LOCAL FRAME (I.E. USER DEFINED CENTER OF MOLECULE)
C
      CR(1) = XY(1)/2.0D0
      CR(2) = XY(1)/2.0D0
      SCAL = 10.0D0/XY(1)
C
      CALL PLOTS (53,0,13)
C
C    CENTER OF PLOT
C
      READ (IVEC,1010) LINE
      LPST = 8
      DO 100 I = 1,3
        IF (NUMBER(LINE,LPST,NUM,CX(I)) .GT. 0) GOTO 2000
100   CONTINUE
C
C    WALKING PARAMETERS
C
      READ (IVEC,*) WORD,R1,R2,R3,ENDPT,IPROJ,IINC
C
C    DEFINE PLANE
C
      READ (IVEC,1010) LINE
      LPST = 8
      IF (NUMBER(LINE,LPST,IEG,DNUM) .GT. 0) GOTO 2070
      IF (IEG .EQ. 0) THEN
        J = 1
110     IF (NUMBER(LINE,LPST,IR(J),DNUM) .GT. 0) GOTO 120
        J = J + 1
        GOTO 110
C
C    DECIDE IF DUMMY ATOMS WERE USED AND READ IN THEIR COORDINATES
C    TO XS(3,N)
C
120     DO 130 I = 1,J
          IF (IR(I) .LT. 0) THEN
            LPST = 8
            READ (IVEC,1010) LINE
            DO 140 K = 1,3
              IF (NUMBER(LINE,LPST,NUM,XS(K,ABS(IR(I)))) .GT. 0) 
     +        GOTO 2010
140         CONTINUE
          END IF
130     CONTINUE
C
C    GENERATE ROTATION MATRIX 
C
      DO 170 I = 1,NCENT
        C(1,I) = CO(IXC+I)
        C(2,I) = CO(IYC+I)
        C(3,I) = CO(IZC+I)
170   CONTINUE
C
      CALL GROCKLE(NCENT,C,IR,XS,A)
C
      ELSE IF (IEG .EQ. 1) THEN
        DO 125 I = 1,3
          IF (NUMBER(LINE,LPST,NUM,EUL(I)) .GT. 0) GOTO 2080
125     CONTINUE
C
        CALL EULER (EUL,A)
C
      END IF
C
C   NUMBER OF ORIGINS
C
      READ (IVEC,1010) LINE
      LPST = 8
      IF (NUMBER(LINE,LPST,NORG,DNUM) .GT. 0) GOTO 2020
C
C    GET ORIGINS
C
      DO 150 I = 1,NORG
        READ (IVEC,*) (ORG(J,I),J=1,3),IAR(I),NUP(I),NDN(I)
150     CONTINUE
C
C    TRANSFORM CENTER OF PLOT TO CENTER OF PLOTTING FRAME
C
      CT(1) = A(1,1)*CX(1) + A(2,1)*CX(2) + A(3,1)*CX(3)
      CT(2) = A(1,2)*CX(1) + A(2,2)*CX(2) + A(3,2)*CX(3)
      CT(3) = A(1,3)*CX(1) + A(2,3)*CX(2) + A(3,3)*CX(3)
C
C    PLOT ATOM POSITIONS
C
      DO 180 I = 1,NCENT
        XX = A(1,1)*CO(IXC+I)+A(2,1)*CO(IYC+I)+A(3,1)*CO(IZC+I)-CT(1)
        YY = A(1,2)*CO(IXC+I)+A(2,2)*CO(IYC+I)+A(3,2)*CO(IZC+I)-CT(2)
        ZZ = A(1,3)*CO(IXC+I)+A(2,3)*CO(IYC+I)+A(3,3)*CO(IZC+I)-CT(3)
        CALL NUCLEI(XX,YY,ZZ)
180   CONTINUE
C
C    LOOP OVER EACH ORIGIN
C
      DO 190 I = 1,NORG
C
C    FIND ORIGIN'S COORDINATES IN PLOT'S FRAME OF REFERENCE
C
          X = A(1,1)*ORG(1,I)+A(2,1)*ORG(2,I)+A(3,1)*ORG(3,I)-CT(1)
          Y = A(1,2)*ORG(1,I)+A(2,2)*ORG(2,I)+A(3,2)*ORG(3,I)-CT(2)
          Z = A(1,3)*ORG(1,I)+A(2,3)*ORG(2,I)+A(3,3)*ORG(3,I)-CT(3)
C
C    THIS ORIGIN IS A (2,-2) ATTRACTOR IN THE PLANE
C
        IF (IAR(I) .EQ. 0) THEN
C
          PATHA = 360.0D0/DFLOAT(NDN(I))
          DO 200 J = 1,NDN(I)
C
C    FIND STARTING POINT FOR EACH GRADIENT PATH
C
            ANG = DFLOAT((J-1))*PATHA
            X0 = X + R1*DCOSD(ANG)
            Y0 = Y + R1*DSIND(ANG)
            Z0 = 0.0D0
            CALL TRUDGE(X0,Y0,Z0,0,IINC,ENDPT)
200       CONTINUE
C
C    THIS ORIGIN IS A (2,0) REPELLOR IN THIS PLANE
C
        ELSE IF (IAR(I) .EQ. 1) THEN
C
C    CALL TO GET HESSIAN AND GRADIENT AT ORIGIN
C
          X0 = ORG(1,I)
          Y0 = ORG(2,I)
          Z0 = ORG(3,I)
          CALL HESS(X0,Y0,Z0,F,H)
C
C    CALL TO GET DESIRED SHIFT VECTOR FROM HESSIAN
C    TAKING INTO ACCOUNT WHETHER OR NOT PROJECTIONS
C    OUT OF THE PLANE ARE TO BE ALLOWED
C
          IF (NUP(I) .GT. 0) THEN
            CALL SHFTVC(A,F,H,1,SV)
            X0 = X+A(1,1)*SV(1)*R2+A(2,1)*SV(2)*R2+A(3,1)*SV(3)*R2
            Y0 = Y+A(1,2)*SV(1)*R2+A(2,2)*SV(2)*R2+A(3,2)*SV(3)*R2
            IF (IPROJ .EQ. 0) THEN 
              Z0 = 0.0D0
            ELSE
              Z0 = Z+A(1,3)*SV(1)*R2+A(2,3)*SV(2)*R2+A(3,3)*SV(3)*R2 
            END IF
            CALL TRUDGE(X0,Y0,Z0,1,IINC,ENDPT)
            X0 = X-A(1,1)*SV(1)*R2-A(2,1)*SV(2)*R2-A(3,1)*SV(3)*R2
            Y0 = Y-A(1,2)*SV(1)*R2-A(2,2)*SV(2)*R2-A(3,2)*SV(3)*R2
            IF (IPROJ .EQ. 0) THEN 
              Z0 = 0.0D0
            ELSE
              Z0 = Z-A(1,3)*SV(1)*R2-A(2,3)*SV(2)*R2-A(3,3)*SV(3)*R2 
            END IF
            CALL TRUDGE(X0,Y0,Z0,1,IINC,ENDPT)
          END IF 
          IF (NDN(I) .GT. 0) THEN
            CALL SHFTVC(A,F,H,0,SV)
            X0 = X+A(1,1)*SV(1)*R3+A(2,1)*SV(2)*R3+A(3,1)*SV(3)*R3
            Y0 = Y+A(1,2)*SV(1)*R3+A(2,2)*SV(2)*R3+A(3,2)*SV(3)*R3
            IF (IPROJ .EQ. 0) THEN 
              Z0 = 0.0D0
            ELSE
              Z0 = Z+A(1,3)*SV(1)*R3+A(2,3)*SV(2)*R3+A(3,3)*SV(3)*R3
            END IF
            CALL TRUDGE(X0,Y0,Z0,0,IINC,ENDPT)
            X0 = X-A(1,1)*SV(1)*R3-A(2,1)*SV(2)*R3-A(3,1)*SV(3)*R3
            Y0 = Y-A(1,2)*SV(1)*R3-A(2,2)*SV(2)*R3-A(3,2)*SV(3)*R3
            IF (IPROJ .EQ. 0) THEN 
              Z0 = 0.0D0
            ELSE
              Z0 = Z-A(1,3)*SV(1)*R3-A(2,3)*SV(2)*R3-A(3,3)*SV(3)*R3
            END IF
            CALL TRUDGE(X0,Y0,Z0,0,IINC,ENDPT)
          END IF
        END IF
C
190   CONTINUE
C
      CALL PLOT (0.,0.,999)
      GOTO 4999
C
C    FORMATS
C
1000  FORMAT(7X,A70)
1010  FORMAT(A80)
1990  WRITE (IOUT,1995)
1995  FORMAT(' ERROR IN PLOT CARD ')
      GOTO 4999
2000  WRITE (IOUT,2005)
2005  FORMAT(' ERROR IN CENTER DEFINITION CARD ')
      GOTO 4999
2010  WRITE (IOUT,2015)
2015  FORMAT(' ERROR IN DUMMY ATOM COORDINATES ')
      GOTO 4999
2020  WRITE (IOUT,2025)
2025  FORMAT(' ERROR IN CARD GIVING NUMBER OF ORIGINS ')
      GOTO 4999
2040  WRITE (IOUT,2045)
2045  FORMAT(' ERROR IN ORIGIN TYPE DECLARATION ')
      GOTO 4999
2050  WRITE (IOUT,2055)
2055  FORMAT(' ERROR IN NUMBER OF ASCENDING GRADIENT VECTORS ')
      GOTO 4999
2060  WRITE (IOUT,2065)
2065  FORMAT(' ERROR IN NUMBER OF DESENDING GRADIENT VECTORS ')
      GOTO 4999
2070  WRITE (IOUT,2075)
2075  FORMAT(' ERROR IN ORIENTATION METHOD ')
      GOTO 4999
2080  WRITE (IOUT,2085)
2085  FORMAT( 'ERROR IN EULER ANGLES CARD ')
      GOTO 4999
3000  WRITE (6,3005) 
3005  FORMAT(' ERROR IN ORIGIN COORDINATES ')
      GOTO 4999
3010  WRITE (6,3015)
3015  FORMAT(' ERROR IN ATTRACTOR TYPE ')
      GOTO 4999
3020  WRITE (6,3025)
3025  FORMAT(' ERROR IN NUMBER OF ASCENDING PATHS ')
      GOTO 4999
3030  WRITE (6,3035)
3035  FORMAT(' ERROR IN NUMBER OF ASCENDING PATHS ')
      GOTO 4999
3040  WRITE (6,3045)
3045  FORMAT(' ERROR IN PARAMETER CARD ')
4999  END
      SUBROUTINE EULER (E,A)
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
C   CALCULATE THE ROTATION MATRICES USING
C   EULER ANGLES
C
      DIMENSION E(3), A(9)
C
      DATA RADIAN /0.01745329252D0/
C
      E1 = E(1)*RADIAN
      E2 = E(2)*RADIAN
      E3 = E(3)*RADIAN
C
      SA = DSIN(E1)
      SB = DSIN(E2)
      SC = DSIN(E3)
      CA = DCOS(E1)
      CB = DCOS(E2)
      CC = DCOS(E3)
C
C     A(1) =  CC*CB*CA - SC*SA
C     A(2) = -SC*CB*CA - CC*SA
C     A(3) =  SB*CA
C     A(4) =  CC*CB*SA + SC*CA
C     A(5) = -SC*CB*SA + CC*CA
C     A(6) =  SB*SA
C     A(7) = -CC*SB 
C     A(8) =  SC*SB
C     A(9) =  CB
C
      A(1) =  CC*CA - CB*SA*SC
      A(2) =  CC*SA + CB*CA*SC
      A(3) =  SC*SB
      A(4) = -SC*CA - CB*SA*CC
      A(5) = -SC*SA + CB*CA*CC
      A(6) =  CC*SB
      A(7) =  SB*SA
      A(8) = -SB*CA
      A(9) =  CB
C
      RETURN
      END
      SUBROUTINE GAUS2
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C                                                                       
C    FOR GAUSSIAN WAVEFUNCTIONS ONLY.  CALCULATES AT A GIVEN POINT      
C    THE VALUE OF EACH MOLECULAR ORBITAL AND THE MO GRADIENT VECTOR     
C    COMPONENTS.                                                        
C
      COMMON CO(30000),IC(30000),MODE,NCENT,NMO,NPRIMS
C                                                                       
       COMMON /OFFSET/ ITYPE,ICENT,KATOM,IEORB,IE,IMO,ICHARG,IXC,
     + IYC,IZC,IXX,IYY,IZZ,IXS,IYS,IZS,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,
     + ID2,IGXX,IGXY,IGXZ,IGYY,IGYZ,IGZZ
C                                                                       
      DATA ZERO /0.0D0/, ONE /1.0D0/, TWO /2.0D0/, THREE /3.0D0/ 
C
      DO 100 J = 1,NMO
        CO(IPSI+J) = ZERO
        CO(IGX+J) = ZERO
        CO(IGY+J) = ZERO
        CO(IGZ+J) = ZERO
100   CONTINUE
C
      DO 200 I = 1,NPRIMS
        EE = CO(IE+I)*TWO
        X = CO(IXX+IC(ICENT+I))
        Y = CO(IYY+IC(ICENT+I))
        Z = CO(IZZ+IC(ICENT+I))
        IF (CO(IE+I)*CO(IR2+IC(ICENT+I)) .GT. 170.0D0) THEN
          EXPON = ZERO
        ELSE
          EXPON = DEXP(-CO(IE+I)*CO(IR2+IC(ICENT+I)))
        ENDIF
        EXPEE = EXPON*EE
        IF (EXPEE.LT.1.D-25) EXPEE = ZERO
        IF (EXPON.LT.1.D-25) EXPON = ZERO
C                                                                       
C    S                                                                  
C                                                                       
        IF (IC(ITYPE+I) .EQ. 1) THEN 
          BF =  EXPON
          GX = -EXPEE*X
          GY = -EXPEE*Y
          GZ = -EXPEE*Z
C                                                                       
C    X                                                                  
C
        ELSE IF (IC(ITYPE+I) .EQ. 2) THEN
          BF =  EXPON*X
          GX =  EXPON*(ONE-EE*X*X)
          GY = -EXPEE*X*Y
          GZ = -EXPEE*X*Z
C                                                                       
C    Y                                                                  
C                                                                       
        ELSE IF (IC(ITYPE+I) .EQ. 3) THEN
          BF =  EXPON*Y
          GX = -EXPEE*Y*X
          GY =  EXPON*(ONE-EE*Y*Y)
          GZ = -EXPEE*Y*Z
C                                                                       
C    Z                                                                  
C                                                                       
        ELSE IF (IC(ITYPE+I) .EQ. 4) THEN
          BF =  EXPON*Z
          GX = -EXPEE*Z*X
          GY = -EXPEE*Z*Y
          GZ =  EXPON*(ONE-EE*Z*Z)
C                                                                       
C    XX                                                                
C                                                                       
        ELSE IF (IC(ITYPE+I) .EQ. 5) THEN
          BF =  EXPON*X*X
          GX = -EXPEE*X*X*X + TWO*EXPON*X
          GY = -EXPEE*X*X*Y
          GZ = -EXPEE*X*X*Z
C                                                                      
C    YY                                                               
C                                                                    
        ELSE IF (IC(ITYPE+I) .EQ. 6) THEN
          BF =  EXPON*Y*Y
          GX = -EXPEE*Y*Y*X
          GY = -EXPEE*Y*Y*Y + TWO*EXPON*Y
          GZ = -EXPEE*Y*Y*Z
C                                                                  
C    ZZ                                                             
C                                                                 
        ELSE IF (IC(ITYPE+I) .EQ. 7) THEN
          BF =  EXPON*Z*Z
          GX = -EXPEE*Z*Z*X
          GY = -EXPEE*Z*Z*Y
          GZ = -EXPEE*Z*Z*Z + TWO*EXPON*Z
C                                                                
C    XY                                                         
C                                                              
        ELSE IF (IC(ITYPE+I) .EQ. 8) THEN
          BF =  EXPON*X*Y
          GX = -EXPEE*X*Y*X + EXPON*Y
          GY = -EXPEE*X*Y*Y + EXPON*X
          GZ = -EXPEE*X*Y*Z
C                                                             
C    XZ                                                      
C                                                           
        ELSE IF (IC(ITYPE+I) .EQ. 9) THEN
          BF =  EXPON*X*Z
          GX = -EXPEE*X*Z*X + EXPON*Z
          GY = -EXPEE*X*Z*Y
          GZ = -EXPEE*X*Z*Z + EXPON*X
C                                                          
C    YZ                                                   
C                                                        
        ELSE IF (IC(ITYPE+I) .EQ. 10) THEN
          BF =  EXPON*Y*Z
          GX = -EXPEE*Y*Z*X
          GY = -EXPEE*Y*Z*Y + EXPON*Z
          GZ = -EXPEE*Y*Z*Z + EXPON*Y
C
C     XXX
C
        ELSE IF (IC(ITYPE+I) .EQ. 11) THEN
          BF =  EXPON*X*X*X
          GX = -EXPON*X*X*EE*X*X + THREE*EXPON*X*X 
          GY = -EXPON*X*X*EE*X*Y
          GZ = -EXPON*X*X*EE*X*Z
C
C     YYY 
C
        ELSE IF (IC(ITYPE+I) .EQ. 12) THEN
          BF =  EXPON*Y*Y*Y
          GX = -EXPON*Y*Y*EE*Y*X  
          GY = -EXPON*Y*Y*EE*Y*Y + THREE*EXPON*Y*Y
          GZ = -EXPON*Y*Y*EE*Y*Z
C
C     ZZZ
C
        ELSE IF (IC(ITYPE+I) .EQ. 13) THEN
          BF =  EXPON*Z*Z*Z
          GX = -EXPON*Z*Z*EE*Z*X 
          GY = -EXPON*Z*Z*EE*Z*Y
          GZ = -EXPON*Z*Z*EE*Z*Z + THREE*EXPON*Z*Z
C
C     XXY
C
        ELSE IF (IC(ITYPE+I) .EQ. 14) THEN
          BF =  EXPON*X*X*Y
          GX = -EXPON*X*Y*EE*X*X + TWO*EXPON*X*Y
          GY = -EXPON*X*X*EE*Y*Y + ONE*EXPON*X*X
          GZ = -EXPON*X*X*EE*Y*Z
C
C     XXZ 
C
        ELSE IF (IC(ITYPE+I) .EQ. 15) THEN
          BF =  EXPON*X*X*Z
          GX = -EXPON*X*Z*EE*X*X + TWO*EXPON*X*Z
          GY = -EXPON*X*X*EE*Y*Z 
          GZ = -EXPON*X*X*EE*Z*Z + ONE*EXPON*X*X
C
C     YYZ 
C
        ELSE IF (IC(ITYPE+I) .EQ. 16) THEN
          BF =  EXPON*Y*Y*Z
          GX = -EXPON*Y*Y*EE*X*Z 
          GY = -EXPON*Y*Z*EE*Y*Y + TWO*EXPON*Y*Z
          GZ = -EXPON*Y*Y*EE*Z*Z + ONE*EXPON*Y*Y
C
C     XYY 
C
        ELSE IF (IC(ITYPE+I) .EQ. 17) THEN
          BF =  EXPON*X*Y*Y
          GX = -EXPON*Y*Y*EE*X*X + ONE*EXPON*Y*Y
          GY = -EXPON*X*Y*EE*Y*Y + TWO*EXPON*X*Y
          GZ = -EXPON*Y*Y*EE*X*Z
C
C     XZZ 
C
        ELSE IF (IC(ITYPE+I) .EQ. 18) THEN
          BF =  EXPON*X*Z*Z
          GX = -EXPON*Z*Z*EE*X*X + ONE*EXPON*Z*Z
          GY = -EXPON*Z*Z*EE*X*Y 
          GZ = -EXPON*X*Z*EE*Z*Z + TWO*EXPON*X*Z
C
C     YZZ 
C
        ELSE IF (IC(ITYPE+I) .EQ. 19) THEN
          BF =  EXPON*Y*Z*Z
          GX = -EXPON*Z*Z*EE*X*Y
          GY = -EXPON*Z*Z*EE*Y*Y + ONE*EXPON*Z*Z
          GZ = -EXPON*Y*Z*EE*Z*Z + TWO*EXPON*Y*Z
C
C     XYZ
C
        ELSE IF (IC(ITYPE+I) .EQ. 20) THEN
          BF =  EXPON*X*Y*Z
          GX = -EXPON*Y*Z*EE*X*X + ONE*EXPON*Y*Z
          GY = -EXPON*X*Z*EE*Y*Y + ONE*EXPON*X*Z
          GZ = -EXPON*X*Y*EE*Z*Z + ONE*EXPON*X*Y
C
        END IF
C
C     END F-FUNCTIONS
C
C    CALCULATE DENSITY AND GRADIENT COMPONENTS FOR EACH MO.            
C
        DO 200 J = 1,NMO
          CO(IPSI+J) = CO(IPSI+J) + CO(IMO+NPRIMS*(J-1)+I)*BF
          CO(IGX+J)  = CO(IGX+J)  + CO(IMO+NPRIMS*(J-1)+I)*GX
          CO(IGY+J)  = CO(IGY+J)  + CO(IMO+NPRIMS*(J-1)+I)*GY
          CO(IGZ+J)  = CO(IGZ+J)  + CO(IMO+NPRIMS*(J-1)+I)*GZ
200   CONTINUE
      RETURN
      END
      SUBROUTINE GAUS4
C+++
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C+++
C
C    FOR GAUSSIAN WAVEFUNCTIONS ONLY.  CALCULATES AT A GIVEN POINT
C    THE VALUE OF EACH MOLECULAR ORBITAL AND THE MO GRADIENT VECTOR
C    COMPONENTS.
C
      COMMON /OFFSET/ ITYPE,ICENT,KATOM,IEORB,IE,IMO,ICHARG,IXC,IYC,IZC,
     +  IXX, IYY, IZZ,IXS,IYS,IZS,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,ID2,
     +  IGXX,IGXY,IGXZ,IGYY,IGYZ,IGZZ
      COMMON CO(30000),IC(30000),MODE,NCENT,NMO,NPRIMS
C
      DO 310 J=1,NMO
      CO(IPSI+J)=0.0
      CO(IGX+J)=0.0
      CO(IGY+J)=0.0
      CO(IGZ+J)=0.0
      CO(IGXX+J)=0.0
      CO(IGXY+J)=0.0
      CO(IGXZ+J)=0.0
      CO(IGYY+J)=0.0
      CO(IGYZ+J)=0.0
      CO(IGZZ+J)=0.
  310 CONTINUE
      DO 360 I=1,NPRIMS
      K=IC(ICENT+I)
      EE=CO(IE+I)*2.0
      X=CO(IXX+K)
      Y=CO(IYY+K)
      Z=CO(IZZ+K)
      CCOO=-CO(IE+I)*CO(IR2+K)
      EXPON=0.D0
      IF(CCOO.GT.-170.D0) EXPON=DEXP(CCOO)
      XX=X*X
      XY=X*Y
      XZ=X*Z
      YY=Y*Y
      YZ=Y*Z
      ZZ=Z*Z
      EXPEE=EXPON*EE
      IF (EXPON.LT.1.D-25) EXPON = 0.D0
      IF (EXPEE.LT.1.D-25) EXPEE = 0.D0
      GO TO (361,362,363,364,365,366,367,368,369,370) IC(ITYPE+I)
C
C    S
C
  361 CONTINUE
      BF=EXPON
      GX=-EXPEE*X
      GY=-EXPEE*Y
      GZ=-EXPEE*Z
      GXX=-GX*X*EE-EXPEE
      GXY=-GX*EE*Y
      GXZ=-GX*EE*Z
      GYY=-GY*Y*EE-EXPEE
      GYZ=-GY*EE*Z
      GZZ=-GZ*Z*EE-EXPEE
      GO TO 340
C
C    X
C
  362 CONTINUE
      BF=EXPON*X
      Q2=-EXPEE*X
      GX=EXPON*(1.0-EE*XX)
       GY=Q2*Y
      GZ=Q2*Z
      Q2=EE*X
      Q1=Q2*X
      GXX=EXPEE*X*(Q1-3.)
      GXY=EXPEE*Y*(Q1-1.)
      GXZ=EXPEE*Z*(Q1-1.)
      GYY=EXPEE*(Q2*YY-X)
      GYZ=EXPEE*Q2*YZ
      GZZ=EXPEE*(Q2*ZZ-X)
      GO TO 340
C
C    Y
C
  363 CONTINUE
      BF=EXPON*Y
      Q2=-EXPEE*Y
      GX=Q2*X
      GY=EXPON*(1.0-EE*YY)
      GZ=Q2*Z
      Q2=EE*Y
      Q1=Q2*Y
      GXX=EXPEE*(Q2*XX-Y)
      GXY=EXPEE*X*(Q1-1.)
      GXZ=EXPEE*Q2*XZ
      GYY=EXPEE*Y*(Q1-3.)
      GYZ=EXPEE*Z*(Q1-1.)
      GZZ=EXPEE*(Q2*ZZ-Y)
      GO TO 340
C
C    Z
C
  364 CONTINUE
      BF=EXPON*Z
      Q2=-EXPEE*Z
      GX=Q2*X
      GY=Q2*Y
      GZ=EXPON*(1.0-EE*ZZ)
      Q2=EE*Z
      Q1=Q2*Z
      GXX=EXPEE*(Q2*XX-Z)
      GXY=EXPEE*Q2*XY
      GXZ=EXPEE*X*(Q1-1.)
      GYY=EXPEE*(Q2*YY-Z)
      GYZ=EXPEE*Y*(Q1-1.)
      GZZ=EXPEE*Z*(Q1-3.)
      GO TO 340
C
C    XX
C
  365 CONTINUE
      BF=EXPON*XX
      Q2=-BF*EE
      GX=   X*(Q2+2.*EXPON)
      GY=   Y*Q2
      GZ=   Z*Q2
      Q1=EE*XX
      GXX=(2.+(Q1-5.)*Q1)*EXPON
      GXY=EXPEE*XY*(Q1-2.)
      GXZ=EXPEE*XZ*(Q1-2.)
      GYY=EXPON*Q1*(EE*YY-1.)
      GYZ=EXPEE*Q1*YZ
      GZZ=EXPON*Q1*(EE*ZZ-1.)
      GO TO 340
C
C    YY
C
  366 CONTINUE
      BF=EXPON*YY
      Q2=-BF*EE
      GX=   X*Q2
      GY=   Y*(Q2+2.*EXPON)
      GZ=   Z*Q2
      Q1=EE*YY
      GXX=EXPON*Q1*(EE*XX-1.)
      GXY=EXPEE*XY*(Q1-2.)
      GXZ=EXPEE*XZ*Q1
      GYY=EXPON*(2.+(Q1-5.)*Q1)
      GYZ=EXPEE*YZ*(Q1-2.)
      GZZ=EXPON*Q1*(EE*ZZ-1.)
      GO TO 340
C
C    ZZ
C
  367 CONTINUE
      BF=EXPON*ZZ
      Q2=-BF*EE
      GX=   X*Q2
      GY=   Y*Q2
      GZ=   Z*(Q2+2.*EXPON)
      Q1=EE*ZZ
      GXX=EXPON*Q1*(EE*XX-1.)
      GXY=EXPEE*Q1*XY
      GXZ=EXPEE*XZ*(Q1-2.)
      GYY=EXPON*Q1*(EE*YY-1.)
      GYZ=EXPEE*YZ*(Q1-2.)
      GZZ=EXPON*(2.+(Q1-5.)*Q1)
      GO TO 340
C
C    XY
C
  368 Q1=   EXPON*X
      BF=Q1*Y
      Q3=-BF*EE
      GX=   Q3*X + EXPON*Y
      GY=   Q3*Y + Q1
      GZ=   Q3*Z
      Q1=EE*XX
      Q2=EE*YY
      GXX=EXPEE*XY*(Q1-3.)
      GXY=(Q1-1.)*(Q2-1.)*EXPON
      GXZ=EXPEE*YZ*(Q1-1.)
      GYY=EXPEE*XY*(Q2-3.)
      GYZ=EXPEE*XZ*(Q2-1.)
      GZZ=EXPEE*XY*(EE*ZZ-1.)
      GO TO 340
C
C    XZ
C
  369 Q1=   EXPON*Z
      BF=Q1*X
      Q3=-BF*EE
      GX=   Q3*X + Q1
      GY=   Q3*Y
      GZ=   Q3*Z + EXPON*X
      Q1=EE*XX
      Q2=EE*ZZ
      GXX=EXPEE*XZ*(Q1-3.)
      GXY=EXPEE*YZ*(Q1-1.)
      GXZ=EXPON*(Q1-1.)*(Q2-1.)
      GYY=EXPEE*XZ*(EE*YY-1.)
      GYZ=EXPEE*XY*(Q2-1.)
      GZZ=EXPEE*XZ*(Q2-3.)
      GO TO 340
C
C    YZ
C
  370 Q1=   EXPON*Y
      BF=Q1*Z
      Q3=-BF*EE
      GX=   Q3*X
      GY=   Q3*Y + EXPON*Z
      GZ=   Q3*Z + Q1
      Q1=EE*YY
      Q2=EE*ZZ
      GXX=EXPEE*YZ*(EE*XX-1.)
      GXY=EXPEE*XZ*(Q1-1.)
      GXZ=EXPEE*XY*(Q2-1.)
      GYY=EXPEE*YZ*(Q1-3.)
      GYZ=(Q1-1.)*(Q2-1.)*EXPON
      GZZ=EXPEE*YZ*(Q2-3.)
  340 CONTINUE
C
C    CALCULATE DENSITY AND GRADIENT COMPONENTS FOR EACH MO.
C
      DO 360 J=1,NMO
      CIJ=CO(IMO+NPRIMS*(J-1)+I)
      CO(IPSI+J)=CO(IPSI+J)+CIJ*BF
      CO(IGX+J)=CO(IGX+J)+CIJ*GX
      CO(IGY+J)=CO(IGY+J)+CIJ*GY
      CO(IGZ+J)=CO(IGZ+J)+CIJ*GZ
      CO(IGXX+J)=CO(IGXX+J)+CIJ*GXX
      CO(IGXY+J)=CO(IGXY+J)+CIJ*GXY
      CO(IGXZ+J)=CO(IGXZ+J)+CIJ*GXZ
      CO(IGYY+J)=CO(IGYY+J)+CIJ*GYY
      CO(IGYZ+J)=CO(IGYZ+J)+CIJ*GYZ
      CO(IGZZ+J)=CO(IGZZ+J)+CIJ*GZZ
  360 CONTINUE
      RETURN
      END
      SUBROUTINE GROCKLE (N, X, IR, S, E)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(3,N), S(3,N), C(3), E(3,3), EV(3), R(3,3)
      INTEGER IR(N)
      PARAMETER (ZERO = 0.0D0)
C
C    ZERO OUT CENTROID AND EIGENVECTORS
C
      C(1) = ZERO
      C(2) = ZERO
      C(3) = ZERO
      E(1,1) = ZERO
      E(1,2) = ZERO
      E(1,3) = ZERO
      E(2,2) = ZERO
      E(2,3) = ZERO
      E(3,3) = ZERO
C
C    CENTROID OF FRAGMENT
C
      M = ZERO
      DO 100 I = 1, N
      IF (IR(I) .GT. 0) THEN
        C(1) = C(1) + X(1,IR(I))
        C(2) = C(2) + X(2,IR(I))
        C(3) = C(3) + X(3,IR(I))
	M = M + 1
      ELSE IF (IR(I) .LT. 0) THEN
        C(1) = C(1) + S(1,-IR(I))
        C(2) = C(2) + S(2,-IR(I))
        C(3) = C(3) + S(3,-IR(I))
	M = M + 1
      END IF
100     CONTINUE
C
      DD = 1.0D0/DFLOAT(M)
      C(1) = DD*C(1)
      C(2) = DD*C(2)
      C(3) = DD*C(3)
C
C    CALCULATE INERTIAL MATRIX.
C
      DO 200 I = 1, N
      IF (IR(I) .GT. 0) THEN
        X1 = X(1,IR(I)) - C(1)
        X2 = X(2,IR(I)) - C(2)
        X3 = X(3,IR(I)) - C(3)
      ELSE IF (IR(I) .LT. 0) THEN
        X1 = S(1,-IR(I)) - C(1)
        X2 = S(2,-IR(I)) - C(2)
        X3 = S(3,-IR(I)) - C(3)
      END IF
      E(1,1) = E(1,1) + X2**2 + X3**2
      E(2,2) = E(2,2) + X1**2 + X3**2
      E(3,3) = E(3,3) + X1**2 + X2**2
      E(1,2) = E(1,2) - X1*X2
      E(1,3) = E(1,3) - X1*X3
      E(2,3) = E(2,3) - X2*X3
      X1 = ZERO
      X2 = ZERO
      X3 = ZERO
200     CONTINUE
C
      E(2,1) = E(1,2)
      E(3,1) = E(1,3)
      E(3,2) = E(2,3)
C
C    GENERATES EIGENVALUES AND EIGENVECTORS OF THE INERTIAL MATRIX.
C
      CALL TRACE(E, EV, C, 3, IFAIL)
C
C    CHECK FOR RIGHT HAND CONVENTION FOR EIGEN-AXES
C
      DET = E(1,1)*(E(2,2)*E(3,3) - E(3,2)*E(2,3)) +
     +      E(1,2)*(E(3,1)*E(2,3) - E(2,1)*E(3,3)) +
     +      E(1,3)*(E(2,1)*E(3,2) - E(3,1)*E(2,2))
C
      IF (DET .LT. ZERO) THEN
        E(1,2) = -E(1,2)
        E(2,2) = -E(2,2)
        E(3,2) = -E(3,2)
      END IF
C
      RETURN
      END
      SUBROUTINE HESS(X,Y,Z,W,H)
C
C    HESS DETERMINES THE HESSIAN OF RHO AT THE SPECIFIED POINT AND
C    RETURNS THE EIGENVALUES (W) AND THE EIGENVECTORS (H) OF THE
C    HESSIAN.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION F(3),H(3,3),W(3)
C
      COMMON CO(30000),IC(30000),MODE,NCENT,NMO,NPRIMS
C
      COMMON /OFFSET/ ITYPE,ICENT,KATOM,IEORB,IE,IMO,ICHARG,IXC,IYC,IZC,
     + IXX,IYY,IZZ,IXS,IYS,IZS,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,ID2,IGXX,
     + IGXY,IGXZ,IGYY,IGYZ,IGZZ
C
      DO 100 I = 1,NCENT
        CO(IXX+I) = X - CO(IXC+I)
        CO(IYY+I) = Y - CO(IYC+I)
        CO(IZZ+I) = Z - CO(IZC+I)
        CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +
     +              CO(IYY+I)*CO(IYY+I) +
     +              CO(IZZ+I)*CO(IZZ+I)
        CO(IRR+I) = DSQRT(CO(IR2+I))
100   CONTINUE
C
      CALL GAUS4
C
      H(1,1) = 0.D0
      H(1,2) = 0.D0
      H(1,3) = 0.D0
      H(2,2) = 0.D0
      H(2,3) = 0.D0
      H(3,3) = 0.D0
C
C     BUILD UP THE GRADIENT AND THE HESSIAN ELEMENTS.
C
      DO 110 I = 1,NMO
        PROD = 2.*CO(IP+I)*CO(IPSI+I)
        H(1,1) = H(1,1)+PROD*CO(IGXX+I)+2.*CO(IP+I)*CO(IGX+I)*CO(IGX+I)
        H(1,2) = H(1,2)+PROD*CO(IGXY+I)+2.*CO(IP+I)*CO(IGX+I)*CO(IGY+I)
        H(1,3) = H(1,3)+PROD*CO(IGXZ+I)+2.*CO(IP+I)*CO(IGX+I)*CO(IGZ+I)
        H(2,2) = H(2,2)+PROD*CO(IGYY+I)+2.*CO(IP+I)*CO(IGY+I)*CO(IGY+I)
        H(2,3) = H(2,3)+PROD*CO(IGYZ+I)+2.*CO(IP+I)*CO(IGY+I)*CO(IGZ+I)
        H(3,3) = H(3,3)+PROD*CO(IGZZ+I)+2.*CO(IP+I)*CO(IGZ+I)*CO(IGZ+I)
110   CONTINUE
      DO 120 I = 1,3
        DO 120 J = 1,I
          H(I,J) = H(J,I)
120   CONTINUE
C
C    DIAGONALIZE THE HESSIAN
C
      CALL TRACE(H,W,F,3,IERR)
C
      RETURN
      END
      SUBROUTINE MAKNAME(I,STRING,L,EXT)
      CHARACTER*(*) STRING,EXT
      INTEGER I,J
      CALL GETARG(I,STRING)
      J = LEN(STRING)
      DO 10 N = 1,J
        IF(STRING(N:N) .EQ. ' ') THEN
          L = N - 1
          STRING = STRING(1:L)//EXT
          RETURN
        ENDIF
10    CONTINUE
      STOP ' FAILED TO MAKE A FILE NAME '
      RETURN
      END
      SUBROUTINE NUCLEI(XN,YN,ZN)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /TRANS/ A(3,3),CR(2),CX(3),CT(3),SCAL,IPROJ
C
      DATA SN1,SN2 /.06,.04/
C
      XN = (XN + CR(1))*SCAL
      YN = (YN + CR(2))*SCAL
C
      CALL PLOT (SNGL(XN-SCAL*SN1),SNGL(YN),3)
C
      IF (DABS(ZN) .GE. 1.0D-01) THEN
        CALL PLOT (SNGL(XN-SCAL*SN2),SNGL(YN),2)
        CALL PLOT (SNGL(XN+SCAL*SN2),SNGL(YN),3)
        CALL PLOT (SNGL(XN+SCAL*SN1),SNGL(YN),2)
        CALL PLOT (SNGL(XN),SNGL(YN-SCAL*SN1),3)
        CALL PLOT (SNGL(XN),SNGL(YN-SCAL*SN2),2)
        CALL PLOT (SNGL(XN),SNGL(YN+SCAL*SN2),3)
        CALL PLOT (SNGL(XN),SNGL(YN+SCAL*SN1),2)
      ELSE 
        CALL PLOT (SNGL(XN+SCAL*SN1),SNGL(YN),2)
        CALL PLOT (SNGL(XN),SNGL(YN-SCAL*SN1),3)
        CALL PLOT (SNGL(XN),SNGL(YN+SCAL*SN1),2)
      END IF
      RETURN
      END
C SKK ================================================================== SKK
C
        FUNCTION        NUMBER	(LINE, LPST, NUM, DEC)
C
C CONVERTS A CHARACTER STRING OF NUMBERS INTO ACTUAL NUMBERS EITHER
C INTEGERS OR DECIMAL MAY BE READ.
C NUMBER = 1 IF ALL THE REMAINING CHARACTERS ARE BLANK
C        = 2 IF CHARACTERS ARE NOT RECOGNISED AS A NUMBER, LPST IS RESET
C SKK ================================================================== SKK

        DOUBLE PRECISION DEC, TEN
        CHARACTER*(*)   LINE
        CHARACTER       BLANK, COMMA, DOT, MINUS, L
        CHARACTER       CTEN(0:9)
        DATA    CTEN    /'0','1','2','3','4','5','6','7','8','9'/
        PARAMETER       (BLANK = ' ', COMMA = ',')
        PARAMETER       (DOT   = '.', MINUS = '-')
        INTEGER         ITEN
        PARAMETER       (ITEN = 10, TEN = 10.0D0)
        NUM     = 0
        DEC     = 0.0D0
        NP      = 0
        ND      = 0
        MS      = 0
        NUMBER	= 0
        LPEND   = LEN (LINE)
5       IF (LINE(LPST:LPST) .EQ. BLANK) THEN
                LPST    = LPST + 1
                IF (LPST .GT. LPEND) THEN
                        NUMBER	= 1
                        RETURN
                END IF
                GOTO 5
        END IF
        LBEFOR  = LPST

        DO 1 I  = LBEFOR, LPEND
        LPST    = I
        L       = LINE(I:I)
        IF (L .EQ. BLANK .OR. L .EQ. COMMA) THEN
                GOTO 2
        ELSE IF (L .EQ. MINUS) THEN
                MS      = 1
                GOTO 1
        ELSE IF (L .EQ. DOT) THEN
                NP      = 1
                GOTO 1
        END IF
        DO 3 J  = 0, 9
        IF (L .EQ. CTEN(J)) THEN
                N       = J
                GOTO 4
        END IF
3       CONTINUE
        NUMBER	= 2
        LPST    = LBEFOR
        RETURN

4       IF (NP .EQ. 1) THEN
                ND      = ND + 1
                DEC     = DEC + DFLOAT(N)/TEN**ND
        ELSE
                NUM     = NUM*ITEN + N
        END IF
1       CONTINUE

2       DEC     = DFLOAT(NUM) + DEC
        IF (MS .EQ. 0) RETURN
        DEC     = -DEC
        NUM     = -NUM
        RETURN
        END

      SUBROUTINE RDWFN
C+++
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C+++
C
C  READS WAVEFUNCTION.  INPUT (PREPARED BY PROGRAM PSI) CONSISTS OF
C     LABEL  - A LABEL FOR THE WAVEFUNCTION
C     MODE   - WAVEFUNCTION TYPE (SLATER OR GAUSSIAN)
C     NMO    - NO. OF MOLECULAR ORBITALS
C     NPRIMS -NO. OF (PRIMITIVE) BASIS FUNCTIONS
C     NCENT  - NO. OF NUCLEI
C
C  THEN FOR EACH NUCLEUS,
C     KATOM      - NAME
C     X/Y/ZCENTR - COORDINATES
C     CHARGE     - ATOMIC NUMBER
C
C  AND FOR EACH BASIS FUNCTION
C     ICENT  - THE NO. OF THE NUCLEUS ON WHICH IT IS CENTERED
C     ITYPE  - FUNCTION TYPE (SEE ARRAYS LABELS AND LABELG)
C     E      - EXPONENT
C
C  AND FOR EACH MOLECULAR ORBITAL,
C     MOLAB - A LABEL
C     P     - OCCUPATION NUMBER
C     EORB  - ORBITAL ENERGY
C     CO    - COEFFICIENTS OF PRIMITIVE BASIS FUNCTIONS.  THESE
C             INCLUDE ALL NORMALIZATION AND CONTRACTION COEFFICIENTS.
C
      CHARACTER*80 WFNTTL,JOBTTL
      CHARACTER*8 ATNAM
      DIMENSION LABELS(31),LABELG(20)
C
      COMMON /DATA/ MODES,MODEG
C
      COMMON /UNITS/  ISRF,INPT,IOUT,IWFN,IDBG
C
      COMMON /C10/ GC
C
      COMMON /STRING/ WFNTTL,JOBTTL,ATNAM(30),NAT
C
      COMMON /C29/ THRESH1,THRESH2,GAMMA,TOTE
C
      DIMENSION MOLAB(2)
C
      COMMON /OFFSET/ ITYPE,ICENT,KATOM,IEORB,IE,IMO,ICHARG,IXC,IYC,IZC,
     +  IXX, IYY, IZZ,IXS,IYS,IZS,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,ID2,IGXX,
     +  IGXY,IGXZ,IGYY,IGYZ,IGZZ
      COMMON CO(30000),IC(30000),MODE,NCENT,NMO,NPRIMS
C
      REAL*8  LABELS
C                             TO HOLD 8 CHARACTERS PER PARCEL
C
      DATA MODEG,  MODES / 4HSIAN, 4HER   /
C
      DATA LABELS /2H1S,2H2S,4H2P X,4H2P Y,4H2P Z,2H3S,4H3P X,4H3P Y,
     14H3P Z,5H3D Z2,8H3D X2-Y2,5H3D XY,5H3D XZ,5H3D YZ,4H4P X,4H4P Y,
     24H4P Z,5H4D Z2,5H4D XY,5H4D XZ,5H4D YZ,5H4F Z3,6H4F Z2X,6H4F Z2Y,
     36H4F Z2Z,2H4S,8H4D X2-Y2,2H5S,4H5P X,4H5P Y,4H5P Z/
C
C
      DATA LABELG /1HS,1HX,1HY,1HZ,2HXX,2HYY,2HZZ,2HXY,2HXZ,2HYZ,3HXXX,
     13HYYY,3HZZZ,3HXXY,3HXXZ,3HXYY,3HYYZ,3HXZZ,3HYZZ,3HXYZ/
C
      DATA JP/1/,JG/1/,JD/1/,JHESS/1/,JPL/0/
C
      DATA ENDATA /8HEND DATA/
C
      READ(IWFN,101) WFNTTL
      READ(IWFN,102) MODE,NMO,NPRIMS,NCENT
C
      CONTINUE
      GCO=0.
      GC=0.
C
C   THE QUANTITIES BELOW ARE OFFSETS IN VECTORS CO AND IO
C
      IXS=0
      IYS=NCENT*JPL
      IZS=IYS+NCENT*JPL
C                             FUNCTION TYPES...
      ITYPE=IZS+NCENT*JPL
C                             FUNCTION CENTRES
      ICENT=ITYPE+NPRIMS
C                             CENTRE NUMBER OF EACH FUNCTION
      KATOM=ICENT+NPRIMS
C                             ORBITAL ENERGIES OF EACH M.O.
      IEORB=KATOM+NCENT
C                             EXPONENTS OF THE FUNCTIONS
      IE=IEORB+NMO
C                             NUCLEAR CHARGE OF EACH CENTRE
      ICHARG=IE+NPRIMS
C
C                             XYZ COORDS OF THE CENTRES IN
C                             ORIGINAL CARTESIAN SYSTEM
      IXC=ICHARG+NCENT
      IYC=IXC+NCENT
      IZC=IYC+NCENT
C
C                             XYZ COORDS OF CENTRES RELATIVE TO A
C                             CURRENT TEST POINT IN THE INTEGRATION
      IXX=IZC+NCENT
      IYY=IXX+NCENT
      IZZ=IYY+NCENT
C
C                             SQUARE OF DISTANCE FROM CENTRES TO
C                             CURRENT POINT
      IRR=IZZ+NCENT
C                             DISTANCE FROM CENTRES TO CURRENT POINT
      IR2=IRR+NCENT
C                             NUMBER OF ELECTRONS IN THE M.O.
      IP=IR2+NCENT
C                             PSI VALUES FOR EACH M.O.
      IPSI=IP+NMO
C
C                              GRAD VALUES (XYZ) FOR EACH M.O.
      IGX=IPSI+NMO*JP
      IGY=IGX+NMO*JG
      IGZ=IGY+NMO*JG
C                             DEL2RHO FOR EACH M.O.
      ID2=IGZ+NMO*JG
C
C                             OFFSET FOR STORAGE OF M.O.S
      IMO=ID2+NMO*JD
C
C     OFFSET FOR STORAGE OF THE HESSIAN
C
      IGXX=IMO+NMO*NPRIMS
      IGXY=IGXX+NMO*JHESS
      IGXZ=IGXY+NMO*JHESS
      IGYY=IGXZ+NMO*JHESS
      IGYZ=IGYY+NMO*JHESS
      IGZZ=IGYZ+NMO*JHESS
C
C    SET CORE TO REQUIRED LENGTH
C
C     J=LOCF(CO(IGZZ+NMO*JHESS+4))
C                             PROVE THAT WE HAVE ENOUGH MEMORY
      COTST=CO(IGZZ+NMO*JHESS+2)
      DO 2 I=1,NCENT
        READ(IWFN,103) ATNAM(I),J,CO(IXC+J),CO(IYC+J),CO(IZC+J),
     +  CO(ICHARG+J)
2     CONTINUE
      READ(IWFN,104) (IC(ICENT+I),I=1,NPRIMS)
      READ(IWFN,104) (IC(ITYPE+I),I=1,NPRIMS)
      READ(IWFN,105) (CO(IE+I),I=1,NPRIMS)
      DO 7 I=1,NMO
        READ(IWFN,106) (MOLAB(J),J=1,2),CO(IP+I),CO(IEORB+I)
        GCO=DMAX1(GCO,DABS(CO(IP+I)))
        K=NPRIMS*(I-1)+IMO
        READ(IWFN,107) (CO(K+J),J=1,NPRIMS)
        DO 7 J=1,NPRIMS
          GC=DMAX1(GC,DABS(CO(J+K)))
7     CONTINUE
      GC=DLOG(GC*GC*GCO*40000.)
      READ(IWFN,108) CHECK
C
C    READ IN TOTAL SCF ENERGY AND -V/T
C
      READ(IWFN,109) TOTE,GAMMA
      RETURN
C
101   FORMAT (A80)
102   FORMAT (BZ,4X,A4,12X,3(I3,17X))
103   FORMAT(A8,11X,I3,2X,3F12.8,10X,F5.1)
104   FORMAT (BZ,20X,20I3)
105   FORMAT (BZ,10X,5E14.7)
106   FORMAT(BZ,10X,2A5,15X,F12.8,15X,F12.8)
107   FORMAT (BZ,5E16.8)
108   FORMAT(BZ,A8)
109   FORMAT(17X,F20.12,18X,F13.8)
      END
      SUBROUTINE SHFTVC(A,F,H,IUD,SV)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(3,3),F(3),H(3,3),SV(3),X(3)
C
      ADOT = 1.0D0
      XDOT = 1.0D0
C
C    DETERMINE NORMAL TO PLANE
C
      X(1) = 1.0D0*A(1,3)
      X(2) = 1.0D0*A(2,3)
      X(3) = 1.0D0*A(3,3)
C
      DO 100 K = 1,3
        IF (F(K) .GT. 0.0D0) 
     +  XDOT = H(K,1)*X(1) + H(K,2)*X(2) + H(K,3)*X(3)
        IF (XDOT .LT. ADOT) THEN 
          MDOT = K
          ADOT = XDOT
        END IF
100   CONTINUE
C
      IF (IUD .EQ. 1) THEN
C
C    TAKE POSITIVE EIGENVECTOR WITH SMALLEST DOT PRODUCT
C    TAKEN WITH NORMAL TO PLANE TO GET POSITIVE 
C    EIGENVECTOR IN PLANE
C
        SV(1) = H(1,MDOT)
        SV(2) = H(2,MDOT)
        SV(3) = H(3,MDOT)
      ELSE IF (IUD .EQ. 0) THEN
C
C    TAKE CROSS PRODUCT OF POSITIVE EIGENVECTOR IN PLANE
C    WITH NORMAL TO PLANE TO GET INITIAL STARTING SHIFT
C    VECTOR FOR DECENT IN GRAD RHO
C
        SV(1) =  H(2,MDOT)*X(3) - H(3,MDOT)*X(2)
        SV(2) =  H(3,MDOT)*X(1) - H(1,MDOT)*X(3)
        SV(3) =  H(1,MDOT)*X(2) - H(2,MDOT)*X(1)
      END IF
      RETURN
      END
        SUBROUTINE      TQLGRM	(N, D, E, Z, IERR)
        IMPLICIT        DOUBLE PRECISION (A-H, O-Z)
        DIMENSION       D(*), E(*), Z(N,N)
	PARAMETER (AMACH = 16.0E-13)
	PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
        IERR    = 0
        IF (N .EQ. 1) RETURN
C
	DO 30 I = 2,N
	E(I-1)  = E(I)
30	CONTINUE

        F       = ZERO
        B       = ZERO
        E(N)    = ZERO

	DO 31 L = 1,N
	J       = 0
	H       = AMACH*(DABS(D(L)) + DABS(E(L)))
	IF (B .LT. H) B = H

105		DO 32 M = L,N
		IF (DABS(E(M)) .LE. B) GOTO 120
32		CONTINUE

120	IF (M .EQ. L) GOTO 220

130	IF (J .EQ. 30) THEN
		IERR    = L
		RETURN
	END IF

	J       = J + 1
	L1      = L + 1
	G       = D(L)
	P       = (D(L1) - G)/(2*E(L))
	IF (DABS(P*AMACH) .GT. ONE) THEN
		R       = P
	ELSE
		R       = DSQRT(P*P + 1)
	END IF
	D(L)    = E(L)/(P + DSIGN(R,P))
	H       = G - D(L)

		DO 33 I = L1,N
		D(I)    = D(I) - H
33		CONTINUE

	F       = F + H
	P       = D(M)
	C       = ONE
	S       = ZERO
	MML     = M - L

		DO 34 II = 1,MML
		I       = M - II
		G       = C*E(I)
		H       = C*P
		IF (DABS(P) .GE. DABS(E(I))) THEN
			C       = E(I)/P
			R       = DSQRT(C*C + 1)
			E(I+1)  = S*P*R
			S       = C/R
			C       = ONE/R
		ELSE
			C       = P/E(I)
			R       = DSQRT(C*C + 1)
			E(I+1)  = S*E(I)*R
			S       = 1.D0/R
			C       = C*S
		END IF
		P       = C*D(I) - S*G
		D(I+1)  = H + S*(C*G + S*D(I))

			DO 35 K = 1,N
			H       = Z(K,I+1)
			Z(K,I+1)= S*Z(K,I) + C*H
			Z(K,I)  = C*Z(K,I) - S*H
35			CONTINUE

34		CONTINUE

	E(L)    = S*P
	D(L)    = C*P
	IF (DABS(E(L)) .GT. B) GOTO 130

220	D(L)    = D(L) + F
31	CONTINUE

	DO 300 II = 2,N
	I       = II - 1
	K       = I
	P       = D(I)

		DO 260 J = II,N
		IF (D(J) .GE. P) GOTO 260
		K       = J
		P       = D(J)
260		CONTINUE

	IF (K .EQ. I) GOTO 300
	D(K)    = D(I)
	D(I)    = P

		DO 37 J = 1,N
		P       = Z(J,I)
		Z(J,I)  = Z(J,K)
		Z(J,K)  = P
37		CONTINUE

300	CONTINUE
        RETURN
        END
        SUBROUTINE	TRACE	(H, E, W, N, IERR)
C
C TRACE CALLS TREDIG AND TLQGRM TO DIAGONALIZE A SYMMETRIC REAL MATRIX.
C THE MATRIX IS PASSED DOWN IN H AND IS REPLACED BY THE EIGENVECTORS.
C THE EIGENVALUES IN E ARE STORED SMALLEST FIRST.
C THE WORK STORE W SHOULD BE AT LEAST OF DIMENSION N.
C SKK ==================================================================
C
        IMPLICIT        DOUBLE PRECISION (A-H, O-Z)
        DIMENSION       H(N,N), E(N), W(N)
C
        CALL TREDIG	(N, E, W, H)
        CALL TQLGRM	(N, E, W, H, IERR)
C
        RETURN
        END
        SUBROUTINE      TREDIG	(N, D, E, Z)
        IMPLICIT        DOUBLE PRECISION (A-H, O-Z)
        DIMENSION       D(N), E(N), Z(N,N)
	PARAMETER	(ZERO = 0.0D0, ONE = 1.0D0)

        IF (N .EQ. 1) GOTO 320

	DO 30 II = 2,N
	I       = N + 2 - II
	L       = I - 1
	H       = ZERO
	SCALE   = ZERO

	IF (L .LT. 2) GOTO 130

		DO 31 K = 1,L
		SCALE   = SCALE + DABS(Z(I,K))
31              CONTINUE

	IF (SCALE .NE. ZERO) GOTO 140
130     E(I)    = Z(I,L)
        GOTO 290

140	RSCALE	= ONE/SCALE
		DO 32 K = 1,L
		Z(I,K)  = Z(I,K)*RSCALE
		H       = H + Z(I,K)*Z(I,K)
32		CONTINUE
	F       = Z(I,L)
	G       = -DSIGN(DSQRT(H),F)
	E(I)    = SCALE*G
	H       = H - F*G
	Z(I,L)  = F - G
	F       = ZERO
	RH	= ONE/H
	RHSCALE	= RH*RSCALE	

		DO 33 J = 1,L
		Z(J,I)  = Z(I,J)*RHSCALE
		G       = ZERO

			DO 34 K = 1,J
			G       = G + Z(J,K)*Z(I,K)
34                      CONTINUE

		JP1     = J + 1
		IF (L .LT. JP1) GOTO 220

			DO 35 K = JP1,L
			G       = G + Z(K,J)*Z(I,K)
35                      CONTINUE

220		E(J)    = G*RH
		F       = F + E(J)*Z(I,J)
33		CONTINUE

	HH      = F/(H + H)

		DO 36 J = 1,L
		F       = Z(I,J)
		G       = E(J) - HH*F
		E(J)    = G
			DO 37 K = 1,J
			Z(J,K)  = Z(J,K) - F*E(K) - G*Z(I,K)
37			CONTINUE
36		CONTINUE

		DO 38 K	= 1,L
		Z(I,K)  =  SCALE*Z(I,K)
38		CONTINUE

290	D(I)    = H
30	CONTINUE

320     D(1)    = ZERO
        E(1)    = ZERO

	DO 500 I = 1,N
	L       = I - 1
	IF (D(I) .EQ. ZERO) GOTO 380

		DO 40 J	= 1,L
		G       = ZERO

			DO 41 K	= 1,L
			G       = G + Z(I,K)*Z(K,J)
41			CONTINUE

			DO 42 K = 1,L
			Z(K,J)  = Z(K,J) - G*Z(K,I)
42			CONTINUE

40		CONTINUE

380	D(I)    = Z(I,I)
	Z(I,I)  = ONE
	IF(L .LT. 1) GOTO 500

		DO 43 J	= 1,L
		Z(J,I)  = ZERO
		Z(I,J)  = ZERO
43		CONTINUE

500	CONTINUE
        RETURN
        END
      SUBROUTINE TRUDGE(X0,Y0,Z0,IUP,IINC,END)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON CO(30000),IC(30000),MODE,NCENT,NMO,NPRIMS
      COMMON /OFFSET/ ITYPE,ICENT,KATOM,IEORB,IE,IMO,ICHARG,IXC,IYC,IZC,
     +  IXX, IYY, IZZ,IXS,IYS,IZS,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,ID2,IGXX,
     +  IGXY,IGXZ,IGYY,IGYZ,IGZZ
      COMMON /TRANS/ A(3,3),CR(2),CX(3),CT(3),SCALE,IPROJ
C
      DATA STP /1.0D-03/, STPINC /1.0D-04/
C
      STEP = STP
      ENDPT = END
      DXO = 0.0D0
      DYO = 0.0D0
      DZO = 0.0D0
C
      XX = (X0 + CR(1))*SCALE
      YY = (Y0 + CR(2))*SCALE
C
      CALL PLOT(SNGL(XX),SNGL(YY),3)
C
C    DETERMINE IF THIS WALK IS AN ASCENT OR DECENT
C
      IF (IUP .EQ. 0) DIR =  1.0D0
      IF (IUP .EQ. 1) DIR = -1.0D0

C    TRANSFORM POINT FROM PLOTTING PLANE TO MOLECULAR
C    FRAME OF REFERENCE
C
50      IF (IPROJ .EQ. 0) Z0=0.0D0
      X = A(1,1)*X0 + A(1,2)*Y0 + A(1,3)*Z0+CX(1)
      Y = A(2,1)*X0 + A(2,2)*Y0 + A(2,3)*Z0+CX(2)
      Z = A(3,1)*X0 + A(3,2)*Y0 + A(3,3)*Z0+CX(3)
C
      XX = (X0 + CR(1))*SCALE
      YY = (Y0 + CR(2))*SCALE
C
C    CHECK TO SEE IF PLOT BOUNDARY CROSSED
C
      IF (XX .GT. 10.0D0 .OR. XX .LT. 0.0D0) GOTO 200
      IF (YY .GT. 10.0D0 .OR. YY .LT. 0.0D0) GOTO 200
C
      CALL PLOT(SNGL(XX),SNGL(YY),2)
C
C    DETERMINE THE GRADIENT VECTOR AT (X,Y,Z)
C
      DO 100 I = 1,NCENT
        CO(IXX+I) = X - CO(IXC+I)
        CO(IYY+I) = Y - CO(IYC+I)
        CO(IZZ+I) = Z - CO(IZC+I)
        CO(IR2+I) = CO(IXX+I)**2 + CO(IYY+I)**2 + CO(IZZ+I)**2
        CO(IRR+I) = DSQRT(CO(IR2+I))
100   CONTINUE
C
      CALL GAUS2
C
      DX = 0.0D0
      DY = 0.0D0
      DZ = 0.0D0
C
      DO 120 I = 1,NMO
        DX = DX - DIR*CO(IP+I)*CO(IPSI+I)*CO(IGX+I)
        DY = DY - DIR*CO(IP+I)*CO(IPSI+I)*CO(IGY+I)
        DZ = DZ - DIR*CO(IP+I)*CO(IPSI+I)*CO(IGZ+I)
120   CONTINUE
C
      DSQ = DSQRT(DX**2 + DY**2 + DZ**2)
C
C    HAS THE WALK BEEN COMPLETED?
C
      IF (DSQ .LT. ENDPT) GOTO 200
C
C    IF FIRST PASS THEN EQUATE DXO, DYO, AND DZO WITH
C    DX, DY, AND DZ RESPECTIVELY
C
      IF (DXO .EQ. 0.0D0 .AND. DYO .EQ. 0.0D0 
     +  .AND. DZO .EQ. 0.0D0) THEN
        DXO = DX
        DYO = DY
        DZO = DZ
      END IF
C
C    AVOID OCILLATION BY CHECKING FOR OPPOSING SIGNS OF THE GRADIENT
C    VECTOR WITH THE PREVIOUS STEP
C
      IF ((DX*DXO + DY*DYO + DZ*DZO) .LT. 0.0D0) GOTO 200
C
      X = X + (STEP/DSQ)*DX
      Y = Y + (STEP/DSQ)*DY
      Z = Z + (STEP/DSQ)*DZ
C
C    STORE THE PRESENT GRADIENT VECTOR TO AVOID OSCILLATION ABOUT
C    A CRITICAL POINT
C
      DXO = DX
      DYO = DY
      DZO = DZ
C
C   IF WALKING DOWN HILL INCREASE STEP SIZE EACH ITERATION
C
      IF (IUP .EQ. 0 .AND. IINC .EQ. 1) STEP = STEP + STPINC
C
C    TRANSFORM NEW POINT TO PLOT FRAME AND DECIDE
C    WHETHER TO ALLOW PROJECTIONS OUT OF PLOTTING
C    PLANE
C
      X0 = A(1,1)*X + A(2,1)*Y + A(3,1)*Z - CT(1)
      Y0 = A(1,2)*X + A(2,2)*Y + A(3,2)*Z - CT(2)
      Z0 = A(1,3)*X + A(2,3)*Y + A(3,3)*Z - CT(3)
C
      GOTO 50
C
200   RETURN
      END

