      PROGRAM SCHUSS 
C
C     VERSION 94 - Revision A
C
C     THIS PROGRAM CALCULATES A SET OF TRAJECTORIES OF THE GRADIENT VECTOR 
C     FIELD OF THE ELECTRON DENSITY STARTING IN A USER SPECIFIED PLANE OF A
C     MOLECULE.  
C
C     GRDVEC MODIFIED, VECTORIZED TO "SCHUSS" BY T.A. KEITH 1991
C
C     SCHUSS CAN HANDLE S,P,D AND F-TYPE GAUSSIAN FUNCTIONS. 
C     MODIFIED TO INCLUDE F-FUNCTIONS BY R.G.A. BONE, FEBRUARY, 1993
C
C     CLEANED UP Feb 1994. TAK
C
C     QUESTIONS AND SUGGESTIONS SHOULD BE DIRECTED TO:
C     Richard Bader McMASTER UNIVERSITY, HAMILTON ONTARIO, CANADA 
C     Bitnet Address:  Bader@mcmail.cis.McMaster.CA
C     or
C     Todd A. Keith:   keith@babbage.chemistry.mcmaster.ca
C
C     TO REDIMENSION SCHUSS TO HANDLE MORE MO's, PRIMITIVES AND NUCLEI
C     CHANGE THE VALUES OF MMO, MPRIMS AND MCENT IN THE PARAMETER
C     STATEMENTS.
C
C     THE PARAMETER MPTS IS THE MAXIMUM NUMBER OF PATHS OF A GIVEN TYPE
C     WHICH CAN BE CALCULATED.  MaxOrg is the Maximum Number of 
C     User-Specified Origins For the Gradient Paths.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 LINE,TITLE
      CHARACTER*40 WVEC,WFN,WGVP
      CHARACTER*7 WORD
      CHARACTER*4 FVEC /'.vec'/, FWFN /'.wfn'/, FGVP /'.gvp'/ 
      PARAMETER (MCENT=100, MPTS=300, MaxOrg=50) 
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /UNITS/  IVEC,IOUT,IWFN
      COMMON /VALUES/ THRESH1,THRESH2,GAMMA,TOTE
      DIMENSION A(3,3),C(3,MCENT),CR(2),CT(3),CX(3),IAR(MAxOrg),
     + IR(MCENT),XS(3,MCENT),YY(MCENT),ZZ(MCENT),NUP(MaxOrg),
     + NDN(MaxOrg),ORG(3,MaxOrg),XYZU(3,MPTS),XYZD(3,MPTS),EUL(3),
     + XYZ0(3,MPTS),XX(MCENT),F(3),SV(3),H(3,3)
C
      Save Zero,One,Two,Ten
      DATA IWFN /10/, IVEC /15/, IOUT /13/,Zero/0.0d0/,One/1.0d0/,
     $  Two/2.0d0/,Ten/10.0d0/
C
      CALL MAKNAME (1,WVEC,ILEN,FVEC)
      IF (ILEN .EQ. 0) STOP 'usage: schuss vecfile wfnfile'
      CALL MAKNAME (1,WGVP,ILEN,FGVP)
      IF (ILEN .EQ. 0) STOP 'usage: schuss vecfile wfnfile'
      CALL MAKNAME (2,WFN,ILEN,FWFN)
      IF (ILEN .EQ. 0) STOP 'usage: schuss vecfile wfnfile'
C
      OPEN (IVEC,FILE=WVEC)
      OPEN (IOUT,FILE=WGVP)
      OPEN (IWFN,FILE=WFN)
C
      PI=DACOS(-ONE)
C
      CALL RDPSI 
C
C     TITLE: 
C
      READ (IVEC,1000) TITLE
C
C     PLOT:
C
      READ (IVEC,1010) LINE
      LPST = 8
      IF (NUMBER(LINE,LPST,NUM,XY) .GT. 0) GOTO 1990
C 
C    DETERMINE CENTER OF LOCAL FRAME (I.E. USER DEFINED CENTER OF PLOT)
C
      CR(1) = XY/Two
      CR(2) = XY/Two
      SCAL = Ten/XY
C
      CALL PLOTS (53,0,13)
C
C    CENTER OF PLOT:
C
      READ (IVEC,1010) LINE
      LPST = 8
      DO 100 I = 1,3
        IF (NUMBER(LINE,LPST,NUM,CX(I)) .GT. 0) GOTO 2000
100   CONTINUE
C
C    WALKING PARAMETERS: 
C
      READ (IVEC,*) R1,R2,R3,ENDPT,IPROJ,IINC
C  
C    DEFINE PLANE: 
      READ (IVEC,1010) LINE
      LPST = 8                                                          
      IF (NUMBER(LINE,LPST,IEG,DNUM) .GT. 0) GOTO 2070 
      IF (IEG .EQ. 0) THEN
      L = 1                                                             
120   IF (NUMBER(LINE,LPST,IR(L),DNUM) .GT. 0) GOTO 10                  
      L = L + 1                                                         
      GOTO 120                                                          
C
C    DECIDE IF DUMMY ATOMS WERE USED AND THEN READ IN THEIR
C    COORDINATES TO XS(3,N)
C
10    DO 130 I = 1,L                                                    
      IF (IR(I) .LT. 0) THEN                                            
      LPST = 8                                                          
      READ (IVEC,1010) LINE
      DO 140 K = 1,3                                                    
      IF (NUMBER(LINE,LPST,NUM,XS(K,ABS(IR(I)))) .GT. 0) GOTO 2010       
140   CONTINUE                                                          
      END IF                                                            
130   CONTINUE                                                          
C
C    GENERATE ROTATION MATRIX 
C
      DO 170 I = 1,NCENT
        C(1,I) = XC(I) 
        C(2,I) = YC(I) 
        C(3,I) = ZC(I) 
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
C    GET ORIGIN COORDINATES
C
      DO 150 I = 1,NORG
      READ (IVEC,*)ORG(1,I),ORG(2,I),ORG(3,I),IAR(I),NUP(I),NDN(I)
150     CONTINUE
C
C    TRANSFORM CENTER OF MOLECULAR FRAME TO CENTER OF PLOTTING FRAME
C
      CT(1) = A(1,1)*CX(1) + A(2,1)*CX(2) + A(3,1)*CX(3)
      CT(2) = A(1,2)*CX(1) + A(2,2)*CX(2) + A(3,2)*CX(3)
      CT(3) = A(1,3)*CX(1) + A(2,3)*CX(2) + A(3,3)*CX(3)
C
C    PLOT NUCLEAR POSITIONS
C
      DO 180 I = 1,NCENT
        XX(I) = A(1,1)*XC(I) + A(2,1)*YC(I) + A(3,1)*ZC(I) 
        YY(I) = A(1,2)*XC(I) + A(2,2)*YC(I) + A(3,2)*ZC(I) 
        ZZ(I) = A(1,3)*XC(I) + A(2,3)*YC(I) + A(3,3)*ZC(I)-CT(3)
        CALL NUCLEI(XX(I),YY(I),ZZ(I),CR,CT,SCAL)
180   CONTINUE
C
C    LOOP OVER EACH ORIGIN
C
         K = 1
         L = 0
         M = 0
        DO 190 I = 1,NORG
C
C    FIND ORIGIN'S COORDINATES IN PLOT'S FRAME OF REFERENCE
C
       X=A(1,1)*ORG(1,I)+A(2,1)*ORG(2,I)+A(3,1)*ORG(3,I)-CT(1)
       Y=A(1,2)*ORG(1,I)+A(2,2)*ORG(2,I)+A(3,2)*ORG(3,I)-CT(2)
       Z=A(1,3)*ORG(1,I)+A(2,3)*ORG(2,I)+A(3,3)*ORG(3,I)-CT(3)
C
C    THIS ORIGIN IS A (2,-2) ATTRACTOR IN THE PLANE
C
        IF (IAR(I) .EQ. 0) THEN
C
          DO 200 J = K,(K+NDN(I)-1)
C
C    FIND STARTING POINT FOR EACH GRADIENT PATH FROM THE ITH ORIGIN (2,-2)
C
            ANG = (J-K)*TWO*PI/NDN(I)
            XYZ0(1,J) = X + R1*DCOS(ANG)
            XYZ0(2,J) = Y + R1*DSIN(ANG)
            XYZ0(3,J) = Zero
200    CONTINUE 
            K = K + NDN(I) 
C
C
        ELSE IF (IAR(I) .EQ. 1) THEN 
C    THIS ORIGIN IS A (2,0) REPELLOR IN THIS PLANE
C
C
C    CALL TO GET HESSIAN AND GRADIENT AT ORIGIN
C
          CALL HESS(ORG(1,I),ORG(2,I),ORG(3,I),F,H)
C
C    DEFINE STARTING POINTS OF TWO ASCENDING AND TWO DESCENDING
C    GRADIENT PATHS FROM THE (2,0) CRITICAL POINT BASED ON THE 
C    PROJECTION INTO THE PLANE OF THE EIGENVECTORS OF THE HESSIAN
C    OF RHO AT THE CRITICAL POINT.
C
       IF (NUP(I) .GT. 0) THEN
       CALL SHFTVC(A,F,H,1,SV)
       L = L + 1
       LL = 2*L
       XYZU(1,LL-1)=X+A(1,1)*SV(1)*R2+A(2,1)*SV(2)*R2+A(3,1)*SV(3)*R2 
       XYZU(2,LL-1)=Y+A(1,2)*SV(1)*R2+A(2,2)*SV(2)*R2+A(3,2)*SV(3)*R2 
       IF (IPROJ .EQ. 1) THEN
       XYZU(3,LL-1)=Z+A(1,3)*SV(1)*R2+A(2,3)*SV(2)*R2+A(3,3)*SV(3)*R2 
       ELSE
       XYZU(3,LL-1) = Zero
       END IF
       XYZU(1,LL)=X-A(1,1)*SV(1)*R2-A(2,1)*SV(2)*R2-A(3,1)*SV(3)*R2 
       XYZU(2,LL)=Y-A(1,2)*SV(1)*R2-A(2,2)*SV(2)*R2-A(3,2)*SV(3)*R2 
       IF (IPROJ .EQ. 1) THEN
       XYZU(3,LL)=Z-A(1,3)*SV(1)*R2-A(2,3)*SV(2)*R2-A(3,3)*SV(3)*R2 
       ELSE
       XYZU(3,LL) = Zero
       END IF
       END IF 
       IF (NDN(I) .GT. 0) THEN
       CALL SHFTVC(A,F,H,0,SV)
       M = M + 1
       MM = 2*M
       XYZD(1,MM-1)=X+A(1,1)*SV(1)*R3+A(2,1)*SV(2)*R3+A(3,1)*SV(3)*R3 
       XYZD(2,MM-1)=Y+A(1,2)*SV(1)*R3+A(2,2)*SV(2)*R3+A(3,2)*SV(3)*R3 
       IF (IPROJ .EQ. 1) THEN
       XYZD(3,MM-1)=Z+A(1,3)*SV(1)*R3+A(2,3)*SV(2)*R3+A(3,3)*SV(3)*R3 
       ELSE
       XYZD(3,MM-1) = Zero
       END IF
       XYZD(1,MM)=X-A(1,1)*SV(1)*R3-A(2,1)*SV(2)*R3-A(3,1)*SV(3)*R3 
       XYZD(2,MM)=Y-A(1,2)*SV(1)*R3-A(2,2)*SV(2)*R3-A(3,2)*SV(3)*R3 
       IF (IPROJ .EQ. 1) THEN
       XYZD(3,MM)=Z-A(1,3)*SV(1)*R3-A(2,3)*SV(2)*R3-A(3,3)*SV(3)*R3 
       ELSE
       XYZD(3,MM) = Zero
       END IF
       END IF
      END IF
C
C    THE STARTING POINTS OF ALL OF THE GRADIENT PATHS TO BE 
C    CALCULATED ARE NOW DEFINED.  THE GRADIENT PATHS DESCENDING
C    FROM THE (2,-2) CRITICAL POINTS ARE CALCULATED FIRST.  THEN
C    THE GRADIENT PATHS DESCENDING FROM (2,0) CRITICAL POINTS.
C    THEN THE GRADIENT PATHS ASCENDING FROM (2,0) CRITICAL POINTS. 
C
190   CONTINUE
      IF ((K-1) .GT. MPTS) STOP 'Too many paths'
      IF (LL .GT. MPTS) STOP 'Too many paths'
      IF (MM .GT. MPTS) STOP 'Too many paths'
      IF(K.GT.1)CALL TRUDGE(XYZ0,K-1,0,A,CT,CR,ENDPT,SCAL,IINC,IPROJ,
     +CX)
      IF(M.GT.0)CALL TRUDGE(XYZD,MM,0,A,CT,CR,ENDPT,SCAL,IINC,IPROJ,
     +CX)
      IF(L.GT.0)CALL TRUDGE(XYZU,LL,1,A,CT,CR,ENDPT,SCAL,IINC,IPROJ,
     +CX)
C
      CALL PLOT (0.,0.,999)
      GOTO 4999
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
2070  WRITE (IOUT,2075)
2075  FORMAT(' ERROR IN ORIENTATION METHOD ')
      GOTO 4999
2080  WRITE (IOUT,2085)
2085  FORMAT( 'ERROR IN EULER ANGLES CARD ')
      GOTO 4999
C
4999  Continue
      STOP 'Schuss is Done'
      END
      SUBROUTINE EULER (E,A)
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C 
      DIMENSION E(3), A(3,3)
      Save One,Two
      DATA ONE /1.0d0/,Two/2.0d0/
C
      PI=DACOS(-ONE)
      RADIAN=Two*Pi/360
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
      A(1,1) =  CC*CA - CB*SA*SC
      A(2,1) =  CC*SA + CB*CA*SC
      A(3,1) =  SC*SB
      A(1,2) = -SC*CA - CB*SA*CC
      A(2,2) = -SC*SA + CB*CA*CC
      A(3,2) =  CC*SB
      A(1,3) =  SB*SA
      A(2,3) = -SB*CA
      A(3,3) =  CB
C
      RETURN
      END
      SUBROUTINE GAUS(XYZ,PTS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER PTS
      PARAMETER(MCENT=100,MMO=100,MPRIMS=420,MPTS=300,NTYPE=20)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),NMO
      COMMON /PRIMS/ COO(MPRIMS,MMO),EXX(MPRIMS),ICT(MPRIMS),
     +       SUM(MPRIMS),Div(MPrims),ITP(NTYPE),NPRIMS
       COMMON /ZZZZ/ PSI(MPTS,MMO),GX(MPTS,MMO),GY(MPTS,MMO),
     + GZ(MPTS,MMO)
      DIMENSION R2(MPTS,MCENT),DX(MPTS,MCENT),DY(MPTS,MCENT),
     +  CHI(MPTS,MPRIMS),CHIX(MPTS,MPRIMS),DZ(MPTS,MCENT),
     +  CHIZ(MPTS,MPRIMS),CHIY(MPTS,MPRIMS),XYZ(MPTS,3),ChiMax(MPrims)
      Save Zero,One,Two,Four,Five,Seven,Three,Cutoff
      DATA ZERO /0.0D0/,ONE/1.0D0/,TWO/2.0D0/,FOUR/4.0D0/,FIVE/5.0D0/,
     +  SEVEN/7.0D0/,THREE/3.0D0/,Cutoff/1.0d-10/
C
      DO 110 J = 1,NCENT
       DO 112 I=1,PTS 
        DX(I,J) = XYZ(I,1) - XC(J)
        DY(I,J) = XYZ(I,2) - YC(J)
        DZ(I,J) = XYZ(I,3) - ZC(J)
        R2(I,J)= DX(I,J)*DX(I,J)+DY(I,J)*DY(I,J)+DZ(I,J)*DZ(I,J)
112   CONTINUE
110   CONTINUE
C
C       FOR S-TYPE
C
        DO 120 J = 1,ITP(1) 
        IS=ICT(J)
        DO 122 I=1,PTS 
C
        A=SUM(J)*DEXP(-EXX(J)*R2(I,is))
C
        CHI(I,J)=A*DIV(J)
        CHIX(I,J)=DX(I,is)*A
        CHIY(I,J)=DY(I,is)*A
        CHIZ(I,J)=DZ(I,is)*A
122     CONTINUE
120     CONTINUE
C
C       FOR Px-TYPE
C
        DO 140 J=ITP(1)+1,ITP(2)
        IS=ICT(J)
        DO 142 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DX(I,is)*A*SUM(J)
C
        CHI(I,J)=A*DX(I,is)
        CHIX(I,J)=A+DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
142     CONTINUE
140     CONTINUE
C
C       FOR Py-TYPE
C
        DO 160 J=ITP(2)+1,ITP(3)
        IS=ICT(J)
        DO 162 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DY(I,is)*A*SUM(J)
C
        CHI(I,J)=A*DY(I,is)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=A+DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
162     CONTINUE
160     CONTINUE
C
C       FOR Pz-TYPE
C
        DO 180 J=ITP(3)+1,ITP(4)
        IS=ICT(J)
        DO 182 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DZ(I,is)*A*SUM(J)
C
        CHI(I,J)=A*DZ(I,is)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=A+DZ(I,is)*B
182     CONTINUE
180     CONTINUE
C
C       FOR Dxx-TYPE
C
        DO 220 J=ITP(4)+1,ITP(5)
        IS=ICT(J)
        DO 222 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DX(I,is)*DX(I,is)*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=(TWO*A+B)*DX(I,is)
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
222     CONTINUE
220     CONTINUE
C
C       FOR Dyy-TYPE
C
        DO 240 J=ITP(5)+1,ITP(6)
        IS=ICT(J)
        DO 242 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DY(I,is)*DY(I,is)*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=(TWO*A+B)*DY(I,is)
        CHIZ(I,J)=DZ(I,is)*B
242     CONTINUE
240     CONTINUE
C
C       FOR Dzz-TYPE
C
        DO 260 J=ITP(6)+1,ITP(7)
        IS=ICT(J)
        DO 262 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DZ(I,is)*DZ(I,is)*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=(TWO*A+B)*DZ(I,is)
262     CONTINUE
260     CONTINUE
C
C       FOR Dxy-TYPE
C
        DO 280 J=ITP(7)+1,ITP(8)
        IS=ICT(J)
        DO 282 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DX(I,is)*DY(I,is)*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,is)*B+DY(I,is)*A
        CHIY(I,J)=DY(I,is)*B+DX(I,is)*A
        CHIZ(I,J)=DZ(I,is)*B
282     CONTINUE
280     CONTINUE
C
C       FOR Dxz-TYPE
C
        DO 320 J=ITP(8)+1,ITP(9)
        IS=ICT(J)
        DO 322 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DX(I,is)*DZ(I,is)*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,is)*B+DZ(I,is)*A
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B+DX(I,is)*A
322     CONTINUE
320     CONTINUE
C
C       FOR Dyz-TYPE
C
        DO 340 J=ITP(9)+1,ITP(10)
        IS=ICT(J)
        DO 342 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DY(I,is)*DZ(I,is)*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B+DZ(I,is)*A
        CHIZ(I,J)=DZ(I,is)*B+DY(I,is)*A
342     CONTINUE
340     CONTINUE
C
C       FOR Fxxx-TYPE
C
        DO 501 J=ITP(10)+1,ITP(11)
        IS=ICT(J)
        DO 502 I=1,PTS
 
        A=DEXP(-EXX(J)*R2(I,is))
        B=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=B*DX(I,is)
        CHIX(I,J)=(THREE + SUM(J)*DX(I,is)*DX(I,is))*B
        CHIY(I,J)=SUM(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=SUM(J)*DZ(I,is)*CHI(I,J)
 502    CONTINUE
 501    CONTINUE
C
C       FOR Fyyy-TYPE
C
        DO 511 J=ITP(11)+1,ITP(12)
        IS=ICT(J)
        DO 512 I=1,PTS
        A=DEXP(-EXX(J)*R2(I,is))
        B=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=B*DY(I,is)
        CHIX(I,J)=SUM(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(THREE + SUM(J)*DY(I,is)*DY(I,is))*B
        CHIZ(I,J)=SUM(J)*DZ(I,is)*CHI(I,J)
C
 512    CONTINUE
 511    CONTINUE
C
C       FOR Fzzz-TYPE
C
        DO 521 J=ITP(12)+1,ITP(13)
        IS=ICT(J)
        DO 523 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        B=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=B*DZ(I,is)
        CHIX(I,J)=SUM(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=SUM(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(THREE + SUM(J)*DZ(I,is)*DZ(I,is))*B
C
 523    CONTINUE
 521    CONTINUE
C
C       FOR Fxxy-TYPE
C
        DO 531 J=ITP(13)+1,ITP(14)
        IS=ICT(J)
        DO 532 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BXX=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=BXY*DX(I,is)
        CHIX(I,J)=(TWO + SUM(J)*DX(I,is)*DX(I,is))*BXY
        CHIY(I,J)=(ONE + SUM(J)*DY(I,is)*DY(I,is))*BXX
        CHIZ(I,J)=SUM(J)*DZ(I,is)*CHI(I,J)
C
 532    CONTINUE
 531    CONTINUE
C
C       FOR Fxxz-TYPE
C
        DO 541 J=ITP(14)+1,ITP(15)
        IS=ICT(J)
        DO 543 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BXZ=DX(I,is)*DZ(I,is)*A
        BXX=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=BXZ*DX(I,is)
        CHIX(I,J)=(TWO + SUM(J)*DX(I,is)*DX(I,is))*BXZ
        CHIY(I,J)=SUM(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(ONE + SUM(J)*DZ(I,is)*DZ(I,is))*BXX
C
 543    CONTINUE
 541    CONTINUE
C
C       FOR Fyyz-TYPE
C
        DO 561 J=ITP(15)+1,ITP(16)
        IS=ICT(J)
        DO 563 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BYZ=DZ(I,is)*DY(I,is)*A
        BYY=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=BYZ*DY(I,is)
        CHIX(I,J)=SUM(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(TWO + SUM(J)*DY(I,is)*DY(I,is))*BYZ
        CHIZ(I,J)=(ONE + SUM(J)*DZ(I,is)*DZ(I,is))*BYY
C
 563    CONTINUE
 561    CONTINUE
C
C       FOR Fxyy-TYPE
C
        DO 551 J=ITP(16)+1,ITP(17)
        IS=ICT(J)
        DO 552 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BYY=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=BXY*DY(I,is)
        CHIX(I,J)=(ONE + SUM(J)*DX(I,is)*DX(I,is))*BYY
        CHIY(I,J)=(TWO + SUM(J)*DY(I,is)*DY(I,is))*BXY
        CHIZ(I,J)=SUM(J)*DZ(I,is)*CHI(I,J)
C
 552    CONTINUE
 551    CONTINUE
C
C       FOR Fxzz-TYPE
C
        DO 571 J=ITP(17)+1,ITP(18)
        IS=ICT(J)
        DO 572 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BXZ=DZ(I,is)*DX(I,is)*A
        BZZ=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=BXZ*DZ(I,is)
        CHIX(I,J)=(ONE + SUM(J)*DX(I,is)*DX(I,is))*BZZ
        CHIY(I,J)=SUM(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(TWO + SUM(J)*DZ(I,is)*DZ(I,is))*BXZ
C
 572    CONTINUE
 571    CONTINUE
C
C       FOR Fyzz-TYPE
C
        DO 581 J=ITP(18)+1,ITP(19)
        IS=ICT(J)
        DO 583 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BYZ=DZ(I,is)*DY(I,is)*A
        BZZ=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=BYZ*DZ(I,is)
        CHIX(I,J)=SUM(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(ONE + SUM(J)*DY(I,is)*DY(I,is))*BZZ
        CHIZ(I,J)=(TWO + SUM(J)*DZ(I,is)*DZ(I,is))*BYZ
C
 583    CONTINUE
 581    CONTINUE
C
C       FOR Fxyz-TYPE
C
        DO 591 J=ITP(19)+1,ITP(20)
        IS=ICT(J)
        DO 592 I=1,PTS
C
        A=DEXP(-EXX(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BYZ=DZ(I,is)*DY(I,is)*A
        BXZ=DX(I,is)*DZ(I,is)*A
C
        CHI(I,J)=DX(I,is)*BYZ
        CHIX(I,J)=(ONE + SUM(J)*DX(I,is)*DX(I,is))*BYZ
        CHIY(I,J)=(ONE + SUM(J)*DY(I,is)*DY(I,is))*BXZ
        CHIZ(I,J)=(ONE + SUM(J)*DZ(I,is)*DZ(I,is))*BXY
C
 592    CONTINUE
 591    CONTINUE
C
       DO 428 J=1,NPRIMS
       temp=Zero
       DO 429 I=1,PTS
       Check=Dmax1(Dabs(CHI(I,J)),DABS(CHIX(I,J)),DABS(CHIY(I,J)),
     $             Dabs(CHIZ(I,J)),Temp)
       If(Check.GT.Temp)Temp=Check
429    Continue
       Chimax(j)=Temp
428    Continue
C
      DO 430 L=1,NMO
      DO 434 I=1,PTS
         PSI(I,L) = ZERO 
         GX(I,L) = ZERO
         GY(I,L) = ZERO
         GZ(I,L) = ZERO
434   CONTINUE      
      DO 431 J=1,NPRIMS
      Check=Dabs(coo(J,L)*ChiMax(j))
      IF(Check.Gt.Cutoff)THEN
      TEMP=COO(J,L)
      DO 432 I=1,PTS
         PSI(I,L) = PSI(I,L) + TEMP*CHI(I,J)
         GX(I,L) = GX(I,L) + TEMP*CHIX(I,J)
         GY(I,L) = GY(I,L) + TEMP*CHIY(I,J)
         GZ(I,L) = GZ(I,L) + TEMP*CHIZ(I,J)
432      CONTINUE
         ENDIF 
431      CONTINUE
430      CONTINUE
C
        RETURN
        END
      SUBROUTINE GAUSH(XYZ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION NINE      
      PARAMETER(MCENT=100,MMO=100,MPRIMS=420,MPTS=300,NTYPE=20)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),NMO
      COMMON /PRIMS/ COO(MPRIMS,MMO),EXX(MPRIMS),ICT(MPRIMS),
     +       SUM(MPRIMS),Div(Mprims),ITP(NTYPE),NPRIMS
       COMMON /ZZZZZ/ PSI(MMO),GX(MMO),GY(MMO),GZ(MMO),GXX(MMO),
     + GXY(MMO),GXZ(MMO),GYY(MMO),GYZ(MMO),GZZ(MMO)
      DIMENSION XYZ(3),DX(MCENT),DY(MCENT),DZ(MCENT),R2(MCENT),
     $ CHI(MPRIMS),CHIX(MPRIMS),CHIY(MPRIMS),CHIXX(MPRIMS),
     $ CHIXY(MPRIMS),CHIXZ(MPRIMS),CHIYY(MPRIMS),CHIYZ(MPRIMS),
     $ CHIZZ(MPRIMS),CHIZ(MPRIMS)
      Save Zero,One,Two,Four,Five,Seven,Three,Six,Nine
      DATA ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,FOUR/4.0D0/,FIVE/5.0D0/,
     $  SEVEN/7.0D0/,THREE/3.0D0/,SIX/6.0D0/,NINE/9.0D0/
C
      DO 110 J = 1,NCENT
        DX(J) = XYZ(1) - XC(J)
        DY(J) = XYZ(2) - YC(J)
        DZ(J) = XYZ(3) - ZC(J)
        R2(J)= DX(J)*DX(J)+DY(J)*DY(J)+DZ(J)*DZ(J)
110   CONTINUE
C
C       FOR S-TYPE
C
        DO 120 J = 1,ITP(1) 
C
        A=SUM(J)*DEXP(-EXX(J)*R2(ICT(J)))
C
        CHI(J)=A*DIV(J)
        CHIX(J)=DX(ICT(J))*A
        CHIY(J)=DY(ICT(J))*A
        CHIZ(J)=DZ(ICT(J))*A
        CHIXX(J)=(ONE+SUM(J)*DX(ICT(J))*DX(ICT(J)))*A
        CHIYY(J)=(ONE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*A
        CHIZZ(J)=(ONE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*A
        CHIXY(J)=SUM(J)*DX(ICT(J))*DY(ICT(J))*A
        CHIXZ(J)=SUM(J)*DX(ICT(J))*DZ(ICT(J))*A
        CHIYZ(J)=SUM(J)*DY(ICT(J))*DZ(ICT(J))*A
120     CONTINUE
C
C       FOR Px-TYPE
C
        DO 140 J=ITP(1)+1,ITP(2)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        B=DX(ICT(J))*A*SUM(J)
C
        CHI(J)=A*DX(ICT(J))
        CHIX(J)=A+DX(ICT(J))*B
        CHIY(J)=DY(ICT(J))*B
        CHIZ(J)=DZ(ICT(J))*B
        CHIXX(J)=(THREE+SUM(J)*DX(ICT(J))*DX(ICT(J)))*B
        CHIYY(J)=(ONE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*B
        CHIZZ(J)=(ONE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*B
        CHIXY(J)=(A+DX(ICT(J))*B)*DY(ICT(J))*SUM(J)
        CHIXZ(J)=(A+DX(ICT(J))*B)*DZ(ICT(J))*SUM(J)
        CHIYZ(J)=SUM(J)*DY(ICT(J))*DZ(ICT(J))*B
140     CONTINUE
C
C       FOR Py-TYPE
C
        DO 160 J=ITP(2)+1,ITP(3)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        B=DY(ICT(J))*A*SUM(J)
C
        CHI(J)=A*DY(ICT(J))
        CHIX(J)=DX(ICT(J))*B
        CHIY(J)=A+DY(ICT(J))*B
        CHIZ(J)=DZ(ICT(J))*B
        CHIXX(J)=(ONE+SUM(J)*DX(ICT(J))*DX(ICT(J)))*B
        CHIYY(J)=(THREE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*B
        CHIZZ(J)=(ONE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*B
        CHIXY(J)=(A+DY(ICT(J))*B)*DX(ICT(J))*SUM(J)
        CHIXZ(J)=SUM(J)*DX(ICT(J))*DZ(ICT(J))*B
        CHIYZ(J)=(A+DY(ICT(J))*B)*DZ(ICT(J))*SUM(J)
160     CONTINUE
C
C       FOR Pz-TYPE
C
        DO 180 J=ITP(3)+1,ITP(4)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        B=DZ(ICT(J))*A*SUM(J)
C
        CHI(J)=A*DZ(ICT(J))
        CHIX(J)=DX(ICT(J))*B
        CHIY(J)=DY(ICT(J))*B
        CHIZ(J)=A+DZ(ICT(J))*B
        CHIXX(J)=(ONE+SUM(J)*DX(ICT(J))*DX(ICT(J)))*B
        CHIYY(J)=(ONE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*B
        CHIZZ(J)=(THREE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*B
        CHIXY(J)=SUM(J)*DX(ICT(J))*DY(ICT(J))*B
        CHIXZ(J)=(A+DZ(ICT(J))*B)*DX(ICT(J))*SUM(J)
        CHIYZ(J)=(A+DZ(ICT(J))*B)*DY(ICT(J))*SUM(J)
180     CONTINUE
C
C       FOR Dxx-TYPE
C
        DO 220 J=ITP(4)+1,ITP(5)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        B=DX(ICT(J))*DX(ICT(J))*A*SUM(J)
C
        CHI(J)=B*DIV(J)
        CHIX(J)=(TWO*A+B)*DX(ICT(J))
        CHIY(J)=DY(ICT(J))*B
        CHIZ(J)=DZ(ICT(J))*B
        CHIXX(J)=TWO*A+(FIVE+SUM(J)*DX(ICT(J))*DX(ICT(J)))*B
        CHIYY(J)=(ONE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*B
        CHIZZ(J)=(ONE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*B
        CHIXY(J)=(TWO*A+B)*SUM(J)*DX(ICT(J))*DY(ICT(J))
        CHIXZ(J)=(TWO*A+B)*SUM(J)*DX(ICT(J))*DZ(ICT(J))
        CHIYZ(J)=SUM(J)*DY(ICT(J))*DZ(ICT(J))*B
220     CONTINUE
C
C       FOR Dyy-TYPE
C
        DO 240 J=ITP(5)+1,ITP(6)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        B=DY(ICT(J))*DY(ICT(J))*A*SUM(J)
C
        CHI(J)=B*DIV(J)
        CHIX(J)=DX(ICT(J))*B
        CHIY(J)=(TWO*A+B)*DY(ICT(J))
        CHIZ(J)=DZ(ICT(J))*B
        CHIXX(J)=(ONE+SUM(J)*DX(ICT(J))*DX(ICT(J)))*B
        CHIYY(J)=TWO*A+(FIVE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*B
        CHIZZ(J)=(ONE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*B
        CHIXY(J)=(TWO*A+B)*SUM(J)*DX(ICT(J))*DY(ICT(J))
        CHIXZ(J)=SUM(J)*DX(ICT(J))*DZ(ICT(J))*B
        CHIYZ(J)=(TWO*A+B)*SUM(J)*DZ(ICT(J))*DY(ICT(J))
240     CONTINUE
C
C       FOR Dzz-TYPE
C
        DO 260 J=ITP(6)+1,ITP(7)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        B=DZ(ICT(J))*DZ(ICT(J))*A*SUM(J)
C
        CHI(J)=B*DIV(J)
        CHIX(J)=DX(ICT(J))*B
        CHIY(J)=DY(ICT(J))*B
        CHIZ(J)=(TWO*A+B)*DZ(ICT(J))
        CHIXX(J)=(ONE+SUM(J)*DX(ICT(J))*DX(ICT(J)))*B
        CHIYY(J)=(ONE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*B
        CHIZZ(J)=TWO*A+(FIVE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*B
        CHIXY(J)=SUM(J)*DX(ICT(J))*DY(ICT(J))*B
        CHIXZ(J)=(TWO*A+B)*SUM(J)*DX(ICT(J))*DZ(ICT(J))
        CHIYZ(J)=(TWO*A+B)*SUM(J)*DZ(ICT(J))*DY(ICT(J))
260     CONTINUE
C
C       FOR Dxy-TYPE
C
        DO 280 J=ITP(7)+1,ITP(8)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        B=SUM(J)*DX(ICT(J))*DY(ICT(J))*A
C
        CHI(J)=B*DIV(J)
        CHIX(J)=DX(ICT(J))*B+DY(ICT(J))*A
        CHIY(J)=DY(ICT(J))*B+DX(ICT(J))*A
        CHIZ(J)=DZ(ICT(J))*B
        CHIXX(J)=(THREE+SUM(J)*DX(ICT(J))*DX(ICT(J)))*B
        CHIYY(J)=(THREE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*B
        CHIZZ(J)=(ONE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*B
        CHIXY(J)=(ONE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*(A+
     +           A*SUM(J)*DX(ICT(J))*DX(ICT(J))) 
        CHIXZ(J)=(ONE+SUM(J)*DX(ICT(J))*DX(ICT(J)))*A*
     +           DY(ICT(J))*DZ(ICT(J))*SUM(J)
        CHIYZ(J)=(ONE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*A*
     +           DX(ICT(J))*DZ(ICT(J))*SUM(J)
280     CONTINUE
C
C       FOR Dxz-TYPE
C
        DO 320 J=ITP(8)+1,ITP(9)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        B=DX(ICT(J))*DZ(ICT(J))*A*SUM(J)
C
        CHI(J)=B*DIV(J)
        CHIX(J)=DX(ICT(J))*B+DZ(ICT(J))*A
        CHIY(J)=DY(ICT(J))*B
        CHIZ(J)=DZ(ICT(J))*B+DX(ICT(J))*A
        CHIXX(J)=(THREE+SUM(J)*DX(ICT(J))*DX(ICT(J)))*B
        CHIYY(J)=(ONE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*B
        CHIZZ(J)=(THREE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*B
        CHIXZ(J)=(ONE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*(A+
     +           A*SUM(J)*DX(ICT(J))*DX(ICT(J))) 
        CHIXY(J)=(ONE+SUM(J)*DX(ICT(J))*DX(ICT(J)))*A*
     +           DZ(ICT(J))*DY(ICT(J))*SUM(J)
        CHIYZ(J)=(ONE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*A*
     +           DX(ICT(J))*DY(ICT(J))*SUM(J)
320     CONTINUE
C
C       FOR Dyz-TYPE
C
        DO 340 J=ITP(9)+1,ITP(10)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        B=DY(ICT(J))*DZ(ICT(J))*A*SUM(J)
C
        CHI(J)=B*DIV(J)
        CHIX(J)=DX(ICT(J))*B
        CHIY(J)=DY(ICT(J))*B+DZ(ICT(J))*A
        CHIZ(J)=DZ(ICT(J))*B+DY(ICT(J))*A
        CHIXX(J)=(ONE+SUM(J)*DX(ICT(J))*DX(ICT(J)))*B
        CHIYY(J)=(THREE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*B
        CHIZZ(J)=(THREE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*B
        CHIYZ(J)=(ONE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*(A+
     +           A*SUM(J)*DY(ICT(J))*DY(ICT(J))) 
        CHIXY(J)=(ONE+SUM(J)*DY(ICT(J))*DY(ICT(J)))*A*
     +           DZ(ICT(J))*DX(ICT(J))*SUM(J)
        CHIXZ(J)=(ONE+SUM(J)*DZ(ICT(J))*DZ(ICT(J)))*A*
     +           DX(ICT(J))*DY(ICT(J))*SUM(J)
340     CONTINUE
C
C       FOR Fxxx-TYPE
C
        DO 501 J=ITP(10)+1,ITP(11)
 
        A=DEXP(-EXX(J)*R2(ICT(J)))
        XSQ=DX(ICT(J))*DX(ICT(J))
        YSQ=DY(ICT(J))*DY(ICT(J))
        ZSQ=DZ(ICT(J))*DZ(ICT(J))
        B=XSQ*A
C
        CHI(J)=B*DX(ICT(J))
        CHIX(J)=(THREE + SUM(J)*XSQ)*B
        CHIY(J)=SUM(J)*DY(ICT(J))*CHI(J)
        CHIZ(J)=SUM(J)*DZ(ICT(J))*CHI(J)
        CHIXX(J)=SIX*DX(ICT(J))*A + SUM(J)*(SEVEN+SUM(J)*XSQ)*CHI(J)
        CHIYY(J)=SUM(J)*CHI(J)*(ONE + SUM(J)*YSQ)
        CHIZZ(J)=SUM(J)*CHI(J)*(ONE + SUM(J)*ZSQ)
        CHIXY(J)=SUM(J)*DY(ICT(J))*(THREE+SUM(J)*XSQ)*B
        CHIXZ(J)=SUM(J)*DZ(ICT(J))*(THREE+SUM(J)*XSQ)*B
        CHIYZ(J)=SUM(J)*DY(ICT(J))*CHIZ(J)
C
 501    CONTINUE
C
C       FOR Fyyy-TYPE
C
        DO 511 J=ITP(11)+1,ITP(12)
        A=DEXP(-EXX(J)*R2(ICT(J)))
        XSQ=DX(ICT(J))*DX(ICT(J))
        YSQ=DY(ICT(J))*DY(ICT(J))
        ZSQ=DZ(ICT(J))*DZ(ICT(J))
        B=YSQ*A
C
        CHI(J)=B*DY(ICT(J))
        CHIX(J)=SUM(J)*DX(ICT(J))*CHI(J)
        CHIY(J)=(THREE + SUM(J)*YSQ)*B
        CHIZ(J)=SUM(J)*DZ(ICT(J))*CHI(J)
        CHIXX(J)=SUM(J)*CHI(J)*(ONE + SUM(J)*XSQ)
        CHIYY(J)=SIX*DY(ICT(J))*A + SUM(J)*(SEVEN+SUM(J)*YSQ)*CHI(J)
        CHIZZ(J)=SUM(J)*CHI(J)*(ONE + SUM(J)*ZSQ)
        CHIXY(J)=SUM(J)*DX(ICT(J))*(THREE+SUM(J)*YSQ)*B
        CHIXZ(J)=SUM(J)*DX(ICT(J))*CHIZ(J)
        CHIYZ(J)=SUM(J)*DZ(ICT(J))*(THREE+SUM(J)*YSQ)*B
C
 511    CONTINUE
C
C       FOR Fzzz-TYPE
C
        DO 521 J=ITP(12)+1,ITP(13)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        XSQ=DX(ICT(J))*DX(ICT(J))
        YSQ=DY(ICT(J))*DY(ICT(J))
        ZSQ=DZ(ICT(J))*DZ(ICT(J))
        B=ZSQ*A
C
        CHI(J)=B*DZ(ICT(J))
        CHIX(J)=SUM(J)*DX(ICT(J))*CHI(J)
        CHIY(J)=SUM(J)*DY(ICT(J))*CHI(J)
        CHIZ(J)=(THREE + SUM(J)*ZSQ)*B
        CHIXX(J)=SUM(J)*CHI(J)*(ONE + SUM(J)*XSQ)
        CHIYY(J)=SUM(J)*CHI(J)*(ONE + SUM(J)*YSQ)
        CHIZZ(J)=SIX*DZ(ICT(J))*A + SUM(J)*(SEVEN+SUM(J)*ZSQ)*CHI(J)
        CHIXY(J)=SUM(J)*DX(ICT(J))*CHIY(J)
        CHIXZ(J)=SUM(J)*DX(ICT(J))*(THREE+SUM(J)*ZSQ)*B
        CHIYZ(J)=SUM(J)*DY(ICT(J))*(THREE+SUM(J)*ZSQ)*B
C
 521    CONTINUE
C
C       FOR Fxxy-TYPE
C
        DO 531 J=ITP(13)+1,ITP(14)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        XY=DX(ICT(J))*DY(ICT(J))
        XZ=DX(ICT(J))*DZ(ICT(J))
        XSQ=DX(ICT(J))*DX(ICT(J))
        YSQ=DY(ICT(J))*DY(ICT(J))
        ZSQ=DZ(ICT(J))*DZ(ICT(J))
C
        CHI(J)=XY*DX(ICT(J))*A
        CHIX(J)=(TWO + SUM(J)*XSQ)*XY*A
        CHIY(J)=(ONE + SUM(J)*YSQ)*XSQ*A
        CHIZ(J)=SUM(J)*DZ(ICT(J))*CHI(J)
        CHIXX(J)=DY(ICT(J))*A*(TWO + SUM(J)*XSQ*(FIVE+XSQ*SUM(J)))
        CHIYY(J)=SUM(J)*(THREE + SUM(J)*YSQ)*CHI(J)
        CHIZZ(J)=SUM(J)*(ONE + SUM(J)*ZSQ)*CHI(J)
        CHIXY(J)=DX(ICT(J))*(TWO+SUM(J)*(XSQ+YSQ*(TWO+SUM(J)*XSQ)))*A
        CHIXZ(J)=SUM(J)*XZ*(TWO*DY(ICT(J))*A + SUM(J)*CHI(J))
        CHIYZ(J)=SUM(J)*DZ(ICT(J))*XSQ*(ONE+ SUM(J)*YSQ)*A
C
 531    CONTINUE
C
C       FOR Fxxz-TYPE
C
        DO 541 J=ITP(14)+1,ITP(15)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        XY=DX(ICT(J))*DY(ICT(J))
        XZ=DX(ICT(J))*DZ(ICT(J))
        XSQ=DX(ICT(J))*DX(ICT(J))
        YSQ=DY(ICT(J))*DY(ICT(J))
        ZSQ=DZ(ICT(J))*DZ(ICT(J))
C
        CHI(J)=XZ*DX(ICT(J))*A
        CHIX(J)=(TWO + SUM(J)*XSQ)*XZ*A
        CHIY(J)=SUM(J)*DY(ICT(J))*CHI(J)
        CHIZ(J)=(ONE + SUM(J)*ZSQ)*XSQ*A
        CHIXX(J)=DZ(ICT(J))*A*(TWO + SUM(J)*XSQ*(FIVE+XSQ*SUM(J)))
        CHIYY(J)=SUM(J)*(ONE + SUM(J)*YSQ)*CHI(J)
        CHIZZ(J)=SUM(J)*(THREE + SUM(J)*ZSQ)*CHI(J)
        CHIXY(J)=SUM(J)*XY*(TWO*DZ(ICT(J))*A + SUM(J)*CHI(J))
        CHIXZ(J)=DX(ICT(J))*(TWO+SUM(J)*(XSQ+ZSQ*(TWO+SUM(J)*XSQ)))*A
        CHIYZ(J)=SUM(J)*DY(ICT(J))*XSQ*(ONE+ SUM(J)*ZSQ)*A
C
 541    CONTINUE
C
C       FOR Fyyz-TYPE
C
        DO 561 J=ITP(15)+1,ITP(16)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        XY=DX(ICT(J))*DY(ICT(J))
        YZ=DZ(ICT(J))*DY(ICT(J))
        XSQ=DX(ICT(J))*DX(ICT(J))
        YSQ=DY(ICT(J))*DY(ICT(J))
        ZSQ=DZ(ICT(J))*DZ(ICT(J))
C
        CHI(J)=YZ*DY(ICT(J))*A
        CHIX(J)=SUM(J)*DX(ICT(J))*CHI(J)
        CHIY(J)=(TWO + SUM(J)*YSQ)*YZ*A
        CHIZ(J)=(ONE + SUM(J)*ZSQ)*YSQ*A
        CHIXX(J)=SUM(J)*(ONE + SUM(J)*XSQ)*CHI(J)
        CHIYY(J)=DZ(ICT(J))*A*(TWO + SUM(J)*YSQ*(FIVE+YSQ*SUM(J)))
        CHIZZ(J)=SUM(J)*(THREE + SUM(J)*ZSQ)*CHI(J)
        CHIXY(J)=SUM(J)*XY*(TWO*DZ(ICT(J))*A + SUM(J)*CHI(J))
        CHIXZ(J)=SUM(J)*DX(ICT(J))*YSQ*(ONE+ SUM(J)*ZSQ)*A
        CHIYZ(J)=DY(ICT(J))*(TWO+SUM(J)*(YSQ+ZSQ*(TWO+SUM(J)*YSQ)))*A
C
 561    CONTINUE
C
C       FOR Fxyy-TYPE
C
        DO 551 J=ITP(16)+1,ITP(17)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        XY=DX(ICT(J))*DY(ICT(J))
        YZ=DY(ICT(J))*DZ(ICT(J))
        XSQ=DX(ICT(J))*DX(ICT(J))
        YSQ=DY(ICT(J))*DY(ICT(J))
        ZSQ=DZ(ICT(J))*DZ(ICT(J))
C
        CHI(J)=XY*DY(ICT(J))*A
        CHIX(J)=(ONE + SUM(J)*XSQ)*YSQ*A
        CHIY(J)=(TWO + SUM(J)*YSQ)*XY*A
        CHIZ(J)=SUM(J)*DZ(ICT(J))*CHI(J)
        CHIXX(J)=SUM(J)*(THREE + SUM(J)*XSQ)*CHI(J)
        CHIYY(J)=DX(ICT(J))*A*(TWO + SUM(J)*YSQ*(FIVE+YSQ*SUM(J)))
        CHIZZ(J)=SUM(J)*(ONE + SUM(J)*ZSQ)*CHI(J)
        CHIXY(J)=DY(ICT(J))*(TWO+SUM(J)*(YSQ+XSQ*(TWO+SUM(J)*YSQ)))*A
        CHIYZ(J)=SUM(J)*YZ*(TWO*DX(ICT(J))*A + SUM(J)*CHI(J))
        CHIXZ(J)=SUM(J)*DZ(ICT(J))*YSQ*(ONE+ SUM(J)*XSQ)*A
C
 551    CONTINUE
C
C       FOR Fxzz-TYPE
C
        DO 571 J=ITP(17)+1,ITP(18)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        XZ=DZ(ICT(J))*DX(ICT(J))
        YZ=DZ(ICT(J))*DY(ICT(J))
        XSQ=DX(ICT(J))*DX(ICT(J))
        YSQ=DY(ICT(J))*DY(ICT(J))
        ZSQ=DZ(ICT(J))*DZ(ICT(J))
C
        CHI(J)=XZ*DZ(ICT(J))*A
        CHIX(J)=(ONE + SUM(J)*XSQ)*ZSQ*A
        CHIY(J)=SUM(J)*DY(ICT(J))*CHI(J)
        CHIZ(J)=(TWO + SUM(J)*ZSQ)*XZ*A
        CHIXX(J)=SUM(J)*(THREE + SUM(J)*XSQ)*CHI(J)
        CHIYY(J)=SUM(J)*(ONE + SUM(J)*YSQ)*CHI(J)
        CHIZZ(J)=DX(ICT(J))*A*(TWO + SUM(J)*ZSQ*(FIVE+ZSQ*SUM(J)))
        CHIXY(J)=SUM(J)*DY(ICT(J))*ZSQ*(ONE+ SUM(J)*XSQ)*A
        CHIXZ(J)=DZ(ICT(J))*(TWO+SUM(J)*(ZSQ+XSQ*(TWO+SUM(J)*ZSQ)))*A
        CHIYZ(J)=SUM(J)*YZ*(TWO*DX(ICT(J))*A + SUM(J)*CHI(J))
C
 571    CONTINUE
C
C       FOR Fyzz-TYPE
C
        DO 581 J=ITP(18)+1,ITP(19)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        XZ=DX(ICT(J))*DZ(ICT(J))
        YZ=DZ(ICT(J))*DY(ICT(J))
        XSQ=DX(ICT(J))*DX(ICT(J))
        YSQ=DY(ICT(J))*DY(ICT(J))
        ZSQ=DZ(ICT(J))*DZ(ICT(J))
C
        CHI(J)=YZ*DZ(ICT(J))*A
        CHIX(J)=SUM(J)*DX(ICT(J))*CHI(J)
        CHIY(J)=(ONE + SUM(J)*YSQ)*ZSQ*A
        CHIZ(J)=(TWO + SUM(J)*ZSQ)*YZ*A
        CHIXX(J)=SUM(J)*(ONE + SUM(J)*XSQ)*CHI(J)
        CHIYY(J)=SUM(J)*(THREE + SUM(J)*YSQ)*CHI(J)
        CHIZZ(J)=DY(ICT(J))*A*(TWO + SUM(J)*ZSQ*(FIVE+ZSQ*SUM(J)))
        CHIXY(J)=SUM(J)*DX(ICT(J))*ZSQ*(ONE+ SUM(J)*YSQ)*A
        CHIYZ(J)=DZ(ICT(J))*(TWO+SUM(J)*(ZSQ+YSQ*(TWO+SUM(J)*ZSQ)))*A
        CHIXZ(J)=SUM(J)*XZ*(TWO*DY(ICT(J))*A + SUM(J)*CHI(J))
C
 581    CONTINUE
C
C       FOR Fxyz-TYPE
C
        DO 591 J=ITP(19)+1,ITP(20)
C
        A=DEXP(-EXX(J)*R2(ICT(J)))
        XY=DX(ICT(J))*DY(ICT(J))
        YZ=DZ(ICT(J))*DY(ICT(J))
        XZ=DX(ICT(J))*DZ(ICT(J))
        XSQ=DX(ICT(J))*DX(ICT(J))
        YSQ=DY(ICT(J))*DY(ICT(J))
        XSQ=DZ(ICT(J))*DZ(ICT(J))
C
        CHI(J)=DX(ICT(J))*YZ*A
        CHIX(J)=(ONE + SUM(J)*XSQ)*YZ*A
        CHIY(J)=(ONE + SUM(J)*YSQ)*XZ*A
        CHIZ(J)=(ONE + SUM(J)*ZSQ)*XY*A
        CHIXX(J)=(THREE + XSQ*SUM(J))*SUM(J)*CHI(J)
        CHIYY(J)=(THREE + YSQ*SUM(J))*SUM(J)*CHI(J)
        CHIZZ(J)=(THREE + ZSQ*SUM(J))*SUM(J)*CHI(J)
        CHIXY(J)=DZ(ICT(J))*(ONE+SUM(J)*(XSQ+YSQ+SUM(J)*XY*XY))*A
        CHIXZ(J)=DY(ICT(J))*(ONE+SUM(J)*(XSQ+ZSQ+SUM(J)*XZ*XZ))*A
        CHIYZ(J)=DX(ICT(J))*(ONE+SUM(J)*(YSQ+ZSQ+SUM(J)*YZ*YZ))*A
C
 591    CONTINUE
C
C
C     CALCULATE THE VALUES OF THE MO's AND THEIR FIRST AND SECOND
C     DERIVATIVE COMPONENTS AT THE PTS (XYZ)'s 
C
        DO 128 L = 1,NMO
          PSI(L) =ZERO
          GX(L) =ZERO
          GY(L) =ZERO
          GZ(L) =ZERO
          GXX(L) =ZERO
          GXY(L) =ZERO
          GXZ(L) =ZERO
          GYY(L) = ZERO
          GYZ(L) = ZERO
          GZZ(L) = ZERO
128      CONTINUE
C
        DO 124 L = 1,NMO
        DO 125 J = 1,NPRIMS
          PSI(L) = PSI(L) + COO(J,L)*CHI(J)
          GX(L) = GX(L) + COO(J,L)*CHIX(J)
          GY(L) = GY(L) + COO(J,L)*CHIY(J)
          GZ(L) = GZ(L) + COO(J,L)*CHIZ(J)
          GXX(L) = GXX(L) + COO(J,L)*CHIXX(J)
          GXY(L) = GXY(L) + COO(J,L)*CHIXY(J)
          GXZ(L) = GXZ(L) + COO(J,L)*CHIXZ(J)
          GYY(L) = GYY(L) + COO(J,L)*CHIYY(J)
          GYZ(L) = GYZ(L) + COO(J,L)*CHIYZ(J)
          GZZ(L) = GZZ(L) + COO(J,L)*CHIZZ(J)
125      CONTINUE
124      CONTINUE
C
        RETURN
        END
      SUBROUTINE GROCKLE (N, X, IR, S, E)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION X(3,N),S(3,N),C(3),E(3,3),EV(3),R(3,3),IR(N)
      Save Zero,One
      Data ZERO/0.0d0/,One/1.0d0/
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
      M = 0
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
      DD = One/DFLOAT(M)
      C(1) = DD*C(1)
      C(2) = DD*C(2)
      C(3) = DD*C(3)
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
      CALL TRACE(E, EV, C, 3, IFAIL)
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MCENT=100, MMO=100, MPRIMS=420, MPTS=300)
      COMMON /ORBTL/ EORB(MMO),PO(MMO),NMO
      COMMON /ZZZZZ/ PSI(MMO),GX(MMO),GY(MMO),GZ(MMO),GXX(MMO),
     + GXY(MMO),GXZ(MMO),GYY(MMO),GYZ(MMO),GZZ(MMO)
      DIMENSION W(3),F(3),H(3,3),XYZ(3)
      Save Zero,Two
      Data Zero/0.0d0/,Two/2.0d0/
C
      XYZ(1)= X 
      XYZ(2)= Y 
      XYZ(3)= Z 
C
      CALL GAUSH(XYZ) 
C
      H(1,1) = Zero
      H(1,2) = Zero
      H(1,3) = Zero
      H(2,2) = Zero
      H(2,3) = Zero
      H(3,3) = Zero
C
      DO 110 I = 1,NMO
        H(1,1) = H(1,1)+Two*PO(I)*(PSI(I)*GXX(I)+GX(I)*GX(I)) 
        H(1,2) = H(1,2)+Two*PO(I)*(PSI(I)*GXY(I)+GX(I)*GY(I)) 
        H(1,3) = H(1,3)+Two*PO(I)*(PSI(I)*GXZ(I)+GX(I)*GZ(I)) 
        H(2,2) = H(2,2)+Two*PO(I)*(PSI(I)*GYY(I)+GY(I)*GY(I))
        H(2,3) = H(2,3)+Two*PO(I)*(PSI(I)*GYZ(I)+GY(I)*GZ(I)) 
        H(3,3) = H(3,3)+Two*PO(I)*(PSI(I)*GZZ(I)+GZ(I)*GZ(I))
110   CONTINUE
C
      DO 120 I = 1,3
        DO 120 J = 1,I
          H(I,J) = H(J,I)
120   CONTINUE
C
      CALL TRACE(H,W,F,3,IERR)
C
      RETURN
      END

      SUBROUTINE MAKNAME (I,STRING,L,EXT)
C
      CHARACTER*(*) STRING,EXT
      INTEGER I,J,L
C
      CALL GETARG (I,STRING)
C
      J=LEN(STRING)
C
      DO 10 N=1,J
         IF (STRING(N:N) .EQ. ' ') THEN
            L=N-1
            STRING=STRING(1:L)//EXT
            RETURN
         ENDIF
10    CONTINUE
      STOP 'Failed to make a filename.'
      RETURN 
      END 
      SUBROUTINE NUCLEI(XN,YN,ZN,CR,CT,SCAL)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION CR(2),CT(3)
      Save Sn1,SN2,Pt1
      DATA SN1,SN2,Pt1 /.06d0,.04d0,0.1d0/
C
      XN = (XN - CT(1))*SCAL+CR(1)
      YN = (YN - CT(2))*SCAL+CR(2)
C
      CALL PLOT (SNGL(XN-SCAL*SN1),SNGL(YN),3)
C
      IF (DABS(ZN) .GE. Pt1) THEN
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
      SUBROUTINE RDPSI
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*80 WFNTTL,JOBTTL,LINE
      CHARACTER*8 ATNAM
      PARAMETER(MCENT=100,MMO=100,MPRIMS=420,MPTS=300,NTYPE=20)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),NMO
      COMMON /PRIMS/ COO(MPRIMS,MMO),EXX(MPRIMS),ICT(MPRIMS),
     +       SUM(MPRIMS),Div(Mprims),ITP(NTYPE),NPRIMS
      COMMON /STRING/ WFNTTL,JOBTTL,ATNAM(MCENT),NAT
      COMMON /UNITS/  IVEC,IOUT,IWFN
      COMMON /VALUES/ THRESH1,THRESH2,GAMMA,TOTE
      DIMENSION ITYPE(MPRIMS),CO(MPRIMS,MMO),ICENT(MPRIMS),EX(MPRIMS)
      Save One,Two
      Data One/1.0d0/,Two/2.0d0/
101   FORMAT (A80)
102   FORMAT (4X,A4,12X,3(I3,17X))
103   FORMAT (A8,11X,I3,2X,3F12.8,10X,F5.1)
104   FORMAT (20X,20I3)
105   FORMAT (10X,5E14.7)
106   FORMAT (35X,F12.8,15X,F12.8)
107   FORMAT (5E16.8)
109   FORMAT (17X,F20.12,18X,F13.8)
C
      READ(IWFN,101) WFNTTL
C
      READ(IWFN,102) MODE,NMO,NPRIMS,NCENT
C
      IF (NMO .GT. MMO) STOP 'Too many MOs'
      IF (NPRIMS .GT. MPRIMS) STOP 'Too many primitives'
      IF (NCENT .GT. MCENT) STOP 'Too many nuclei'
C
      DO 100 I = 1,NCENT
      READ(IWFN,103) ATNAM(I),J,XC(J),YC(J),ZC(J),CHARG(J)
100   CONTINUE
C
      READ(IWFN,104) (ICENT(I),I=1,NPRIMS)
      READ(IWFN,104) (ITYPE(I),I=1,NPRIMS)
      READ(IWFN,105) (EX(I),I=1,NPRIMS)
C
      DO 110 I = 1,NMO
        READ(IWFN,106) PO(I),EORB(I)
        READ(IWFN,107) (CO(J,I),J=1,NPRIMS)
110   CONTINUE
C
      READ(IWFN,101) LINE
C
      READ(IWFN,109) TOTE,GAMMA
C
      N=0
      DO 160 K=1,NTYPE
      DO 170 J=1,NPRIMS
      IF(ITYPE(J).EQ.K) THEN
      N=N+1
      EXX(N)=EX(J)
      ICT(N)=ICENT(J)
      SUM(N)=-Two*EX(J)
      DIV(N)=One/SUM(N)
      DO 180 L=1,NMO
      COO(N,L)=CO(J,L)
180   CONTINUE
      ENDIF
170   CONTINUE
      ITP(K)=N
160   CONTINUE 
C
      RETURN
      END
      SUBROUTINE SHFTVC(A,F,H,IUD,SV)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(3,3),F(3),H(3,3),SV(3),X(3)
      Save Zero,One
      Data Zero/0.0d0/,One/1.0d0/
C
      ADOT = One
      XDOT = One
C
      X(1) = One*A(1,3)
      X(2) = One*A(2,3)
      X(3) = One*A(3,3)
C
      DO 100 K = 1,3
        IF (F(K) .GT. Zero) 
     +  XDOT = H(K,1)*X(1) + H(K,2)*X(2) + H(K,3)*X(3)
        IF (XDOT .LT. ADOT) THEN 
          MDOT = K
          ADOT = XDOT
        END IF
100   CONTINUE
C
      IF (IUD .EQ. 1) THEN
C
        SV(1) = H(1,MDOT)
        SV(2) = H(2,MDOT)
        SV(3) = H(3,MDOT)
      ELSE IF (IUD .EQ. 0) THEN
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

		DO 32 M = L,N
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
        DIMENSION       H(N,N), E(*), W(*)
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
      SUBROUTINE TRUDGE(XYZ0,NPTH,IUP,A,CT,CR,ENDPT,SCAL,IINC,IPROJ,CX)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MCENT=100,MMO=100,MPRIMS=420,MPTS=300,LLL=2000)
      COMMON /ORBTL/ EORB(MMO),PO(MMO),NMO
      COMMON /ZZZZ/ PSI(MPTS,MMO),GX(MPTS,MMO),GY(MPTS,MMO),
     + GZ(MPTS,MMO)
      DIMENSION XYZ(MPTS,3),XX(LLL,MPTS),YY(LLL,MPTS),DX(MPTS),
     $ DY(MPTS),DZ(MPTS),DSQ(MPTS),X(MPTS),Y(MPTS),XYZ0(3,MPTS),
     $ L(MPTS),CX(3),M(MPTS),CT(3),CR(2),A(3,3),DX0(MPTS),
     $ DY0(MPTS),DZ0(MPTS),MX(MPTS)
      Save Zero,Step,StepInc,Five,One,Two
      DATA Zero/0.0d0/,STEP/1.0D-3/,STEPINC/1.0D-5/,Five/5.0d0/,
     $ One/1.0d0/,Two/2.0d0/
C
      NPTH0=NPTH
C
      LC=0
      K = 0
      DO 25 I = 1,NPTH
        M(I)=1
        L(I) = 0
        DX0(I) = Zero
        DY0(I) = Zero
        DZ0(I) = Zero
        MX(I)=I
25    CONTINUE
C
      IF (IUP .EQ. 0) DIR =  One
      IF (IUP .EQ. 1) DIR = -One
C
      DO 27 I=1,NPTH
      XYZ(I,1)=A(1,1)*XYZ0(1,I)+A(1,2)*XYZ0(2,I)+A(1,3)*XYZ0(3,I)+CX(1)
      XYZ(I,2)=A(2,1)*XYZ0(1,I)+A(2,2)*XYZ0(2,I)+A(2,3)*XYZ0(3,I)+CX(2) 
      XYZ(I,3)=A(3,1)*XYZ0(1,I)+A(3,2)*XYZ0(2,I)+A(3,3)*XYZ0(3,I)+CX(3) 
      XX(1,I) = XYZ0(1,I)*SCAL + CR(1)
      YY(1,I) = XYZ0(2,I)*SCAL + CR(2)
27    CONTINUE
C
30    DO 28 I=1,NPTH
      M(I) = M(I) + 1
      DX(I) = Zero
      DY(I) = Zero
      DZ(I) = Zero
      L(I) = 0
28    CONTINUE
C
      CALL GAUS(XYZ,NPTH) 
C
        DO 120 I = 1,NMO
        DO 122 J = 1,NPTH
        DX(J) = DX(J) - DIR*GX(J,I)*PSI(J,I)*PO(I)*Two
        DY(J) = DY(J) - DIR*GY(J,I)*PSI(J,I)*PO(I)*Two
        DZ(J) = DZ(J) - DIR*GZ(J,I)*PSI(J,I)*PO(I)*Two
122   CONTINUE
120   CONTINUE
C
        DO 127 I = 1,NPTH
        DSQ(I) = DSQRT(DX(I)**2 + DY(I)**2 + DZ(I)**2)
        IF (DSQ(I) .LT. ENDPT) GOTO 189
        IF(IPROJ.NE.1)THEN
        IF ((XX(M(I)-1,MX(I)) .GT. (FIVE+CR(1)))
     +.OR.(XX(M(I)-1,MX(I)).LT.(CR(1)-FIVE)).OR.(YY(M(I)-1,MX(I)).GT.
     + (Five+CR(2))).OR.(YY(M(I)-1,MX(I)) .LT. (CR(2)-Five))) GOTO 189 
       ENDIF
       APAR=(DX0(I)*DX(I) + DY0(I)*DY(I) + DZ0(I)*DZ(I))/Dsq(i)
       IF(APAR.LT.ZERO) GOTO 189
       IF (M(I)-1 .EQ. LLL-1) GOTO 190
 127    CONTINUE 
C
      NLP=0
      DO 500 J=1,NPTH
      IF(L(J).EQ.0) THEN
      NLP=NLP+1
      XYZ(NLP,1)=XYZ(J,1) 
      XYZ(NLP,2)=XYZ(J,2) 
      XYZ(NLP,3)=XYZ(J,3) 
      DSQ(NLP)=DSQ(J)
      DX(NLP)=DX(J)
      DY(NLP)=DY(J)
      DZ(NLP)=DZ(J)
      M(NLP)=M(J)
      MX(NLP)=MX(J)
      ENDIF
500   CONTINUE
C
      NPTH=NPTH-LC
      LC=0
C
      IF(IPROJ .EQ. 1) THEN
      DO 128 I = 1,NPTH
      Dmult2=One/Dsq(i)
      DMult1=Dmult2*Step
      XYZ(I,1) = XYZ(I,1) + DMult1*DX(I) 
      XYZ(I,2) = XYZ(I,2) + DMult1*DY(I) 
      XYZ(I,3) = XYZ(I,3) + DMult1*DZ(I) 
      DX0(I) = DX(I)*DMult2
      DY0(I) = DY(I)*DMult2
      DZ0(I) = DZ(I)*DMult2
128   CONTINUE
      DO 129 I=1,NPTH 
      XX(M(I),MX(I))=A(1,1)*XYZ(I,1)+A(2,1)*XYZ(I,2)+A(3,1)*XYZ(I,3)-
     + CT(1)*SCAL+CR(1)
      YY(M(I),MX(I))=A(1,2)*XYZ(I,1)+A(2,2)*XYZ(I,2)+A(3,2)*XYZ(I,3)-
     + CT(2)*SCAL+CR(2)
129    CONTINUE
C
       ELSE 
       DO 555 I=1,NPTH
       DMult2=One/DSQ(i)
       DMult1=DMult2*Step
      XYZ(I,1) = XYZ(I,1) + DMult1*DX(I) 
      XYZ(I,2) = XYZ(I,2) + DMult1*DY(I) 
      XYZ(I,3) = XYZ(I,3) + DMult1*DZ(I) 
      DX0(I) = DX(I)*DMult2
      DY0(I) = DY(I)*DMult2
      DZ0(I) = DZ(I)*DMult2
555   CONTINUE
      DO 556 I=1,NPTH
      X(I) = A(1,1)*XYZ(I,1) + A(2,1)*XYZ(I,2) + A(3,1)*XYZ(I,3)-CT(1)
      Y(I) = A(1,2)*XYZ(I,1) + A(2,2)*XYZ(I,2) + A(3,2)*XYZ(I,3)-CT(2)
556   CONTINUE
      DO 557 I=1,NPTH
      XYZ(I,1) = A(1,1)*X(I) + A(1,2)*Y(I)+CX(1)
      XYZ(I,2) = A(2,1)*X(I) + A(2,2)*Y(I)+CX(2)
      XYZ(I,3) = A(3,1)*X(I) + A(3,2)*Y(I)+CX(3)
557   CONTINUE
      DO 558 I=1,NPTH
      XX(M(I),MX(I))=X(I)*SCAL+CR(1)
      YY(M(I),MX(I))=Y(I)*SCAL+CR(2)
558   CONTINUE
      ENDIF
C
      IF ((IUP .EQ. 0) .AND. (IINC .EQ. 1)) STEP = STEP + STEPINC
C
      GOTO 30
C
189     K=K+1
        LC=LC+1
        L(I) = 1
        NM=M(I)-1
        GOTO 191
190     NM=M(I)-1
        M(I)=1
191     CALL PLOT(SNGL(XX(1,MX(I))),SNGL(YY(1,MX(I))),3)
        DO 192 J = 1,NM
        CALL PLOT(SNGL(XX(J,MX(I))),SNGL(YY(J,MX(I))),2) 
192     CONTINUE
        IF (K .EQ. NPTH0) GOTO 200
        GOTO 127
C
200   RETURN
      END
