
      PROGRAM ENVELOPE
C
C    ROUTINE TO CREATE CONSTANT-VALUE ENVELOPES FROM THE
C    CUBES OF DATA PRODUCED BY THE ROUTINE CUBE
C
      CHARACTER*80 TITLE
      CHARACTER*40 WQUB,WENV
      CHARACTER*4  FQUB,FENV
      DIMENSION T(50,50,50),EYE(3),SLAB(52,52)
      DATA IFLAG /7/, FQUB /'.qub'/, FENV /'.env'/
C
      CALL PLOTS (53,0,2)
      CALL PLOT (0.0,0.0,3)
C
C   FOR UNIX IMPLEMENTATION
C
C     CALL MAKNAME(1,WQUB,ILEN,FQUB)
C     CALL MAKNAME(1,WENV,ILEN,FENV)
C     IF (ILEN .EQ. 0) STOP ' usage: envelope qubfile '
C
      OPEN (30,FILE=WQUB)
      OPEN (2,FILE=WENV)
C
      WRITE (6,100)
100   FORMAT(' EYE POSITION ',$)
      READ (5,*) (EYE(I),I=1,3)
      WRITE (6,110) 
110   FORMAT(' OUTER CONTOUR VALUE ',$)
      READ (5,*) TISO
      WRITE (6,130)
130   FORMAT(' LARGER VALUES INSIDE OR OUTSIDE ENVELOPE ',$)
      READ (5,*) INO
C
      IFLAG = IFLAG*INO
C
      OPEN (30)
      READ (30,*) IX,IY,IZ
      DO 1000 I = 1,IX
        DO 1000 J = 1,IY
          READ (30,*) (T(I,J,K),K=1,IZ)
1000  CONTINUE
C
      CALL ISOSRF(T,50,IX,50,IY,IZ,EYE,52,SLAB,TISO,IFLAG)
C
      CALL PLOT(0.,0.,999)
      END      SUBROUTINE DRAWI (IXA,IYA,IXB,IYB)                                ISO02076
C                                                                       ISO02077
C INCLUDED FOR USE BY PWRZ                                              ISO02078
C                                                                       ISO02079
      CALL FRSTC (IXA,IYA,1)                                            ISO02080
      CALL FRSTC (IXB,IYB,2)                                            ISO02081
      RETURN                                                            ISO02082
      END                                                               ISO02083
      SUBROUTINE DRCNTR (Z,L,MM,NN)                                     ISO00933
C                                                                       ISO00934
      DIMENSION       Z(L,NN)                                           ISO00935
C                                                                       ISO00936
C THIS ROUTINE TRACES A CONTOUR LINE WHEN GIVEN THE BEGINNING BY STLINE.ISO00937
C TRANSFORMATIONS CAN BE ADDED BY DELETING THE STATEMENT FUNCTIONS FOR  ISO00938
C FX AND FY IN DRLINE AND MINMAX AND ADDING EXTERNAL FUNCTIONS.         ISO00939
C X=1. AT Z(1,J), X=FLOAT(M) AT Z(M,J). X TAKES ON NON-INTEGER VALUES.  ISO00940
C Y=1. AT Z(I,1), Y=FLOAT(N) AT Z(I,N). Y TAKES ON NON-INTEGER VALUES.  ISO00941
C                                                                       ISO00942
      COMMON /ISOSR6/ IX         ,IY         ,IDX        ,IDY        ,  ISO00943
     1                IS         ,ISS        ,NP         ,CV         ,  ISO00944
     2                INX(8)     ,INY(8)     ,IR(500)    ,NR            ISO00945
      COMMON /ISOSR9/ BIG        ,IXBIT                                 ISO00946
C                                                                       ISO00947
      LOGICAL         IPEN       ,IPENO                                 ISO00948
C                                                                       ISO00949
      DATA IOFFP,SPVAL/0,0./                                            ISO00950
      DATA IPEN,IPENO/.TRUE.,.TRUE./                                    ISO00951
C                                                                       ISO00952
C  PACK X AND Y                                                         ISO00953
C                                                                       ISO00954
      IPXY(I1,J1) = ISHFT(I1,IXBIT)+J1                                  ISO00955
      FX(X1,Y1) = X1                                                    ISO00956
      FY(X1,Y1) = Y1                                                    ISO00957
      C(P11,P21) = (P11-CV)/(P11-P21)                                   ISO00958
C                                                                       ISO00959
      M = MM                                                            ISO00960
      N = NN                                                            ISO00961
      IF (IOFFP .EQ. 0) GO TO  10                                       ISO00962
      ASSIGN 100 TO JUMP1                                               ISO00963
      ASSIGN 150 TO JUMP2                                               ISO00964
      GO TO  20                                                         ISO00965
   10 ASSIGN 120 TO JUMP1                                               ISO00966
      ASSIGN 160 TO JUMP2                                               ISO00967
   20 IX0 = IX                                                          ISO00968
      IY0 = IY                                                          ISO00969
      IS0 = IS                                                          ISO00970
      IF (IOFFP .EQ. 0) GO TO  30                                       ISO00971
      IX2 = IX+INX(IS)                                                  ISO00972
      IY2 = IY+INY(IS)                                                  ISO00973
      IPEN = Z(IX,IY).NE.SPVAL .AND. Z(IX2,IY2).NE.SPVAL                ISO00974
      IPENO = IPEN                                                      ISO00975
   30 IF (IDX .EQ. 0) GO TO  40                                         ISO00976
      Y = IY                                                            ISO00977
      ISUB = IX+IDX                                                     ISO00978
      X = C(Z(IX,IY),Z(ISUB,IY))*FLOAT(IDX)+FLOAT(IX)                   ISO00979
      GO TO  50                                                         ISO00980
   40 X = IX                                                            ISO00981
      ISUB = IY+IDY                                                     ISO00982
      Y = C(Z(IX,IY),Z(IX,ISUB))*FLOAT(IDY)+FLOAT(IY)                   ISO00983
   50 IF (IPEN) CALL FRSTS (FX(X,Y),FY(X,Y),1)                          ISO00984
   60 IS = IS+1                                                         ISO00985
      IF (IS .GT. 8) IS = IS-8                                          ISO00986
      IDX = INX(IS)                                                     ISO00987
      IDY = INY(IS)                                                     ISO00988
      IX2 = IX+IDX                                                      ISO00989
      IY2 = IY+IDY                                                      ISO00990
      IF (ISS .NE. 0) GO TO  70                                         ISO00991
      IF (IX2.GT.M .OR. IY2.GT.N .OR. IX2.LT.1 .OR. IY2.LT.1) GO TO 190 ISO00992
   70 IF (CV-Z(IX2,IY2))  80, 80, 90                                    ISO00993
   80 IS = IS+4                                                         ISO00994
      IX = IX2                                                          ISO00995
      IY = IY2                                                          ISO00996
      GO TO  60                                                         ISO00997
   90 IF (IS/2*2 .EQ. IS) GO TO  60                                     ISO00998
      GO TO JUMP1,(100,120)                                             ISO00999
  100 ISBIG = IS+(8-IS)/6*8                                             ISO01000
      IX3 = IX+INX(ISBIG-1)                                             ISO01001
      IY3 = IY+INY(ISBIG-1)                                             ISO01002
      IX4 = IX+INX(ISBIG-2)                                             ISO01003
      IY4 = IY+INY(ISBIG-2)                                             ISO01004
      IPENO = IPEN                                                      ISO01005
      IF (ISS .NE. 0) GO TO 110                                         ISO01006
      IF (IX3.GT.M .OR. IY3.GT.N .OR. IX3.LT.1 .OR. IY3.LT.1) GO TO 190 ISO01007
      IF (IX4.GT.M .OR. IY4.GT.N .OR. IX4.LT.1 .OR. IY4.LT.1) GO TO 190 ISO01008
  110 IPEN = Z(IX,IY).NE.SPVAL .AND. Z(IX2,IY2).NE.SPVAL .AND.          ISO01009
     1       Z(IX3,IY3).NE.SPVAL .AND. Z(IX4,IY4).NE.SPVAL              ISO01010
  120 IF (IDX .EQ. 0) GO TO 130                                         ISO01011
      Y = IY                                                            ISO01012
      ISUB = IX+IDX                                                     ISO01013
      X = C(Z(IX,IY),Z(ISUB,IY))*FLOAT(IDX)+FLOAT(IX)                   ISO01014
      GO TO 140                                                         ISO01015
  130 X = IX                                                            ISO01016
      ISUB = IY+IDY                                                     ISO01017
      Y = C(Z(IX,IY),Z(IX,ISUB))*FLOAT(IDY)+FLOAT(IY)                   ISO01018
  140 GO TO JUMP2,(150,160)                                             ISO01019
  150 IF (.NOT.IPEN) GO TO 170                                          ISO01020
      IF (IPENO) GO TO 160                                              ISO01021
C                                                                       ISO01022
C END OF LINE SEGMENT                                                   ISO01023
C                                                                       ISO01024
      CALL FRSTS (D1,D2,3)                                              ISO01025
      CALL FRSTS (FX(XOLD,YOLD),FY(XOLD,YOLD),1)                        ISO01026
C                                                                       ISO01027
C CONTINUE LINE SEGMENT                                                 ISO01028
C                                                                       ISO01029
  160 CALL FRSTS (FX(X,Y),FY(X,Y),2)                                    ISO01030
  170 XOLD = X                                                          ISO01031
      YOLD = Y                                                          ISO01032
      IF (IS .NE. 1) GO TO 180                                          ISO01033
      NP = NP+1                                                         ISO01034
      IF (NP .GT. NR) GO TO 190                                         ISO01035
      IR(NP) = IPXY(IX,IY)                                              ISO01036
  180 IF (ISS .EQ. 0) GO TO  60                                         ISO01037
      IF (IX.NE.IX0 .OR. IY.NE.IY0 .OR. IS.NE.IS0) GO TO  60            ISO01038
C                                                                       ISO01039
C END OF LINE                                                           ISO01040
C                                                                       ISO01041
  190 CALL FRSTS (D1,D2,3)                                              ISO01042
      RETURN                                                            ISO01043
      END                                                               ISO01044

      SUBROUTINE FILLIN                                                 ISO01994
C                                                                       ISO01995
      COMMON /ISOSR2/ LX         ,NX         ,NY         ,ISCR(8,128),  ISO01996
     1                ISCA(8,128)                                       ISO01997
      COMMON /ISOSR5/ NBPW       ,MASK(16)   ,GENDON                    ISO01998
      LOGICAL         GENDON                                            ISO01999
      COMMON /ISOSR7/ IENTRY     ,IONES                                 ISO02000
      LOGICAL         YCHANG     ,HBFLAG     ,FIRST      ,IHF           ISO02001
C                                                                       ISO02002
      IF (IENTRY .EQ. 0) RETURN                                         ISO02003
C                                                                       ISO02004
C  THIS IS A SHADING ALGORITHM IT IS USED TO DETERMINE CONTOUR LINES    ISO02005
C  THAT ARE HIDDEN BY THE PRESENT LINE. THE ALGORITHM PROCESSES         ISO02006
C  HORIZONTAL ROWS. IT ASSUMES THAT THE BIT PATTERN PASSED TO IT        ISO02007
C  HAS ONLY BITS SET TO MARK THE START AND END OF SHADING. THE          ISO02008
C  ALGORITHM ALSO ASSUMES THAT WHEN AN ON BIT IS ENCOUNTERED THAT A     ISO02009
C  CORRESPONDING OFF BIT IS INCLUDED IN THE SAME ROW.                   ISO02010
C                                                                       ISO02011
C                                                                       ISO02012
C  PULL OUT ROWS OF THE CONTOUR PATTERN                                 ISO02013
C                                                                       ISO02014
      IBVAL = 0                                                         ISO02015
      DO  80 IYNOW=1,NY                                                 ISO02016
         DO  40 IXNOW=1,LX                                              ISO02017
C                                                                       ISO02018
C  IF NO ACTIVATED BITS BRANCH                                          ISO02019
C                                                                       ISO02020
            ICRWD = ISCR(IXNOW,IYNOW)                                   ISO02021
            IF (ICRWD .EQ. 0) GO TO  30                                 ISO02022
C                                                                       ISO02023
C  ACTIVATED BITS IN WORD SET SHADING FLAG                              ISO02024
C                                                                       ISO02025
C  CHECK BIT BY BIT FOR ON/OFF FLAGS                                    ISO02026
C                                                                       ISO02027
            DO  20 IB=1,NBPW                                            ISO02028
               IBIT = (NBPW+1)-IB                                       ISO02029
C                                                                       ISO02030
C                                                                       ISO02031
C  PULL OUT THE CURRENT GRID POINT VALUE                                ISO02032
C                                                                       ISO02033
               IVAL = IAND(ICRWD,MASK(IBIT))                            ISO02034
C                                                                       ISO02035
C  IF IVAL SET, THIS IS AN ON/OFF FLAG                                  ISO02036
C                                                                       ISO02037
               IF (IVAL .EQ. 0) GO TO  10                               ISO02038
C                                                                       ISO02039
C  FLAG BIT, ALWAYS SET                                                 ISO02040
C                                                                       ISO02041
               IBVAL = MOD(IBVAL+1,2)                                   ISO02042
               GO TO  20                                                ISO02043
C                                                                       ISO02044
C  SHADE THE SCREEN ACCORDING TO THE STATUS OF IBVAL                    ISO02045
C                                                                       ISO02046
   10          IF (IBVAL .NE. 0) ICRWD = IOR(ICRWD,MASK(IBIT))          ISO02047
C                                                                       ISO02048
   20       CONTINUE                                                    ISO02049
C                                                                       ISO02050
C  ZERO OUT THE SCREEN                                                  ISO02051
C                                                                       ISO02052
            ISCR(IXNOW,IYNOW) = 0                                       ISO02053
            ISCA(IXNOW,IYNOW) = IOR(ICRWD,ISCA(IXNOW,IYNOW))            ISO02054
            GO TO  40                                                   ISO02055
C                                                                       ISO02056
   30       IF (IBVAL .NE. 0) ISCA(IXNOW,IYNOW) = IONES                 ISO02057
   40    CONTINUE                                                       ISO02058
C                                                                       ISO02059
C  FIX FOR NONCORRECTABLE RUNAWAYS                                      ISO02060
C                                                                       ISO02061
         IF (IBVAL .EQ. 0) GO TO  80                                    ISO02062
         IBVAL = 0                                                      ISO02063
         DO  70 K=1,LX                                                  ISO02064
            ITEST = 0                                                   ISO02065
            IF (IYNOW .EQ. 1) GO TO  50                                 ISO02066
            ITEST = ISCA(K,IYNOW-1)                                     ISO02067
            IF (IYNOW .EQ. NY) GO TO  60                                ISO02068
   50       ITEST = IOR(ITEST,ISCA(K,IYNOW+1))                          ISO02069
   60       ISCA(K,IYNOW) = ITEST                                       ISO02070
   70    CONTINUE                                                       ISO02071
C                                                                       ISO02072
   80 CONTINUE                                                          ISO02073
      RETURN                                                            ISO02074
      END                                                               ISO02075
      SUBROUTINE FRSTC (MX,MY,IENT)                                     ISO01710
C                                                                       ISO01711
      COMMON /ISOSR2/ LX         ,NX         ,NY         ,ISCR(8,128),  ISO01712
     1                ISCA(8,128)                                       ISO01713
      COMMON /ISOSR4/ RX         ,RY                                    ISO01714
      COMMON /ISOSR5/ NBPW       ,MASK(16)   ,GENDON                    ISO01715
      LOGICAL         GENDON                                            ISO01716
      COMMON /ISOSR8/ NMASK(16)  ,IXOLD      ,IYOLD      ,IBTOLD     ,  ISO01717
     1                HBFLAG     ,IOSLSN     ,LRLX       ,IFSX       ,  ISO01718
     2                IFSY       ,FIRST      ,IYDIR      ,IHX        ,  ISO01719
     3                IHB        ,IHS        ,IHV        ,IVOLD      ,  ISO01720
     4                IVAL       ,IHRX       ,YCHANG     ,ITPD       ,  ISO01721
     5                IHF                                               ISO01722
      LOGICAL         YCHANG     ,HBFLAG     ,FIRST      ,IHF           ISO01723
C                                                                       ISO01724
C                                                                       ISO01725
C  DRAW LINE TO THE POINT MX,MY                                         ISO01726
C                                                                       ISO01727
C  ENTER THE POINT INTO THE CURRENT SCREEN, ISCR, IF THE POINT CONFORMS ISO01728
C  TO THE SHADING ALGORITHM.                                            ISO01729
C  THE POINT IS NOT ENTERED WHEN;                                       ISO01730
C  1. IT IS THE SAME POINT USED IN THE LAST CALL, RESOLUTION PROBLEM    ISO01731
C  2. IT IS PART OF A HORIZONTAL LINE BUT NOT AN END POINT              ISO01732
C  3.  THE ENTIRE CONTOUR RESTS ON A HORIZONTAL PLANE                   ISO01733
C                                                                       ISO01734
C  WHEN DRAWING A HORIZONTAL LINE THREE CONDITIONS EXIST;               ISO01735
C  1. WHEN THE LINE IS A HORIZONTAL STEP ENTER ONLY THE OUTSIDE POINT.  ISO01736
C     A HORIZONTAL STEP IS DEFINED BY THE ENTERING AND EXITING Y        ISO01737
C     DIRECTION THAT IS THE SAME.                                       ISO01738
C  2. ENTER BOTH END POINTS OF A HORIZONTAL TURNING POINT. A HORIZONTAL ISO01739
C     TURNING POINT IS A LINE WITH GREATER THAN 1 HORIZONTAL BITS       ISO01740
C     AND THE ENTERING AND EXITING Y DIRECTION IS DIFFIRENT.            ISO01741
C  3.  WHEN THE ENTIRE CONTOUR IS A HORIZONTAL LINE NO POINTS ARE       ISO01742
C      ENTERED. THIS CONDITION IS DETECTED BY THE STATUS OF YCHANG.     ISO01743
C      IF IT IS TRUE THEN THE CONTOUR IS NOT A SINGLE HORIZONTAL LINE.  ISO01744
C                                                                       ISO01745
C  THE PREVIOUS POINT IS ERASED IF IT IS A VERTICAL TURNING POINT.      ISO01746
C  A VERTICAL TURNING POINT IS A HORIZONTAL LINE WITH ONLY 1 POINT      ISO01747
C  AND THE ENTERING AND EXITING Y DIRECTION DIFFERS.THIS DATA IS        ISO01748
C  IN THE VARIABLES IOSLSN-OLD SLOPE AND ISLSGN-NEW SLOPE.              ISO01749
C  THE CHANGE IN SLOPE MUST BE -1 TO 1 OR 1 TO -1.                      ISO01750
C                                                                       ISO01751
C  OTHERWISE THE POINT IS ENTERED INTO ISCR.                            ISO01752
C                                                                       ISO01753
C  THE TWO ENTRY POINTS ARE REQUIRED BY THE HARDWARE DRAWING ROUTINES.  ISO01754
C  FIRSTC IS USED FOR THE FIRST POINT ON THE CONTOUR. THE REMAINING     ISO01755
C  POINTS ON THE SAME CONTOUR ARE ENTERED VIA VECTC.                    ISO01756
C                                                                       ISO01757
      DATA IONE/1/                                                      ISO01758
      AVE(A,B) = (A+B)*.5                                               ISO01759
C                                                                       ISO01760
C  COMPUTE VISIBILITY OF THIS POINT                                     ISO01761
C                                                                       ISO01762
C  WARNING                                                              ISO01763
C  IF X OR Y PLOTTER MAXIMUM VALUE RANGES FALL BELOW 101 THEN THE       ISO01764
C  FOLLOWING TWO STATEMENTS WHICH SET IX AND IY MUST BE CHANGED.        ISO01765
C  REPLACE THE CONSTANT 1.0 BY 0.5 IN THE STATEMENTS WHERE THE          ISO01766
C  MAXIMUM PLOTTER VALUE IS LESS THAN 101 FOR THAT DIRECTION.  THE      ISO01767
C  PLOTTER CORDINATE RANGES ARE SET IN SET32.                           ISO01768
C                                                                       ISO01769
      IX = FLOAT(MX-1)*RX+1.0                                           ISO01770
      NRLX = IX                                                         ISO01771
      IY = FLOAT(MY-1)*RY+1.0                                           ISO01772
      IBIT = NBPW-MOD(IX,NBPW)                                          ISO01773
      IX = IX/NBPW+1                                                    ISO01774
      IVNOW = IAND(ISHFT(ISCA(IX,IY),1-IBIT),IONE)                      ISO01775
C                                                                       ISO01776
C  DECIDE IF FRSTC OR VECTC CALL                                        ISO01777
C                                                                       ISO01778
      IF (IENT .NE. 1) GO TO  10                                        ISO01779
C                                                                       ISO01780
      XOLD = MX                                                         ISO01781
      YOLD = MY                                                         ISO01782
C                                                                       ISO01783
C                                                                       ISO01784
C  SET INITIAL VALUES                                                   ISO01785
C                                                                       ISO01786
      IHF = .FALSE.                                                     ISO01787
      IYDIR = 0                                                         ISO01788
      ITPD = 0                                                          ISO01789
      IVAL = 0                                                          ISO01790
      IOSLSN = 0                                                        ISO01791
      IFSX = NRLX                                                       ISO01792
      IFSY = IY                                                         ISO01793
      LASTV = IVNOW                                                     ISO01794
      HBFLAG = .FALSE.                                                  ISO01795
      YCHANG = .FALSE.                                                  ISO01796
      CALL PLOTIT (IFIX(XOLD),IFIX(YOLD),0)                             ISO01797
      GO TO 180                                                         ISO01798
C                                                                       ISO01799
C****************************  ENTRY VECTC  ****************************ISO01800
C     ENTRY VECTC (MX,MY)                                               ISO01801
C                                                                       ISO01802
   10 XNOW = MX                                                         ISO01803
      YNOW = MY                                                         ISO01804
      JUMP = IVNOW*2+LASTV+1                                            ISO01805
      GO TO ( 20, 30, 40, 50),JUMP                                      ISO01806
C                                                                       ISO01807
C BOTH VISIBLE                                                          ISO01808
C                                                                       ISO01809
   20 CALL PLOTIT (IFIX(XNOW),IFIX(YNOW),1)                             ISO01810
      GO TO  50                                                         ISO01811
C                                                                       ISO01812
C JUST TURNED VISIBLE                                                   ISO01813
C                                                                       ISO01814
   30 CALL PLOTIT (IFIX(AVE(XNOW,XOLD)),IFIX(AVE(YNOW,YOLD)),0)         ISO01815
      GO TO  50                                                         ISO01816
C                                                                       ISO01817
C JUST TURNED INVISIBLE                                                 ISO01818
C                                                                       ISO01819
   40 CALL PLOTIT (IFIX(AVE(XNOW,XOLD)),IFIX(AVE(YNOW,YOLD)),1)         ISO01820
C                                                                       ISO01821
C BOTH INVISIBLE                                                        ISO01822
C                                                                       ISO01823
   50 XOLD = XNOW                                                       ISO01824
      YOLD = YNOW                                                       ISO01825
      LASTV = IVNOW                                                     ISO01826
C                                                                       ISO01827
C  TEST FOR RESOLUTION PROBLEM                                          ISO01828
C                                                                       ISO01829
      IF (NRLX.EQ.LRLX .AND. IY.EQ.IYOLD) RETURN                        ISO01830
C                                                                       ISO01831
C  TEST FOR HORIZONTAL BITS                                             ISO01832
C                                                                       ISO01833
      IF (IYOLD .NE. IY) GO TO  70                                      ISO01834
C                                                                       ISO01835
C  HORIZONTAL BITS DETECTED. SET FLAG AND EXIT.                         ISO01836
C  THIS AND THE NEXT HORIZONTAL BIT TEST IS NECESSARY FOR ISCR TO       ISO01837
C  CONFORM TO THE SHADING ALGORITHM IN SUBROUTINE FILLIN                ISO01838
C                                                                       ISO01839
C                                                                       ISO01840
C  IF HORIZONTAL LINE PREVIOUSLY DETECTED EXIT                          ISO01841
C                                                                       ISO01842
      IF (.NOT.HBFLAG) GO TO  60                                        ISO01843
C                                                                       ISO01844
C  IF END OF CONTOUR ON A HORIZONTAL LINE BRANCH FOR SPECIAL PROCESSING.ISO01845
C                                                                       ISO01846
      IF (NRLX.EQ.IFSX .AND. IY.EQ.IFSY) GO TO 210                      ISO01847
      GO TO 200                                                         ISO01848
C                                                                       ISO01849
C  SAVE SLOPE PRIOR TO HORIZONTAL LINE                                  ISO01850
C                                                                       ISO01851
   60 IHX = IXOLD                                                       ISO01852
      IHB = IBTOLD                                                      ISO01853
      IHS = IOSLSN                                                      ISO01854
      IOSLSN = 0                                                        ISO01855
      HBFLAG = .TRUE.                                                   ISO01856
      IHRX = LRLX                                                       ISO01857
      IHV = IVOLD                                                       ISO01858
      IF (LRLX.EQ.IFSX .AND. IYOLD.EQ.IFSY) IHF = .TRUE.                ISO01859
C                                                                       ISO01860
C  THIS IS THE SECOND TRAP FOR END OF CONTOUR ON A HORIZONTAL LINE.     ISO01861
C                                                                       ISO01862
      IF (NRLX.EQ.IFSX .AND. IY.EQ.IFSY) GO TO 210                      ISO01863
      GO TO 200                                                         ISO01864
C                                                                       ISO01865
C  COMPUTE THE SLOPE TO THIS POINT                                      ISO01866
C                                                                       ISO01867
   70 IF (IY-IYOLD)  80, 90,100                                         ISO01868
   80 ISLSGN = 1                                                        ISO01869
      GO TO 110                                                         ISO01870
   90 ISLSGN = 0                                                        ISO01871
      GO TO 120                                                         ISO01872
  100 ISLSGN = -1                                                       ISO01873
  110 IF (IYDIR .EQ. 0) IYDIR = ISLSGN                                  ISO01874
  120 CONTINUE                                                          ISO01875
C                                                                       ISO01876
C  IF PROCESS REACHES THIS CODE THE CONTOUR IS NOT CONTAINED ON A SINGLEISO01877
C  HORIZONTAL  PLANE, SO RECORD THIS FACT BY SETTING Y CHANGE FLAG.     ISO01878
C                                                                       ISO01879
      YCHANG = .TRUE.                                                   ISO01880
C                                                                       ISO01881
C  TEST FOR END OF HORIZONTAL LINE                                      ISO01882
C                                                                       ISO01883
      IF (.NOT.HBFLAG) GO TO 160                                        ISO01884
      HBFLAG = .FALSE.                                                  ISO01885
C                                                                       ISO01886
C  HORIZONTAL LINE JUST ENDED                                           ISO01887
C                                                                       ISO01888
C  TEST FOR REDRAW                                                      ISO01889
C                                                                       ISO01890
      ITEMP = IAND(ISCR(IXOLD,IYOLD),MASK(IBTOLD))                      ISO01891
      IF ((IHV .EQ. 0) .AND. (ITEMP .EQ. 0)) GO TO 130                  ISO01892
C                                                                       ISO01893
C  REDRAWING ERASE THIS POINT                                           ISO01894
C                                                                       ISO01895
      ISCR(IXOLD,IYOLD) = IAND(ISCR(IXOLD,IYOLD),NMASK(IBTOLD))         ISO01896
      ISCR(IHX,IYOLD) = IAND(ISCR(IHX,IYOLD),NMASK(IHB))                ISO01897
      GO TO 170                                                         ISO01898
C                                                                       ISO01899
C  TEST FOR STEP PROBLEM                                                ISO01900
C                                                                       ISO01901
  130 IF (IHS .NE. ISLSGN) GO TO 140                                    ISO01902
C                                                                       ISO01903
C  STEP PROBLEM                                                         ISO01904
C                                                                       ISO01905
      GO TO 170                                                         ISO01906
C                                                                       ISO01907
C  TURNING PROBLEM HORIZONTAL LINE IS A TURNING POINT                   ISO01908
C                                                                       ISO01909
  140 CONTINUE                                                          ISO01910
C                                                                       ISO01911
C  ENTER THE TURNING POINT ONLY IF IT IS NOT THE SECOND SUCCEEDING      ISO01912
C  EVENT IN A ROW                                                       ISO01913
C                                                                       ISO01914
      ICTPD = 1                                                         ISO01915
      IF (IHRX .GT. NRLX) ICTPD = -1                                    ISO01916
      IF (ICTPD .NE. ITPD) GO TO 150                                    ISO01917
      ITPD = 0                                                          ISO01918
C                                                                       ISO01919
C  ERASE THE FIRST POINT                                                ISO01920
C                                                                       ISO01921
      ISCR(IHX,IYOLD) = IAND(ISCR(IHX,IYOLD),NMASK(IHB))                ISO01922
      GO TO 170                                                         ISO01923
C                                                                       ISO01924
C  ENTER THE TURNING POINT                                              ISO01925
C                                                                       ISO01926
  150 CONTINUE                                                          ISO01927
      ITPD = ICTPD                                                      ISO01928
C                                                                       ISO01929
C  ENTER THE SECOND POINT                                               ISO01930
C                                                                       ISO01931
      ISCR(IXOLD,IYOLD) = IOR(ISCR(IXOLD,IYOLD),MASK(IBTOLD))           ISO01932
      GO TO 170                                                         ISO01933
C                                                                       ISO01934
C  CHECK IF PREVIOUS ENTRY WAS A VERTICAL TURNING POINT.                ISO01935
C  IF SO ERASE IT.                                                      ISO01936
C                                                                       ISO01937
  160 IF (ISLSGN.EQ.IOSLSN .OR. (IOSLSN.EQ.0 .OR. ISLSGN.EQ.0))         ISO01938
     1    GO TO 170                                                     ISO01939
      ITPD = 0                                                          ISO01940
      ISCR(IXOLD,IYOLD) = IAND(ISCR(IXOLD,IYOLD),NMASK(IBTOLD))         ISO01941
C                                                                       ISO01942
  170 IOSLSN = ISLSGN                                                   ISO01943
C                                                                       ISO01944
C  CHECK IF THIS GRID POINT PREVIOUSLY ACTIVATED                        ISO01945
C                                                                       ISO01946
      IVAL = IAND(ISCR(IX,IY),MASK(IBIT))                               ISO01947
C                                                                       ISO01948
C  IF GRID POINTS ACTIVATED BRANCH                                      ISO01949
C                                                                       ISO01950
      IF (IVAL .NE. 0) GO TO 190                                        ISO01951
C                                                                       ISO01952
C  GRID POINT NOT ACTIVATED   SET AND EXIT                              ISO01953
C                                                                       ISO01954
  180 CONTINUE                                                          ISO01955
      ISCR(IX,IY) = IOR(ISCR(IX,IY),MASK(IBIT))                         ISO01956
      GO TO 200                                                         ISO01957
C                                                                       ISO01958
C  THIS POINT IS BEING REDRAWN SO ERASE IT.                             ISO01959
C  (THIS IS TO CONFORM WITH THE SHADING ALGORITHM, FILLIN.              ISO01960
C  HOWEVER IF BACK TO STARTING POINT DO NOT ERASE                       ISO01961
C                                                                       ISO01962
  190 IF (NRLX.EQ.IFSX .AND. IY.EQ.IFSY) RETURN                         ISO01963
      ISCR(IX,IY) = IAND(ISCR(IX,IY),NMASK(IBIT))                       ISO01964
C                                                                       ISO01965
C                                                                       ISO01966
  200 IXOLD = IX                                                        ISO01967
      LRLX = NRLX                                                       ISO01968
      IYOLD = IY                                                        ISO01969
      IBTOLD = IBIT                                                     ISO01970
      IVOLD = IVAL                                                      ISO01971
      RETURN                                                            ISO01972
C                                                                       ISO01973
C  PERFORM THIS OPERATION WHEN A CONTOUR STARTS OR ENDS ON A HORIZONTAL ISO01974
C  LINE.                                                                ISO01975
C                                                                       ISO01976
  210 CONTINUE                                                          ISO01977
C                                                                       ISO01978
C  ERASE THE FIRST POINT OF A CONTOUR WHEN IT IS PART OF A HORIZONTAL   ISO01979
C  LINE SEGMENT AND IS NOT THE ENDPOINT OF THE SEGMENT                  ISO01980
C                                                                       ISO01981
      IF (.NOT.IHF) GO TO 220                                           ISO01982
      ISCR(IX,IY) = IAND(ISCR(IX,IY),NMASK(IBIT))                       ISO01983
  220 CONTINUE                                                          ISO01984
C                                                                       ISO01985
C  ERASE THE FIRST POINT OF A HORIZONTAL LINE SEGMENT WHEN IT ENDS      ISO01986
C  THE CONTOUR AND IS NOT THE HIGHEST LINE SEG ON THS SIDE.             ISO01987
C                                                                       ISO01988
      IF (.NOT.YCHANG) GO TO 230                                        ISO01989
      IF (IYDIR .NE. IHS) GO TO 200                                     ISO01990
  230 ISCR(IHX,IY) = IAND(ISCR(IHX,IY),NMASK(IHB))                      ISO01991
      GO TO 200                                                         ISO01992
      END                                                               ISO01993
      SUBROUTINE FRSTS (XX,YY,IENT)                                     ISO01067
C                                                                       ISO01068
C THIS IS A SPECIAL VERSION OF THE SMOOTHING DASHED LINE PACKAGE.  LINESISO01069
C ARE SMOOTHED IN THE SAME WAY, BUT NO SOFTFARE DASHED LINES ARE USED.  ISO01070
C CONDITIONAL PLOTTING ROUTINES ARE CALL WHICH DETERMINE THE VISIBILITY ISO01071
C OF A LINE SEGMENT BEFORE PLOTTING.                                    ISO01072
C                                                                       ISO01073
      DIMENSION       XSAVE(70)  ,YSAVE(70)  ,XP(70)     ,YP(70)     ,  ISO01074
     1                TEMP(70)                                          ISO01075
C                                                                       ISO01076
      COMMON /ISOSR7/ IENTRY     ,IONES                                 ISO01077
C                                                                       ISO01078
      DATA NP/150/                                                      ISO01079
      DATA L1/70/                                                       ISO01080
      DATA TENSN/2.5/                                                   ISO01081
      DATA PI/3.14159265358/                                            ISO01082
      DATA SMALL/128./                                                  ISO01083
C                                                                       ISO01084
      AVE(A,B) = .5*(A+B)                                               ISO01085
C                                                                       ISO01086
C  DECIDE IF FRSTS,VECTS,LASTS CALL                                     ISO01087
C                                                                       ISO01088
      GO TO ( 10, 20, 40),IENT                                          ISO01089
   10 DEG = 180./PI                                                     ISO01090
      X = XX                                                            ISO01091
      Y = YY                                                            ISO01092
      LASTFL = 0                                                        ISO01093
      SSLP1 = 0.0                                                       ISO01094
      SSLPN = 0.0                                                       ISO01095
      XSVN = 0.0                                                        ISO01096
      YSVN = 0.0                                                        ISO01097
C                                                                       ISO01098
C INITIALIZE THE POINT AND SEGMENT COUNTER                              ISO01099
C N COUNTS THE NUMBER OF POINTS/SEGMENT                                 ISO01100
C                                                                       ISO01101
      N = 0                                                             ISO01102
C                                                                       ISO01103
C NSEG = 0       FIRST SEGMENT                                          ISO01104
C NSEG = 1       MORE THAN ONE SEGMENT                                  ISO01105
C                                                                       ISO01106
      NSEG = 0                                                          ISO01107
      CALL TR32 (X,Y,MX,MY)                                             ISO01108
C                                                                       ISO01109
C SAVE THE X,Y COORDINATES OF THE FIRST POINT                           ISO01110
C XSV1           CONTAINS THE X COORDINATE OF THE FIRST POINT           ISO01111
C                OF A LINE                                              ISO01112
C YSV1           CONTAINS THE Y COORDINATE OF THE FIRST POINT           ISO01113
C                OF A LINE                                              ISO01114
C                                                                       ISO01115
      XSV1 = MX                                                         ISO01116
      YSV1 = MY                                                         ISO01117
      GO TO  30                                                         ISO01118
C                                                                       ISO01119
C ************************* ENTRY VECTS *************************       ISO01120
C     ENTRY VECTS (XX,YY)                                               ISO01121
C                                                                       ISO01122
   20 X = XX                                                            ISO01123
      Y = YY                                                            ISO01124
C                                                                       ISO01125
C VECTS          SAVES THE X,Y COORDINATES OF THE ACCEPTED              ISO01126
C                POINTS ON A LINE SEGMENT                               ISO01127
C                                                                       ISO01128
      CALL TR32 (X,Y,MX,MY)                                             ISO01129
C                                                                       ISO01130
CIF THE NEW POINT IS TOO CLOSE TO THE PREVIOUS POINT, IGNORE IT         ISO01131
C                                                                       ISO01132
      IF (ABS(FLOAT(IFIX(XSVN)-MX))+ABS(FLOAT(IFIX(YSVN)-MY)) .LT.      ISO01133
     1    SMALL) RETURN                                                 ISO01134
      IFLAG = 0                                                         ISO01135
   30 N = N+1                                                           ISO01136
C                                                                       ISO01137
C SAVE THE X,Y COORDINATES OF EACH POINT OF THE SEGMENT                 ISO01138
C XSAVE          THE ARRAY OF X COORDINATES OF LINE SEGMENT             ISO01139
C YSAVE          THE ARRAY OF Y COORDINATES OF LINE SEGMENT             ISO01140
C                                                                       ISO01141
      XSAVE(N) = MX                                                     ISO01142
      YSAVE(N) = MY                                                     ISO01143
      XSVN = XSAVE(N)                                                   ISO01144
      YSVN = YSAVE(N)                                                   ISO01145
      IF (N .GE. L1-1) GO TO  50                                        ISO01146
      RETURN                                                            ISO01147
C                                                                       ISO01148
C ************************* ENTRY LASTS *************************       ISO01149
C     ENTRY LASTS                                                       ISO01150
C                                                                       ISO01151
   40 LASTFL = 1                                                        ISO01152
C                                                                       ISO01153
C LASTS          CHECKS FOR PERIODIC LINES AND SETS UP                  ISO01154
C                  THE CALLS TO KURV1S AND KURV2S                       ISO01155
C                                                                       ISO01156
C IFLAG = 0      OK TO CALL LASTS DIRECTLY                              ISO01157
C IFLAG = 1      LASTS WAS JUST CALLED FROM BY VECTS                    ISO01158
C                IGNORE CALL TO LASTS                                   ISO01159
C                                                                       ISO01160
      IF (IFLAG .EQ. 1) RETURN                                          ISO01161
C                                                                       ISO01162
C COMPARE THE LAST POINT OF SEGMENT WITH FIRST POINT OF LINE            ISO01163
C                                                                       ISO01164
   50 IFLAG = 1                                                         ISO01165
C                                                                       ISO01166
C IPRD = 0       PERIODIC LINE                                          ISO01167
C IPRD = 1       NON-PERIODIC LINE                                      ISO01168
C                                                                       ISO01169
      IPRD = 1                                                          ISO01170
      IF (ABS(XSV1-XSVN)+ABS(YSV1-YSVN) .LT. SMALL) IPRD = 0            ISO01171
C                                                                       ISO01172
C TAKE CARE OF THE CASE OF ONLY TWO DISTINCT P0INTS ON A LINE           ISO01173
C                                                                       ISO01174
      IF (NSEG .GE. 1) GO TO  70                                        ISO01175
      IF (N-2) 160,150, 60                                              ISO01176
   60 IF (N .GE. 4) GO TO  70                                           ISO01177
      DX = XSAVE(2)-XSAVE(1)                                            ISO01178
      DY = YSAVE(2)-YSAVE(1)                                            ISO01179
      SLOPE = ATAN2(DY,DX)*DEG+90.                                      ISO01180
      IF (SLOPE .GE. 360.) SLOPE = SLOPE-360.                           ISO01181
      IF (SLOPE .LE. 0.) SLOPE = SLOPE+360.                             ISO01182
      SLP1 = SLOPE                                                      ISO01183
      SLPN = SLOPE                                                      ISO01184
      ISLPSW = 0                                                        ISO01185
      SIGMA = TENSN                                                     ISO01186
      GO TO 110                                                         ISO01187
   70 SIGMA = TENSN                                                     ISO01188
      IF (IPRD .GE. 1) GO TO  90                                        ISO01189
      IF (NSEG .GE. 1) GO TO  80                                        ISO01190
C                                                                       ISO01191
C SET UP FLAGS FOR A  1  SEGMENT, PERIODIC LINE                         ISO01192
C                                                                       ISO01193
      ISLPSW = 4                                                        ISO01194
      XSAVE(N) = XSV1                                                   ISO01195
      YSAVE(N) = YSV1                                                   ISO01196
      GO TO 110                                                         ISO01197
C                                                                       ISO01198
C SET UP FLAGS FOR AN N-SEGMENT, PERIODIC LINE                          ISO01199
C                                                                       ISO01200
   80 SLP1 = SSLPN                                                      ISO01201
      SLPN = SSLP1                                                      ISO01202
      ISLPSW = 0                                                        ISO01203
      GO TO 110                                                         ISO01204
   90 IF (NSEG .GE. 1) GO TO 100                                        ISO01205
C                                                                       ISO01206
C SET UP FLAGS FOR THE 1ST SEGMENT OF A NON-PERIODIC LINE               ISO01207
C                                                                       ISO01208
      ISLPSW = 3                                                        ISO01209
      GO TO 110                                                         ISO01210
C                                                                       ISO01211
C SET UP FLAGS FOR THE NTH SEGMENT OF A NON-PERIODIC LINE               ISO01212
C                                                                       ISO01213
  100 SLP1 = SSLPN                                                      ISO01214
      ISLPSW = 1                                                        ISO01215
C                                                                       ISO01216
C CALL THE SMOOTHING ROUTINES                                           ISO01217
C                                                                       ISO01218
  110 CALL KURV1S (N,XSAVE,YSAVE,SLP1,SLPN,XP,YP,TEMP,S,SIGMA,ISLPSW)   ISO01219
      IF (IPRD.EQ.0 .AND. NSEG.EQ.0 .AND. S.LT.70.) GO TO 170           ISO01220
      IENTRY = 1                                                        ISO01221
C                                                                       ISO01222
C DETERMINE THE NUMBER OF POINTS TO INTERPOLATE FOR EACH SEGMENT        ISO01223
C                                                                       ISO01224
      IF (NSEG.GE.1 .AND. N.LT.L1-1) GO TO 120                          ISO01225
      NPRIME = FLOAT(NP)-(S*FLOAT(NP))/(2.*32768.)                      ISO01226
      IF (S .GE. 32768.) NPRIME = .5*FLOAT(NP)                          ISO01227
      NPL = FLOAT(NPRIME)*S/32768.                                      ISO01228
      IF (NPL .LT. 2) NPL = 2                                           ISO01229
  120 DT = 1./FLOAT(NPL)                                                ISO01230
      IF (NSEG .LE. 0) CALL FRSTC (IFIX(XSAVE(1)),IFIX(YSAVE(1)),1)     ISO01231
      T = 0.0                                                           ISO01232
      NSLPSW = 1                                                        ISO01233
      IF (NSEG .GE. 1) NSLPSW = 0                                       ISO01234
      NSEG = 1                                                          ISO01235
      CALL KURV2S (T,XS,YS,N,XSAVE,YSAVE,XP,YP,S,SIGMA,NSLPSW,SLP)      ISO01236
C                                                                       ISO01237
C SAVE SLOPE AT THE FIRST POINT OF THE LINE                             ISO01238
C                                                                       ISO01239
      IF (NSLPSW .GE. 1) SSLP1 = SLP                                    ISO01240
      NSLPSW = 0                                                        ISO01241
      XSOLD = XSAVE(1)                                                  ISO01242
      YSOLD = YSAVE(1)                                                  ISO01243
      DO 130 I=1,NPL                                                    ISO01244
         T = T+DT                                                       ISO01245
         TT = -T                                                        ISO01246
         IF (I .EQ. NPL) NSLPSW = 1                                     ISO01247
         CALL KURV2S (TT,XS,YS,N,XSAVE,YSAVE,XP,YP,S,SIGMA,NSLPSW,SLP)  ISO01248
C                                                                       ISO01249
C SAVE THE LAST SLOPE OF THIS LINE SEGMENT                              ISO01250
C                                                                       ISO01251
         IF (NSLPSW .GE. 1) SSLPN = SLP                                 ISO01252
C                                                                       ISO01253
C DRAW EACH PART OF THE LINE SEGMENT                                    ISO01254
C                                                                       ISO01255
         CALL FRSTC (IFIX(AVE(XSOLD,XS)),IFIX(AVE(YSOLD,YS)),2)         ISO01256
         CALL FRSTC (IFIX(XS),IFIX(YS),2)                               ISO01257
         XSOLD = XS                                                     ISO01258
         YSOLD = YS                                                     ISO01259
  130 CONTINUE                                                          ISO01260
      IF (IPRD .NE. 0) GO TO 140                                        ISO01261
C                                                                       ISO01262
C CONNECT THE LAST POINT WITH THE FIRST POINT OF A PERIODIC LINE        ISO01263
C                                                                       ISO01264
      CALL FRSTC (IFIX(AVE(XSOLD,XS)),IFIX(AVE(YSOLD,YS)),2)            ISO01265
      CALL FRSTC (IFIX(XSV1),IFIX(YSV1),2)                              ISO01266
C                                                                       ISO01267
C BEGIN THE NEXT LINE SEGMENT WITH THE LAST POINT OF THIS SEGMENT       ISO01268
C                                                                       ISO01269
  140 XSAVE(1) = XS                                                     ISO01270
      YSAVE(1) = YS                                                     ISO01271
      N = 1                                                             ISO01272
  150 CONTINUE                                                          ISO01273
  160 RETURN                                                            ISO01274
  170 N = 0                                                             ISO01275
      RETURN                                                            ISO01276
      END                                                               ISO01277
      BLOCKDATA ISOSRB                                                  ISO02084
C                                                                       ISO02085
C     BLOCK DATA                                                        ISO02086
C                                                                       ISO02087
      COMMON /ISOSR2/ LX         ,NX         ,NY         ,ISCR(8,128),  ISO02088
     1                ISCA(8,128)                                       ISO02089
      COMMON /ISOSR4/ RX         ,RY                                    ISO02090
      COMMON /ISOSR5/ NBPW       ,MASK(16)   ,GENDON                    ISO02091
      LOGICAL         GENDON                                            ISO02092
      COMMON /ISOSR6/ IX         ,IY         ,IDX        ,IDY        ,  ISO02093
     1                IS         ,ISS        ,NP         ,CV         ,  ISO02094
     2                INX(8)     ,INY(8)     ,IR(500)    ,NR            ISO02095
      COMMON /ISOSR7/ IENTRY     ,IONES                                 ISO02096
      COMMON /ISOSR8/ NMASK(16)  ,IXOLD      ,IYOLD      ,IBTOLD     ,  ISO02097
     1                HBFLAG     ,IOSLSN     ,LRLX       ,IFSX       ,  ISO02098
     2                IFSY       ,FIRST      ,IYDIR      ,IHX        ,  ISO02099
     3                IHB        ,IHS        ,IHV        ,IVOLD      ,  ISO02100
     4                IVAL       ,IHRX       ,YCHANG     ,ITPD       ,  ISO02101
     5                IHF                                               ISO02102
      COMMON /ISOSR9/ BIG        ,IXBIT                                 ISO02103
      COMMON /TEMPR/  RZERO                                             ISO02104
      LOGICAL         YCHANG     ,HBFLAG     ,FIRST      ,IHF           ISO02105
C                                                                       ISO02106
      DATA LX,NX,NY/8,128,128/                                          ISO02107
      DATA INX(1),INX(2),INX(3),INX(4),INX(5),INX(6),INX(7),INX(8)/     ISO02108
     1        -1 ,   -1 ,    0 ,    1 ,    1 ,    1 ,    0 ,   -1 /     ISO02109
      DATA INY(1),INY(2),INY(3),INY(4),INY(5),INY(6),INY(7),INY(8)/     ISO02110
     1         0 ,    1 ,    1 ,    1 ,    0 ,   -1 ,   -1 ,   -1 /     ISO02111
      DATA NR/500/                                                      ISO02112
      DATA NBPW/16/                                                     ISO02113
      DATA IHF/.FALSE./                                                 ISO02114
C                                                                       ISO02115
      DATA GENDON /.FALSE./                                             ISO02116
      DATA RZERO/0./                                                    ISO02117
C                                                                       ISO02118
C                                                                       ISO02119
C RX = (NX-1)/SCREEN WIDTH FROM TRN32I                                  ISO02120
C RY = (NY-1)/SCREEN HEIGHT FROM TRN32I                                 ISO02121
C                                                                       ISO02122
      DATA RX,RY/.00389,.00389/                                         ISO02123
C                                                                       ISO02124
      END                                                               ISO02125
      SUBROUTINE ISOSRF (T,LU,MU,LV,MV,MW,EYE,MUVWP2,SLAB,TISO,IFLAG)   ISO00001
C                                                                       ISO00002
C                                                                       ISO00003
C DIMENSION OF           T(LU,LV,MW),EYE(3),SLAB(MUVWP2,MUVWP2)         ISO00004
C ARGUMENTS                                                             ISO00005
C                                                                       ISO00006
C                       JUNE 1980                                       ISO00007
C                                                                       ISO00008
C PURPOSE                ISOSRF DRAWS AN APPROXIMATION OF AN ISO-VALUED ISO00009
C                        SURFACE FROM A THREE-DIMENSIONAL ARRAY WITH    ISO00010
C                        HIDDEN LINES REMOVED.                          ISO00011
C                                                                       ISO00012
C ACCESS CARDS           *FORTRAN,S=ULIB,N=ISOSRF                       ISO00013
C                                                                       ISO00014
C                                                                       ISO00015
C USAGE                  IF THE FOLLOWING ASSUMPTIONS ARE MET, USE      ISO00016
C                            CALL EZISOS (T,MU,MV,MW,EYE,SLAB,TISO)     ISO00017
C                               ASSUMPTIONS:                            ISO00018
C                                  ALL OF THE T ARRAY IS TO BE USED.    ISO00019
C                                  IFLAG IS CHOSEN INTERNALLY.          ISO00020
C                                  FRAME IS CALLED BY EZISOS.           ISO00021
C                        IF THE ASSUMPTIONS ARE NOT MET, USE            ISO00022
C                            CALL ISOSRF (T,LU,MU,LV,MV,MW,EYE,MUVWP2,  ISO00023
C                                         SLAB,TISO,IFLAG)              ISO00024
C                                                                       ISO00025
C ARGUMENTS                                                             ISO00026
C                                                                       ISO00027
C ON INPUT               T                                              ISO00028
C                          THREE DIMENSIONAL ARRAY OF DATA THAT IS USED ISO00029
C                          TO DETERMINE THE ISO-VALUED SURFACE.         ISO00030
C                                                                       ISO00031
C                        LU                                             ISO00032
C                          FIRST DIMENSION OF T IN THE CALLING PROGRAM. ISO00033
C                                                                       ISO00034
C                        MU                                             ISO00035
C                          THE NUMBER OF DATA VALUES OF T TO BE         ISO00036
C                          PROCESSED IN THE U DIRECTION (THE FIRST      ISO00037
C                          SUBSCRIPT DIRECTION).  WHEN PROCESSING THE   ISO00038
C                          ENTIRE ARRAY, LU = MU (AND LV = MV).  SEE    ISO00039
C                          APPENDIX 1 OF THE GRAPHICS CHAPTER FOR AN    ISO00040
C                          EXPLANATION OF USING THIS ARGUMENT LIST TO   ISO00041
C                          PROCESS ANY PART OF AN ARRAY.                ISO00042
C                                                                       ISO00043
C                        LV                                             ISO00044
C                          SECOND DIMENSION OF T IN THE CALLING PROGRAM.ISO00045
C                                                                       ISO00046
C                        MV                                             ISO00047
C                          THE NUMBER OF DATA VALUES OF T TO BE         ISO00048
C                          PROCESSED IN THE V DIRECTION (THE SECOND     ISO00049
C                          SUBSCRIPT DIRECTION).                        ISO00050
C                                                                       ISO00051
C                        MW                                             ISO00052
C                          THE NUMBER OF DATA VALUES OF T TO BE         ISO00053
C                          PROCESSED IN THE W DIRECTION (THE THIRD      ISO00054
C                          SUBSCRIPT DIRECTION).                        ISO00055
C                                                                       ISO00056
C                        EYE                                            ISO00057
C                          THE POSITION OF THE EYE IN THREE-SPACE.  T ISISO00058
C                          CONSIDERED TO BE IN A BOX WITH OPPOSITE      ISO00059
C                          CORNERS (1,1,1) AND (MU,MV,MW).  THE EYE IS  ISO00060
C                          AT (EYE(1),EYE(2),EYE(3)), WHICH MUST BE     ISO00061
C                          OUTSIDE THE BOX THAT T IS IN.  WHILE GAINING ISO00062
C                          EXPERIENCE WITH THE ROUTINE, A GOOD CHOICE   ISO00063
C                          FOR EYE MIGHT BE (5.0*MU,3.5*MV,2.0*MW).     ISO00064
C                                                                       ISO00065
C                        MUVWP2                                         ISO00066
C                          THE MAXIMUM OF (MU,MV,MW)+2                  ISO00067
C                          (MUVWP2 = MAX0(MU,MV,MW)+2).                 ISO00068
C                                                                       ISO00069
C                        SLAB                                           ISO00070
C                          A WORK SPACE USED FOR INTERNAL STORAGE.  SLABISO00071
C                          MUST BE AT LEAST MUVWP2*MUVWP2 WORDS LONG.   ISO00072
C                                                                       ISO00073
C                        TISO                                           ISO00074
C                          THE ISO-VALUE USED TO DEFINE THE SURFACE; THEISO00075
C                          SURFACE DRAWN WILL SEPARATE VOLUMES OF T THATISO00076
C                          HAVE VALUE GREATER THAN TISO FROM VOLUMES OF ISO00077
C                          T THAT HAVE VALUE LESS THAN TISO.            ISO00078
C                                                                       ISO00079
C                        IFLAG                                          ISO00080
C                          A FLAG WHICH SERVES TWO PURPOSES.            ISO00081
C                          .  FIRST, THE ABSOLUTE VALUE OF IFLAG        ISO00082
C                             DETERMINES WHICH TYPES OF LINES ARE DRAWN ISO00083
C                             TO APPROXIMATE THE SURFACE.  THREE TYPES  ISO00084
C                             OF LINES ARE CONSIDERED:  LINES OF        ISO00085
C                             CONSTANT U, LINES OF CONSTANT V AND LINES ISO00086
C                             OF CONSTANT W.  THE FOLLOWING TABLE LISTS ISO00087
C                             THE TYPES OF LINES DRAWN.                 ISO00088
C                                                                       ISO00089
C                                                LINES OF CONSTANT      ISO00090
C                                                -----------------      ISO00091
C                             IABS(IFLAG)        U    V    W            ISO00092
C                                  1             NO   NO   YES          ISO00093
C                                  2             NO   YES  NO           ISO00094
C                                  3             NO   YES  YES          ISO00095
C                                  4             YES  NO   NO           ISO00096
C                                  5             YES  NO   YES          ISO00097
C                                  6             YES  YES  NO           ISO00098
C                             0, 7 OR MORE       YES  YES  YES          ISO00099
C                                                                       ISO00100
C                          .  SECOND, THE SIGN OF IFLAG DETERMINES WHAT ISO00101
C                             IS INSIDE AND WHAT IS OUTSIDE, HENCE,     ISO00102
C                             WHICH LINES ARE VISIBLE AND WHAT IS DONE  ISO00103
C                             AT THE BOUNDARY OF T.  FOR IFLAG:         ISO00104
C                                                                       ISO00105
C                             POSITIVE   T VALUES GREATER THAN TISO ARE ISO00106
C                                        ASSUMED TO BE INSIDE THE SOLID ISO00107
C                                        FORMED BY THE DRAWN SURFACE.   ISO00108
C                             NEGATIVE   T VALUES LESS THAN TISO ARE    ISO00109
C                                        ASSUMED TO BE INSIDE THE SOLID ISO00110
C                                        FORMED BY THE DRAWN SURFACE.   ISO00111
C                             IF THE ALGORITHM DRAWS A CUBE, REVERSE THEISO00112
C                             SIGN OF IFLAG.                            ISO00113
C                                                                       ISO00114
C ON OUTPUT              T,LU,MU,LV,MV,MW,EYE,MUVWP2,TISO AND IFLAG ARE ISO00115
C                        UNCHANGED.  SLAB HAS BEEN WRITTEN IN.          ISO00116
C                                                                       ISO00117
C NOTE                   .   THIS ROUTINE IS FOR LOWER RESOLUTION ARRAYSISO00118
C                            THAN ISOSRFHR.  40 BY 40 BY 30 IS A        ISO00119
C                            PRACTICAL MAXIMUM.                         ISO00120
C                        .   TRANSFORMATIONS CAN BE ACHIEVED BY         ISO00121
C                            ADJUSTING SCALING STATEMENT FUNCTIONS IN   ISO00122
C                            ISOSRF, SET3D AND TR32.                    ISO00123
C                        .   THE HIDDEN-LINE ALGORITHM IS NOT EXACT, SO ISO00124
C                            VISIBILITY ERRORS CAN OCCUR.               ISO00125
C                                                                       ISO00126
C ENTRY POINTS           ISOSRF, EZISOS, SET3D, TRN32I, ZEROSC,         ISO00127
C                        STCNTR, DRCNTR, TR32, FRSTS, KURV1S, KURV2S,   ISO00128
C                        FRSTC, FILLIN, DRAWI, ISOSRB, MMASK            ISO00129
C                                                                       ISO00130
C COMMON BLOCKS             NAME     LENGTH                             ISO00131
C                          -----     ------                             ISO00132
C                          ISOSR1         4                             ISO00133
C                          ISOSR2      4003 (OCTAL)                     ISO00134
C                          ISOSR3         7                             ISO00135
C                          ISOSR4         2                             ISO00136
C                          ISOSR5        22 (OCTAL)                     ISO00137
C                          ISOSR6      1015 (OCTAL)                     ISO00138
C                          ISOSR7          2                            ISO00139
C                          ISOSR8        44 (OCTAL)                     ISO00140
C                          ISOSR9          2                            ISO00141
C                          TEMPR           1                            ISO00142
C                          PWRZ1I        12 (OCTAL)                     ISO00143
C                                                                       ISO00144
C I/O                    PLOTS SURFACE                                  ISO00145
C                                                                       ISO00146
C PRECISION              SINGLE                                         ISO00147
C                                                                       ISO00148
C REQUIRED ULIB          NONE                                           ISO00149
C ROUTINES                                                              ISO00150
C                                                                       ISO00151
C                                                                       ISO00152
C LANGUAGE               FORTRAN                                        ISO00153
C                                                                       ISO00154
C HISTORY                DEVELOPED FOR USERS OF ISOSRFHR WITH SMALLER   ISO00155
C                        ARRAYS.                                        ISO00156
C                                                                       ISO00157
C ALGORITHM              CUTS THROUGH THE THREE-DIMENSIONAL ARRAY ARE   ISO00158
C                        CONTOURED WITH A SMOOTHING CONTOURER WHICH ALSOISO00159
C                        MARKS A MODEL OF THE PLOTTING PLANE.  INTERIORSISO00160
C                        OF BOUNDARIES ARE FILLED IN AND THE RESULT IS  ISO00161
C                        .OR.ED INTO ANOTHER MODEL OF THE PLOTTING PLANEISO00162
C                        WHICH IS USED TO TEST SUBSEQUENT CONTOUR LINES ISO00163
C                        FOR VISIBILITY.                                ISO00164
C                                                                       ISO00165
C SPACE REQUIRED         ABOUT 11000 (OCTAL) NOT INCLUDING THE SYSTEM   ISO00166
C                        PLOT PACKAGE.                                  ISO00167
C                                                                       ISO00168
C TIMING                 VARIES WIDELY WITH SIZE OF T AND THE VOLUME OF ISO00169
C                        THE SPACE ENCLOSED BY THE SURFACE DRAWN.  THE  ISO00170
C                        SAMPLE PICTURE TOOK ABOUT 1 SECOND OF 7600     ISO00171
C                        TIME.                                          ISO00172
C                                                                       ISO00173
C IMPLEMENTATION                                                        ISO00174
C               THE IMPLEMENTATION OF ISOSRF REQUIRES THE CODING OF     ISO00175
C               SOME LOCAL ROUTINES LISTED BELOW:                       ISO00176
C                                                                       ISO00177
C               SUBROUTINE PLOTIT(IX,IY,IPEN)                           ISO00178
C                       MOVE TO INTEGER IX,IY COORDINATES BASED ON A    ISO00179
C                       0 TO 32767 COORDINATE GRID.                     ISO00180
C                       IPEN=0 MOVE ONLY (PEN UP)                       ISO00181
C                       IPEN=1 DRAW (PEN DOWN)                          ISO00182
C                                                                       ISO00183
C               FUNCTION ISHFT(IWORD,N)                                 ISO00184
C                       IWORD IS SHIFTED N BITS                         ISO00185
C                       IF N.GT.0 LEFT CIRCULAR SHIFT                   ISO00186
C                       IF N.LT.0 RIGHT END OFF SHIFT                   ISO00187
C                                 (LOGICAL OR SIGN EXTEND)              ISO00188
C                       IF N=0 NO ACTION                                ISO00189
C                       RETURN SHIFTED RESULT AS FUNCTION VALUE         ISO00190
C                                                                       ISO00191
C               FUNCTION IAND(I1,I2)                                    ISO00192
C                       LOGICAL AND OF I1 AND I2                        ISO00193
C                       RETURN RESULT AS FUNCTION VALUE                 ISO00194
C                                                                       ISO00195
C               FUNCTION IOR(I1,I2)                                     ISO00196
C                                                                       ISO00197
C                       LOGICAL OR OF I1 TO I2                          ISO00198
C                       RETURN RESULT AS FUNCTION VALUE                 ISO00199
C                                                                       ISO00200
C               FUNCTION R1MACH(IARG)                                   ISO00201
C                       THIS IS FROM THE PORT LIBRARY (A BELL LABS      ISO00202
C                       PRODUCT).  IF YOU DO NOT HAVE THE PORT ROUTINES ISO00203
C                       THE FUNCTION OF R1MACH IN ISOSRF IS TO RETURN   ISO00204
C                       THE LARGEST FLOATING POINT NUMBER USEABLE ON THEISO00205
C                       MACHINE RUNNING ISOSRF.  RETURN THIS FLOATING   ISO00206
C                       POINT VALUE AS THE FUNCTION VALUE AND IGNORE    ISO00207
C                       THE INPUT PARAMETER.                            ISO00208
C                                                                       ISO00209
C               SUBROUTINE ULIBER(IERR,MESS,LMESS)                      ISO00210
C                       PRINTS ERROR MESSAGE OR PRINTS ERROR NUMBER AND ISO00211
C                       ERROR MESSAGE.                                  ISO00212
C                       IERR=ERROR NUMBER (NOT PRINTED IF 0)            ISO00213
C                       MESS=MESSAGE TO BE PRINTED INCLUDING CARRIAGE   ISO00214
C                            CONTROL                                    ISO00215
C                       LMESS=NUMBER OF CHARACTERS IN MESS( .LE. 130)   ISO00216
C                                                                       ISO00217
C  **NOTE**                  SPACE REQUIREMENTS CAN BE REDUCED BY       ISO00218
C                            CHANGING THE SIZE OF THE ARRAYS ISCR, ISCA ISO00219
C                            (FOUND IN COMMON ISOSR2), MASK(FOUND IN    ISO00220
C                            COMMON ISOSR5) AND THE VARIABLE NBPW       ISO00221
C                            (COMMON ISOSR5).                           ISO00222
C                            ISCR AND ISCA NEED 128X128 BITS.  SO ON A  ISO00223
C                            64 BIT MACHINE ISCR, ISCA CAN BE           ISO00224
C                            DIMENSIONED TO (2,128). NBPW SET IN        ISO00225
C                            SUBROUTINE MMASK SHOULD CONTAIN THE        ISO00226
C                            NUMBER OF BITS PER WORD YOU WISH TO        ISO00227
C                            UTILIZE.                                   ISO00228
C                            THE DIMENSION OF MASK AND NMASK SHOULD     ISO00229
C                            EQUAL THE VALUE OF NBPW.                   ISO00230
C                            LX SHOULD BE SET TO THE FIRST DIMENSION    ISO00231
C                            OF ISCA AND ISCR.                          ISO00232
C                            EXAMPLE                                    ISO00233
C                                 SUPPOSE YOU HAVE A 60 BIT MACHINE     ISO00234
C                                 DIMENSION ISCA, ISCR TO 4,128         ISO00235
C                                 NBPW = 32                             ISO00236
C                                 DIMENSION MASK TO 32                  ISO00237
C                                 SUPPOSE YOU HAVE A 64 BIT MACHINE     ISO00238
C                                 DIMENSION ISCA, ISCR TO 2,128         ISO00239
C                                 NBPW = 64                             ISO00240
C                                 DIMENSION MASK TO 64                  ISO00241
C                                                                       ISO00242
C                                                                       ISO00243
C PLOTTING ROUTINES      PLOTIT, FRAME                                  ISO00244
C USED                                                                  ISO00245
C                                                                       ISO00246
C REQUIRED RESIDENT      SQRT, ACOS, SIN, COS, ATAN2, EXP               ISO00247
C ROUTINES                                                              ISO00248
C                                                                       ISO00249
C INTERNAL PARAMETERS    NAME  DEFAULT           FUNCTION               ISO00250
C                        ----  -------           --------               ISO00251
C                        IREF     1     FLAG TO CONTROL DRAWING OF AXES.ISO00252
C                                       .IREF'0  MEANS DRAW AXES.       ISO00253
C                                       .IREF=0  MEANS DO NOT.          ISO00254
C                                                                       ISO00255
C                                                                       ISO00256
C                                                                       ISO00257
C                                                                       ISO00258
      DIMENSION T(LU,LV,MW),EYE(3),SLAB(MUVWP2,MUVWP2)                  ISO00260
      COMMON /ISOSR1/ ISLBT      ,U          ,V          ,W             ISO00261
      COMMON /ISOSR2/ LX         ,NX         ,NY         ,ISCR(8,128),  ISO00262
     1                ISCA(8,128)                                       ISO00263
      COMMON /ISOSR3/ ISCALE     ,XMIN       ,XMAX       ,YMIN       ,  ISO00264
     1                YMAX       ,BIGD       ,R0                        ISO00265
      COMMON /ISOSR4/ RX         ,RY                                    ISO00266
      COMMON /ISOSR5/ NBPW       ,MASK(16)   ,GENDON                    ISO00267
      COMMON /ISOSR6/ IX         ,IY         ,IDX        ,IDY        ,  ISO00268
     1                IS         ,ISS        ,NP         ,CV         ,  ISO00269
     2                INX(8)     ,INY(8)     ,IR(500)    ,NR            ISO00270
      COMMON /ISOSR7/ IENTRY     ,IONES                                 ISO00271
      COMMON /ISOSR9/ BIG        ,IXBIT                                 ISO00272
      EXTERNAL        ISOSRB                                            ISO00273
C                                                                       ISO00274
      LOGICAL         GENDON                                            ISO00275
      DATA IREF/1/                                                      ISO00276
C                                                                       ISO00277
      AVE(A,B) = (A+B)*.5                                               ISO00278
C                                                                       ISO00279
C A.S.F. FOR SCALING                                                    ISO00280
C                                                                       ISO00281
      SU(UTEMP) = UTEMP                                                 ISO00282
      SV(VTEMP) = VTEMP                                                 ISO00283
      SW(WTEMP) = WTEMP                                                 ISO00284
C                                                                       ISO00285
C 3-SPACE  U,V,W,IU,IV,IW,ETC                                           ISO00290
C 2-SPACE  X,Y,IX,IY,ETC                                                ISO00291
C                                                                       ISO00292
C  INITIALIZE MASKS                                                     ISO00293
C                                                                       ISO00294
      IF (.NOT.GENDON) CALL MMASK                                       ISO00295
C                                                                       ISO00296
C  SET SHIFT VALUE FOR X,Y PACKING                                      ISO00297
C                                                                       ISO00298
C  IF YOUR MACHINE HAS MORE THAN 16 BITS PER WORD THIS CHECK MAY BE     ISO00299
C  MODIFIED                                                             ISO00300
C                                                                       ISO00301
      IF (LU .LE. 256) GO TO  10                                        ISO00302
      CALL ULIBER (0,29HDIMENSION OF CUBE EXCEEDS 256,29)               ISO00303
      RETURN                                                            ISO00304
   10 DO  20 J=1,30                                                     ISO00305
         IF (LU .LE. 2**(J-1)) GO TO  30                                ISO00306
   20 CONTINUE                                                          ISO00307
   30 IXBIT = J                                                         ISO00308
      NU = MU                                                           ISO00309
      NUP2 = NU+2                                                       ISO00310
      NV = MV                                                           ISO00311
      NVP2 = NV+2                                                       ISO00312
      NW = MW                                                           ISO00313
      NWP2 = NW+2                                                       ISO00314
      FNU = NU                                                          ISO00315
      FNV = NV                                                          ISO00316
      FNW = NW                                                          ISO00317
      SU1 = SU(1.)                                                      ISO00318
      SV1 = SV(1.)                                                      ISO00319
      SW1 = SW(1.)                                                      ISO00320
      SUNU = SU(FNU)                                                    ISO00321
      SVNV = SV(FNV)                                                    ISO00322
      SWNW = SW(FNW)                                                    ISO00323
      AVEU = AVE(SU1,SUNU)                                              ISO00324
      AVEV = AVE(SV1,SVNV)                                              ISO00325
      AVEW = AVE(SW1,SWNW)                                              ISO00326
      EYEU = EYE(1)                                                     ISO00327
      EYEV = EYE(2)                                                     ISO00328
      EYEW = EYE(3)                                                     ISO00329
      NUVWP2 = MUVWP2                                                   ISO00330
      TVAL = TISO                                                       ISO00331
      NFLAG = IABS(IFLAG)                                               ISO00332
      IF (NFLAG.EQ.0 .OR. NFLAG.GE.8) NFLAG = 7                         ISO00333
C                                                                       ISO00334
C SET UP SCALING                                                        ISO00335
C                                                                       ISO00336
      FACT = -ISIGN(1,IFLAG)                                            ISO00337
      CALL SET3D (EYE,1.,FNU,1.,FNV,1.,FNW)                             ISO00338
C                                                                       ISO00339
C BOUND LOWER AND LEFT EDGE OF SLAB                                     ISO00340
C                                                                       ISO00341
      EDGE = SIGN(BIG,FACT)                                             ISO00342
      DO  40 IUVW=1,NUVWP2                                              ISO00343
         SLAB(IUVW,1) = EDGE                                            ISO00344
         SLAB(1,IUVW) = EDGE                                            ISO00345
   40 CONTINUE                                                          ISO00346
C                                                                       ISO00347
C SLICES PERPENDICULAR TO U. THAT IS, V W SLICES. T OF CONSTANT U.      ISO00348
C                                                                       ISO00349
      IF (NFLAG .LT. 4) GO TO 100                                       ISO00350
      CALL ZEROSC                                                       ISO00351
      ISLBT = -1                                                        ISO00352
C                                                                       ISO00353
C BOUND UPPER AND RIGHT EDGE OF SLAB.                                   ISO00354
C                                                                       ISO00355
      DO  50 IV=2,NVP2                                                  ISO00356
         SLAB(IV,NWP2) = EDGE                                           ISO00357
   50 CONTINUE                                                          ISO00358
      DO  60 IW=2,NWP2                                                  ISO00359
         SLAB(NVP2,IW) = EDGE                                           ISO00360
   60 CONTINUE                                                          ISO00361
C                                                                       ISO00362
C GO THRU 3-D ARRAY IN U DIRECTION.  IUEW=IU EITHER WAY.                ISO00363
C PICK IU BASED ON EYEU.                                                ISO00364
C                                                                       ISO00365
      DO  90 IUEW=1,NU                                                  ISO00366
         IU = IUEW                                                      ISO00367
         IF (EYEU .GT. AVEU) IU = NU+1-IUEW                             ISO00368
         U = IU                                                         ISO00369
C                                                                       ISO00370
C LOAD THIS SLICE OF T INTO SLAB.                                       ISO00371
C                                                                       ISO00372
         DO  80 IV=1,NV                                                 ISO00373
            DO  70 IW=1,NW                                              ISO00374
               SLAB(IV+1,IW+1) = T(IU,IV,IW)                            ISO00375
   70       CONTINUE                                                    ISO00376
   80    CONTINUE                                                       ISO00377
C                                                                       ISO00378
C CONTOUR THIS SLAB.                                                    ISO00379
C                                                                       ISO00380
         CALL STCNTR (SLAB,NUVWP2,NVP2,NWP2,TVAL)                       ISO00381
C                                                                       ISO00382
C CONSTRUCT VISIBILITY ARRAY.                                           ISO00383
C                                                                       ISO00384
         CALL FILLIN                                                    ISO00385
   90 CONTINUE                                                          ISO00386
C                                                                       ISO00387
C SLICES PERPENDICULAR TO V. U W SLICES. T OF CONSTANT V.               ISO00388
C                                                                       ISO00389
  100 IF (MOD(NFLAG/2,2) .EQ. 0) GO TO 160                              ISO00390
      CALL ZEROSC                                                       ISO00391
      ISLBT = 0                                                         ISO00392
C                                                                       ISO00393
C BOUND UPPER AND RIGHT EDGE OF SLAB.                                   ISO00394
C                                                                       ISO00395
      DO 110 IU=2,NUP2                                                  ISO00396
         SLAB(IU,NWP2) = EDGE                                           ISO00397
  110 CONTINUE                                                          ISO00398
      DO 120 IW=2,NWP2                                                  ISO00399
         SLAB(NUP2,IW) = EDGE                                           ISO00400
  120 CONTINUE                                                          ISO00401
C                                                                       ISO00402
C GO THRU T IN V DIRECTION.  IVEW=IV EITHER WAY.                        ISO00403
C                                                                       ISO00404
      DO 150 IVEW=1,NV                                                  ISO00405
         IV = IVEW                                                      ISO00406
         IF (EYEV .GT. AVEV) IV = NV+1-IVEW                             ISO00407
         V = IV                                                         ISO00408
C                                                                       ISO00409
C LOAD THIS SLICE OF T INTO SLAB.                                       ISO00410
C                                                                       ISO00411
         DO 140 IU=1,NU                                                 ISO00412
            DO 130 IW=1,NW                                              ISO00413
               SLAB(IU+1,IW+1) = T(IU,IV,IW)                            ISO00414
  130       CONTINUE                                                    ISO00415
  140    CONTINUE                                                       ISO00416
C                                                                       ISO00417
C CONTOUR THIS SLAB.                                                    ISO00418
C                                                                       ISO00419
         CALL STCNTR (SLAB,NUVWP2,NUP2,NWP2,TVAL)                       ISO00420
C                                                                       ISO00421
C CONSTRUCT VISIBILITY ARRAY.                                           ISO00422
C                                                                       ISO00423
         CALL FILLIN                                                    ISO00424
  150 CONTINUE                                                          ISO00425
C                                                                       ISO00426
C SLICES PERPENDICULAR TO W. U V SLICES. T OF CONSTANT W.               ISO00427
C                                                                       ISO00428
  160 IF (MOD(NFLAG,2) .EQ. 0) GO TO 220                                ISO00429
      CALL ZEROSC                                                       ISO00430
C                                                                       ISO00431
      ISLBT = 1                                                         ISO00432
C                                                                       ISO00433
C BOUND UPPER AND RIGHT EDGE OF SLAB.                                   ISO00434
C                                                                       ISO00435
      DO 170 IU=2,NUP2                                                  ISO00436
         SLAB(IU,NVP2) = EDGE                                           ISO00437
  170 CONTINUE                                                          ISO00438
      DO 180 IV=2,NVP2                                                  ISO00439
         SLAB(NUP2,IV) = EDGE                                           ISO00440
  180 CONTINUE                                                          ISO00441
C                                                                       ISO00442
C GO THRU T IN W DIRECTION.                                             ISO00443
C                                                                       ISO00444
      DO 210 IWEW=1,NW                                                  ISO00445
         IW = IWEW                                                      ISO00446
         IF (EYEW .GT. AVEW) IW = NW+1-IWEW                             ISO00447
         W = IW                                                         ISO00448
C                                                                       ISO00449
C LOAD THIS SLICE OF T INTO SLAB.                                       ISO00450
C                                                                       ISO00451
         DO 200 IU=1,NU                                                 ISO00452
            DO 190 IV=1,NV                                              ISO00453
               SLAB(IU+1,IV+1) = T(IU,IV,IW)                            ISO00454
  190       CONTINUE                                                    ISO00455
  200    CONTINUE                                                       ISO00456
C                                                                       ISO00457
C CONTOUR THIS SLAB.                                                    ISO00458
C                                                                       ISO00459
         CALL STCNTR (SLAB,NUVWP2,NUP2,NVP2,TVAL)                       ISO00460
C                                                                       ISO00461
C CONSTRUCT VISIBILITY ARRAY.                                           ISO00462
C                                                                       ISO00463
         CALL FILLIN                                                    ISO00464
  210 CONTINUE                                                          ISO00465
C                                                                       ISO00466
C DRAW REFERENCE PLANE EDGES AND W AXIS.                                ISO00467
C                                                                       ISO00468
  220 IF (IREF .EQ. 0) RETURN                                           ISO00469
      CALL TRN32I (SU1,SV1,SW1,XT,YT,DUM,2)                             ISO00470
      IF (EYEV .LT. SV1) GO TO 240                                      ISO00471
      CALL FRSTC (IFIX(XT),IFIX(YT),1)                                  ISO00472
      DO 230 IU=2,NU                                                    ISO00473
         CALL TRN32I (SU(FLOAT(IU)),SV1,SW1,XT,YT,DUM,2)                ISO00474
         CALL FRSTC (IFIX(XT),IFIX(YT),2)                               ISO00475
  230 CONTINUE                                                          ISO00476
      GO TO 250                                                         ISO00477
  240 CALL PLOTIT (IFIX(XT),IFIX(YT),0)                                 ISO00478
      CALL TRN32I (SUNU,SV1,SW1,XT,YT,DUM,2)                            ISO00479
      CALL PLOTIT (IFIX(XT),IFIX(YT),1)                                 ISO00480
  250 IF (EYEU .GT. SUNU) GO TO 270                                     ISO00481
      CALL FRSTC (IFIX(XT),IFIX(YT),1)                                  ISO00482
      DO 260 IV=2,NV                                                    ISO00483
         CALL TRN32I (SUNU,SV(FLOAT(IV)),SW1,XT,YT,DUM,2)               ISO00484
         CALL FRSTC (IFIX(XT),IFIX(YT),2)                               ISO00485
  260 CONTINUE                                                          ISO00486
      GO TO 280                                                         ISO00487
  270 CALL PLOTIT (IFIX(XT),IFIX(YT),0)                                 ISO00488
      CALL TRN32I (SUNU,SVNV,SW1,XT,YT,DUM,2)                           ISO00489
      CALL PLOTIT (IFIX(XT),IFIX(YT),1)                                 ISO00490
  280 IF (EYEV .GT. SVNV) GO TO 300                                     ISO00491
      CALL FRSTC (IFIX(XT),IFIX(YT),1)                                  ISO00492
      DO 290 IUOW=2,NU                                                  ISO00493
         CALL TRN32I (SU(FLOAT(NU-IUOW+1)),SVNV,SW1,XT,YT,DUM,2)        ISO00494
         CALL FRSTC (IFIX(XT),IFIX(YT),2)                               ISO00495
  290 CONTINUE                                                          ISO00496
      GO TO 310                                                         ISO00497
  300 CALL PLOTIT (IFIX(XT),IFIX(YT),0)                                 ISO00498
      CALL TRN32I (SU1,SVNV,SW1,XT,YT,DUM,2)                            ISO00499
      CALL PLOTIT (IFIX(XT),IFIX(YT),1)                                 ISO00500
  310 IF (EYEU .LT. SU1) GO TO 330                                      ISO00501
      CALL FRSTC (IFIX(XT),IFIX(YT),1)                                  ISO00502
      DO 320 IVOW=2,NV                                                  ISO00503
         CALL TRN32I (SU1,SV(FLOAT(NV-IVOW+1)),SW1,XT,YT,DUM,2)         ISO00504
         CALL FRSTC (IFIX(XT),IFIX(YT),2)                               ISO00505
  320 CONTINUE                                                          ISO00506
      GO TO 340                                                         ISO00507
  330 CALL PLOTIT (IFIX(XT),IFIX(YT),0)                                 ISO00508
      CALL TRN32I (SU1,SV1,SW1,XT,YT,DUM,2)                             ISO00509
      CALL PLOTIT (IFIX(XT),IFIX(YT),1)                                 ISO00510
  340 IF (EYEU.LE.SU1 .OR. EYEV.LE.SV1) GO TO 360                       ISO00511
      CALL FRSTC (IFIX(XT),IFIX(YT),1)                                  ISO00512
      DO 350 IW=2,NW                                                    ISO00513
         CALL TRN32I (SU1,SV1,SW(FLOAT(IW)),XT,YT,DUM,2)                ISO00514
         CALL FRSTC (IFIX(XT),IFIX(YT),2)                               ISO00515
  350 CONTINUE                                                          ISO00516
      RETURN                                                            ISO00517
  360 CALL PLOTIT (IFIX(XT),IFIX(YT),0)                                 ISO00518
      CALL TRN32I (SU1,SV1,SWNW,XT,YT,DUM,2)                            ISO00519
      CALL PLOTIT (IFIX(XT),IFIX(YT),1)                                 ISO00520
      RETURN                                                            ISO00521
      END                                                               ISO00522
      SUBROUTINE KURV1S (N,X,Y,SLOP1,SLOPN,XP,YP,TEMP,S,SIGMA,ISLPSW)   ISO01278
C                                                                       ISO01279
C                                                                       ISO01280
C DIMENSION OF           X(N),Y(N),XP(N),YP(N),TEMP(N)                  ISO01281
C ARGUMENTS                                                             ISO01282
C                                                                       ISO01283
C LATEST REVISION        JUNE 1979                                      ISO01284
C                                                                       ISO01285
C PURPOSE                KURV1S DETERMINES THE PARAMETERS NECESSARY TO  ISO01286
C                        COMPUTE A SPLINE UNDER TENSION PASSING THROUGH ISO01287
C                        A SEQUENCE OF PAIRS                            ISO01288
C                        (X(1),Y(1)),...,(X(N),Y(N)) IN THE PLANE.      ISO01289
C                        THE SLOPES AT THE TWO ENDS OF THE CURVE MAY BE ISO01290
C                        SPECIFIED OR OMITTED.  FOR ACTUAL COMPUTATION  ISO01291
C                        OF POINTS ON THE CURVE IT IS NECESSARY TO CALL ISO01292
C                        THE SUBROUTINE KURV2.                          ISO01293
C                                                                       ISO01294
C USAGE                  CALL KURV1S(N,X,Y,SLP1,SLPN,XP,YP,TEMP,S,SIGMA)ISO01295
C                                                                       ISO01296
C ARGUMENTS                                                             ISO01297
C                                                                       ISO01298
C ON INPUT               N                                              ISO01299
C                          IS THE NUMBER OF POINTS TO BE INTERPOLATED   ISO01300
C                          (N .GE. 2).                                  ISO01301
C                                                                       ISO01302
C                        X                                              ISO01303
C                          IS AN ARRAY CONTAINING THE N X-COORDINATES   ISO01304
C                          OF THE POINTS.                               ISO01305
C                                                                       ISO01306
C                        Y                                              ISO01307
C                          IS AN ARRAY CONTAINING THE N Y-COORDINATES   ISO01308
C                          OF THE POINTS.                               ISO01309
C                                                                       ISO01310
C                        SLOP1 AND SLOPN                                ISO01311
C                          CONTAIN THE DESIRED VALUES FOR THE SLOPE OF  ISO01312
C                          THE CURVE AT (X(1),Y(1)) AND (X(N),Y(N)),    ISO01313
C                          RESPECTIVELY.  THESE QUANTITIES ARE IN       ISO01314
C                          DEGREES AND MEASURED COUNTER-CLOCKWISE       ISO01315
C                          FROM THE POSITIVE X-AXIS.  THE POSITIVE      ISO01316
C                          FROM THE POSITIVE X-AXIS.  IF ISLPSW IS NON- ISO01317
C                          ZERO, ONE OR BOTH OF SLP1 AND SLPN MAY BE    ISO01318
C                          DETERMINED INTERNALLY BY KURV1S.             ISO01319
C                                                                       ISO01320
C                        XP AND YP                                      ISO01321
C                          ARE ARRAYS OF LENGTH AT LEAST N.             ISO01322
C                                                                       ISO01323
C                        TEMP                                           ISO01324
C                          IS AN ARRAY OF LENGTH AT LEAST N WHICH IS    ISO01325
C                          USED FOR SCRATCH STORAGE.                    ISO01326
C                                                                       ISO01327
C                        SIGMA                                          ISO01328
C                          CONTAINS THE TENSION FACTOR.  THIS IS        ISO01329
C                          NON-ZERO AND INDICATES THE CURVINESS DESIRED.ISO01330
C                          IF ABS(SIGMA) IS VERY LARGE (E.G., 50.) THE  ISO01331
C                          RESULTING CURVE IS VERY NEARLY A POLYGONAL   ISO01332
C                          LINE.  A STANDARD VALUE FOR SIGMA IS ABOUT 2.ISO01333
C                                                                       ISO01334
C                        ISLPSW                                         ISO01335
C                          IS AN INTEGER INDICATING WHICH END SLOPES    ISO01336
C                          HAVE BEEN USER PROVIDED AND WHICH MUST BE    ISO01337
C                          COMPUTED BY KURV1S.  FOR ISLPSW              ISO01338
C                            = 0  INDICATES BOTH SLOPES ARE PROVIDED,   ISO01339
C                            = 1  ONLY SLOP1 IS PROVIDED,               ISO01340
C                            = 2  ONLY SLOPN IS PROVIDED,               ISO01341
C                            = 3  NEITHER SLOP1 NOR SLOPN IS PROVIDED.  ISO01342
C                            = 4  NEITHER SLOP1 NOR SLOPN IS PROVIDED,  ISO01343
C                                 BUT SLOP1=SLOPN.  IN THIS CASE X(1)=  ISO01344
C                                 X(N), Y(1)=Y(N) AND N.GE.3.           ISO01345
C ON OUTPUT              XP AND YP                                      ISO01346
C                          CONTAIN INFORMATION ABOUT THE CURVATURE OF   ISO01347
C                          THE CURVE AT THE GIVEN NODES.                ISO01348
C                                                                       ISO01349
C                        S                                              ISO01350
C                          CONTAINS THE POLYGONAL ARCLENGTH OF THE      ISO01351
C                          CURVE.                                       ISO01352
C                                                                       ISO01353
C                        N, X, Y, SLP1, SLPN, SIGMA AND ISLPSW ARE      ISO01354
C                        UNCHANGED.                                     ISO01355
C                                                                       ISO01356
C ENTRY POINTS           KURV1S                                         ISO01357
C                                                                       ISO01358
C SPECIAL CONDITIONS     NONE                                           ISO01359
C                                                                       ISO01360
C COMMON BLOCKS          NONE                                           ISO01361
C                                                                       ISO01362
C I/O                    NONE                                           ISO01363
C                                                                       ISO01364
C PRECISION              SINGLE                                         ISO01365
C                                                                       ISO01366
C REQUIRED ULIB          NONE                                           ISO01367
C ROUTINES                                                              ISO01368
C                                                                       ISO01369
C SPECIALIST             RUSSELL K. REW, NCAR, BOULDER, COLORADO  80302 ISO01370
C                                                                       ISO01371
C LANGUAGE               FORTRAN                                        ISO01372
C                                                                       ISO01373
C HISTORY                ORIGINALLY WRITTEN BY A. K. CLINE, MARCH 1972. ISO01374
C                                                                       ISO01375
C                                                                       ISO01376
C                                                                       ISO01377
C                                                                       ISO01378
      INTEGER         N                                                 ISO01379
      REAL            X(N)       ,Y(N)       ,XP(N)      ,YP(N)      ,  ISO01380
     1                TEMP(N)    ,S          ,SIGMA                     ISO01381
C                                                                       ISO01382
      NN = N                                                            ISO01383
      JSLPSW = ISLPSW                                                   ISO01384
      SLP1 = SLOP1                                                      ISO01385
      SLPN = SLOPN                                                      ISO01386
      DEGRAD = 3.1415926535897932/180.                                  ISO01387
      NM1 = NN-1                                                        ISO01388
      NP1 = NN+1                                                        ISO01389
      DELX1 = X(2)-X(1)                                                 ISO01390
      DELY1 = Y(2)-Y(1)                                                 ISO01391
      DELS1 = SQRT(DELX1*DELX1+DELY1*DELY1)                             ISO01392
      DX1 = DELX1/DELS1                                                 ISO01393
      DY1 = DELY1/DELS1                                                 ISO01394
C                                                                       ISO01395
C DETERMINE SLOPES IF NECESSARY                                         ISO01396
C                                                                       ISO01397
      IF (JSLPSW .NE. 0) GO TO  70                                      ISO01398
   10 SLPP1 = SLP1*DEGRAD                                               ISO01399
      SLPPN = SLPN*DEGRAD                                               ISO01400
C                                                                       ISO01401
C SET UP RIGHT HAND SIDES OF TRIDIAGONAL LINEAR SYSTEM FOR XP           ISO01402
C AND YP                                                                ISO01403
C                                                                       ISO01404
      XP(1) = DX1-COS(SLPP1)                                            ISO01405
      YP(1) = DY1-SIN(SLPP1)                                            ISO01406
      TEMP(1) = DELS1                                                   ISO01407
      SS = DELS1                                                        ISO01408
      IF (NN .EQ. 2) GO TO  30                                          ISO01409
      DO  20 I=2,NM1                                                    ISO01410
         DELX2 = X(I+1)-X(I)                                            ISO01411
         DELY2 = Y(I+1)-Y(I)                                            ISO01412
         DELS2 = SQRT(DELX2*DELX2+DELY2*DELY2)                          ISO01413
         DX2 = DELX2/DELS2                                              ISO01414
         DY2 = DELY2/DELS2                                              ISO01415
         XP(I) = DX2-DX1                                                ISO01416
         YP(I) = DY2-DY1                                                ISO01417
         TEMP(I) = DELS2                                                ISO01418
         DELX1 = DELX2                                                  ISO01419
         DELY1 = DELY2                                                  ISO01420
         DELS1 = DELS2                                                  ISO01421
         DX1 = DX2                                                      ISO01422
         DY1 = DY2                                                      ISO01423
C                                                                       ISO01424
C ACCUMULATE POLYGONAL ARCLENGTH                                        ISO01425
C                                                                       ISO01426
         SS = SS+DELS1                                                  ISO01427
   20 CONTINUE                                                          ISO01428
   30 XP(NN) = COS(SLPPN)-DX1                                           ISO01429
      YP(NN) = SIN(SLPPN)-DY1                                           ISO01430
C                                                                       ISO01431
C DENORMALIZE TENSION FACTOR                                            ISO01432
C                                                                       ISO01433
      SIGMAP = ABS(SIGMA)*FLOAT(NN-1)/SS                                ISO01434
C                                                                       ISO01435
C PERFORM FORWARD ELIMINATION ON TRIDIAGONAL SYSTEM                     ISO01436
C                                                                       ISO01437
      S = SS                                                            ISO01438
      DELS = SIGMAP*TEMP(1)                                             ISO01439
      EXPS = EXP(DELS)                                                  ISO01440
      SINHS = .5*(EXPS-1./EXPS)                                         ISO01441
      SINHIN = 1./(TEMP(1)*SINHS)                                       ISO01442
      DIAG1 = SINHIN*(DELS*.5*(EXPS+1./EXPS)-SINHS)                     ISO01443
      DIAGIN = 1./DIAG1                                                 ISO01444
      XP(1) = DIAGIN*XP(1)                                              ISO01445
      YP(1) = DIAGIN*YP(1)                                              ISO01446
      SPDIAG = SINHIN*(SINHS-DELS)                                      ISO01447
      TEMP(1) = DIAGIN*SPDIAG                                           ISO01448
      IF (NN .EQ. 2) GO TO  50                                          ISO01449
      DO  40 I=2,NM1                                                    ISO01450
         DELS = SIGMAP*TEMP(I)                                          ISO01451
         EXPS = EXP(DELS)                                               ISO01452
         SINHS = .5*(EXPS-1./EXPS)                                      ISO01453
         SINHIN = 1./(TEMP(I)*SINHS)                                    ISO01454
         DIAG2 = SINHIN*(DELS*(.5*(EXPS+1./EXPS))-SINHS)                ISO01455
         DIAGIN = 1./(DIAG1+DIAG2-SPDIAG*TEMP(I-1))                     ISO01456
         XP(I) = DIAGIN*(XP(I)-SPDIAG*XP(I-1))                          ISO01457
         YP(I) = DIAGIN*(YP(I)-SPDIAG*YP(I-1))                          ISO01458
         SPDIAG = SINHIN*(SINHS-DELS)                                   ISO01459
         TEMP(I) = DIAGIN*SPDIAG                                        ISO01460
         DIAG1 = DIAG2                                                  ISO01461
   40 CONTINUE                                                          ISO01462
   50 DIAGIN = 1./(DIAG1-SPDIAG*TEMP(NM1))                              ISO01463
      XP(NN) = DIAGIN*(XP(NN)-SPDIAG*XP(NM1))                           ISO01464
      YP(NN) = DIAGIN*(YP(NN)-SPDIAG*YP(NM1))                           ISO01465
C                                                                       ISO01466
C PERFORM BACK SUBSTITUTION                                             ISO01467
C                                                                       ISO01468
      DO  60 I=2,NN                                                     ISO01469
         IBAK = NP1-I                                                   ISO01470
         XP(IBAK) = XP(IBAK)-TEMP(IBAK)*XP(IBAK+1)                      ISO01471
         YP(IBAK) = YP(IBAK)-TEMP(IBAK)*YP(IBAK+1)                      ISO01472
   60 CONTINUE                                                          ISO01473
      RETURN                                                            ISO01474
   70 IF (NN .EQ. 2) GO TO 100                                          ISO01475
C                                                                       ISO01476
C IF NO SLOPES ARE GIVEN, USE SECOND ORDER INTERPOLATION ON             ISO01477
C INPUT DATA FOR SLOPES AT ENDPOINTS                                    ISO01478
C                                                                       ISO01479
      IF (JSLPSW .EQ. 4) GO TO  90                                      ISO01480
      IF (JSLPSW .EQ. 2) GO TO  80                                      ISO01481
      DELNM1 = SQRT((X(NN-2)-X(NM1))**2+(Y(NN-2)-Y(NM1))**2)            ISO01482
      DELN = SQRT((X(NM1)-X(NN))**2+(Y(NM1)-Y(NN))**2)                  ISO01483
      DELNN = DELNM1+DELN                                               ISO01484
      C1 = (DELNN+DELN)/DELNN/DELN                                      ISO01485
      C2 = -DELNN/DELN/DELNM1                                           ISO01486
      C3 = DELN/DELNN/DELNM1                                            ISO01487
      SX = C3*X(NN-2)+C2*X(NM1)+C1*X(NN)                                ISO01488
      SY = C3*Y(NN-2)+C2*Y(NM1)+C1*Y(NN)                                ISO01489
C                                                                       ISO01490
      SLPN = ATAN2(SY,SX)/DEGRAD                                        ISO01491
   80 IF (JSLPSW .EQ. 1) GO TO  10                                      ISO01492
      DELS2 = SQRT((X(3)-X(2))**2+(Y(3)-Y(2))**2)                       ISO01493
      DELS12 = DELS1+DELS2                                              ISO01494
      C1 = -(DELS12+DELS1)/DELS12/DELS1                                 ISO01495
      C2 = DELS12/DELS1/DELS2                                           ISO01496
      C3 = -DELS1/DELS12/DELS2                                          ISO01497
      SX = C1*X(1)+C2*X(2)+C3*X(3)                                      ISO01498
      SY = C1*Y(1)+C2*Y(2)+C3*Y(3)                                      ISO01499
C                                                                       ISO01500
      SLP1 = ATAN2(SY,SX)/DEGRAD                                        ISO01501
      GO TO  10                                                         ISO01502
   90 DELN = SQRT((X(NM1)-X(NN))**2+(Y(NM1)-Y(NN))**2)                  ISO01503
      DELNN = DELS1+DELN                                                ISO01504
      C1 = -DELS1/DELN/DELNN                                            ISO01505
      C2 = (DELS1-DELN)/DELS1/DELN                                      ISO01506
      C3 = DELN/DELNN/DELS1                                             ISO01507
      SX = C1*X(NM1)+C2*X(1)+C3*X(2)                                    ISO01508
      SY = C1*Y(NM1)+C2*Y(1)+C3*Y(2)                                    ISO01509
      IF (SX.EQ.0. .AND. SY.EQ.0.) SX = 1.                              ISO01510
      SLP1 = ATAN2(SY,SX)/DEGRAD                                        ISO01511
      SLPN = SLP1                                                       ISO01512
      GO TO  10                                                         ISO01513
C                                                                       ISO01514
C IF ONLY TWO POINTS AND NO SLOPES ARE GIVEN, USE STRAIGHT              ISO01515
C LINE SEGMENT FOR CURVE                                                ISO01516
C                                                                       ISO01517
  100 IF (JSLPSW .NE. 3) GO TO 110                                      ISO01518
      XP(1) = 0.                                                        ISO01519
      XP(2) = 0.                                                        ISO01520
      YP(1) = 0.                                                        ISO01521
      YP(2) = 0.                                                        ISO01522
C                                                                       ISO01523
      SLP1 = ATAN2(Y(2)-Y(1),X(2)-X(1))/DEGRAD                          ISO01524
      SLPN = SLP1                                                       ISO01525
      RETURN                                                            ISO01526
C                                                                       ISO01527
  110 IF (JSLPSW .EQ. 2)                                                ISO01528
     1    SLP1 = ATAN2(Y(2)-Y(1)-SLPN*(X(2)-X(1)),                      ISO01529
     2                                X(2)-X(1)-SLPN*(Y(2)-Y(1)))/DEGRADISO01530
C                                                                       ISO01531
      IF (JSLPSW .EQ. 1)                                                ISO01532
     1    SLPN = ATAN2(Y(2)-Y(1)-SLP1*(X(2)-X(1)),                      ISO01533
     2                                X(2)-X(1)-SLP1*(Y(2)-Y(1)))/DEGRADISO01534
      GO TO  10                                                         ISO01535
      END                                                               ISO01536
      SUBROUTINE KURV2S (T,XS,YS,N,X,Y,XP,YP,S,SIGMA,NSLPSW,SLP)        ISO01537
C                                                                       ISO01538
C                                                                       ISO01539
C                                                                       ISO01540
C DIMENSION OF           X(N),Y(N),XP(N),YP(N)                          ISO01541
C ARGUMENTS                                                             ISO01542
C                                                                       ISO01543
C LATEST REVISION        JUNE 1979                                      ISO01544
C                                                                       ISO01545
C PURPOSE                KURV2 PERFORMS THE MAPPING OF POINTS IN THE    ISO01546
C                        INTERVAL (0.,1.) ONTO A CURVE IN THE PLANE.    ISO01547
C                        THE SUBROUTINE KURV1 SHOULD BE CALLED EARLIER  ISO01548
C                        TO DETERMINE CERTAIN NECESSARY PARAMETERS.     ISO01549
C                        THE RESULTING CURVE HAS A PARAMETRIC           ISO01550
C                        REPRESENTATION BOTH OF WHOSE COMPONENTS ARE    ISO01551
C                        SPLINES UNDER TENSION AND FUNCTIONS OF THE     ISO01552
C                        POLYGONAL ARCLENGTH PARAMETER.                 ISO01553
C                                                                       ISO01554
C ACCESS CARDS           *FORTRAN,S=ULIB,N=KURV                         ISO01555
C                                                                       ISO01556
C                                                                       ISO01557
C USAGE                  CALL KURV2S (T,XS,YS,N,X,Y,XP,YP,S,SIGMA)      ISO01558
C                                                                       ISO01559
C ARGUMENTS                                                             ISO01560
C                                                                       ISO01561
C ON INPUT               T                                              ISO01562
C                          CONTAINS A REAL VALUE OF ABSOLUTE VALUE LESS ISO01563
C                          THAN OR EQUAL TO 1. TO BE MAPPED TO A POINT  ISO01564
C                          ON THE CURVE.  THE SIGN OF T IS IGNORED AND  ISO01565
C                          THE INTERVAL (0.,1.) IS MAPPED ONTO THE      ISO01566
C                          ENTIRE CURVE.  IF T IS NEGATIVE, THIS        ISO01567
C                          INDICATES THAT THE SUBROUTINE HAS BEEN CALLEDISO01568
C                          PREVIOUSLY (WITH ALL OTHER INPUT VARIABLES   ISO01569
C                          UNALTERED) AND THAT THIS VALUE OF T EXCEEDS  ISO01570
C                          THE PREVIOUS VALUE IN ABSOLUTE VALUE.  WITH  ISO01571
C                          SUCH INFORMATION THE SUBROUTINE IS ABLE TO   ISO01572
C                          MAP THE POINT MUCH MORE RAPIDLY.  THUS IF THEISO01573
C                          USER SEEKS TO MAP A SEQUENCE OF POINTS ONTO  ISO01574
C                          THE SAME CURVE, EFFICIENCY IS GAINED BY      ISO01575
C                          ORDERING THE VALUES INCREASING IN MAGNITUDE  ISO01576
C                          AND SETTING THE SIGNS OF ALL BUT THE FIRST   ISO01577
C                          NEGATIVE.                                    ISO01578
C                                                                       ISO01579
C                        N                                              ISO01580
C                          CONTAINS THE NUMBER OF POINTS WHICH WERE     ISO01581
C                          INTERPOLATED TO DETERMINE THE CURVE.         ISO01582
C                                                                       ISO01583
C                        X AND Y                                        ISO01584
C                          ARRAYS CONTAINING THE X- AND Y-COORDINATES   ISO01585
C                          OF THE INTERPOLATED POINTS.                  ISO01586
C                                                                       ISO01587
C                        XP AND YP                                      ISO01588
C                          ARE THE ARRAYS OUTPUT FROM KURV1 CONTAINING  ISO01589
C                          CURVATURE INFORMATION.                       ISO01590
C                                                                       ISO01591
C                        S                                              ISO01592
C                          CONTAINS THE POLYGONAL ARCLENGTH OF THE      ISO01593
C                          CURVE.                                       ISO01594
C                                                                       ISO01595
C                        SIGMA                                          ISO01596
C                          CONTAINS THE TENSION FACTOR (ITS SIGN IS     ISO01597
C                          IGNORED).                                    ISO01598
C                                                                       ISO01599
C                        NSLPSW                                         ISO01600
C                          IS AN INTEGER SWITCH WHICH TURNS ON OR OFF   ISO01601
C                          THE CALCULATION OF SLP                       ISO01602
C                          NSLPSW                                       ISO01603
C                                 = 0 INDICATES THAT SLP WILL NOT BE    ISO01604
C                                     CALCULATED                        ISO01605
C                                 = 1 SLP WILL BE CALCULATED            ISO01606
C                                                                       ISO01607
C                        THE PARAMETERS N, X, Y, XP, YP, S AND SIGMA    ISO01608
C                        SHOULD BE INPUT UNALTERED FROM THE OUTPUT OF   ISO01609
C                        KURV1S.                                        ISO01610
C                                                                       ISO01611
C ON OUTPUT              XS AND YS                                      ISO01612
C                          CONTAIN THE X- AND Y-COORDINATES OF THE IMAGEISO01613
C                          POINT ON THE CURVE.                          ISO01614
C                                                                       ISO01615
C                        SLP                                            ISO01616
C                          CONTAINS THE SLOPE OF THE CURVE IN DEGREES ATISO01617
C                          THIS POINT.                                  ISO01618
C                                                                       ISO01619
C                        T, N, X, Y, XP, YP, S AND SIGMA ARE UNALTERED. ISO01620
C                                                                       ISO01621
C ENTRY POINTS           KURV2S                                         ISO01622
C                                                                       ISO01623
C SPECIAL CONDITIONS     NONE                                           ISO01624
C                                                                       ISO01625
C COMMON BLOCKS          NONE                                           ISO01626
C                                                                       ISO01627
C I/O                    NONE                                           ISO01628
C                                                                       ISO01629
C PRECISION              SINGLE                                         ISO01630
C                                                                       ISO01631
C REQUIRED ULIB          NONE                                           ISO01632
C ROUTINES                                                              ISO01633
C                                                                       ISO01634
C SPECIALIST             RUSSELL K. REW, NCAR, BOULDER, COLORADO  80302 ISO01635
C                                                                       ISO01636
C LANGUAGE               FORTRAN                                        ISO01637
C                                                                       ISO01638
C HISTORY                ORIGINALLY WRITTEN BY A. K. CLINE, MARCH 1972. ISO01639
C                                                                       ISO01640
C                                                                       ISO01641
C                                                                       ISO01642
C                                                                       ISO01643
      INTEGER         N                                                 ISO01644
      REAL            T          ,XS         ,YS         ,X(N)       ,  ISO01645
     1                Y(N)       ,XP(N)      ,YP(N)      ,S          ,  ISO01646
     2                SIGMA      ,SLP                                   ISO01647
C                                                                       ISO01648
C DENORMALIZE SIGMA                                                     ISO01649
C                                                                       ISO01650
      SIGMAP = ABS(SIGMA)*FLOAT(N-1)/S                                  ISO01651
C                                                                       ISO01652
C STRETCH UNIT INTERVAL INTO ARCLENGTH DISTANCE                         ISO01653
C                                                                       ISO01654
      TN = ABS(T*S)                                                     ISO01655
C                                                                       ISO01656
C FOR NEGATIVE T START SEARCH WHERE PREVIOUSLY TERMINATED,              ISO01657
C OTHERWISE START FROM BEGINNING                                        ISO01658
C                                                                       ISO01659
      IF (T .LT. 0.) GO TO  10                                          ISO01660
      DEGRAD = 3.1415926535897932/180.                                  ISO01661
      I1 = 2                                                            ISO01662
      XS = X(1)                                                         ISO01663
      YS = Y(1)                                                         ISO01664
      SUM = 0.                                                          ISO01665
      IF (T .LT. 0.) RETURN                                             ISO01666
C                                                                       ISO01667
C DETERMINE INTO WHICH SEGMENT TN IS MAPPED                             ISO01668
C                                                                       ISO01669
   10 DO  30 I=I1,N                                                     ISO01670
         DELX = X(I)-X(I-1)                                             ISO01671
         DELY = Y(I)-Y(I-1)                                             ISO01672
         DELS = SQRT(DELX*DELX+DELY*DELY)                               ISO01673
         IF (SUM+DELS-TN)  20, 40, 40                                   ISO01674
   20    SUM = SUM+DELS                                                 ISO01675
   30 CONTINUE                                                          ISO01676
C                                                                       ISO01677
C IF ABS(T) IS GREATER THAN 1., RETURN TERMINAL POINT ON                ISO01678
C CURVE                                                                 ISO01679
C                                                                       ISO01680
      XS = X(N)                                                         ISO01681
      YS = Y(N)                                                         ISO01682
      RETURN                                                            ISO01683
C                                                                       ISO01684
C SET UP AND PERFORM INTERPOLATION                                      ISO01685
C                                                                       ISO01686
   40 DEL1 = TN-SUM                                                     ISO01687
      DEL2 = DELS-DEL1                                                  ISO01688
      EXPS1 = EXP(SIGMAP*DEL1)                                          ISO01689
      SINHD1 = .5*(EXPS1-1./EXPS1)                                      ISO01690
      EXPS2 = EXP(SIGMAP*DEL2)                                          ISO01691
      SINHD2 = .5*(EXPS2-1./EXPS2)                                      ISO01692
      EXPS = EXPS1*EXPS2                                                ISO01693
      SINHS = .5*(EXPS-1./EXPS)                                         ISO01694
      XS = (XP(I)*SINHD1+XP(I-1)*SINHD2)/SINHS+                         ISO01695
     1     ((X(I)-XP(I))*DEL1+(X(I-1)-XP(I-1))*DEL2)/DELS               ISO01696
      YS = (YP(I)*SINHD1+YP(I-1)*SINHD2)/SINHS+                         ISO01697
     1     ((Y(I)-YP(I))*DEL1+(Y(I-1)-YP(I-1))*DEL2)/DELS               ISO01698
      I1 = I                                                            ISO01699
      IF (NSLPSW .EQ. 0) RETURN                                         ISO01700
      COSHD1 = .5*(EXPS1+1./EXPS1)*SIGMAP                               ISO01701
      COSHD2 = .5*(EXPS2+1./EXPS2)*SIGMAP                               ISO01702
      XT = (XP(I)*COSHD1-XP(I-1)*COSHD2)/SINHS+                         ISO01703
     1     ((X(I)-XP(I))-(X(I-1)-XP(I-1)))/DELS                         ISO01704
      YT = (YP(I)*COSHD1-YP(I-1)*COSHD2)/SINHS+                         ISO01705
     1     ((Y(I)-YP(I))-(Y(I-1)-YP(I-1)))/DELS                         ISO01706
      SLP = ATAN2(YT,XT)/DEGRAD                                         ISO01707
      RETURN                                                            ISO01708
      END                                                               ISO01709
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
      SUBROUTINE MMASK                                                  ISO02126
C                                                                       ISO02127
C  MAKE THE MACHINE DEPENDENT MASKS USED IN THE CONTOUR DRAWING         ISO02128
C  AND SHADING ALGORITHMS                                               ISO02129
C                                                                       ISO02130
      COMMON /ISOSR5/ NBPW       ,MASK(16)   ,GENDON                    ISO02131
      LOGICAL         GENDON                                            ISO02132
      COMMON /ISOSR7/ IENTRY     ,IONES                                 ISO02133
      COMMON /ISOSR8/ NMASK(16)  ,IXOLD      ,IYOLD      ,IBTOLD     ,  ISO02134
     1                HBFLAG     ,IOSLSN     ,LRLX       ,IFSX       ,  ISO02135
     2                IFSY       ,FIRST      ,IYDIR      ,IHX        ,  ISO02136
     3                IHB        ,IHS        ,IHV        ,IVOLD      ,  ISO02137
     4                IVAL       ,IHRX       ,YCHANG     ,ITPD       ,  ISO02138
     5                IHF                                               ISO02139
      COMMON /ISOSR9/ BIG        ,IXBIT                                 ISO02140
      LOGICAL         YCHANG     ,HBFLAG     ,FIRST      ,IHF           ISO02141
      GENDON = .TRUE.                                                   ISO02142
      NBPW = 16                                                         ISO02143
C                                                                       ISO02144
C  BIGGEST REAL NUMBER FOR STELLAR                                      ISO02145
C                                                                       ISO02146
      BIG = 2.E0**(31)                                                  ISO02147
C                                                                       ISO02148
C  MASKS TO SELECT A SPECIFIC BIT                                       ISO02149
C                                                                       ISO02150
      DO  10 K=1,NBPW                                                   ISO02151
         MASK(K) = ISHFT(1,K-1)                                         ISO02152
   10 CONTINUE                                                          ISO02153
C                                                                       ISO02154
C  GENERATE  THE BIT PATTERN 177777 OCTAL                               ISO02155
C                                                                       ISO02156
      ITEMP1 = 0                                                        ISO02157
      ITEMP = MASK(NBPW)                                                ISO02158
      IST = NBPW-1                                                      ISO02159
      DO  20 K=1,IST                                                    ISO02160
         ITEMP1 = IOR(ITEMP,ISHFT(ITEMP1,-1))                           ISO02161
   20 CONTINUE                                                          ISO02162
      MFIX = IOR(ITEMP1,1)                                              ISO02163
C                                                                       ISO02164
C  MASKS TO CLEAR A SPECIFIC BIT                                        ISO02165
C                                                                       ISO02166
      DO  30 K=1,NBPW                                                   ISO02167
         NMASK(K) = IAND(ITEMP1,MFIX)                                   ISO02168
         ITEMP1 = IOR(ISHFT(ITEMP1,1),1)                                ISO02169
   30 CONTINUE                                                          ISO02170
      IONES = MFIX                                                      ISO02171
      RETURN                                                            ISO02172
C                                                                       ISO02173
C REVISION HISTORY---                                                   ISO02174
C                                                                       ISO02175
C JANUARY 1978     DELETED REFERENCES TO THE  *COSY  CARDS AND          ISO02176
C                  ADDED REVISION HISTORY                               ISO02177
C JANUARY 1979     NEW SHADING ALGORITHM                                ISO02178
C MARCH 1979       MADE CODE MACHINE INDEPENDENT AND CONFORM            ISO02179
C                  TO 66 FORTRAN STANDARD                               ISO02180
C JUNE 1979        THIS VERSION PLACED ON ULIB.                         ISO02181
C SEPTEMBER 1979   FIXED PROBLEM IN EZISOS DEALING WITH                 ISO02182
C                  DETERMINATION OF VISIBILITY OF W PLANE.              ISO02183
C DECEMBER 1979    FIXED PROBLEM WITH PEN DOWN ON CONTOUR               ISO02184
C                  INITIALIZATION IN SUBROUTINE FRSTC                   ISO02185
C MARCH            CHANGED ROUTINE NAMES  TRN32I  AND  DRAW  TO         ISO02186
C                  TRN32I  AND  DRAWI  TO BE CONSISTENT WITH THE        ISO02187
C                  USAGE OF THE NEW ROUTINE  PWRZI.                     ISO02188
C JUNE  1980       FIXED PROBLEM WITH ZERO INDEX COMPUTATION IN         ISO02189
C                  SUBROUTINE FRSTC.  ADDED INPUT PARAMETER             ISO02190
C                  DIMENSION STATEMENT MISSING IN EZISOS.               ISO02191
C                  FIXED ERROR IN COMPUTATION OF ARCCOSINE              ISO02192
C                  IN EZISOS AND TRN32I.                                ISO02193
C-----------------------------------------------------------------------ISO02194
C                                                                       ISO02195
      END                                                               ISO02196
      SUBROUTINE PLOTIT(IX,IY,IPEN)
C
C    IPEN = 0 PEN UP
C    IPEN = 1 PEN DOWN
C
      DATA CONV /32767.0/
C
      IF (IPEN .EQ. 0) IP = 3
      IF (IPEN .EQ. 1) IP = 2
C
      X = 10.0*FLOAT(IX)/CONV
      Y = 10.0*FLOAT(IY)/CONV
C
      CALL PLOT(X,Y,IP)
C
      RETURN
      END
      SUBROUTINE SET3D (EYE,ULO,UHI,VLO,VHI,WLO,WHI)                    ISO00613
      COMMON /TEMPR/  RZERO                                             ISO00614
C                                                                       ISO00615
      DIMENSION       EYE(3)                                            ISO00616
C                                                                       ISO00617
      COMMON /ISOSR3/ ISCALE     ,XMIN       ,XMAX       ,YMIN       ,  ISO00618
     1                YMAX       ,BIGD       ,R0                        ISO00619
      COMMON /PWRZ1I/ UUMIN      ,UUMAX      ,VVMIN      ,VVMAX      ,  ISO00620
     1                WWMIN      ,WWMAX      ,DELCRT     ,EYEU       ,  ISO00621
     2                EYEV       ,EYEW                                  ISO00622
C                                                                       ISO00623
C                                                                       ISO00624
      AVE(A,B) = (A+B)*.5                                               ISO00625
C                                                                       ISO00626
C A.S.F. FOR SCALING                                                    ISO00627
C                                                                       ISO00628
      SU(UTEMP) = UTEMP                                                 ISO00629
      SV(VTEMP) = VTEMP                                                 ISO00630
      SW(WTEMP) = WTEMP                                                 ISO00631
C                                                                       ISO00632
C CONSTANTS FOR PWRZ                                                    ISO00633
C                                                                       ISO00634
      UUMIN = ULO                                                       ISO00635
      UUMAX = UHI                                                       ISO00636
      VVMIN = VLO                                                       ISO00637
      VVMAX = VHI                                                       ISO00638
      WWMIN = WLO                                                       ISO00639
      WWMAX = WHI                                                       ISO00640
      EYEU = EYE(1)                                                     ISO00641
      EYEV = EYE(2)                                                     ISO00642
      EYEW = EYE(3)                                                     ISO00643
C                                                                       ISO00644
C FIND CORNERS IN 2-SPACE FOR 3-SPACE BOX CONTAINING OBJECT             ISO00645
C                                                                       ISO00646
      ISCALE = 0                                                        ISO00647
      ATU = AVE(SU(UUMIN),SU(UUMAX))                                    ISO00648
      ATV = AVE(SV(VVMIN),SV(VVMAX))                                    ISO00649
      ATW = AVE(SW(WWMIN),SW(WWMAX))                                    ISO00650
      BIGD = 0.                                                         ISO00651
      IF (RZERO .LE. 0.) GO TO  10                                      ISO00652
C                                                                       ISO00653
C RELETIVE SIZE FEATURE IN USE.                                         ISO00654
C GENERATE EYE POSITION THAT MAKES BOX HAVE MAXIMUM PROJECTED SIZE.     ISO00655
C                                                                       ISO00656
      ALPHA = -(VVMIN-ATV)/(UUMIN-ATU)                                  ISO00657
      VVEYE = -RZERO/SQRT(1.+ALPHA*ALPHA)                               ISO00658
      UUEYE = VVEYE*ALPHA                                               ISO00659
      VVEYE = VVEYE+ATV                                                 ISO00660
      UUEYE = UUEYE+ATU                                                 ISO00661
      WWEYE = ATW                                                       ISO00662
      CALL TRN32I (ATU,ATV,ATW,UUEYE,VVEYE,WWEYE,1)                     ISO00663
      CALL TRN32I (UUMIN,VVMIN,ATW,XMIN,DUMM,DUMM,2)                    ISO00664
      CALL TRN32I (UUMAX,VVMIN,WWMIN,DUMM,YMIN,DUMM,2)                  ISO00665
      CALL TRN32I (UUMAX,VVMAX,ATW,XMAX,DUMM,DUMM,2)                    ISO00666
      CALL TRN32I (UUMAX,VVMIN,WWMAX,DUMM,YMAX,DUMM,2)                  ISO00667
      BIGD = SQRT((UUMAX-UUMIN)**2+(VVMAX-VVMIN)**2+(WWMAX-WWMIN)**2)*.5ISO00668
      R0 = RZERO                                                        ISO00669
      GO TO  20                                                         ISO00670
   10 CALL TRN32I (ATU,ATV,ATW,EYE(1),EYE(2),EYE(3),1)                  ISO00671
      CALL TRN32I (SU(UUMIN),SV(VVMIN),SW(WWMIN),X1,Y1,DUM,2)           ISO00672
      CALL TRN32I (SU(UUMIN),SV(VVMIN),SW(WWMAX),X2,Y2,DUM,2)           ISO00673
      CALL TRN32I (SU(UUMIN),SV(VVMAX),SW(WWMIN),X3,Y3,DUM,2)           ISO00674
      CALL TRN32I (SU(UUMIN),SV(VVMAX),SW(WWMAX),X4,Y4,DUM,2)           ISO00675
      CALL TRN32I (SU(UUMAX),SV(VVMIN),SW(WWMIN),X5,Y5,DUM,2)           ISO00676
      CALL TRN32I (SU(UUMAX),SV(VVMIN),SW(WWMAX),X6,Y6,DUM,2)           ISO00677
      CALL TRN32I (SU(UUMAX),SV(VVMAX),SW(WWMIN),X7,Y7,DUM,2)           ISO00678
      CALL TRN32I (SU(UUMAX),SV(VVMAX),SW(WWMAX),X8,Y8,DUM,2)           ISO00679
      XMIN = AMIN1(X1,X2,X3,X4,X5,X6,X7,X8)                             ISO00680
      XMAX = AMAX1(X1,X2,X3,X4,X5,X6,X7,X8)                             ISO00681
      YMIN = AMIN1(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)                             ISO00682
      YMAX = AMAX1(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)                             ISO00683
C                                                                       ISO00684
C ADD RIGHT AMOUNT TO KEEP PICTURE SQUARE                               ISO00685
C                                                                       ISO00686
   20 WIDTH = XMAX-XMIN                                                 ISO00687
      HIGHT = YMAX-YMIN                                                 ISO00688
      DIF = .5*(WIDTH-HIGHT)                                            ISO00689
      IF (DIF)  30, 50, 40                                              ISO00690
   30 XMIN = XMIN+DIF                                                   ISO00691
      XMAX = XMAX-DIF                                                   ISO00692
      GO TO  50                                                         ISO00693
   40 YMIN = YMIN-DIF                                                   ISO00694
      YMAX = YMAX+DIF                                                   ISO00695
   50 ISCALE = 1                                                        ISO00696
      CALL TRN32I (ATU,ATV,ATW,EYE(1),EYE(2),EYE(3),1)                  ISO00697
      RETURN                                                            ISO00698
      END                                                               ISO00699
      SUBROUTINE STCNTR (Z,L,M,N,CONV)                                  ISO00845
C                                                                       ISO00846
      DIMENSION       Z(L,N)                                            ISO00847
C                                                                       ISO00848
C THIS ROUTINE FINDS THE BEGINNINGS OF ALL CONTOUR LINES AT LEVEL CONV. ISO00849
C FIRST THE EDGES ARE SEARCHED FOR LINES INTERSECTING THE EDGE (OPEN    ISO00850
C LINES) THEN THE INTERIOR IS SEARCHED FOR LINES WHICH DO NOT INTERSECT ISO00851
C THE EDGE (CLOSED LINES).  BEGINNINGS ARE STORED IN IR TO PREVENT RE-  ISO00852
C TRACING OF LINES.  IF IR IS FILLED, THE SEARCH IS STOPPED FOR THIS    ISO00853
C CONV.                                                                 ISO00854
C                                                                       ISO00855
      COMMON /ISOSR6/ IX         ,IY         ,IDX        ,IDY        ,  ISO00856
     1                IS         ,ISS        ,NP         ,CV         ,  ISO00857
     2                INX(8)     ,INY(8)     ,IR(500)    ,NR            ISO00858
      COMMON /ISOSR7/ IENTRY     ,IONES                                 ISO00859
      COMMON /ISOSR9/ BIG        ,IXBIT                                 ISO00860
C                                                                       ISO00861
C  PACK X AND Y                                                         ISO00862
C                                                                       ISO00863
      IPXY(I1,J1) = ISHFT(I1,IXBIT)+J1                                  ISO00864
C                                                                       ISO00865
      IENTRY = 0                                                        ISO00866
      NP = 0                                                            ISO00867
      CV = CONV                                                         ISO00868
C                                                                       ISO00869
C THE FOLLOWING CODE SHOULD BE RE-ENABLED IF THIS ROUTINE IS USED FOR   ISO00870
C GENERAL CONTOURING                                                    ISO00871
C                                                                       ISO00872
C     ISS=0                                                             ISO00873
C     DO 2 IP1=2,M                                                      ISO00874
C     I=IP1-1                                                           ISO00875
C     IF(Z(I,1).GE.CV.OR.Z(IP1,1).LT.CV) GO TO 1                        ISO00876
C     IX=IP1                                                            ISO00877
C     IY=1                                                              ISO00878
C     IDX=-1                                                            ISO00879
C     IDY=0                                                             ISO00880
C     IS=1                                                              ISO00881
C     CALL DRLINE(Z,L,M,N)                                              ISO00882
C   1 IF(Z(IP1,N).GE.CV.OR.Z(I,N).LT.CV) GO TO 2                        ISO00883
C     IX=I                                                              ISO00884
C     IY=N                                                              ISO00885
C     IDX=1                                                             ISO00886
C     IDY=0                                                             ISO00887
C     IS=5                                                              ISO00888
C     CALL DRLINE(Z,L,M,N)                                              ISO00889
C   2 CONTINUE                                                          ISO00890
C     DO 4 JP1=2,N                                                      ISO00891
C     J=JP1-1                                                           ISO00892
C     IF(Z(M,J).GE.CV.OR.Z(M,JP1).LT.CV) GO TO 3                        ISO00893
C     IX=M                                                              ISO00894
C     IY=JP1                                                            ISO00895
C     IDX=0                                                             ISO00896
C     IDY=-1                                                            ISO00897
C     IS=7                                                              ISO00898
C     CALL DRLINE(Z,L,M,N)                                              ISO00899
C   3 IF(Z(1,JP1).GE.CV.OR.Z(1,J).LT.CV) GO TO 4                        ISO00900
C     IX=1                                                              ISO00901
C     IY=J                                                              ISO00902
C     IDX=0                                                             ISO00903
C     IDY=1                                                             ISO00904
C     IS=3                                                              ISO00905
C     CALL DRLINE(Z,L,M,N)                                              ISO00906
C   4 CONTINUE                                                          ISO00907
C                                                                       ISO00908
      ISS = 1                                                           ISO00909
      DO  40 JP1=3,N                                                    ISO00910
         J = JP1-1                                                      ISO00911
      DO  30 IP1=2,M                                                    ISO00912
            I = IP1-1                                                   ISO00913
            IF (Z(I,J).GE.CV .OR. Z(IP1,J).LT.CV) GO TO  30             ISO00914
            IXY = IPXY(IP1,J)                                           ISO00915
            IF (NP .EQ. 0) GO TO  20                                    ISO00916
            DO  10 K=1,NP                                               ISO00917
               IF (IR(K) .EQ. IXY) GO TO  30                            ISO00918
   10       CONTINUE                                                    ISO00919
   20       NP = NP+1                                                   ISO00920
            IF (NP .GT. NR) RETURN                                      ISO00921
            IR(NP) = IXY                                                ISO00922
            IX = IP1                                                    ISO00923
            IY = J                                                      ISO00924
            IDX = -1                                                    ISO00925
            IDY = 0                                                     ISO00926
            IS = 1                                                      ISO00927
            CALL DRCNTR (Z,L,M,N)                                       ISO00928
   30    CONTINUE                                                       ISO00929
   40 CONTINUE                                                          ISO00930
      RETURN                                                            ISO00931
      END                                                               ISO00932
      SUBROUTINE TR32 (X,Y,MX,MY)                                       ISO01045
C                                                                       ISO01046
      COMMON /ISOSR1/ ISLBT      ,U          ,V          ,W             ISO01047
C                                                                       ISO01048
C A.S.F. FOR SCALING                                                    ISO01049
C                                                                       ISO01050
      SU(UTEMP) = UTEMP                                                 ISO01051
      SV(VTEMP) = VTEMP                                                 ISO01052
      SW(WTEMP) = WTEMP                                                 ISO01053
C                                                                       ISO01054
      XX = X                                                            ISO01055
      YY = Y                                                            ISO01056
      IF (ISLBT)  10, 20, 30                                            ISO01057
   10 CALL TRN32I (SU(U),SV(XX-1.),SW(YY-1.),XT,YT,DUM,2)               ISO01058
      GO TO  40                                                         ISO01059
   20 CALL TRN32I (SU(XX-1.),SV(V),SW(YY-1.),XT,YT,DUM,2)               ISO01060
      GO TO  40                                                         ISO01061
   30 CALL TRN32I (SU(XX-1.),SV(YY-1.),SW(W),XT,YT,DUM,2)               ISO01062
   40 MX = XT                                                           ISO01063
      MY = YT                                                           ISO01064
      RETURN                                                            ISO01065
      END                                                               ISO01066
      SUBROUTINE TRN32I (U,V,W,XT,YT,ZT,IENT)                           ISO00700
C                                                                       ISO00701
C THIS ROUTINE IMPLEMENTS THE 3-SPACE TO 2-SPACE TRANSFOR-              ISO00702
C MATION BY KUBER, SZABO AND GIULIERI, THE PERSPECTIVE                  ISO00703
C REPRESENTATION OF FUNCTIONS OF TWO VARIABLES. J. ACM 15,              ISO00704
C 2, 193-204,1968.                                                      ISO00705
C ARGUMENTS FOR SET                                                     ISO00706
C U,V,W    ARE THE 3-SPACE COORDINATES OF THE INTERSECTION              ISO00707
C          OF THE LINE OF SIGHT AND THE IMAGE PLANE.  THIS              ISO00708
C          POINT CAN BE THOUGHT OF AS THE POINT LOOKED AT.              ISO00709
C XT,YT,ZT ARE THE 3-SPACE COORDINATES OF THE EYE POSITION.             ISO00710
C                                                                       ISO00711
C TRN32 ARGUMENTS                                                       ISO00712
C U,V,W    ARE THE 3-SPACE COORDINATES OF A POINT TO BE                 ISO00713
C          TRANSFORMED.                                                 ISO00714
C XT,YT    THE RESULTS OF THE 3-SPACE TO 2-SPACE TRANSFOR-              ISO00715
C          MATION.  WHEN ISCALE=0, XT AND YT ANR IN THE SAME            ISO00716
C          UNITS AS U,V, AND W.  WHEN ISCALE'0, XT AND YT               ISO00717
C          ARE IN PLOTTER COORDINATES.                                  ISO00718
C ZT       NOT USED.                                                    ISO00719
C                                                                       ISO00720
      COMMON /PWRZ1I/ UUMIN      ,UUMAX      ,VVMIN      ,VVMAX      ,  ISO00721
     1                WWMIN      ,WWMAX      ,DELCRT     ,EYEU       ,  ISO00722
     2                EYEV       ,EYEW                                  ISO00723
      COMMON /ISOSR3/ ISCALE     ,XMIN       ,XMAX       ,YMIN       ,  ISO00724
     1                YMAX       ,BIGD       ,R0                        ISO00725
C                                                                       ISO00726
C RANGE OF PLOTTER COORDINATES                                          ISO00727
C                                                                       ISO00728
C                                                                       ISO00729
C  WARNING                                                              ISO00730
C  IF PLOTTER MAXIMUM VALUE RANGES (IN X OR Y DIRECTION) FALL BELOW     ISO00731
C  101, THEN CHANGES MUST BE MADE IN SUBROUTINE FRSTC. THE REQUIRED     ISO00732
C  CHANGES ARE MARKED BY WARNING COMMENTS IN FRSTC.                     ISO00733
      DATA NLX,NBY,NRX,NTY/10,10,32760,32760/                           ISO00734
        DATA PI/3.141592/                                               ISO00735
C                                                                       ISO00736
C STORE THE PARAMETERS OF THE SET CALL FOR USE                          ISO00737
C WITH THE TRANSLATE CALL                                               ISO00738
C                                                                       ISO00739
C DECIDE IF SET OR TRANSLATE CALL                                       ISO00740
C                                                                       ISO00741
      IF (IENT .NE. 1) GO TO  50                                        ISO00742
      AU = U                                                            ISO00743
      AV = V                                                            ISO00744
      AW = W                                                            ISO00745
      EU = XT                                                           ISO00746
      EV = YT                                                           ISO00747
      EW = ZT                                                           ISO00748
C                                                                       ISO00749
C                                                                       ISO00750
C                                                                       ISO00751
C                                                                       ISO00752
C                                                                       ISO00753
      DU = AU-EU                                                        ISO00754
      DV = AV-EV                                                        ISO00755
      DW = AW-EW                                                        ISO00756
      D = SQRT(DU*DU+DV*DV+DW*DW)                                       ISO00757
      COSAL = DU/D                                                      ISO00758
      COSBE = DV/D                                                      ISO00759
      COSGA = DW/D                                                      ISO00760
C                                                                       ISO00761
C  COMPUTE THE ARCCOSINE                                                ISO00762
C                                                                       ISO00763
C     AL = ATAN(ABS(SQRT(1.-COSAL*COSAL)/COSAL))                        ISO00764
C     IF (COSAL .LE. 0.) AL = PI-AL                                     ISO00765
      IF(ABS(COSBE).LT.1.E-5)THEN
      BE=PI*0.5
      ELSE
      BE = ATAN(ABS(SQRT(1.-COSBE*COSBE)/COSBE))                        ISO00766
      IF (COSBE .LE. 0.) BE = PI-BE                                     ISO00767
      ENDIF
      IF(ABS(COSGA).LE.1.E-5)THEN
      GA=PI*0.5
      ELSE
      GA = ATAN(ABS(SQRT(1.-COSGA*COSGA)/COSGA))                        ISO00768
      IF (COSGA .LE. 0.) GA = PI-GA                                     ISO00769
      ENDIF
      SINGA = SIN(GA)                                                   ISO00770
C                                                                       ISO00771
C THE 3-SPACE POINT LOOKED AT IS TRANSFORMED INTO (0,0) OF              ISO00772
C THE 2-SPACE.  THE 3-SPACE W AXIS IS TRANSFORMED INTO THE              ISO00773
C 2-SPACE Y AXIS.  IF THE LINE OF SIGHT IS CLOSE TO PARALLEL            ISO00774
C TO THE 3-SPACE W AXIS, THE 3-SPACE V AXIS IS CHOSEN (IN-              ISO00775
C STEAD OF THE 3-SPACE W AXIS) TO BE TRANSFORMED INTO THE               ISO00776
C 2-SPACE Y AXIS.                                                       ISO00777
C                                                                       ISO00778
      ASSIGN  90 TO JDONE                                               ISO00779
      IF (ISCALE)  10, 30, 10                                           ISO00780
   10 X0 = XMIN                                                         ISO00781
      Y0 = YMIN                                                         ISO00782
      X1 = NLX                                                          ISO00783
      Y1 = NBY                                                          ISO00784
      X2 = NRX-NLX                                                      ISO00785
      Y2 = NTY-NBY                                                      ISO00786
      X3 = X2/(XMAX-XMIN)                                               ISO00787
      Y3 = Y2/(YMAX-YMIN)                                               ISO00788
      X4 = NRX                                                          ISO00789
      Y4 = NTY                                                          ISO00790
      FACT = 1.                                                         ISO00791
      IF (BIGD .LE. 0.) GO TO  20                                       ISO00792
      X0 = -BIGD                                                        ISO00793
      Y0 = -BIGD                                                        ISO00794
      X3 = X2/(2.*BIGD)                                                 ISO00795
      Y3 = Y2/(2.*BIGD)                                                 ISO00796
      FACT = R0/D                                                       ISO00797
   20 DELCRT = X2                                                       ISO00798
      ASSIGN  80 TO JDONE                                               ISO00799
   30 IF (SINGA .LT. 0.0001) GO TO  40                                  ISO00800
      R = 1./SINGA                                                      ISO00801
      ASSIGN  70 TO JUMP                                                ISO00802
      RETURN                                                            ISO00803
   40 SINBE = SIN(BE)                                                   ISO00804
      R = 1./SINBE                                                      ISO00805
      ASSIGN  60 TO JUMP                                                ISO00806
      RETURN                                                            ISO00807
C                                                                       ISO00808
C********************  ENTRY TRN32  ************************            ISO00809
C     ENTRY TRN32 (U,V,W,XT,YT,ZT)                                      ISO00810
C                                                                       ISO00811
   50 UU = U                                                            ISO00812
      VV = V                                                            ISO00813
      WW = W                                                            ISO00814
      Q = D/((UU-EU)*COSAL+(VV-EV)*COSBE+(WW-EW)*COSGA)                 ISO00815
      GO TO JUMP,( 60, 70)                                              ISO00816
   60 UU = ((EW+Q*(WW-EW)-AW)*COSAL-(EU+Q*(UU-EU)-AU)*COSGA)*R          ISO00817
      VV = (EV+Q*(VV-EV)-AV)*R                                          ISO00818
      GO TO JDONE,( 80, 90)                                             ISO00819
   70 UU = ((EU+Q*(UU-EU)-AU)*COSBE-(EV+Q*(VV-EV)-AV)*COSAL)*R          ISO00820
      VV = (EW+Q*(WW-EW)-AW)*R                                          ISO00821
      GO TO JDONE,( 80, 90)                                             ISO00822
   80 XT = AMIN1(X4,AMAX1(X1,X1+X3*(FACT*UU-X0)))                       ISO00823
      YT = AMIN1(Y4,AMAX1(Y1,Y1+Y3*(FACT*VV-Y0)))                       ISO00824
      RETURN                                                            ISO00825
   90 XT = UU                                                           ISO00826
      YT = VV                                                           ISO00827
      RETURN                                                            ISO00828
      END                                                               ISO00829
      SUBROUTINE ULIBER (IERR,MESSG,NMESSG)                             ULI00001
C                                                                       ULI00002
C DIMENSION OF           MESSG(1+(NMESSG-1)/NCPWD), WHERE NCPWD IS THE  ULI00003
C ARGUMENTS              NUMBER OF CHARACTERS THAT CAN BE STORED IN ONE ULI00004
C                        INTEGER WORD.                                  ULI00005
C                                                                       ULI00006
C LATEST REVISION        APRIL 1977.                                    ULI00007
C                                                                       ULI00008
C PURPOSE                PRINTS AN ERROR MESSAGE.                       ULI00009
C                                                                       ULI00010
C USAGE                  CALL ULIBER (IERR,MESSG,NMESSG)                ULI00011
C                                                                       ULI00012
C ARGUMENTS                                                             ULI00013
C                                                                       ULI00014
C ON INPUT               IERR                                           ULI00015
C                          ERROR NUMBER.  IF IERR .LT. 32, THE ERROR IS ULI00016
C                          CONSIDERED TO BE NON-FATAL AND ULIBER RETURNSULI00017
C                          AFTER PRINTING THE ERROR MESSAGE.  IF IERR   ULI00018
C                          .GE. 32 THE ERROR IS CONSIDERED FATAL, AND   ULI00019
C                          ULIBER STOPS AFTER PRINTING THE MESSAGE.     ULI00020
C                        MESSG                                          ULI00021
C                          THE ERROR MESSAGE TO BE PRINTED.             ULI00022
C                        NMESSG                                         ULI00023
C                          THE LENGTH OF THE MESSAGE, IN CHARACTERS.    ULI00024
C                                                                       ULI00025
      DIMENSION MESSG(1)                                                ULI00026
      INTEGER ERUNIT,PRUNIT                                             ULI00027
C                                                                       ULI00028
C ERUNIT SHOULD BE SET TO THE LOCAL UNIT NUMBER FOR ERROR MESSAGES.     ULI00029
C                                                                       ULI00030
      DATA ERUNIT/0/,PRUNIT/6/                                          ULI00031
C                                                                       ULI00032
      IF (ERUNIT .EQ. 0) WRITE(PRUNIT,1000)                             ULI00033
C                                                                       ULI00034
 1000 FORMAT(50H1ULIBER, A SUBROUTINE TO PRINT ERROR MESSAGES, HAS/     ULI00035
     1       50H BEEN CALLED, BUT NO LOCAL IMPLEMENTATION OF      /     ULI00036
     2       50H ULIBER HAS BEEN PROVIDED.  TO PROPERLY IMPLEMENT /     ULI00037
     3       50H THIS SUBROUTINE FOR THE LOCAL MACHINE AND        /     ULI00038
     4       50H ENVIRONMENT, PLEASE SEE THE COMMENTS IN ULIBER.  /)    ULI00039
C                                                                       ULI00040
C REPLACE THE IF STATEMENT AND FORMAT STATEMENT ABOVE WITH THE FOLLOWINGULI00041
C CODE WHERE $NCPWD SHOULD BE REPLACED WITH THE NUMBER OF CHARACTERS    ULI00042
C THAT CAN BE STORED IN ONE INTEGER WORD, AND $N SHOULD BE REPLACED BY  ULI00043
C THE VALUE OF 1+(79/$NCPWD).                                           ULI00044
C                                                                       ULI00045
C     MM=1+(NMESSG-1)/$NCPWD                                            ULI00046
C     WRITE(ERUNIT,1000) IERR,(MESSG(I),I=1,MM)                         ULI00047
C1000 FORMAT(18H ****ERROR NUMBER ,I5,25H,ERROR MESSAGE FOLLOWS.../     ULI00048
C    1       ($N A $NCPWD ))                                            ULI00049
C                                                                       ULI00050
   10 IF (IERR .GE. 32) STOP                                            ULI00051
      RETURN                                                            ULI00052
      END                                                               ULI00053
      SUBROUTINE ZEROSC                                                 ISO00830
C                                                                       ISO00831
      COMMON /ISOSR2/ LX         ,NX         ,NY         ,ISCR(8,128),  ISO00832
     1                ISCA(8,128)                                       ISO00833
C                                                                       ISO00834
C ZERO BOTH SCREEN MODELS.                                              ISO00835
C                                                                       ISO00836
      DO  20 I=1,LX                                                     ISO00837
         DO  10 J=1,NY                                                  ISO00838
            ISCR(I,J) = 0                                               ISO00839
            ISCA(I,J) = 0                                               ISO00840
   10    CONTINUE                                                       ISO00841
   20 CONTINUE                                                          ISO00842
      RETURN                                                            ISO00843
      END                                                               ISO00844

