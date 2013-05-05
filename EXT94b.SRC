      PROGRAM EXTREME
C
C    VERSION 94 - Revision B
C    
C    EXTREME SEARCHES 3-D MOLECULAR SPACE FOR CRITICAL POINTS OF THE 
C    FOLLOWING SCALAR FIELDS:
C
C    1) THE ELECTRON DENSITY DISTRIBUTION (RHO).
C    2) THE LAPLACIAN OF THE ELECTRON DENSITY DISTRIBUTION (DEL**2(RHO)).
C    3) THE LAGRANGIAN KINETIC ENERGY DENSITY DISTRIBUTION (KEG).
C    4) THE HAMILTONIAN KINETIC ENERGY DENSITY DISTRIBUTION (KEK).
C    5) THE NUCLEAR POTENTIAL DISTRIBUTION (Vnuc).
C    6) THE POTENTIAL ENERGY DENSITY - TRACE OF THE STRESS TENSOR. (V).
C
C    AS INPUT, AN AIMPAC WAVEFUNCTION FILE IS REQUIRED.  THE REST OF THE
C    INPUT IS SUPPLIED BY THE USER INTERACTIVELY.
C
C    BASICALLY, THE CRITICAL POINTS ARE FOUND USING A NEWTON-RAPHSON
C    SEARCH, WITH INITIAL GUESSES SUPPLIED BY THE USER, OR AUTOMATICALLY.
C
C    EXTREME CAN HANDLE S, P(3), D(6) and F(10) TYPE GAUSSIAN BASIS 
C    FUNCTIONS.
C
C    ORIGINAL EXTREME (SADDLE) DEVELOPED BY VARIOUS MEMBERS OF
C    RICHARD BADER'S RESEARCH GROUP - MOST EXTENSIVELY BY KEITH E
C    LAIDIG APRIL 1989.
C
C    HEAVILY MODIFIED, EXTENDED AND CLEANED UP:  ADDITION OF THE FIELDS
C    DEL**2(RHO), G, K, NUCLEAR POTENTIAL AND POTENTIAL ENERGY DENSITY.
C    INCORPORATION OF AUTOMATED SEARCHING OPTIONS, AND MANY OTHER 
C    NICETIES.
C    TAK 8/94
C
C    Questions and Suggestions Should be Directed to:
C    R.F.W. Bader:  bader@mcmail.cis.mcmaster.ca
C    or
C    T.A. Keith:    keith@babbage.chemistry.mcmaster.ca
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*8  ATNAM
      CHARACTER*80 JOBTTL,WFNTTL
      CHARACTER*4 FWFN /'.wfn'/, FWLP /'.crt'/
      CHARACTER*8 BLANK /'        '/
      CHARACTER*40 WFN,WLP
      Parameter(MaxAtm=100,MaxOff=200000,MaxCrt=500)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /ANG/ ANGLE(3,MaxAtm,MaxAtm)
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
      COMMON /STRING/ WFNTTL,JOBTTL,ATNAM(MaxAtm)
      COMMON /UNITS/  INPT,IOUT,IWFN,IWLP
      COMMON /OPTIONS/ ICut,IPrint,Eps,EpsNuc,Dmp,DmpNuc
      DIMENSION JI1(MaxCrt),JI2(MaxCrt),X(MaxCrt),Y(MaxCrt),Z(MaxCrt),
     $  XYZ(3),EVSAVE(MaxCrt,3)
      DATA Zero/0.d0/,One/1.d0/,Two/2.d0/,Three/3.d0/,Ten/10.0d0/
     $  Four/4.d0/,pt5/0.5d0/,degen/1.d-8/
1000  FORMAT(' ENTER THE TITLE OF THIS RUN ')
1010  FORMAT(A80)
1020  FORMAT(' EXTREME ')
1030  Format(' NUCLEAR COORDINATES')
1040  FORMAT(6X,A8,4X,3F15.8)
1050  FORMAT(' WHICH SCALAR FIELD TO ANALYZE:  ',/,' ELECTRON',
     +' DENSITY, RHO (1)',/,' LAPLACIAN OF RHO, DEL**2(RHO) (2)',/,
     +  ' LAGRANGIAN KINETIC ENERGY DENSITY, G (3)',/,' HAMILTONIAN ',
     +  'KINETIC ENERGY DENSITY, K (4)',/,
     +  ' NUCLEAR POTENTIAL, Vnuc (5)',/,
     + ' POTENTIAL ENERGY DENSITY, V (6)',/,' Enter (1-6):  ', $)
1060  FORMAT(/,' TRY AGAIN ')
1070  FORMAT(' ANALYZING RHO, THE ELECTRON DENSITY')
1080  FORMAT(' ANALYZING DEL**2(RHO), THE LAPLACIAN OF THE ELECTRON',
     +' DENSITY')
1090  FORMAT(' ANALYZING G, THE LAGRANGIAN KINETIC ENERGY DENSITY')
1100  FORMAT(' ANALYZING K, THE HAMILTONIAN KINETIC ENERGY DENSITY')
1110  FORMAT(' ANALYZING Vnuc, THE NUCLEAR POTENTIAL FIELD')
1115  FORMAT(' ANALYZING V, THE POTENTIAL ENERGY DENSITY ')
1120  FORMAT(/,' COORDS or NUCLEAR (0), BOND (1), RING (2),',
     + ' CAGE (3), ANGLE (4), POINT (5),',/,
     +' MEGA (6), OPTIONS (7), LIST CURRENT CPs (8), STOP (9)',/,
     +' Enter (1-9)  ',$)
1121  Format(/,' Nuclear Coordinates (1) or Arbitrary',
     +' Coordinates (0) ? ',$)
1122  Format(/' Enter Nucleus Number for Nuclear Critical Point',
     +' Search  ',$)
1130  FORMAT(/,' INPUT STARTING COORDINATES FOR NEWTON-RAPHSON',
     $' SEARCH ',/,' ',$)
1140  FORMAT(/,' STARTING COORDINATES FOR NEWTON-RAPHSON SEARCH ',
     $ /,1X,3F15.8)
1141  FORMAT(/,' NUCLEUS ',I4,' STARTING COORDINATES FOR',
     $' NEWTON-RAPHSON SEARCH ',
     $ /,1X,3F15.8)
1150  FORMAT(/,' INPUT NUMBERS OF TWO ATOMS TO SEARCH BETWEEN: ',$)
1160  FORMAT(/,' ONE OR MORE OF THE ATOMS DOES NOT EXIST:',
     + ' INPUT NUMBERS AGAIN ',$)
1170  FORMAT(/,' START HALFWAY BETWEEN ATOMS (1) OR ANOTHER FRACTION',
     $' OF DISTANCE (0) ? ',$)
1180  FORMAT(/,' FRACTION (eg. 0.5) OF DISTANCE FROM ATOM ',I4,
     $' TO ATOM ', I4,' TO START SEARCH ? ',/,$)
1190  FORMAT(/,' SEARCHING BETWEEN ATOMS ',2I4)
1200  FORMAT(/,' INPUT NUMBERS OF THREE ATOMS TO SEARCH BETWEEN: ',$)
1210  FORMAT(/,' SEARCHING BETWEEN ATOMS ',3I4)
1220  FORMAT(/,' INPUT NUMBERS OF FOUR ATOMS TO SEARCH BETWEEN: ',$)
1230  FORMAT(/,' SEARCHING BETWEEN ATOMS ',4I4)
1240  FORMAT(/,' INPUT COORDINATES FOR PROPERTY EVALUATION ',/,' ')
1250  FORMAT(/,' COORDINATES FOR PROPERTY EVALUATION ',/,' ',3F9.6)
1260  FORMAT(/,' Using Avg of Corresponding Nuclear Coordinates as',
     +' Initial Guess,',
     +/,' Search:  For Nuclear Critical Points (1); Between',
     +' Atom Pairs (2);',/,' Between Atom Triads (3);',
     +' Between Atom Quads (4); All of the Above (5);',/,
     +' or None of the Above (0) ? ',$)
1265  FORMAT(/,' Search Along Whole Line Between Two Atoms (1)',/,
     +' Search Another Line (2)',/,
     +' None (0)',/,' Enter (0-2)  ',$)
1266  Format(/,' Enter Two Atom Numbers and the Number of Starting',
     + ' Points  ',$)
1267  Format(/,' Enter Coordinates of One End of Desired Line',
     + /,$)
1268  Format(/,' Enter Coordinates of Other End of ',
     + ' Line and Number of Starting Points',/,$)
1270  FORMAT(/,' Search Over one or all atom-centered spheres ?',
     +' (1=yes/0=no) ',$)
1280  FORMAT(/,' What Atom Number ?  Type 0 for all atoms ',$)
1290  FORMAT(/,' Input Spherical Grid Search Parameters: ',/,
     $' NPhi,NTheta,RadMin,RadMax and NradPt (eg. 6 6 0.5 2.0 10) ',
     $ /,$)
1295  FORMAT(' Change Accuracy to be Used in Evaluating',
     $' Functions (1)',/,
     $' Change Print Options (2)',/,' Change Gradient Threshold',
     $' For Calling',
     $' a Point a Critical Point (3) ',/,
     $' Change Damping Factor for Newton-Raphson Search (4)',/,
     $' Finished With Options (0)',/,
     $' Enter (0-4)  ',$)
1296  FORMAT(' Enter Integer N (Default is 20) to Change the Function',
     $/,' Accuracy to 10**(-N)  ',$)
1297  FORMAT(' Print Output to Screen ? (1=yes/0=no(default))  ',$)
1298  FORMAT(' Function Accuracy Set to ',1PE16.8)
1299  FORMAT(' Print Option Set to ',I3)
1300  FORMAT(/,' CRITICAL POINTS: ')
1301  FORMAT(' Enter Integer N (Default is 10) to Change the',/,
     $' Gradient Threshhold Value to 10**(-N) for Non-Nuclear',/,
     $' Critical Points and 10**(-(N-2)) for Nuclear Critical',
     $' Points  ',$)
1302  FORMAT(' Gradient Threshold Set to ',1PE16.8,
     $' for Non-Nuclear',/, ' Critical Points and ',1PE16.8,
     $' for Nuclear Critical Points')
1303  FORMAT(' Enter Damping Factor F (Default is 0.25) for',
     $' Newton-Raphson Search.',/, ' The Value of F for',
     $' Nuclear Searches Will be (F/10.0)  ',$)
1304  FORMAT(' Newton-Raphson Damping Factor Set to ',1PE16.8,
     $' for Non-Nuclear',/, ' Critical Points and ',1PE16.8,
     $' for Nuclear Critical Points')
1310  FORMAT(I4,1X,3(1PE16.8),2A8,2X,'(',I1,',',I2,')')
C
      CALL MAKNAME(1,WFN,ILEN,FWFN)
      IF (ILEN .EQ. 0) STOP 'usage:  extreme wfnfile '
      CALL MAKNAME(1,WLP,ILEN,FWLP)
      IF (ILEN .EQ. 0) STOP 'usage:  extreme wfnfile '
C
      OPEN (IWFN,FILE=WFN,status='unknown')
      OPEN (IWLP,FILE=WLP,status='unknown')
C
      INP=0
C
      CALL RDWFN
C
      WRITE (IOUT,1010) WFNTTL
      WRITE (IOUT,1000)
      READ (INPT,1010) JOBTTL
      WRITE (IWLP,1020)
      WRITE (IWLP,1010) WFNTTL
      WRITE (IWLP,1010) JOBTTL
C
      Write(IWLP,1030)
      DO 10 I = 1,NCENT
        WRITE(IWLP,1040) ATNAM(I),CO(IXC+I),CO(IYC+I),CO(IZC+I)
10    CONTINUE
C
20    WRITE(IOUT,1050)
      READ(INPT,*)IFUNC
      If(IFUNC.lt.1.or.IFUNC.gt.6)Then
      WRITE(Iout,1060)
      GOTO 20
      ENDIF
C
      WRITE(IOUT,*) 
      WRITE(IWLP,*)
      IF(IFUNC.EQ.1) WRITE(IWLP,1070)
      IF(IFUNC.EQ.2) WRITE(IWLP,1080)
      IF(IFUNC.EQ.3) WRITE(IWLP,1090)
      IF(IFUNC.EQ.4) WRITE(IWLP,1100)
      IF(IFUNC.EQ.5) WRITE(IWLP,1110)
      IF(IFUNC.EQ.6) WRITE(IWLP,1115)
      IF(IFUNC.EQ.1) WRITE(IOUT,1070)
      IF(IFUNC.EQ.2) WRITE(IOUT,1080)
      IF(IFUNC.EQ.3) WRITE(IOUT,1090)
      IF(IFUNC.EQ.4) WRITE(IOUT,1100)
      IF(IFUNC.EQ.5) WRITE(IOUT,1110)
      IF(IFUNC.EQ.6) WRITE(IOUT,1115)
      WRITE(IOUT,*) 
C
30    WRITE (IOUT,1120)
      READ (INPT,*) IST
C
      IF (IST .EQ. 0) THEN
C
        write(iout,1121)
        Read(Inpt,*) nucarb
        If(nucarb.eq.0)Then
        WRITE (IOUT,1130)
        READ (INPT,*) XYZ(1),XYZ(2),XYZ(3)
        If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
        WRITE (IOUT,*)
        WRITE (IWLP,*)
        WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
        J1 = 0
        J2 = 0
        Inuc=0
        ElseIf(nucarb.eq.1)Then
        J2=0
        Inuc=1
        write(iout,1122)
        Read(Inpt,*)J1
        XYZ(1)=CO(IXC+J1)
        XYZ(2)=CO(IYC+J1)
        XYZ(3)=CO(IZC+J1)
        If(Iprint.eq.1)WRITE (IOUT,1141) XYZ(1),XYZ(2),XYZ(3)
        WRITE (IOUT,*)
        WRITE (IWLP,*)
        WRITE (IWLP,1141) XYZ(1),XYZ(2),XYZ(3)
        Endif
        IWHOLE=0
        Ifail=0
        CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,Inuc)
        If(Ifail.eq.0)Then
        Icrit=1
        CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
        If(Icancl.eq.0)Then
        Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
        Endif
        Endif
        Goto 30
C
      ELSE IF (IST .EQ. 1) THEN
C
        WRITE (IOUT,1150)
40      READ (INPT,*) J1,J2
        IF(J1.GT.NCENT.OR.J2.GT.NCENT) THEN
        WRITE (IOUT, 1160)
        GOTO 40
        ENDIF
        Write(iout,1170) 
        Read(inpt,*) IHALF
        Fract=pt5
        If(IHALF.eq.0)Then
        Write(Iout,1180) J1, J2
        Read(Inpt,*)Fract
        Endif
        WRITE (IOUT,1190) J1,J2
        WRITE (IWLP,1190) J1,J2
        XYZ(1) = CO(IXC+J1)*fract+CO(IXC+J2)*(one-fract)
        XYZ(2) = CO(IYC+J1)*fract+CO(IYC+J2)*(one-fract)
        XYZ(3) = CO(IZC+J1)*fract+CO(IZC+J2)*(one-fract)
        IWHOLE=0
        Ifail=0
        If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
        WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
        CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
        If(Ifail.eq.0)Then
        Icrit=1
        CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
        If(Icancl.eq.0)Then
        Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
        Endif
        Endif
        Goto 30
C
      ELSE IF (IST .EQ. 2) THEN
C
        WRITE (IOUT,1200)
50      READ (INPT,*) J1,J2,J3
        IF(J1.GT.NCENT.OR.J2.GT.NCENT.OR.J3.GT.NCENT) THEN
        WRITE (IOUT, 1160)
        GOTO 50
        ENDIF
        WRITE (IOUT,1210) J1,J2,J3
        WRITE (IWLP,1210) J1,J2,J3
        XYZ(1) = (CO(IXC+J1)+CO(IXC+J2)+CO(IXC+J3))/THREE
        XYZ(2) = (CO(IYC+J1)+CO(IYC+J2)+CO(IYC+J3))/THREE
        XYZ(3) = (CO(IZC+J1)+CO(IZC+J2)+CO(IZC+J3))/THREE
        J1=0
        J2=0
        IWHOLE=0
        IFAIL=0
        If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
        WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
        CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
        If(IFAIL.eq.0)Then
        Icrit=1
        CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
        If(Icancl.eq.0)Then
        Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
        Endif
        Endif
        Goto 30
C
      ELSE IF (IST .EQ. 3) THEN
C
        WRITE (IOUT,1220)
60      READ (INPT,*) J1,J2,J3,J4
        IF(J1.GT.NCENT.OR.J2.GT.NCENT.OR.J3.GT.NCENT
     $.OR.J4.GT.NCENT) THEN
        WRITE (IOUT, 1160)
        GOTO 60
        ENDIF
        WRITE (IOUT,1230) J1,J2,J3,J4
        WRITE (IWLP,1230) J1,J2,J3,J4
        XYZ(1) = (CO(IXC+J1)+CO(IXC+J2)+CO(IXC+J3)+CO(IXC+J4))/FOUR
        XYZ(2) = (CO(IYC+J1)+CO(IYC+J2)+CO(IYC+J3)+CO(IYC+J4))/FOUR
        XYZ(3) = (CO(IZC+J1)+CO(IZC+J2)+CO(IZC+J3)+CO(IZC+J4))/FOUR
        J1=0 
        J2=0
        IWHOLE=0
        IFAIL=0
        If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
        WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
        CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
        If(IFAIL.eq.0)Then
        Icrit=1
        CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
        If(Icancl.eq.0)Then
        Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
        Endif
        Endif
        Goto 30
C
      ELSE IF (IST .EQ. 4) THEN
C
        CALL ALPHA
        GOTO 30
C
      ELSE IF (IST .EQ. 5) THEN
C
        WRITE (IOUT,1240)
        READ (INPT,*) (XYZ(II),II=1,3)
        IF(IPRINT.eq.1)WRITE (IOUT,1250) (XYZ(II),II=1,3)
        WRITE (IWLP,1250) (XYZ(II),II=1,3)
        IWHOLE=0
        Icrit=0
        J1=0
        J2=0
        Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
        GOTO 30
C
        ELSEIF (IST.EQ.6)Then
C
71      WRITE(IOUT,1260)
        READ(INPT,*) IWHOLE
C
        If(Iwhole.eq.1.or.Iwhole.eq.5)Then
        Do 75 IAtom=1,Ncent
        XYZ(1) = CO(IXC+Iatom)
        XYZ(2) = CO(IYC+Iatom)
        XYZ(3) = CO(IZC+Iatom)
        Ifail=0
        J1=Iatom
        J2=0
        If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
        If(Iprint.eq.1)WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
        CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,1)
        If(IFAIL.eq.0)Then
        Icrit=1
        CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
        If(Icancl.eq.0)Then
        Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
        Endif
        Endif
75      Continue
        If(Iwhole.eq.1)Goto 71
        Endif
C
        If(Iwhole.eq.2.or.Iwhole.eq.5)Then
        Do 80 IAtom=1,Ncent
        Do 90 Katom=1,Iatom-1
        J1 = IATOM
        J2 = KATOM
        WRITE (IWLP,1190) J1,J2
        XYZ(1) = (CO(IXC+J1)+CO(IXC+J2))/TWO
        XYZ(2) = (CO(IYC+J1)+CO(IYC+J2))/TWO
        XYZ(3) = (CO(IZC+J1)+CO(IZC+J2))/TWO
        Ifail=0
        If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
        If(Iprint.eq.1)WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
        CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
        If(IFAIL.eq.0)Then
        Icrit=1
        CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
        If(Icancl.eq.0)Then
        Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
        Endif
        Endif
90      Continue
80      Continue
        If(Iwhole.eq.2)Goto 71
        Endif
C
        If(IWhole.eq.3.or.Iwhole.eq.5)Then
        Do 100 IAtom=1,Ncent
        Do 110 Jatom=1,Iatom-1
        Do 120 Katom=1,Jatom-1
        J1 = Iatom
        J2 = Jatom
        J3 = Katom
        WRITE (IWLP,1210) J1,J2,J3
        XYZ(1) = (CO(IXC+J1)+CO(IXC+J2)+CO(IXC+J3))/THREE
        XYZ(2) = (CO(IYC+J1)+CO(IYC+J2)+CO(IYC+J3))/THREE
        XYZ(3) = (CO(IZC+J1)+CO(IZC+J2)+CO(IZC+J3))/THREE
        J1=0
        J2=0
        Ifail=0
        If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
        If(Iprint.eq.1)WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
        CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
        If(IFAIL.eq.0)Then
        Icrit=1
        CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
        If(Icancl.eq.0)Then
        Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
        Endif
        Endif
120     Continue
110     Continue
100     Continue
        If(Iwhole.eq.3)Goto 71
        Endif
C
        If(Iwhole.eq.4.or.iwhole.eq.5)Then
        Do 130 IAtom=1,Ncent
        Do 140 Jatom=1,Iatom-1
        Do 150 Katom=1,Jatom-1
        Do 160 Latom=1,Katom-1
        J1 = Iatom
        J2 = Jatom
        J3 = Katom
        J4 = Latom
        WRITE (IWLP,1230) J1,J2,J3,J4
        XYZ(1) = (CO(IXC+J1)+CO(IXC+J2)+CO(IXC+J3)+CO(IXC+J4))/FOUR
        XYZ(2) = (CO(IYC+J1)+CO(IYC+J2)+CO(IYC+J3)+CO(IYC+J4))/FOUR
        XYZ(3) = (CO(IZC+J1)+CO(IZC+J2)+CO(IZC+J3)+CO(IYC+J4))/FOUR
        J1=0
        J2=0
        Ifail=0
        If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
        If(Iprint.eq.1)WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
        CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
        If(IFAIL.eq.0)Then
        Icrit=1
        CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
        If(Icancl.eq.0)Then
        Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
        Endif
        Endif
160     Continue
150     Continue
140     Continue
130     Continue
        If(Iwhole.eq.4)Goto 71
        Endif
C
72     Write(Iout,1265)
       Read(Inpt,*) ILine
       Iwhole=6
       If(Iline.eq.1.or.iline.eq.2)Then
       If(ILine.eq.1)Then
       Write(iout,1266)
       Read(Inpt,*) J1,J2,Npts
       Dist=dsqrt((CO(IXC+J1)-CO(IXC+J2))**2+
     $(CO(IYC+J1)-CO(IYC+J2))**2+
     $(CO(IZC+J1)-CO(IZC+J2))**2)
       Unitx = (CO(IXC+J2)-CO(IXC+J1))/Dist
       Unity = (CO(IYC+J2)-CO(IYC+J1))/Dist
       Unitz = (CO(IZC+J2)-CO(IZC+J1))/Dist
       X0=CO(IXC+J1)
       Y0=CO(IYC+J1)
       Z0=CO(IZC+J1)
       ElseIf(Iline.eq.2)Then
       Write(Iout,1267)
       Read(Inpt,*) X1,Y1,Z1
       Write(Iout,1268)
       Read(Inpt,*) X2,Y2,Z2,NPts
       Dist=Dsqrt((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
       Unitx=(X2-X1)/Dist
       Unity=(Y2-Y1)/Dist
       Unitz=(Z2-Z1)/Dist
       X0=X1
       Y0=Y1
       Z0=Z1
       Endif
       Dinc=Dist/NPts
       Do 162 I=1,Npts
       XYZ(1)= X0 + I*Dinc*UNITX
       XYZ(2)= Y0 + I*Dinc*UNITY
       XYZ(3)= Z0 + I*Dinc*UNITZ
        J1=0
        J2=0
        Ifail=0
        If(Iprint.eq.1)WRITE (IOUT,1140) XYZ(1),XYZ(2),XYZ(3)
        If(Iprint.eq.1)WRITE (IWLP,1140) XYZ(1),XYZ(2),XYZ(3)
        CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
        If(IFAIL.eq.0)Then
        Icrit=1
        CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
        If(Icancl.eq.0)Then
        Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
        Endif
        Endif
162     Continue
        Goto 72
        Endif
C 
73     Write(IOUT,1270)
       READ(INPT,*) Isphere
       If(ISphere.ne.0)Then
       Iwhole=7
       Pi=Dacos(-one)
       Write(Iout,1280)
       Read(inpt,*)Iall
       Write(Iout,1290)
       Read(inpt,*)Iphi,ITheta,RadMin,Radmax,Nrad
       If(iall.ne.0)Then
       iatom=iall
       natoms=1
       Else
       iatom=1
       natoms=ncent
       Endif
       Phinc=Two*pi/Iphi
       HPhinc=Pt5*Phinc
       Thinc=Two*pi/ITheta
       HThinc=Pt5*Thinc
       Radinc=Radmax/Nrad
       Do 170 jatom=1,natoms
       Do 180 I=1,Iphi
       Phi=Phinc*I-HPhinc
       Do 190 J=1,Itheta
       Theta=THinc*J-HThinc
       unitx=dsin(theta)*dcos(phi)
       unity=dsin(theta)*dsin(phi)
       unitz=dcos(theta)
       Do 200 K=1,Nrad
       XYZ(1)=CO(IXC+jatom+iatom-1)+(RadMin+RADINC*K)*UNITX
       XYZ(2)=CO(IYC+jatom+iatom-1)+(RadMin+RADINC*K)*UNITY
       XYZ(3)=CO(IZC+jatom+iatom-1)+(RadMin+RADINC*K)*UNITZ
       J1=0
       J2=0
       Ifail=0
       CALL NEWTON(XYZ,IFAIL,IFUNC,IWHOLE,NITER,0)
       If(IFAIL.eq.0)Then
       Icrit=1
       CALL STORE(XYZ,X,Y,Z,J1,J2,ICancl,INP,Ji1,Ji2,NITER)
        If(Icancl.eq.0)Then
        Call Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,EVSAVE)
        Endif
       Endif
200    Continue
190    Continue
180    Continue
170    Continue
       Goto 73
       Endif 
       Goto 30
C
      ELSEIF(IST.EQ.7)Then
C
220   Write(iout,*)
      Write(iout,1295)
      Read(inpt,*) Iopt
      If(Iopt.eq.1)Then
      Write(Iwlp,*)
      Write(iout,*)
      Write(iout,1296)
      Read(inpt,*)ICut
      Cut=Ten**(-Icut)
      Write(iwlp,1298)cut
      Goto 220
      ElseIf(Iopt.eq.2)Then
      Write(Iwlp,*)
      Write(iout,*)
      Write(iout,1297)
      Read(Inpt,*) IPrint
      Write(iwlp,1299)IPrint
      Goto 220
      ElseIf(Iopt.eq.3)Then
      Write(Iwlp,*)
      Write(iout,*)
      Write(iout,1301)
      Read(Inpt,*) Ieps
      eps=Ten**(-ieps)
      iepnuc=ieps-2
      epsnuc=Ten**(-iepnuc)
      Write(iwlp,1302)eps,epsnuc
      Goto 220
      ElseIf(Iopt.eq.4)Then
      Write(Iwlp,*)
      Write(iout,*)
      Write(iout,1303)
      Read(Inpt,*) Dmp
      DmpNuc=Dmp/Ten
      Write(iwlp,1304)Dmp,DmpNuc
      Goto 220
      Endif
      Write(iout,*)
      Goto 30
C
      ELSE IF ((IST .EQ. 8).or.(Ist.eq.9)) THEN
C
        WRITE (IOUT,1300)
        If(IST.eq.9)WRITE (IWLP,1300)
        DO 210 I = 1,INP
        ev1=evsave(i,1)
        ev2=evsave(i,2)
        ev3=evsave(i,3)
        aev1=dabs(ev1)
        aev2=dabs(ev2)
        aev3=dabs(ev3)
        ineg=0
        ipos=0
        If(ev1.lt.zero.and.aev1.gt.degen)ineg=ineg+1
        If(ev2.lt.zero.and.aev2.gt.degen)ineg=ineg+1
        If(ev3.lt.zero.and.aev3.gt.degen)ineg=ineg+1
        If(ev1.gt.zero.and.aev1.gt.degen)ipos=ipos+1
        If(ev2.gt.zero.and.aev2.gt.degen)ipos=ipos+1
        If(ev3.gt.zero.and.aev3.gt.degen)ipos=ipos+1
        irank=ineg+ipos
        isig=ipos-ineg 
        IF(JI1(I).NE.0.AND.JI2(I).NE.0) THEN
        WRITE (IOUT,1310) I,X(I),Y(I),Z(I),ATNAM(JI1(I)),
     $  ATNAM(JI2(I)),irank,isig
        If(Ist.eq.9)WRITE (IWLP,1310) I,X(I),Y(I),Z(I),ATNAM(JI1(I)),
     $  ATNAM(JI2(I)),irank,isig
        ELSEIF(JI1(I).NE.0.AND.JI2(I).EQ.0) THEN
        WRITE (IOUT,1310) I,X(I),Y(I),Z(I),ATNAM(JI1(I)),
     $  BLANK,irank,isig
        If(Ist.eq.9)WRITE (IWLP,1310) I,X(I),Y(I),Z(I),ATNAM(JI1(I)),
     $  BLANK,irank,isig
        ELSE
        WRITE (IOUT,1310) I,X(I),Y(I),Z(I),BLANK,BLANK,irank,isig
        If(Ist.eq.9)WRITE (IWLP,1310) I,X(I),Y(I),Z(I),BLANK,BLANK,
     $  irank,isig
        ENDIF
210     CONTINUE
        If(Ist.eq.8)GOTO 30
C
        STOP ' EXTREME IS DONE '
C
      ELSE
        WRITE (IOUT,1100)
        GOTO 30
      END IF
C
      END
      SUBROUTINE ALPHA
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter (MaxOff=200000,MaxAtm=100)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /ANG/ ANGLE(3,MaxAtm,MaxAtm)
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICOFMx
      COMMON /UNITS/ INPT, IOUT, IWFN, IWLP
      Save One,ThSxty
      Data One/1.d0/,ThSxty/3.6d2/ 
1000  FORMAT(' HERE IS A LIST OF BOND PATHS THAT HAVE BEEN TRACED ')
1010  FORMAT(/,' INPUT ATOMS FOR WHICH ANGLE IS DESIRED (EG.1,2,3) ',$)
1020  FORMAT(/,' BOND PATH ANGLE = ',F8.4)
1030  FORMAT(/,' GEOMETRIC BOND ANGLE = ',F8.4)
1040  FORMAT(/,' DEVIATION OF FIRST BOND PATH VECTOR FROM FIRST',
     + ' GEOMETRIC BOND VECTOR   ',F8.4)
1050  FORMAT(/,' DEVIATION OF SECOND BOND PATH VECTOR FROM SECOND',
     + ' GEOMETRIC BOND VECTOR ',F8.4)
1060  FORMAT(/,' GEOMETRIC BOND ANGLE MINUS BOND PATH ANGLE ',F8.4)
1070  FORMAT(/,' WOULD YOU LIKE TO DO ANOTHER COMPARISION ? (0/1) ',$)
1080  FORMAT(/,' COMPARISION FOR ATOMS ',3I4)
C
      Pi=DaCos(-One)
      Degree=thsxty/(Two*Pi)
C
100   WRITE (IOUT,1010)
      READ (INPT,*) I,J,K
      WRITE (IOUT,1080) I,J,K
      WRITE (IWLP,1080) I,J,K
C
C    GET BOND PATH ANGLE
C
      S = ANGLE(1,J,K)*ANGLE(1,J,I) +
     +    ANGLE(2,J,K)*ANGLE(2,J,I) +
     +    ANGLE(3,J,K)*ANGLE(3,J,I)
      BANGLE = (DACOS(S))*Degree
      WRITE (IOUT,1020) BANGLE
      WRITE (IWLP,1020) BANGLE
C
C    USE UNIT VECTORS ALONG THE BONDS TO CALCULATE
C    GEOMETRIC BOND ANGLE
C
      C1 = CO(IXC+I) - CO(IXC+J)
      D1 = CO(IXC+K) - CO(IXC+J)
      C2 = CO(IYC+I) - CO(IYC+J)
      D2 = CO(IYC+K) - CO(IYC+J)
      C3 = CO(IZC+I) - CO(IZC+J)
      D3 = CO(IZC+K) - CO(IZC+J)
      S1 = SQRT(C1**2 + C2**2 + C3**2)
      S2 = SQRT(D1**2 + D2**2 + D3**2)
      C1 = C1/S1
      C2 = C2/S1
      C3 = C3/S1
      D1 = D1/S2
      D2 = D2/S2
      D3 = D3/S2
C
C      GET ANGLES BETWEEN BOND PATHS AND BOND DIRECTIONS
C
      S1 = ANGLE(1,J,I)*C1 + ANGLE(2,J,I)*C2 + ANGLE(3,J,I)*C3
      S2 = ANGLE(1,J,K)*D1 + ANGLE(2,J,K)*D2 + ANGLE(3,J,K)*D3
      S3 = C1*D1 + C2*D2 + C3*D3
      ANGLE1 = (DACOS(S1))*Degree
      ANGLE2 = (DACOS(S2))*Degree
      ANGLE3 = (DACOS(S3))*Degree
      WRITE (IOUT,1030) ANGLE3
      WRITE (IWLP,1030) ANGLE3
      WRITE (IOUT,1040) ANGLE1
      WRITE (IWLP,1040) ANGLE1
      WRITE (IOUT,1050) ANGLE2
      WRITE (IWLP,1050) ANGLE2
C
C    DETERMINE THE DIFFERENCE BETWEEN THE BOND PATH ANGLE AND THE
C    GEOMETRIC BOND ANGLE
C
      DANGLE = ANGLE3 - BANGLE
      WRITE (IOUT,1060) DANGLE
      WRITE (IWLP,1060) DANGLE
C
      WRITE (IOUT,1070)
      READ (INPT,*) IGO
      IF (IGO .EQ. 1) GOTO 100
C
      RETURN
C
C    FORMATS
C
      END
      Subroutine AZero
C
      Implicit Double Precision (A-H,O-Z)
C
      Common /AFG/ A0,AX,AY,AZ,AXX,AXY,AXZ,AYY,AYZ,AZZ,AXXX,
     $  AXXY,AXXZ,AXYY,AXZZ,AXYZ,AYYY,AYYZ,AYZZ,AZZZ,AXXXX,AXXXY,
     $  AXXXZ,AXXYY,AXXZZ,AXYYY,AXZZZ,AXXYZ,AXYYZ,AXYZZ,AYYYY,
     $  AYYYZ,AYYZZ,AYZZZ,AZZZZ,F0,FX,FY,FZ,FXX,FXY,FXZ,FYY,FYZ,
     $  FZZ,FXXX,FXXY,FXXZ,FXYY,FXZZ,FXYZ,FYYY,FYYZ,FYZZ,FZZZ,
     $  FXXXX,FXXXY,FXXXZ,FXXYY,FXXZZ,FXYYY,FXZZZ,FXXYZ,FXYYZ,
     $  FXYZZ,FYYYY,FYYYZ,FYYZZ,FYZZZ,FZZZZ,G0,GX,GY,GZ,GXX,GXY,
     $  GXZ,GYY,GYZ,GZZ,GXXX,GXXY,GXXZ,GXYY,GXZZ,GXYZ,GYYY,GYYZ,
     $  GYZZ,GZZZ,GXXXX,GXXXY,GXXXZ,GXXYY,GXXZZ,GXYYY,GXZZZ,
     $  GXXYZ,GXYYZ,GXYZZ,GYYYY,GYYYZ,GYYZZ,GYZZZ,GZZZZ
      Save Zero
      Data Zero/0.0d0/
C
      A0=Zero
      AX=Zero
      AY=Zero
      AZ=Zero
      AXX=Zero
      AXY=Zero
      AXZ=Zero
      AYY=Zero
      AYZ=Zero
      AZZ=Zero
      AXXX=Zero
      AXXY=Zero
      AXXZ=Zero
      AXYY=Zero
      AXZZ=Zero
      AXYZ=Zero
      AYYY=Zero
      AYYZ=Zero
      AYZZ=Zero
      AZZZ=Zero
      AXXXX=Zero
      AXXXY=Zero
      AXXXZ=Zero
      AXXYY=Zero
      AXXZZ=Zero
      AXYYY=Zero
      AXZZZ=Zero
      AXXYZ=Zero
      AXYYZ=Zero
      AXYZZ=Zero
      Ayyyy=Zero
      Ayyyz=Zero
      Ayyzz=Zero
      Ayzzz=Zero
      Azzzz=Zero
C
      Return
      End
      BLOCK DATA
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /C2/ EPSD,AB(9,9),AM(9,10)
      COMMON /UNITS/  INPT,IOUT,IWFN,IWLP
      COMMON /OPTIONS/ ICUT,IPrint,Eps,EpsNuc,Dmp,DmpNuc
      DATA INPT /5/, IOUT /6/, IWFN /10/, IWLP /12/, ICUT /20/,
     $IPrint/0/,Eps/1.d-10/,EpsNuc/1.d-8/,Dmp/0.25d0/,DmpNuc/0.025d0/
      DATA AB/1.D0,1.5D0,1.9166666666667D0,2.2916666666667D0,
     *2.6402777777778D0,
     *2.9701388888889D0,3.2857308201058D0,3.5899553571428D0,
     *3.8848233575837D0,0.D0,
     *-.5D0,-1.3333333333333D0,-2.4583333333333D0,-3.8527777777778D0,
     *-5.5020833333333D0,-7.3956349206348D0,-9.5252066798940D0,
     *-11.884150683421D0,0D0,0D0,.4166666666667D0,1.5416666666667D0,
     *3.6333333333333D0,6.9319444444444D0,11.6658234126982D0,
     *18.054538690476D0,26.310842702821D0,0D0,0D0,0D0,-.375D0,
     *-1.7694444444444D0,-5.0680555555556D0,-11.379894179894D0,
     *-22.027752976190D0,-38.540361000881D0,0D0,0D0,0D0,0D0,
     *.3486111111111D0,1.9979166666667D0,6.7317956349205D0,
     *17.379654431217D0,38.020414462081D0,0D0,0D0,0D0,0D0,0D0,
     *-.3298611111111D0,-2.2234126984127D0,-8.6121279761903D0,
     *-25.124736000882D0,0D0,0D0,0D0,0D0,0D0,0D0,.3155919312169D0,
     *2.4451636904761D0,10.701467702822D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,
     *-.3042245370370D0,-2.6631685405644D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,
     *0D0,.2948680004409D0/
      DATA AM/.5D0,.4166666666667D0,.375D0,.3486111111111D0,
     *.3298611111111D0,.3155919312169D0,.3042245370370D0,
     *.2948680004409D0,.2869754464286D0,.5D0,.6666666666667D0,
     *.7916666666667D0,.8972222222222D0,.9909722222222D0,
     *1.0765873015873D0,1.1561590608466D0,1.2310113536155D0,
     *1.3020443397266D0,0D0,-.0833333333333D0,-.2083333333333D0,
     *-.3666666666667D0,-.5541666666667D0,-.7682043650794D0,
     *-1.0069196428571D0,-1.2689026675485D0,-1.55303461199292D0,0D0,0D0,
     *.0416666666667D0,.1472222222222D0,.3347222222222D0,
     *.6201058201058D0,1.0179646164021D0,1.5419306657848D0,
     *2.2049052028218D0,0D0,0D0,0D0,-.0263888888889D0,-.1201388888889D0,
     *-.3341765873016D0,-.7320353835979D0,
     *-1.3869929453262D0,-2.3814547508818D0,0D0,0D0,0D0,0D0,.01875D0,
     *.1043650793651D0,.3430803571429D0,.8670464065255D0,
     *1.8615082120811D0,0D0,0D0,0D0,0D0,0D0,
     *-.0142691798942D0,-.0938409391534D0,-.3558239638448D0,
     *-1.0187985008818D0,0D0,0D0,0D0,0D0,0D0,0D0,
     *.0113673941799D0,.0862196869489D0,.3703516313933D0,
     *0D0,0D0,0D0,0D0,0D0,0D0,0D0,-.0093565365961D0,-.0803895227072D0,
     *0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,.0078925540124D0/
C
      END
      SUBROUTINE BOND(X,Y,Z,BL,NCAL,PN,IFunc,Iwhole)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter(MaxOff=200000,MaxAtm=100,MaxStp=141,MinStp=6)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
      Common /What/ IWhat
      DIMENSION R1(3,MaxStp),PN(4),PR(3,MaxStp),H(3,3),SG(3,3),
     $  PPN(3),GRHO(3),tmpxyz(3)
      Save Zero,dxyz,pt15,fifty,pt2,one,two,three,fvhund,
     $  pt01,pt001,small
      Data Zero/0.0d0/,dxyz/1.d-3/,pt15/0.15d0/,fifty/50.0d0/,
     $  pt2/0.2d0/,One/1.0d0/,Two/2.0d0/,Three/3.0d0/,
     $  fvhund/500.0d0/,pt01/0.01d0/,pt001/0.001d0/,small/1.d-9/,
     $  Hund/1.d2/
C
      Iwhat=0
      Pi=Dacos(-one)
      BL = Zero
      NCAL = 0
      R1(1,1) = X
      R1(2,1) = Y
      R1(3,1) = Z
      tmpxyz(1)=R1(1,1)
      tmpxyz(2)=R1(2,1)
      tmpxyz(3)=R1(3,1)
      CALL GRDRHO(0,tmpxyz,RHO,GRHO,GRAD,H,SG)
      PR(1,1)=GRHO(1)
      PR(2,1)=GRHO(2)
      PR(3,1)=GRHO(3)
      NCAL = NCAL+1
C
110   NN = INT((RMIN-Pt15)*Fifty)
      IF (NN.LT.MinStp) NN = MinStp
      IF (NN.GT.(Maxstp-1)) NN = MaxStp-1
      DS = (RMIN-Pt15)/DFLOAT(NN)
      CALL DES(R1,PR,DS,NN)
      NCAL = NCAL+NN+NN
      BL = BL+DFLOAT(NN)*DS
      DO 100 I = 1,3
        R1(I,1) = R1(I,NN+1)
        PR(I,1) = PR(I,NN+1)
100   CONTINUE
      IF (RMIN.GT.Pt2) GOTO 110
      PPN(1) = CO(IXC+MINR)
      PPN(2) = CO(IYC+MINR)
      PPN(3) = CO(IZC+MINR)
      BL = BL+RMIN+DXYZ
      IFAIL = 1
      CALL NEWTON(PPN,IFAIL,IFUNC,Iwhole,NITER,1)
      PN(1) = PPN(1)
      PN(2) = PPN(2)
      PN(3) = PPN(3)
      PN(4) = RMIN
130   GR = DSQRT((PN(1)-R1(1,1))**2+(PN(2)-R1(2,1))**2+
     +           (PN(3)-R1(3,1))**2)
      IF (GR .GT. pt01) THEN
        NN = INT((GR-dxyz)*fvhund)
        IF (NN.LT.minstp) NN = minstp
        IF (NN.GT.(maxstp-1)) NN = MaxStp-1
        DS = (GR-dxyz)/DFLOAT(NN)
        CALL DES(R1,PR,DS,NN)
        NCAL = NCAL+NN+NN
        DO 120 I = 1,3
          R1(I,1) = R1(I,NN+1)
          PR(I,1) = PR(I,NN+1)
120     CONTINUE
        GOTO 130
      END IF
C
      BANGLE(4) = one
      BANGLE(5) = DACOS(BANGLE(3))
      IF (DABS(BANGLE(1)) .LT. Small) THEN
        BANGLE(6) = Pi/Two
        IF (BANGLE(2).LT.Zero) BANGLE(6) = Three*BANGLE(6)
      ELSE
        BANGLE(6) = DATAN(BANGLE(2)/BANGLE(1))
      ENDIF
      IF (BANGLE(1).LT.Zero) BANGLE(6) = BANGLE(6)+Pi
      IF(BANGLE(1) .GE. Zero .AND. BANGLE(2) .LT. Zero) BANGLE(6)=
     1 BANGLE(6)+Two*Pi
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DTEMP
      INTEGER I,INCX,M,MP1,N,NINCX
C
      DASUM = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DTEMP = DTEMP + DABS(DX(I))
   10 CONTINUE
      DASUM = DTEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I + 1)) + DABS(DX(I + 2))
     *  + DABS(DX(I + 3)) + DABS(DX(I + 4)) + DABS(DX(I + 5))
   50 CONTINUE
   60 DASUM = DTEMP
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,IXIY,M,MP1,N
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
      SUBROUTINE DES(R1,PR,DS,NN)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter (MaxStp=141)
      COMMON /C2/ EPSD,AB(9,9),AM(9,10)
      DIMENSION R1(3,MaxStp),W1(3,9),R(3),W(3),PR(3,MaxStp),
     $  H(3,3), SG(3,3),tmpxyz(3),grho(3)
      Save Zero,One
      Data Zero/0.0d0/,One/1.0d0/
C
C     IF IT IS THE FIRST STEP, START WITH A SINGLE STEP METHOD.
C
      DO 100 I = 1,5
        DO 110 L = 1,3
        W(L) = R1(L,I)
        DO 110 J = 1,I
          G = DSQRT(PR(1,I-J+1)**2+PR(2,I-J+1)**2+PR(3,I-J+1)**2)
          W1(L,J) = PR(L,I-J+1)/G
110     CONTINUE
      CALL MULTI1(W,R1(1,I+1),W1,I,DS,DFLOAT(I)*DS)
      tmpxyz(1)=R1(1,I+1)
      tmpxyz(2)=R1(2,I+1)
      tmpxyz(3)=R1(3,I+1)
      CALL GRDRHO(0,tmpxyz,RHO,GRHO,GRAD,H,SG)
      PR(1,I+1)=GRHO(1)
      PR(2,I+1)=GRHO(2)
      PR(3,I+1)=GRHO(3)
100   CONTINUE
C
C     PERFORM THE REMAINING STEPS WITH A 6 STEP METHOD OF ORDER 7.
C
      DO 130 I = 6,NN
        DO 140 L = 1,3
        W(L) = R1(L,I)
        DO 140 J = 1,6
          G = DSQRT(PR(1,I-J+1)**2+PR(2,I-J+1)**2+PR(3,I-J+1)**2)
          W1(L,J) = PR(L,I-J+1)/G
140     CONTINUE
        CALL MULTI1(W,R1(1,I+1),W1,6,DS,DFLOAT(I)*DS)
      tmpxyz(1)=R1(1,I+1)
      tmpxyz(2)=R1(2,I+1)
      tmpxyz(3)=R1(3,I+1)
        CALL GRDRHO(0,tmpxyz,RHO,GRHO,GRAD,H,SG)
      PR(1,I+1)=GRHO(1)
      PR(2,I+1)=GRHO(2)
      PR(3,I+1)=GRHO(3)
        IF (I .NE. (I/20)*20) GOTO 130
        DO 150 L = 1,3
          W(L) = R1(L,I)
150     CONTINUE
        CALL MULTI1(W,R,W1,5,DS,DFLOAT(I)*DS)
C
C     ESTIMATE THE ERROR.
C
        EPS1 = Zero
        DO 160 L = 1,3
          EPS1 = EPS1+(R(L)-R1(L,I+1))**2
160     CONTINUE
        EPS1 = DSQRT(EPS1)
        RH = DMAX1(One,DSQRT(R1(1,I+1)**2+R1(2,I+1)**2+R1(3,I+1)**2))
        EPSD = DMAX1(EPSD,EPS1/RH)
130   CONTINUE
      RETURN
      END
      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)
      INTEGER LDA,N,IPVT(1)
      DOUBLE PRECISION A(LDA,1),Z(1)
      DOUBLE PRECISION RCOND
C
C     DGECO FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DGEFA
C     BLAS DAXPY,DDOT,DSCAL,DASUM
C     FORTRAN DABS,DMAX1,DSIGN
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
C
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL DGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30
            S = DABS(A(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + DABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150
            S = DABS(A(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
      SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      INTEGER LDA,N,IPVT(1),JOB
      DOUBLE PRECISION A(LDA,1),DET(2),WORK(1)
C
C     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        WORK    DOUBLE PRECISION(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     DOUBLE PRECISION(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,DSWAP
C     FORTRAN DABS,MOD
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      DOUBLE PRECISION TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (DABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(1),INFO
      DOUBLE PRECISION A(LDA,1)
C
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.
C
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,IDAMAX
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
      SUBROUTINE DIVSTR(dsig,dsigm)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter(MaxOff=200000)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICOFMx
      DIMENSION dsig(3)
      Save Zero,Pt5
      Data Zero/0.d0/,Pt5/0.5d0/
C
      DO 100 I = 1,3
100     Dsig(I) = Dsig(I)
C
      DO 110 J = 1,NMO
        Temp = -(CO(IGXX+J)+CO(IGYY+J)+CO(IGZZ+J))
        DSig(1) = DSig(1)+CO(IP+J)*(CO(IPSI+J)*(CO(IGXXX+J)+
     $            CO(IGXYY+J)+CO(IGXZZ+J)) + temp*CO(IGX+J))
        DSig(2) = DSig(2)+CO(IP+J)*(CO(IPSI+J)*(CO(IGXXY+J)+
     $            CO(IGYYY+J)+CO(IGYZZ+J)) + temp*CO(IGY+J))
110     DSig(3) = DSig(3)+CO(IP+J)*(CO(IPSI+J)*(CO(IGXXZ+J)+
     $            CO(IGYYZ+J)+CO(IGZZZ+J)) + temp*CO(IGZ+J))
C
      DO 120 I = 1,3
120     Dsig(I) = pt5*Dsig(I)
C
      Dsigm = dsqrt(dsig(1)**2+dsig(2)**2+dsig(3)**2)

      RETURN
      END
      SUBROUTINE  DSCAL(N,DA,DX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DA,DX(1)
      INTEGER I,INCX,M,MP1,N,NINCX
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
   50 CONTINUE
      RETURN
      END
      SUBROUTINE GAUS4
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter (MaxOff=200000)
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICOFMX
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      Common /AFG/ A0,AX,AY,AZ,AXX,AXY,AXZ,AYY,AYZ,AZZ,AXXX,
     $  AXXY,AXXZ,AXYY,AXZZ,AXYZ,AYYY,AYYZ,AYZZ,AZZZ,AXXXX,AXXXY,
     $  AXXXZ,AXXYY,AXXZZ,AXYYY,AXZZZ,AXXYZ,AXYYZ,AXYZZ,AYYYY,
     $  AYYYZ,AYYZZ,AYZZZ,AZZZZ,F0,FX,FY,FZ,FXX,FXY,FXZ,FYY,FYZ,
     $  FZZ,FXXX,FXXY,FXXZ,FXYY,FXZZ,FXYZ,FYYY,FYYZ,FYZZ,FZZZ,
     $  FXXXX,FXXXY,FXXXZ,FXXYY,FXXZZ,FXYYY,FXZZZ,FXXYZ,FXYYZ,
     $  FXYZZ,FYYYY,FYYYZ,FYYZZ,FYZZZ,FZZZZ,G0,GX,GY,GZ,GXX,GXY,
     $  GXZ,GYY,GYZ,GZZ,GXXX,GXXY,GXXZ,GXYY,GXZZ,GXYZ,GYYY,GYYZ,
     $  GYZZ,GZZZ,GXXXX,GXXXY,GXXXZ,GXXYY,GXXZZ,GXYYY,GXZZZ,
     $  GXXYZ,GXYYZ,GXYZZ,GYYYY,GYYYZ,GYYZZ,GYZZZ,GZZZZ
      Common /What/ Iwhat
      Common /Options/ Icut, Iprint,Eps,EpsNuc,Dmp,DmpNuc
      Save Zero,Ifill,IEmpty,Ten
      Data Zero/0.0d0/,Ifill/1/,IEmpty/0/,Ten/10.0d0/
C
      IPlac=35
      Cutoff=Ten**(-Icut)
C
      DO 10 J=1,Iplac*NMO
      CO(IPSI+J)=ZERO
10    CONTINUE
C
      DO 20 I=1,NPRIMS
      K=IC(ICENT+I)
      ENT=CO(IE+I)
      X=CO(IXX+K)
      Y=CO(IYY+K)
      Z=CO(IZZ+K)
      EXPON=DEXP(-ENT*CO(IR2+K))
C
      Call PrimeF(Ent,Expon,X,Y,Z)
C
      IWhich=IC(ITYPE+I)
C
      Call PrimeA(IWhich,Ifill,X,Y,Z)
C
      Call PrimeG
C
      Call PrimeA(IWhich,IEmpty,X,Y,Z)
C
      Check=Dmax1(dabs(G0),dabs(Gx),dabs(Gy),dabs(Gz))
      If(check*CO(Icofmx+I).gt.cutoff)Then
      DO 30 J=1,NMO
      CIJ=CO(IMO+NPRIMS*(J-1)+I)
      CO(IPSI+J)=CO(IPSI+J)+CIJ*G0
      CO(IGX+J)=CO(IGX+J)+CIJ*GX
      CO(IGY+J)=CO(IGY+J)+CIJ*GY
30    CO(IGZ+J)=CO(IGZ+J)+CIJ*GZ
      Endif
      If(Iwhat.gt.0)Then
      Check=Dmax1(dabs(GXX),dabs(GXY),dabs(GXZ),dabs(GYY),
     $  dabs(GYZ),dabs(GZZ))
      If(check*CO(Icofmx+I).gt.cutoff)Then
      DO 40 J=1,NMO
      CIJ=CO(IMO+NPRIMS*(J-1)+I)
      CO(IGXX+J)=CO(IGXX+J)+CIJ*GXX
      CO(IGXY+J)=CO(IGXY+J)+CIJ*GXY
      CO(IGXZ+J)=CO(IGXZ+J)+CIJ*GXZ
      CO(IGYY+J)=CO(IGYY+J)+CIJ*GYY
      CO(IGYZ+J)=CO(IGYZ+J)+CIJ*GYZ
40    CO(IGZZ+J)=CO(IGZZ+J)+CIJ*GZZ
      Endif
      Endif
      If(Iwhat.gt.1)Then
      Check=Dmax1(dabs(GXXX),dabs(GXXY),dabs(GXXZ),dabs(GXYY),
     $  dabs(GXZZ),dabs(GXYZ),dabs(GYYY),dabs(GYYZ),dabs(GYZZ),
     $  dabs(GZZZ))
      If(check*CO(Icofmx+I).gt.cutoff)Then
      DO 50 J=1,NMO
      CIJ=CO(IMO+NPRIMS*(J-1)+I)
      CO(IGXXX+J)=CO(IGXXX+J)+CIJ*GXXX
      CO(IGXYY+J)=CO(IGXYY+J)+CIJ*GXYY
      CO(IGXZZ+J)=CO(IGXZZ+J)+CIJ*GXZZ
      CO(IGXXY+J)=CO(IGXXY+J)+CIJ*GXXY
      CO(IGXXZ+J)=CO(IGXXZ+J)+CIJ*GXXZ
      CO(IGXYZ+J)=CO(IGXYZ+J)+CIJ*GXYZ
      CO(IGYYY+J)=CO(IGYYY+J)+CIJ*GYYY
      CO(IGYYZ+J)=CO(IGYYZ+J)+CIJ*GYYZ
      CO(IGYZZ+J)=CO(IGYZZ+J)+CIJ*GYZZ
50    CO(IGZZZ+J)=CO(IGZZZ+J)+CIJ*GZZZ
      Endif
      Endif
      If(Iwhat.gt.2)Then
      Check=Dmax1(dabs(GXXXX),dabs(GXXXY),dabs(GXXXZ),
     $  dabs(GXXYY),dabs(GXXZZ),dabs(GXYYY),dabs(GXZZZ),
     $  dabs(GXYYZ),dabs(GXYZZ),dabs(GYYYY),dabs(GYYYZ),
     $  dabs(GYYZZ),dabs(GYZZZ),dabs(GZZZZ))
      If(check*CO(Icofmx+I).gt.cutoff)Then
      DO 60 J=1,NMO
      CIJ=CO(IMO+NPRIMS*(J-1)+I)
      CO(IGXXXX+J)=CO(IGXXXX+J)+CIJ*GXXXX
      CO(IGXXXY+J)=CO(IGXXXY+J)+CIJ*GXXXY
      CO(IGXXXZ+J)=CO(IGXXXZ+J)+CIJ*GXXXZ
      CO(IGXXYY+J)=CO(IGXXYY+J)+CIJ*GXXYY
      CO(IGXXZZ+J)=CO(IGXXZZ+J)+CIJ*GXXZZ
      CO(IGXYYY+J)=CO(IGXYYY+J)+CIJ*GXYYY
      CO(IGXZZZ+J)=CO(IGXZZZ+J)+CIJ*GXZZZ
      CO(IGXXYZ+J)=CO(IGXXYZ+J)+CIJ*GXXYZ
      CO(IGXYYZ+J)=CO(IGXYYZ+J)+CIJ*GXYYZ
      CO(IGXYZZ+J)=CO(IGXYZZ+J)+CIJ*GXYZZ
      CO(IGYYYY+J)=CO(IGYYYY+J)+CIJ*GYYYY
      CO(IGYYYZ+J)=CO(IGYYYZ+J)+CIJ*GYYYZ
      CO(IGYYZZ+J)=CO(IGYYZZ+J)+CIJ*GYYZZ
      CO(IGYZZZ+J)=CO(IGYZZZ+J)+CIJ*GYZZZ
60    CO(IGZZZZ+J)=CO(IGZZZZ+J)+CIJ*GZZZZ
      Endif
      Endif
20    CONTINUE
C
      RETURN
      END
      SUBROUTINE GEOM (N,XYZ,RN,AYZ,AXZ,AXY)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter(MaxOff=200000)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
      Dimension XYZ(3)
      Save Zero,One,Two,thsxty
      Data Zero/0.d0/,One/1.0d0/,Two/2.0d0/,Thsxty/360.0d0/
C
      Pi=Dacos(-one)
      Degree=Thsxty/(Two*pi)
C
      XN2 = (XYZ(1)-CO(IXC+N))*(XYZ(1)-CO(IXC+N))
      YN2 = (XYZ(2)-CO(IYC+N))*(XYZ(2)-CO(IYC+N))
      ZN2 = (XYZ(3)-CO(IZC+N))*(XYZ(3)-CO(IZC+N))
      RN = DSQRT(XN2+YN2+ZN2)
      IF (RN.EQ.Zero) THEN
        AXY = Zero
        AXZ = Zero
        AYZ = Zero
        RETURN
      ELSE
        RYZ = DSQRT(YN2+ZN2)
        RXZ = DSQRT(XN2+ZN2)
        RXY = DSQRT(XN2+YN2)
        AYZ = DACOS(RYZ/RN)*DEGREE
        AXZ = DACOS(RXZ/RN)*DEGREE
        AXY = DACOS(RXY/RN)*DEGREE
        RETURN
      ENDIF
      END
      SUBROUTINE GRDD2R(Iopt,R,VALUE,GRADD2,GRADD,HDEL2)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      Parameter(MaxOff=200000)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
      DIMENSION GRADD2(3),R(3),HDEL2(3,3)
      Save Zero,Small,Two,Thou
      Data Zero/0.d0/,Small/1.d-13/,Two/2.d0/,Thou/1.d3/
C
      If(Iopt.eq.1)Goto 30
C
      RMIN = Thou
      RMAX = Zero
C
      DO 10 I = 1,NCENT
        CO(IXX+I) = R(1)-CO(IXC+I)
        CO(IYY+I) = R(2)-CO(IYC+I)
        CO(IZZ+I) = R(3)-CO(IZC+I)
        CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +
     +              CO(IYY+I)*CO(IYY+I) +
     +              CO(IZZ+I)*CO(IZZ+I)
        CO(IRR+I) = DSQRT(CO(IR2+I))
        IF(CO(IRR+I) .LT. RMIN) THEN
          RMIN = CO(IRR+I)
          MINR = I
          BANGLE(1) = CO(IXX+I)/(CO(IRR+I)+Small)
          BANGLE(2) = CO(IYY+I)/(CO(IRR+I)+Small)
          BANGLE(3) = CO(IZZ+I)/(CO(IRR+I)+Small)
        IF (CO(IRR+I) .GT. RMAX) THEN
          RMAX = CO(IRR+I)
        END IF
        END IF
   10 CONTINUE
C
      CALL GAUS4
C
      DO 20 I = 1,3
        GRADD2(I) = Zero
        DO 20 J = 1,3
          HDEL2(I,J) = Zero
20    CONTINUE
C
30    VALUE=Zero
      DO 40 I = 1,NMO
40      VALUE=VALUE+Two*CO(IP+I)*(CO(IPSI+I)*
     $(CO(IGXX+I)+CO(IGYY+I)+CO(IGZZ+I))+
     $CO(IGX+I)**2+CO(IGY+I)**2+CO(IGZ+I)**2)
C
      If(Iopt.eq.1)Goto 70
C
      DO 50 I = 1,NMO
        PSI=CO(IPSI+I)
        COMM1=Two*CO(IP+I)
        COMM2=CO(IGXX+I)+CO(IGYY+I)+CO(IGZZ+I)
        COMM3=CO(IGXXX+I)+CO(IGXYY+I)+CO(IGXZZ+I)
        COMM4=CO(IGXXY+I)+CO(IGYYY+I)+CO(IGYZZ+I)
        COMM5=CO(IGXXZ+I)+CO(IGYYZ+I)+CO(IGZZZ+I)
C
        GRADD2(1) = GRADD2(1) + COMM1*(Two*(CO(IGX+I)*CO(IGXX+I)+
     1  CO(IGY+I)*CO(IGXY+I)+CO(IGZ+I)*CO(IGXZ+I))+
     2  CO(IGX+I)*COMM2+PSI*COMM3)
        GRADD2(2)=GRADD2(2)+COMM1*(Two*(CO(IGX+I)*CO(IGXY+I)+
     1  CO(IGY+I)*CO(IGYY+I)+CO(IGZ+I)*CO(IGYZ+I))+
     2  CO(IGY+I)*COMM2+PSI*COMM4)
        GRADD2(3)=GRADD2(3)+COMM1*(Two*(CO(IGX+I)*CO(IGXZ+I)+
     1  CO(IGY+I)*CO(IGYZ+I)+CO(IGZ+I)*CO(IGZZ+I))+
     2  CO(IGZ+I)*COMM2+PSI*COMM5)
C
      HDEL2(1,1)=HDEL2(1,1)+COMM1*(Two*(CO(IGXX+I)*CO(IGXX+I)+
     1 CO(IGXY+I)*CO(IGXY+I)+CO(IGXZ+I)*CO(IGXZ+I)+
     2 CO(IGX+I)*CO(IGXXX+I)+CO(IGY+I)*CO(IGXXY+I)+
     3 CO(IGZ+I)*CO(IGXXZ+I)+CO(IGX+I)*COMM3)+CO(IGXX+I)*COMM2+
     4 PSI*(CO(IGXXXX+I)+CO(IGXXYY+I)+CO(IGXXZZ+I)))
      HDEL2(1,2)=HDEL2(1,2)+COMM1*(Two*(CO(IGXY+I)*
     1 (CO(IGXX+I)+CO(IGYY+I))+CO(IGXZ+I)*CO(IGYZ+I)+
     2 CO(IGX+I)*CO(IGXXY+I)+CO(IGY+I)*CO(IGXYY+I)+
     3 CO(IGZ+I)*CO(IGXYZ+I))+CO(IGXY+I)*COMM2+
     4 CO(IGX+I)*COMM4+CO(IGY+I)*COMM3+
     5 PSI*(CO(IGXXXY+I)+CO(IGXYYY+I)+CO(IGXYZZ+I)))
      HDEL2(1,3)=HDEL2(1,3)+COMM1*(Two*(CO(IGXZ+I)*
     1 (CO(IGXX+I)+CO(IGZZ+I))+CO(IGXY+I)*CO(IGYZ+I)+
     2 CO(IGX+I)*CO(IGXXZ+I)+CO(IGY+I)*CO(IGXYZ+I)+
     3 CO(IGZ+I)*CO(IGXZZ+I))+CO(IGXZ+I)*COMM2+
     4 CO(IGX+I)*COMM5+CO(IGZ+I)*COMM3+
     5 PSI*(CO(IGXXXZ+I)+CO(IGXYYZ+I)+CO(IGXZZZ+I)))
      HDEL2(2,2)=HDEL2(2,2)+COMM1*(Two*(CO(IGYY+I)*CO(IGYY+I)+
     1 CO(IGXY+I)*CO(IGXY+I)+CO(IGYZ+I)*CO(IGYZ+I)+
     2 CO(IGX+I)*CO(IGXYY+I)+CO(IGY+I)*CO(IGYYY+I)+
     3 CO(IGZ+I)*CO(IGYYZ+I)+CO(IGY+I)*COMM4)+CO(IGYY+I)*COMM2+
     4 PSI*(CO(IGYYYY+I)+CO(IGXXYY+I)+CO(IGYYZZ+I)))
      HDEL2(2,3)=HDEL2(2,3)+COMM1*(Two*(CO(IGYZ+I)*
     1 (CO(IGYY+I)+CO(IGZZ+I))+CO(IGXY+I)*CO(IGXZ+I)+
     2 CO(IGX+I)*CO(IGXYZ+I)+CO(IGY+I)*CO(IGYYZ+I)+
     3 CO(IGZ+I)*CO(IGYZZ+I))+CO(IGYZ+I)*COMM2+
     4 CO(IGY+I)*COMM5+CO(IGZ+I)*COMM4+
     5 PSI*(CO(IGXXYZ+I)+CO(IGYYYZ+I)+CO(IGYZZZ+I)))
      HDEL2(3,3)=HDEL2(3,3)+COMM1*(Two*(CO(IGZZ+I)*CO(IGZZ+I)+
     1 CO(IGXZ+I)*CO(IGXZ+I)+CO(IGYZ+I)*CO(IGYZ+I)+
     2 CO(IGX+I)*CO(IGXZZ+I)+CO(IGY+I)*CO(IGYZZ+I)+
     3 CO(IGZ+I)*CO(IGZZZ+I)+CO(IGZ+I)*COMM5)+CO(IGZZ+I)*COMM2+
     4 PSI*(CO(IGZZZZ+I)+CO(IGXXZZ+I)+CO(IGYYZZ+I)))
C
50    CONTINUE
C
      HDEL2(2,1)=HDEL2(1,2)
      HDEL2(3,1)=HDEL2(1,3)
      HDEL2(3,2)=HDEL2(2,3)
      GRADD=Zero
      DO 60 I=1,3
60    GRADD=GRADD+GRADD2(I)*GRADD2(I)
      GRADD=DSQRT(GRADD)
C
70    Continue
C
      RETURN
      END
      SUBROUTINE GRDKEG(IOpt,R,VALUE,W,GRAD,H)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter(MaxOff=200000)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
      DIMENSION W(3),R(3),H(3,3)
      Save Zero,Small,Pt5,Thou
      Data Zero/0.d0/,Small/1.d-13/,Pt5/0.5d0/,Thou/1.d3/
C
      If(Iopt.eq.1)Goto 30
C
      RMIN = Thou
      RMAX = Zero
C
      DO 10 I=1,NCENT
        CO(IXX+I) = R(1)-CO(IXC+I)
        CO(IYY+I) = R(2)-CO(IYC+I)
        CO(IZZ+I) = R(3)-CO(IZC+I)
        CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +
     1              CO(IYY+I)*CO(IYY+I) +
     2              CO(IZZ+I)*CO(IZZ+I)
        CO(IRR+I) = DSQRT(CO(IR2+I))
        IF(CO(IRR+I) .LT. RMIN) THEN
          RMIN = CO(IRR+I)
          MINR = I
          BANGLE(1) = CO(IXX+I)/(CO(IRR+I)+Small)
          BANGLE(2) = CO(IYY+I)/(CO(IRR+I)+Small)
          BANGLE(3) = CO(IZZ+I)/(CO(IRR+I)+Small)
        ELSE IF (CO(IRR+I) .GT. RMAX) THEN
          RMAX = CO(IRR+I)
        END IF
10    CONTINUE
C
      CALL GAUS4
C
      DO 20 I = 1,3
        W(I) = Zero
        DO 20 J = I,3
          H(I,J) = Zero
20    CONTINUE
C
30    Value=Zero
      DO 40 I = 1,NMO
40      Value=Value+CO(IP+I)*
     $        (CO(IGX+I)**2+CO(IGY+I)**2+CO(IGZ+I)**2)
      VALUE=Pt5*VALUE
C
      If(Iopt.eq.1)Goto 70
C
      DO 50 I = 1,NMO
        P0 = CO(IP+I)
        W(1) = W(1) + P0*(CO(IGX+I)*CO(IGXX+I)+
     +                    CO(IGY+I)*CO(IGXY+I)+
     +                    CO(IGZ+I)*CO(IGXZ+I))
        W(2) = W(2) + P0*(CO(IGY+I)*CO(IGYY+I)+
     +                    CO(IGX+I)*CO(IGXY+I)+
     +                    CO(IGZ+I)*CO(IGYZ+I))
        W(3) = W(3) + P0*(CO(IGZ+I)*CO(IGZZ+I)+
     +                    CO(IGY+I)*CO(IGYZ+I)+
     +                    CO(IGX+I)*CO(IGXZ+I))
        H(1,1) = H(1,1) + P0*(CO(IGXX+I)**2+
     +           CO(IGX+I)*CO(IGXXX+I)+CO(IGXY+I)**2+
     +           CO(IGY+I)*CO(IGXXY+I)+
     +           CO(IGXZ+I)**2+CO(IGZ+I)*CO(IGXXZ+I))
        H(1,2) = H(1,2) + P0*(CO(IGXX+I)*CO(IGXY+I)+
     +           CO(IGX+I)*CO(IGXXY+I)+CO(IGXY+I)*CO(IGYY+I)+
     +           CO(IGY+I)*CO(IGXYY+I)+
     +           CO(IGXZ+I)*CO(IGYZ+I)+CO(IGZ+I)*CO(IGXYZ+I))
        H(1,3) = H(1,3) + P0*(CO(IGXX+I)*CO(IGXZ+I)+
     +           CO(IGX+I)*CO(IGXXZ+I)+CO(IGXY+I)*CO(IGYZ+I)+
     +           CO(IGY+I)*CO(IGXYZ+I)+
     +           CO(IGXZ+I)*CO(IGZZ+I)+CO(IGZ+I)*CO(IGXZZ+I))
        H(2,2) = H(2,2) + P0*(CO(IGYY+I)**2+
     +           CO(IGY+I)*CO(IGYYY+I)+CO(IGXY+I)**2+
     +           CO(IGX+I)*CO(IGXYY+I)+
     +           CO(IGYZ+I)**2+CO(IGZ+I)*CO(IGYYZ+I))
        H(2,3) = H(2,3) + P0*(CO(IGXY+I)*CO(IGXZ+I)+
     +           CO(IGX+I)*CO(IGXYZ+I)+CO(IGYY+I)*CO(IGYZ+I)+
     +           CO(IGY+I)*CO(IGYYZ+I)+
     +           CO(IGYZ+I)*CO(IGZZ+I)+CO(IGZ+I)*CO(IGYZZ+I))
        H(3,3) = H(3,3) + P0*(CO(IGZZ+I)**2+
     +           CO(IGZ+I)*CO(IGZZZ+I)+CO(IGYZ+I)**2+
     +           CO(IGY+I)*CO(IGYZZ+I)+
     +           CO(IGXZ+I)**2+CO(IGX+I)*CO(IGXZZ+I))
50    CONTINUE
C
      GRAD=Zero
C
      DO 60 I = 1,3
        GRAD = GRAD+W(I)*W(I)
        DO 60 J = I,3
          H(J,I) = H(I,J)
60    CONTINUE
C
      GRAD = DSQRT(GRAD)
C
70    Continue
C
      RETURN
      END
      SUBROUTINE GRDKEK(Iopt,R,VALUE,W,GRADK,H)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      Parameter(MaxOff=200000)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
      DIMENSION RS(3),W(3),R(3),H(3,3)
      Save Zero,Small,Pt5,Two,Thou
      Data Zero/0.d0/,Small/1.d-13/,Pt5/0.5d0/,Two/2.d0/,Thou/1.d3/
C
      If(Iopt.eq.1)Goto 30
C
      RMIN = Thou
      RMAX = Zero
C
      DO 10 I = 1,NCENT
        CO(IXX+I) = R(1)-CO(IXC+I)
        CO(IYY+I) = R(2)-CO(IYC+I)
        CO(IZZ+I) = R(3)-CO(IZC+I)
        CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +
     +              CO(IYY+I)*CO(IYY+I) +
     +              CO(IZZ+I)*CO(IZZ+I)
        CO(IRR+I) = DSQRT(CO(IR2+I))
        IF(CO(IRR+I) .LT. RMIN) THEN
          RMIN = CO(IRR+I)
          MINR = I
          BANGLE(1) = CO(IXX+I)/(CO(IRR+I)+Small)
          BANGLE(2) = CO(IYY+I)/(CO(IRR+I)+Small)
          BANGLE(3) = CO(IZZ+I)/(CO(IRR+I)+Small)
        ELSE IF (CO(IRR+I) .GT. RMAX) THEN
          RMAX = CO(IRR+I)
        END IF
   10 CONTINUE
C
      CALL GAUS4
C
      DO 20 I = 1,3
        W(I) = Zero
        DO 20 J = 1,3
          H(I,J) = Zero
20    CONTINUE
C
30    Value=Zero
      DO 40 I = 1,NMO
40      Value=Value+CO(IP+I)*CO(IPSI+I)*(CO(IGXX+I)+
     $              CO(IGYY+I)+CO(IGZZ+I))
      Value=-Pt5*Value
C
      If(Iopt.eq.1)Goto 70
C
      DO 50 I = 1,NMO
        PSI=CO(IPSI+I)
        COMM1=-CO(IP+I)/Two
        COMM2=CO(IGXX+I)+CO(IGYY+I)+CO(IGZZ+I)
        COMM3=CO(IGXXX+I)+CO(IGXYY+I)+CO(IGXZZ+I)
        COMM4=CO(IGXXY+I)+CO(IGYYY+I)+CO(IGYZZ+I)
        COMM5=CO(IGXXZ+I)+CO(IGYYZ+I)+CO(IGZZZ+I)
C
        W(1) = W(1) + COMM1*(CO(IGX+I)*COMM2+PSI*COMM3)
        W(2) = W(2) + COMM1*(CO(IGY+I)*COMM2+PSI*COMM4)
        W(3) = W(3) + COMM1*(CO(IGZ+I)*COMM2+PSI*COMM5)
C
      H(1,1)=H(1,1)+COMM1*(TWO*CO(IGX+I)*COMM3+CO(IGXX+I)*COMM2+
     $ PSI*(CO(IGXXXX+I)+CO(IGXXYY+I)+CO(IGXXZZ+I)))
C
      H(1,2)=H(1,2)+COMM1*(CO(IGXY+I)*COMM2+CO(IGX+I)*COMM4+
     $ CO(IGY+I)*COMM3+PSI*(CO(IGXXXY+I)+CO(IGXYYY+I)+CO(IGXYZZ+I)))
C
      H(1,3)=H(1,3)+COMM1*(CO(IGXZ+I)*COMM2+CO(IGX+I)*COMM5+
     $ CO(IGZ+I)*COMM3+PSI*(CO(IGXXXZ+I)+CO(IGXYYZ+I)+CO(IGXZZZ+I)))
C
      H(2,2)=H(2,2)+COMM1*(TWO*CO(IGY+I)*COMM4+CO(IGYY+I)*COMM2+
     $ PSI*(CO(IGXXYY+I)+CO(IGYYYY+I)+CO(IGYYZZ+I)))
C
      H(2,3)=H(2,3)+COMM1*(CO(IGYZ+I)*COMM2+CO(IGY+I)*COMM5+
     $ CO(IGZ+I)*COMM4+PSI*(CO(IGXXYZ+I)+CO(IGYYYZ+I)+CO(IGYZZZ+I)))
C
      H(3,3)=H(3,3)+COMM1*(TWO*CO(IGZ+I)*COMM5+CO(IGZZ+I)*COMM2+
     $ PSI*(CO(IGXXZZ+I)+CO(IGYYZZ+I)+CO(IGZZZZ+I)))
C
50    CONTINUE
C
      H(2,1)=H(1,2)
      H(3,1)=H(1,3)
      H(3,2)=H(2,3)
      GRADK=Zero
      DO 60 I=1,3
60    GRADK=GRADK+W(I)*W(I)
      GRADK=DSQRT(GRADK)
C
70    Continue
C
      RETURN
      END
      SUBROUTINE GRDRHO(Iopt,R,VALUE,W,GRAD,H,SG)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter(MaxOff=200000)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
      DIMENSION W(3),R(3),H(3,3),SG(3,3)
      Save Zero,Small,pt5,Two,Thou
      Data Zero/0.d0/,Small/1.d-13/,pt5/0.5d0/,Two/2.d0/,Thou/1.d3/
C
      If(Iopt.eq.2)Goto 30
C
      RMIN = Thou
      RMAX = Zero
C
      DO 10 I=1,NCENT
        CO(IXX+I) = R(1)-CO(IXC+I)
        CO(IYY+I) = R(2)-CO(IYC+I)
        CO(IZZ+I) = R(3)-CO(IZC+I)
        CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +
     1              CO(IYY+I)*CO(IYY+I) +
     2              CO(IZZ+I)*CO(IZZ+I)
        CO(IRR+I) = DSQRT(CO(IR2+I))
        IF(CO(IRR+I) .LT. RMIN) THEN
          RMIN = CO(IRR+I)
          MINR = I
          BANGLE(1) = CO(IXX+I)/(CO(IRR+I)+Small)
          BANGLE(2) = CO(IYY+I)/(CO(IRR+I)+Small)
          BANGLE(3) = CO(IZZ+I)/(CO(IRR+I)+Small)
        ELSE IF (CO(IRR+I) .GT. RMAX) THEN
          RMAX = CO(IRR+I)
        END IF
10    CONTINUE
C
      CALL GAUS4
C
      DO 20 I = 1,3
        W(I) = Zero
        DO 20 J = I,3
          H(I,J) = Zero
          SG(I,J) = Zero
20    CONTINUE
C
30    Value=Zero
      DO 40 I = 1,NMO
        P0 = CO(IP+I)
        P1 = P0*CO(IPSI+I)
        Value=Value+P1*CO(IPSI+I)
        W(1) = W(1) + P1*CO(IGX+I)
        W(2) = W(2) + P1*CO(IGY+I)
40      W(3) = W(3) + P1*CO(IGZ+I)
C
        If(Iopt.eq.2)Goto 70
C
      If(IOpt.eq.1)Then
      DO 50 I = 1,NMO
        P0 = CO(IP+I)
        P1 = P0*CO(IPSI+I)
        H(1,1) = H(1,1) + P0*CO(IGX+I)**2
        H(1,2) = H(1,2) + P0*CO(IGX+I)*CO(IGY+I)
        H(1,3) = H(1,3) + P0*CO(IGX+I)*CO(IGZ+I)
        H(2,2) = H(2,2) + P0*CO(IGY+I)**2
        H(2,3) = H(2,3) + P0*CO(IGY+I)*CO(IGZ+I)
        H(3,3) = H(3,3) + P0*CO(IGZ+I)**2
        SG(1,1) = SG(1,1) + P1*CO(IGXX+I)
        SG(1,2) = SG(1,2) + P1*CO(IGXY+I)
        SG(1,3) = SG(1,3) + P1*CO(IGXZ+I)
        SG(2,2) = SG(2,2) + P1*CO(IGYY+I)
        SG(2,3) = SG(2,3) + P1*CO(IGYZ+I)
50      SG(3,3) = SG(3,3) + P1*CO(IGZZ+I)
      Endif 
C
      Do 60 I=1,3
        DO 60 J = I,3
          DM = H(I,J)
          H(I,J) = (DM+SG(I,J))*Two
          SG(I,J) = (-DM+SG(I,J))*pt5
          H(J,I) = H(I,J)
          SG(J,I) = SG(I,J)
60    CONTINUE
C
70    GRAD=Zero
      DO 80  I = 1,3
      W(I) = W(I)*Two
80      GRAD = GRAD+W(I)*W(I)
      GRAD = DSQRT(GRAD)
C
      RETURN
      END
      SUBROUTINE GRDV(Iopt,R,VALUE,GRADD2,GRADD,HDEL2)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      Parameter(MaxOff=200000)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
      DIMENSION GRADD2(3),R(3),HDEL2(3,3)
      Save Zero,Small,Pt5,Two,Thou
      Data Zero/0.d0/,Small/1.d-13/,Pt5/0.5d0/,Two/2.d0/,Thou/1.d3/
C
      If(Iopt.eq.1)Goto 30
C
      RMIN = Thou
      RMAX = Zero
C
      DO 10 I = 1,NCENT
        CO(IXX+I) = R(1)-CO(IXC+I)
        CO(IYY+I) = R(2)-CO(IYC+I)
        CO(IZZ+I) = R(3)-CO(IZC+I)
        CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +
     +              CO(IYY+I)*CO(IYY+I) +
     +              CO(IZZ+I)*CO(IZZ+I)
        CO(IRR+I) = DSQRT(CO(IR2+I))
        IF(CO(IRR+I) .LT. RMIN) THEN
          RMIN = CO(IRR+I)
          MINR = I
          BANGLE(1) = CO(IXX+I)/(CO(IRR+I)+Small)
          BANGLE(2) = CO(IYY+I)/(CO(IRR+I)+Small)
          BANGLE(3) = CO(IZZ+I)/(CO(IRR+I)+Small)
        IF (CO(IRR+I) .GT. RMAX) THEN
          RMAX = CO(IRR+I)
        END IF
        END IF
   10 CONTINUE
C
      CALL GAUS4
C
      DO 20 I = 1,3
        GRADD2(I) = Zero
        DO 20 J = 1,3
          HDEL2(I,J) = Zero
20    CONTINUE
C
30    VALUE=Zero
      DO 40 I = 1,NMO
40      VALUE=VALUE+Pt5*CO(IP+I)*(CO(IPSI+I)*
     $(CO(IGXX+I)+CO(IGYY+I)+CO(IGZZ+I))-
     $(CO(IGX+I)**2+CO(IGY+I)**2+CO(IGZ+I)**2))
C
      If(Iopt.eq.1)Goto 70
C
      DO 50 I = 1,NMO
        PSI=CO(IPSI+I)
        COMM1=Pt5*CO(IP+I)
        COMM2=CO(IGXX+I)+CO(IGYY+I)+CO(IGZZ+I)
        COMM3=CO(IGXXX+I)+CO(IGXYY+I)+CO(IGXZZ+I)
        COMM4=CO(IGXXY+I)+CO(IGYYY+I)+CO(IGYZZ+I)
        COMM5=CO(IGXXZ+I)+CO(IGYYZ+I)+CO(IGZZZ+I)
C
        GRADD2(1) = GRADD2(1) + COMM1*(-Two*(CO(IGX+I)*CO(IGXX+I)+
     1  CO(IGY+I)*CO(IGXY+I)+CO(IGZ+I)*CO(IGXZ+I))+
     2  CO(IGX+I)*COMM2+PSI*COMM3)
        GRADD2(2)=GRADD2(2)+COMM1*(-Two*(CO(IGX+I)*CO(IGXY+I)+
     1  CO(IGY+I)*CO(IGYY+I)+CO(IGZ+I)*CO(IGYZ+I))+
     2  CO(IGY+I)*COMM2+PSI*COMM4)
        GRADD2(3)=GRADD2(3)+COMM1*(-Two*(CO(IGX+I)*CO(IGXZ+I)+
     1  CO(IGY+I)*CO(IGYZ+I)+CO(IGZ+I)*CO(IGZZ+I))+
     2  CO(IGZ+I)*COMM2+PSI*COMM5)
C
      HDEL2(1,1)=HDEL2(1,1)+COMM1*(-Two*(CO(IGXX+I)*CO(IGXX+I)+
     1 CO(IGXY+I)*CO(IGXY+I)+CO(IGXZ+I)*CO(IGXZ+I)+
     2 CO(IGX+I)*CO(IGXXX+I)+CO(IGY+I)*CO(IGXXY+I)+
     3 CO(IGZ+I)*CO(IGXXZ+I)-CO(IGX+I)*COMM3)+CO(IGXX+I)*COMM2+
     4 PSI*(CO(IGXXXX+I)+CO(IGXXYY+I)+CO(IGXXZZ+I)))
      HDEL2(1,2)=HDEL2(1,2)+COMM1*(-Two*(CO(IGXY+I)*
     1 (CO(IGXX+I)+CO(IGYY+I))+CO(IGXZ+I)*CO(IGYZ+I)+
     2 CO(IGX+I)*CO(IGXXY+I)+CO(IGY+I)*CO(IGXYY+I)+
     3 CO(IGZ+I)*CO(IGXYZ+I))+CO(IGXY+I)*COMM2+
     4 CO(IGX+I)*COMM4+CO(IGY+I)*COMM3+
     5 PSI*(CO(IGXXXY+I)+CO(IGXYYY+I)+CO(IGXYZZ+I)))
      HDEL2(1,3)=HDEL2(1,3)+COMM1*(-Two*(CO(IGXZ+I)*
     1 (CO(IGXX+I)+CO(IGZZ+I))+CO(IGXY+I)*CO(IGYZ+I)+
     2 CO(IGX+I)*CO(IGXXZ+I)+CO(IGY+I)*CO(IGXYZ+I)+
     3 CO(IGZ+I)*CO(IGXZZ+I))+CO(IGXZ+I)*COMM2+
     4 CO(IGX+I)*COMM5+CO(IGZ+I)*COMM3+
     5 PSI*(CO(IGXXXZ+I)+CO(IGXYYZ+I)+CO(IGXZZZ+I)))
      HDEL2(2,2)=HDEL2(2,2)+COMM1*(-Two*(CO(IGYY+I)*CO(IGYY+I)+
     1 CO(IGXY+I)*CO(IGXY+I)+CO(IGYZ+I)*CO(IGYZ+I)+
     2 CO(IGX+I)*CO(IGXYY+I)+CO(IGY+I)*CO(IGYYY+I)+
     3 CO(IGZ+I)*CO(IGYYZ+I)-CO(IGY+I)*COMM4)+CO(IGYY+I)*COMM2+
     4 PSI*(CO(IGYYYY+I)+CO(IGXXYY+I)+CO(IGYYZZ+I)))
      HDEL2(2,3)=HDEL2(2,3)+COMM1*(-Two*(CO(IGYZ+I)*
     1 (CO(IGYY+I)+CO(IGZZ+I))+CO(IGXY+I)*CO(IGXZ+I)+
     2 CO(IGX+I)*CO(IGXYZ+I)+CO(IGY+I)*CO(IGYYZ+I)+
     3 CO(IGZ+I)*CO(IGYZZ+I))+CO(IGYZ+I)*COMM2+
     4 CO(IGY+I)*COMM5+CO(IGZ+I)*COMM4+
     5 PSI*(CO(IGXXYZ+I)+CO(IGYYYZ+I)+CO(IGYZZZ+I)))
      HDEL2(3,3)=HDEL2(3,3)+COMM1*(-Two*(CO(IGZZ+I)*CO(IGZZ+I)+
     1 CO(IGXZ+I)*CO(IGXZ+I)+CO(IGYZ+I)*CO(IGYZ+I)+
     2 CO(IGX+I)*CO(IGXZZ+I)+CO(IGY+I)*CO(IGYZZ+I)+
     3 CO(IGZ+I)*CO(IGZZZ+I)-CO(IGZ+I)*COMM5)+CO(IGZZ+I)*COMM2+
     4 PSI*(CO(IGZZZZ+I)+CO(IGXXZZ+I)+CO(IGYYZZ+I)))
C
50    CONTINUE
C
      HDEL2(2,1)=HDEL2(1,2)
      HDEL2(3,1)=HDEL2(1,3)
      HDEL2(3,2)=HDEL2(2,3)
      GRADD=Zero
      DO 60 I=1,3
60    GRADD=GRADD+GRADD2(I)*GRADD2(I)
      GRADD=DSQRT(GRADD)
C
70    Continue
C
      RETURN
      END
      SUBROUTINE GRDVNE(Iopt,R,VALUE,W,GRAD,H)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter(MaxOff=200000)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
      DIMENSION W(3),R(3),H(3,3)
      Save Zero,Small,One,Three,Thou
      Data Zero/0.d0/,Small/1.d-13/,One/1.d0/,Three/3.d0/,Thou/1.d3/
C
      If(Iopt.eq.1)Goto 30
C
      RMIN = Thou
      RMAX = Zero
C
      DO 10 I=1,NCENT
        CO(IXX+I) = R(1)-CO(IXC+I)
        CO(IYY+I) = R(2)-CO(IYC+I)
        CO(IZZ+I) = R(3)-CO(IZC+I)
        CO(IR2+I) = CO(IXX+I)*CO(IXX+I) +
     1              CO(IYY+I)*CO(IYY+I) +
     2              CO(IZZ+I)*CO(IZZ+I)
        CO(IRR+I) = DSQRT(CO(IR2+I))
        IF(CO(IRR+I) .LT. RMIN) THEN
          RMIN = CO(IRR+I)
          MINR = I
          BANGLE(1) = CO(IXX+I)/(CO(IRR+I)+Small)
          BANGLE(2) = CO(IYY+I)/(CO(IRR+I)+Small)
          BANGLE(3) = CO(IZZ+I)/(CO(IRR+I)+Small)
        ELSE IF (CO(IRR+I) .GT. RMAX) THEN
          RMAX = CO(IRR+I)
        END IF
10    CONTINUE
C
      DO 20 I = 1,3
        W(I) = Zero
        DO 20 J = I,3
          H(I,J) = Zero
20    CONTINUE
C
30    Value=Zero
      DO 40 I=1,NCENT
      OVERR=One/CO(IRR+I)
40    Value=Value+CO(ICHARG+I)*OVERR
C
      If(Iopt.eq.1)Goto 70
C
      DO 50 I=1,NCENT
      OVERR=One/CO(IRR+I)
      OVERR3=OVERR**3
      OVERR5=OVERR**5
      W(1) = W(1) - CO(ICHARG+I)*CO(IXX+I)*OVERR3
      W(2) = W(2) - CO(ICHARG+I)*CO(IYY+I)*OVERR3
      W(3) = W(3) - CO(ICHARG+I)*CO(IZZ+I)*OVERR3
      H(1,1) = H(1,1) - CO(ICHARG+I)*(OVERR3-
     +                  Three*CO(IXX+I)**2*OVERR5)
      H(2,2) = H(2,2) - CO(ICHARG+I)*(OVERR3-
     +                  Three*CO(IYY+I)**2*OVERR5)
      H(3,3) = H(3,3) - CO(ICHARG+I)*(OVERR3-
     +                  Three*CO(IZZ+I)**2*OVERR5)
      H(1,2) = H(1,2) + Three*CO(ICHARG+I)*CO(IXX+I)*CO(IYY+I)*
     +                  OVERR5
      H(1,3) = H(1,3) + Three*CO(ICHARG+I)*CO(IXX+I)*CO(IZZ+I)*
     +                  OVERR5
50    H(2,3) = H(2,3) + Three*CO(ICHARG+I)*CO(IYY+I)*CO(IZZ+I)*
     +                  OVERR5
C
      GRAD=Zero
C
      DO 60 I = 1,3
        GRAD = GRAD+W(I)*W(I)
        DO 60 J = I,3
          H(J,I) = H(I,J)
60    CONTINUE
C
      GRAD = DSQRT(GRAD)
C
70    Continue
C
      RETURN
      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
      DOUBLE PRECISION DX(1),DMAX
      INTEGER I,INCX,IX,N
      IDAMAX = 0
      IF( N .LT. 1 ) RETURN
      IDAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5
         IDAMAX = I
         DMAX = DABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
         IF(DABS(DX(I)).LE.DMAX) GO TO 30
         IDAMAX = I
         DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END
      SUBROUTINE MAKNAME(I,STRING,L,EXT)
      CHARACTER*(*) STRING,EXT
      INTEGER I,J,L
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
      END
      SUBROUTINE MULTI1(W,R,W1,I,DS,S)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /C2/ EPSD,AB(9,9),AM(9,10)
      DIMENSION W1(3,I),W(3),R(3),R1(3),HESS(3,3),SG(3,3)
C
C     PERFORM A PREDICTOR STEP OF ORDER I.
C
      DO 100 L = 1,3
        H = 0.0D0
        DO 110 J = 1,I
          H = H+AB(I,J)*W1(L,J)
110     CONTINUE
        R(L) = W(L)+H*DS
100   CONTINUE
C
C     PERFORM A CORRECTOR STEP OF ORDER I+1.
C
      CALL GRDRHO(0,R,RHO,R1,Gnorm,HESS,SG)
      DO 120 L = 1,3
        H = AM(I,1)*R1(L)/Gnorm
        DO 130 J = 1,I
          H = H+AM(I,J+1)*W1(L,J)
130     CONTINUE
        R(L) = W(L)+H*DS
120   CONTINUE
      RETURN
      END
      SUBROUTINE NEWTON(R,IFAIL,IFUNC,IWHOLE,NITER,INuc)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter(MaxSch=20)
      COMMON /UNITS/ INPT, IOUT, IWFN, IWLP
      Common /What/ IWhat
      Common /Options/ ICut,IPrint,Eps,EpsNuc,Dmp,Dmpnuc
      DIMENSION W(3),R1(3),IW(3),XX(2),R(3),SG(3,3),H(3,3)
      Save Zero,pt1,pt05,pt025,One
      Data Zero/0.0d0/,pt1/0.1d0/,pt05/0.05d0/,pt025/0.025d0/,
     $  One/1.0d0/,pt25/0.25d0/
1000  FORMAT(' STEP',3X,'CURRENT COORDINATES ',40X,'|GRAD(Rho)|',/)
1001  FORMAT(' STEP',3X,'CURRENT COORDINATES ',40X,
     $  '|GRAD(DEL**2(Rho))|',/)
1002  FORMAT(' STEP',3X,'CURRENT COORDINATES ',40X,
     $  '|GRAD(G)|',/)
1003  FORMAT(' STEP',3X,'CURRENT COORDINATES ',40X,
     $  '|GRAD(K)|',/)
1004  FORMAT(' STEP',3X,'CURRENT COORDINATES ',40X,
     $  '|GRAD(Vnuc)|',/)
1005  FORMAT(' STEP',3X,'CURRENT COORDINATES ',40X,
     $  '|GRAD(V)|',/)
1010  FORMAT(1X,I4,3X,1P4E18.10)
1020  FORMAT(' NUMBER OF NEWTON-RAPHSON ITERATIONS : ',I4)
1030  FORMAT(' CRITICAL POINT NOT YET FOUND.  CONTINUE SEARCH',
     $' ? (0=no/1=yes) ',$)
1040  FORMAT(' TRY AGAIN ')
C
      Cutoff=eps
      If(Inuc.eq.1)cutoff=epsnuc
      Damp=Dmp
      If(Inuc.eq.1)Damp=Dmpnuc
C
      If(Ifunc.eq.1)IWhat=1
      If(Ifunc.eq.2)IWhat=3
      If(Ifunc.eq.3)IWhat=2
      If(Ifunc.eq.4)IWhat=3
      If(Ifunc.eq.5)IWhat=0
      If(Ifunc.eq.6)IWhat=3
C
      IF(IFUNC.EQ.1)CALL GRDRHO(1,R,RHO,W,GRAD,H,SG)
      IF(IFUNC.EQ.2)CALL GRDD2R(0,R,D2RHO,W,GRAD,H)
      IF(IFUNC.EQ.3)CALL GRDKEG(0,R,GEE,W,GRAD,H)
      IF(IFUNC.EQ.4)CALL GRDKEK(0,R,QUAY,W,GRAD,H)
      IF(IFUNC.EQ.5)CALL GRDVNE(0,R,VNE,W,GRAD,H)
      IF(IFUNC.EQ.6)CALL GRDV(0,R,VEE,W,GRAD,H)
C
      IS = 0
      If(Iprint.eq.1)Then
      IF (IFAIL .EQ. 0) Then
      If(Iwhole.eq.0)Then
      If(Ifunc.eq.1) WRITE (IOUT,1000)
      If(Ifunc.eq.2) WRITE (IOUT,1001)
      If(Ifunc.eq.3) WRITE (IOUT,1002)
      If(Ifunc.eq.4) WRITE (IOUT,1003)
      If(Ifunc.eq.5) WRITE (IOUT,1004)
      If(Ifunc.eq.6) WRITE (IOUT,1005)
      WRITE (IOUT,1010) IS,R(1),R(2),R(3),GRAD
      Endif
      Endif
      Endif
C
C    BEGIN NEWTON RAPHSON SEARCH
C
99    DO 100 I = 1,MaxSch
        IS = IS + 1
        CALL DGECO(H,3,3,IW,ER,R1)
        CALL DGEDI(H,3,3,IW,XX,R1,1)
        DO 110 J = 1,3
          R1(J) = Zero
          DO 110 K = 1,3
            R1(J) = R1(J) + H(J,K)*W(K)
            DJ = DABS(R1(J))
            IF (DJ.GT.Damp) R1(J) = Damp*R1(J)/DJ
110     CONTINUE
C
C    GENERATE SHIFT VECTOR
C
        DO 120 J = 1,3
          R(J) = R(J) - R1(J)
120     CONTINUE
C
      IF(IFUNC.EQ.1)CALL GRDRHO(1,R,RHO,W,GRAD,H,SG)
      IF(IFUNC.EQ.2)CALL GRDD2R(0,R,D2RHO,W,GRAD,H)
      IF(IFUNC.EQ.3)CALL GRDKEG(0,R,GEE,W,GRAD,H)
      IF(IFUNC.EQ.4)CALL GRDKEK(0,R,QUAY,W,GRAD,H)
      IF(IFUNC.EQ.5)CALL GRDVNE(0,R,VNE,W,GRAD,H)
      IF(IFUNC.EQ.6)CALL GRDV(0,R,VEE,W,GRAD,H)
C
        IF (IFAIL .EQ. 0.and.Iwhole.eq.0.and.iprint.eq.1) Then
        WRITE (IOUT,1010) IS,R(1),R(2),R(3),GRAD
        Endif
C
        IF (GRAD .LE. Cutoff) THEN
          NITER=I
          Goto 150
        END IF
100   CONTINUE
C
C    SHALL WE CONTINUE ?
C
      IF (IFAIL .EQ. 1) Goto 150
C
130   CONTINUE
      IYN=0
      IF(IWHOLE.EQ.0)Then
      WRITE (IOUT,1030)
      READ (INPT,*) IYN
      Endif
      IF (IYN .EQ. 0) THEN
        IFAIL = 1
        Goto 150
      ELSE IF (IYN .EQ. 1) THEN
        GOTO 99
      ELSE
        WRITE (IOUT,1040)
        GOTO 130
      END IF
C
150   Continue
C
      RETURN
      END
      Subroutine PrimeA(IWhich,IMode,X,Y,Z)
C
      Implicit Double Precision (A-H,O-Z)
C
      Common /AFG/ A0,AX,AY,AZ,AXX,AXY,AXZ,AYY,AYZ,AZZ,AXXX,
     $  AXXY,AXXZ,AXYY,AXZZ,AXYZ,AYYY,AYYZ,AYZZ,AZZZ,AXXXX,AXXXY,
     $  AXXXZ,AXXYY,AXXZZ,AXYYY,AXZZZ,AXXYZ,AXYYZ,AXYZZ,AYYYY,
     $  AYYYZ,AYYZZ,AYZZZ,AZZZZ,F0,FX,FY,FZ,FXX,FXY,FXZ,FYY,FYZ,
     $  FZZ,FXXX,FXXY,FXXZ,FXYY,FXZZ,FXYZ,FYYY,FYYZ,FYZZ,FZZZ,
     $  FXXXX,FXXXY,FXXXZ,FXXYY,FXXZZ,FXYYY,FXZZZ,FXXYZ,FXYYZ,
     $  FXYZZ,FYYYY,FYYYZ,FYYZZ,FYZZZ,FZZZZ,G0,GX,GY,GZ,GXX,GXY,
     $  GXZ,GYY,GYZ,GZZ,GXXX,GXXY,GXXZ,GXYY,GXZZ,GXYZ,GYYY,GYYZ,
     $  GYZZ,GZZZ,GXXXX,GXXXY,GXXXZ,GXXYY,GXXZZ,GXYYY,GXZZZ,
     $  GXXYZ,GXYYZ,GXYZZ,GYYYY,GYYYZ,GYYZZ,GYZZZ,GZZZZ
      Save Zero,One,Two,Three,Six
      Data Zero/0.0d0/,One/1.0d0/,Two/2.0d0/,Three/3.0d0/,Six/6.0d0/
C
      If(IWhich.eq.1)Then
        If(IMode.eq.1)Then
          A0=One
        ElseIf(Imode.eq.0)Then
          A0=Zero
        Endif
C  
      ElseIf(Iwhich.eq.2)Then
        If(IMode.eq.1)Then
          A0=X
          Ax=One
        ElseIf(Imode.eq.0)Then
          A0=Zero
          Ax=Zero
        Endif
C
      ElseIf(Iwhich.eq.3)Then
        If(IMode.eq.1)Then
          A0=Y
          Ay=One
        ElseIf(Imode.eq.0)Then
          A0=Zero
          Ay=Zero
        Endif
C
      ElseIf(Iwhich.eq.4)Then
        If(IMode.eq.1)Then
          A0=Z
          AZ=One
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AZ=Zero
        Endif
C
      ElseIf(Iwhich.eq.5)Then
        If(IMode.eq.1)Then
          A0=X**2
          AX=Two*X
          AXX=Two
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AX=Zero
          AXX=Zero
        Endif
C
      ElseIf(Iwhich.eq.6)Then
        If(IMode.eq.1)Then
          A0=Y**2
          AY=Two*Y
          AYY=Two
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AY=Zero
          AYY=Zero
        Endif
C
      ElseIf(Iwhich.eq.7)Then
        If(IMode.eq.1)Then
          A0=Z**2
          AZ=Two*Z
          AZZ=Two
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AZ=Zero
          AZZ=Zero
        Endif
C
      ElseIf(Iwhich.eq.8)Then
        If(IMode.eq.1)Then
          A0=X*Y
          AX=Y
          AY=X
          AXY=One
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AX=Zero
          AY=Zero
          AXY=Zero
        Endif
C
      ElseIf(Iwhich.eq.9)Then
        If(IMode.eq.1)Then
          A0=X*Z
          AX=Z
          AZ=X
          AXZ=One
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AX=Zero
          AZ=Zero
          AXZ=Zero
        Endif
C
      ElseIf(Iwhich.eq.10)Then
        If(IMode.eq.1)Then
          A0=Y*Z
          AY=Z
          AZ=Y
          AYZ=One
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AY=Zero
          AZ=Zero
          AYZ=Zero
        Endif
C
      ElseIf(Iwhich.eq.11)Then
        If(IMode.eq.1)Then
          A0=X**3
          AX=Three*X**2
          AXX=Six*X
          AXXX=Six
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AX=Zero
          AXX=Zero
          AXXX=Zero
        Endif
C
      ElseIf(Iwhich.eq.12)Then
        If(IMode.eq.1)Then
          A0=Y**3
          AY=Three*Y**2
          AYY=Six*Y
          AYYY=Six
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AY=Zero
          AYY=Zero
          AYYY=Zero
        Endif
C
      ElseIf(Iwhich.eq.13)Then
        If(IMode.eq.1)Then
          A0=Z**3
          AZ=Three*Z**2
          AZZ=Six*Z
          AZZZ=Six
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AZ=Zero
          AZZ=Zero
          AZZZ=Zero
        Endif
C
      ElseIf(Iwhich.eq.14)Then
        If(IMode.eq.1)Then
          A0=X**2*Y
          AX=Two*X*Y
          AY=X**2
          AXX=Two*Y
          AXY=Two*X
          AXXY=Two
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AX=Zero
          AY=Zero
          AXX=Zero
          AXY=Zero
          AXXY=Zero
        Endif
C
      ElseIf(Iwhich.eq.15)Then
        If(IMode.eq.1)Then
          A0=X**2*Z
          AX=Two*X*Z
          AZ=X**2
          AXX=Two*Z
          AXZ=Two*X
          AXXZ=Two
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AX=Zero
          AZ=Zero
          AXX=Zero
          AXZ=Zero
          AXXZ=Zero
        Endif
C
      ElseIf(Iwhich.eq.16)Then
        If(IMode.eq.1)Then
          A0=Y**2*Z
          AY=Two*Y*Z
          AZ=Y**2
          AYY=Two*Z
          AYZ=Two*Y
          AYYZ=Two
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AY=Zero
          AZ=Zero
          AYY=Zero
          AYZ=Zero
          AYYZ=Zero
        Endif
C
      ElseIf(Iwhich.eq.17)Then
        If(IMode.eq.1)Then
          A0=X*Y**2
          AX=Y**2
          AY=Two*Y*X
          AYY=Two*X
          AXY=Two*Y
          AXYY=Two
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AX=Zero
          AY=Zero
          AYY=Zero
          AXY=Zero
          AXYY=Zero
        Endif
C
      ElseIf(Iwhich.eq.18)Then
        If(IMode.eq.1)Then
          A0=X*Z**2
          AX=Z**2
          AZ=Two*Z*X
          AZZ=Two*X
          AXZ=Two*Z
          AXZZ=Two
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AX=Zero
          AZ=Zero
          AZZ=Zero
          AXZ=Zero
          AXZZ=Zero
        Endif
C
      ElseIf(Iwhich.eq.19)Then
        If(IMode.eq.1)Then
          A0=Y*Z**2
          AY=Z**2
          AZ=Two*Z*Y
          AZZ=Two*Y
          AYZ=Two*Z
          AYZZ=Two
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AY=Zero
          AZ=Zero
          AZZ=Zero
          AYZ=Zero
          AYZZ=Zero
        Endif
C
      ElseIf(Iwhich.eq.20)Then
        If(IMode.eq.1)Then
          A0=X*Y*Z
          AX=Y*Z
          AY=X*Z
          AZ=X*Y
          AXY=Z
          AXZ=Y
          AYZ=X
          AXYZ=One
        ElseIf(Imode.eq.0)Then
          A0=Zero
          AX=Zero
          AY=Zero
          AZ=Zero
          AXY=Zero
          AXZ=Zero
          AYZ=Zero
          AXYZ=Zero
        Endif
      Endif
C
      Return
      End
      Subroutine PrimeF(Ent,EXPON,X,Y,Z)
C
      Implicit Double Precision (A-H,O-Z)
C
      Common /AFG/ A0,AX,AY,AZ,AXX,AXY,AXZ,AYY,AYZ,AZZ,AXXX,
     $  AXXY,AXXZ,AXYY,AXZZ,AXYZ,AYYY,AYYZ,AYZZ,AZZZ,AXXXX,AXXXY,
     $  AXXXZ,AXXYY,AXXZZ,AXYYY,AXZZZ,AXXYZ,AXYYZ,AXYZZ,AYYYY,
     $  AYYYZ,AYYZZ,AYZZZ,AZZZZ,F0,FX,FY,FZ,FXX,FXY,FXZ,FYY,FYZ,
     $  FZZ,FXXX,FXXY,FXXZ,FXYY,FXZZ,FXYZ,FYYY,FYYZ,FYZZ,FZZZ,
     $  FXXXX,FXXXY,FXXXZ,FXXYY,FXXZZ,FXYYY,FXZZZ,FXXYZ,FXYYZ,
     $  FXYZZ,FYYYY,FYYYZ,FYYZZ,FYZZZ,FZZZZ,G0,GX,GY,GZ,GXX,GXY,
     $  GXZ,GYY,GYZ,GZZ,GXXX,GXXY,GXXZ,GXYY,GXZZ,GXYZ,GYYY,GYYZ,
     $  GYZZ,GZZZ,GXXXX,GXXXY,GXXXZ,GXXYY,GXXZZ,GXYYY,GXZZZ,
     $  GXXYZ,GXYYZ,GXYZZ,GYYYY,GYYYZ,GYYZZ,GYZZZ,GZZZZ
      Common /What/ IWhat
C
      Save Two,Three
      Data Two/2.0d0/,Three/3.d0/
C
      toalp=-Two*Ent
      toalpe=toalp*expon
C
      F0=EXPON
      FX=toalpe*X
      FY=toalpe*Y
      FZ=toalpe*Z
      If(Iwhat.gt.0)Then
      FXX=toalp*X*FX+toalpe
      FXY=toalp*Y*FX
      FXZ=Toalp*Z*FX
      FYY=toalp*Y*FY+toalpe
      FYZ=Toalp*Z*FY
      FZZ=toalp*Z*FZ+toalpe
      Endif
      If(Iwhat.gt.1)Then
      FXXX=Toalp*(Two*FX+X*FXX)
      FYYY=Toalp*(Two*FY+Y*FYY)
      FZZZ=Toalp*(Two*FZ+Z*FZZ)
      FXXY=toalp*Y*FXX
      FXXZ=toalp*Z*FXX
      FXYY=toalp*X*FYY
      FXZZ=toalp*X*FZZ
      FXYZ=toalp*Z*FXY
      FYYZ=toalp*Z*FYY
      FYZZ=toalp*Y*FZZ
      Endif
      If(Iwhat.gt.2)Then
      FXXXX=Toalp*(X*FXXX+Three*FXX)
      FXXXY=FXXX*Toalp*Y
      FXXXZ=FXXX*Toalp*Z
      FXXYY=Toalp*(FXX+Y*FXXY)
      FXXZZ=Toalp*(FXX+Z*FXXZ)
      FXYYY=FYYY*ToAlp*X
      FXZZZ=FZZZ*ToAlp*X
      FXXYZ=FXXY*ToAlp*Z
      FXYYZ=FXYY*ToAlp*Z
      FXYZZ=FXZZ*ToAlp*Y
      FYYYY=Toalp*(Y*FYYY+Three*FYY)
      FYYYZ=FYYY*ToAlp*Z
      FYYZZ=ToAlp*(FYY+Z*FYYZ)
      FYZZZ=FZZZ*ToAlp*Y
      FZZZZ=Toalp*(Z*FZZZ+Three*FZZ)
      Endif
C
      Return
      End
      Subroutine PrimeG
C
      Implicit Double Precision (A-H,O-Z)
C
      Common /AFG/ A0,AX,AY,AZ,AXX,AXY,AXZ,AYY,AYZ,AZZ,AXXX,
     $  AXXY,AXXZ,AXYY,AXZZ,AXYZ,AYYY,AYYZ,AYZZ,AZZZ,AXXXX,AXXXY,
     $  AXXXZ,AXXYY,AXXZZ,AXYYY,AXZZZ,AXXYZ,AXYYZ,AXYZZ,AYYYY,
     $  AYYYZ,AYYZZ,AYZZZ,AZZZZ,F0,FX,FY,FZ,FXX,FXY,FXZ,FYY,FYZ,
     $  FZZ,FXXX,FXXY,FXXZ,FXYY,FXZZ,FXYZ,FYYY,FYYZ,FYZZ,FZZZ,
     $  FXXXX,FXXXY,FXXXZ,FXXYY,FXXZZ,FXYYY,FXZZZ,FXXYZ,FXYYZ,
     $  FXYZZ,FYYYY,FYYYZ,FYYZZ,FYZZZ,FZZZZ,G0,GX,GY,GZ,GXX,GXY,
     $  GXZ,GYY,GYZ,GZZ,GXXX,GXXY,GXXZ,GXYY,GXZZ,GXYZ,GYYY,GYYZ,
     $  GYZZ,GZZZ,GXXXX,GXXXY,GXXXZ,GXXYY,GXXZZ,GXYYY,GXZZZ,
     $  GXXYZ,GXYYZ,GXYZZ,GYYYY,GYYYZ,GYYZZ,GYZZZ,GZZZZ
      Common /What/ IWhat
      Save Two,Three,Four,Six
      Data Two/2.0d0/,Three/3.0d0/,Four/4.0d0/,Six/6.0d0/
C
      G0=A0*F0
      GX=AX*F0+A0*FX
      GY=AY*F0+A0*FY
      GZ=AZ*F0+A0*FZ
      If(IWhat.gt.0)Then
      GXX=AXX*F0+TWO*AX*FX+A0*FXX
      GXY=AXY*F0+AX*FY+AY*FX+A0*FXY
      GXZ=AXZ*F0+AX*FZ+AZ*FX+A0*FXZ
      GYY=AYY*F0+TWO*AY*FY+A0*FYY
      GYZ=AYZ*F0+AY*FZ+AZ*FY+A0*FYZ
      GZZ=AZZ*F0+TWO*AZ*FZ+A0*FZZ
      Endif
      If(IWhat.gt.1)Then
      GXXX=AXXX*F0+THREE*AXX*FX+THREE*AX*FXX+A0*FXXX
      GXXY=AXXY*F0+AXX*FY+TWO*AXY*FX+TWO*AX*FXY+AY*FXX+A0*FXXY
      GXXZ=AXXZ*F0+AXX*FZ+TWO*AXZ*FX+TWO*AX*FXZ+AZ*FXX+A0*FXXZ
      GYYY=AYYY*F0+THREE*AYY*FY+THREE*AY*FYY+A0*FYYY
      GXYY=AXYY*F0+AYY*FX+TWO*AXY*FY+TWO*AY*FXY+AX*FYY+A0*FXYY
      GYYZ=AYYZ*F0+AYY*FZ+TWO*AYZ*FY+TWO*AY*FYZ+AZ*FYY+A0*FYYZ
      GZZZ=AZZZ*F0+THREE*AZZ*FZ+THREE*AZ*FZZ+A0*FZZZ
      GXZZ=AXZZ*F0+AZZ*FX+TWO*AXZ*FZ+TWO*AZ*FXZ+AX*FZZ+A0*FXZZ
      GYZZ=AYZZ*F0+AZZ*FY+TWO*AYZ*FZ+TWO*AZ*FYZ+AY*FZZ+A0*FYZZ
      GXYZ=AXYZ*F0+AXY*FZ+AYZ*FX+AY*FXZ+AXZ*FY+AX*FYZ+AZ*FXY+
     +     A0*FXYZ
      Endif
      If(IWhat.gt.2)Then
      GXXXX=AXXXX*F0+FOUR*AXXX*FX+SIX*AXX*FXX+FOUR*AX*FXXX+A0*FXXXX
      GXXXY=AXXXY*F0+AXXX*FY+THREE*AXXY*FX+THREE*AXX*FXY+THREE*AXY*FXX+
     $      THREE*AX*FXXY+AY*FXXX+A0*FXXXY
      GXXXZ=AXXXZ*F0+AXXX*FZ+THREE*AXXZ*FX+THREE*AXX*FXZ+THREE*AXZ*FXX+
     $      THREE*AX*FXXZ+AZ*FXXX+A0*FXXXZ
      GXXYY=AXXYY*F0+TWO*AXXY*FY+AXX*FYY+TWO*AXYY*FX+FOUR*AXY*FXY+
     $      TWO*AX*FXYY+AYY*FXX+TWO*AY*FXXY+A0*FXXYY
      GXXZZ=AXXZZ*F0+TWO*AXXZ*FZ+AXX*FZZ+TWO*AXZZ*FX+FOUR*AXZ*FXZ+
     $      TWO*AX*FXZZ+AZZ*FXX+TWO*AZ*FXXZ+A0*FXXZZ
      GXYYY=AXYYY*F0+AYYY*FX+THREE*AXYY*FY+THREE*AYY*FXY+THREE*AXY*FYY+
     $      THREE*AY*FXYY+AX*FYYY+A0*FXYYY
      GXZZZ=AXZZZ*F0+AZZZ*FX+THREE*AXZZ*FZ+THREE*AZZ*FXZ+THREE*AXZ*FZZ+
     $      THREE*AZ*FXZZ+AX*FZZZ+A0*FXZZZ
      GXXYZ=AXXYZ*F0+AXXY*FZ+AXXZ*FY+AXX*FYZ+TWO*AXYZ*FX+TWO*AXY*FXZ+
     $      TWO*AXZ*FXY+TWO*AX*FXYZ+AYZ*FXX+AY*FXXZ+AZ*FXXY+A0*FXXYZ
      GXYYZ=AXYYZ*F0+AYYZ*FX+AXYY*FZ+AYY*FXZ+TWO*AXYZ*FY+TWO*AYZ*FXY+
     $      TWO*AXY*FYZ+TWO*AY*FXYZ+AXZ*FYY+AZ*FXYY+AX*FYYZ+A0*FXYYZ
      GXYZZ=AXYZZ*F0+AYZZ*FX+AXZZ*FY+AZZ*FXY+TWO*AXYZ*FZ+TWO*AYZ*FXZ+
     $      TWO*AXZ*FYZ+TWO*AZ*FXYZ+AXY*FZZ+AY*FXZZ+AX*FYZZ+A0*FXYZZ
      GYYYY=AYYYY*F0+FOUR*AYYY*FY+SIX*AYY*FYY+FOUR*AY*FYYY+A0*FYYYY
      GYYYZ=AYYYZ*F0+AYYY*FZ+THREE*AYYZ*FY+THREE*AYY*FYZ+THREE*AYZ*FYY+
     $      THREE*AY*FYYZ+AZ*FYYY+A0*FYYYZ
      GYYZZ=AYYZZ*F0+TWO*AYYZ*FZ+AYY*FZZ+TWO*AYZZ*FY+FOUR*AYZ*FYZ+
     $      TWO*AY*FYZZ+AZZ*FYY+TWO*AZ*FYYZ+A0*FYYZZ
      GYZZZ=AYZZZ*F0+AZZZ*FY+THREE*AYZZ*FZ+THREE*AZZ*FYZ+THREE*AYZ*FZZ+
     $      THREE*AZ*FYZZ+AY*FZZZ+A0*FYZZZ
      GZZZZ=AZZZZ*F0+FOUR*AZZZ*FZ+SIX*AZZ*FZZ+FOUR*AZ*FZZZ+A0*FZZZZ
      Endif
C
      Return
      End 
      Subroutine Props(XYZ,IFUNC,ICrit,Iwhole,J1,J2,INP,Evsave)
C
      Implicit Double Precision (A-H,O-Z)
C
      Character*8 Atnam
      Character*80 WFNTTL,JOBTTL
      Parameter(MaxAtm=100,MaxOff=200000,MaxCrt=500)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /ANG/ ANGLE(3,MaxAtm,MaxAtm)
      COMMON /DIST/ BANGLE(6),RMIN,RMAX,MINR
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICofMx
      COMMON /STRING/ WFNTTL,JOBTTL,ATNAM(MaxAtm)
      COMMON /UNITS/  INPT,IOUT,IWFN,IWLP
      COMMON /WHAT/ Iwhat
      Common /Options/ Icut,Iprint,Eps,Epsnuc,Dmp,DmpNuc
      Dimension XYZ(3),HRHO(3,3),HD2RHO(3,3),HGEE(3,3),HQUAY(3,3),
     $  HVNE(3,3),WRHO(3),WD2RHO(3),WGEE(3),WQUAY(3),WVNE(3),
     $  SG(3,3),EV(3),EU(3),DSIG(3),Evsave(maxCrt,3),SV(3),PN(4),
     $  WORK(3,3),HVEE(3,3),WVEE(3)
      Save Zero,dxyz,One,Four
      Data Zero/0.d0/,dxyz/1.d-3/,One/1.d0/,Four/4.0d0/
1000  FORMAT(/,' COORDINATES OF CRITICAL POINT AND',
     $' DISTANCE FROM MOLECULAR ORIGIN')
1010  FORMAT(/,' COORDINATES OF POINT AND DISTANCE',
     $' FROM MOLECULAR ORIGIN')
1020  FORMAT(10X,'X = ',1PE16.8,/,10X,'Y = ',1PE16.8,/,
     $       10X,'Z = ',1PE16.8,/,10X,'R = ',1PE16.8)
1030  FORMAT(/,' VECTORS AND DISTANCES FROM NUCLEI TO CRITICAL POINT')
1040  FORMAT(/,' VECTORS AND DISTANCES FROM NUCLEI TO POINT')
1050  FORMAT(/,' NUCLEUS',8X,'   X   ',8X,'   Y   ',8X,'   Z   ',
     $8X,'  Dist  ')
1060  FORMAT(1X,A8,1X,1P4E16.8)
1070  FORMAT(/,' ANGLES OF VECTORS WITH RESPECT TO YZ/XZ/XY PLANES',
     $' OF MCS')
1080  FORMAT(/,' NUCLEUS',8X,'YZ ANGLE',8X,'XZ ANGLE',8X,'XY ANGLE')
1090  FORMAT(1X,A8,1X,1P3E16.8)
1100  FORMAT(/,' EIGENVALUES OF THE HESSIAN ')
1110  FORMAT(1X,1P3E18.8)
1120  FORMAT(/,' THE ELLIPTICITY IS ',1PE16.8)
1130  FORMAT(/,' EIGENVECTORS OF THE HESSIAN ')
1140  FORMAT(/,' EIGENVALUES OF THE STRESS TENSOR ')
1150  FORMAT(/,' THE TRACE OF THE STRESS TENSOR IS ',1PE16.8)
1160  FORMAT(/,' EIGENVECTORS OF THE STRESS TENSOR ')
1170  FORMAT(/,' COMPONENTS OF THE DIVERGENCE OF THE STRESS TENSOR ')
1180  FORMAT(/,' MAGNITUDE OF THE DIVERGENCE OF THE STRESS TENSOR ',
     $1PE16.8)
1190  FORMAT(/,' VALUES ',/,' Rho(r)',18X,1PE17.10,/,
     $' |GRAD(Rho(r))|',10X,1PE17.10,/,' GRAD(Rho(r))x',
     $11X,1PE17.10,/,
     $' GRAD(Rho(r))y',11X,1PE17.10,/,' GRAD(Rho(r))z',11X,1PE17.10,/,
     $' DEL**2(Rho(r))',10X,1PE17.10,/,' G(r)',20X,1PE17.10,/,
     $' K(r)',20X,1PE17.10,/,' L(r)',20X,1PE17.10,/,
     $' Vnuc(r)',17X,1PE17.10,/,' V(r)',20X,1PE17.10)
1200  FORMAT(/,' DO YOU WANT TO TRACE THE BOND PATH ? (0=no/1=yes) ',$)
1210  FORMAT(/,' BOND PATH LINKED TO ',A8,' (',I3,' CALLS)',/' EPS=',
     1 1PE16.8,'   LENGTH=',1PE16.8)
1220  FORMAT(/,' CHARGE DENSITY MAXIMUM OCCURS AT ',/,
     1 12X,' X = ',1PE16.8,/,12X,' Y = ',1PE16.8,/,12X,' Z = ',1PE16.8,
     2 //,'  DISPLACEMENT = ',1PE16.8)
1230  FORMAT(/,' UNIT VECTOR FROM NUCLEUS IN BOND PATH DIRECTION ',
     1 /,'    CARTESIAN COORD.    SPHERICAL COORD.')
1240  FORMAT(1H ,3(2X,F16.8,4X,F16.8,/,1X))
1250  FORMAT(/,' TOTAL BOND PATH LENGTH = ',1PE16.8)
1260  FORMAT(' GEOMETRIC BOND LENGTH = ',1PE16.8)
1270  FORMAT(' BOND PATH LENGTH MINUS GEOMETRICAL BOND LENGTH = ',
     $  1PE16.8)
1280  FORMAT(/,' VALUES ',/,' Del**2*(Rho(r))',14X,1PE17.10,/,
     $' |GRAD(Del**2(Rho(r)))|',7X,1PE17.10,/,
     $' GRAD(Del**2(Rho(r)))x',8X,1PE17.10,/,
     $' GRAD(Del**2(Rho(r)))y',8X,1PE17.10,/,
     $' GRAD(Del**2(Rho(r)))z',8X,1PE17.10,/,
     $' DEL**2(DEL**2(Rho(r)))',7X,1PE17.10,/,
     $' Rho(r)',23X,1PE17.10,/,' G(r)',25X,1PE17.10,/,
     $' K(r)',25X,1PE17.10,/,' L(r)',25X,1PE17.10,/,
     $' Vnuc(r)',22X,1PE17.10,/,' V(r)',25X,1PE17.10)
1290  FORMAT(/,' VALUES ',/,' G(r)',20X,1PE17.10,/,
     $' |GRAD(G(r)|',13X,1PE17.10,/,' GRAD(G(r))x',
     $13X,1PE17.10,/,' GRAD(G(r))y',13X,1PE17.10,/,
     $' GRAD(G(r))z',13X,1PE17.10,/,' DEL**2(G(r))',12X,1PE17.10,/,
     $' Rho(r)',18X,1PE17.10,/,' DEL**2(Rho(r))',10X,1PE17.10,/,
     $' K(r)',20X,1PE17.10,/,' L(r)',20X,1PE17.10,/,' Vnuc(r)',
     $17X,1PE17.10,/,' V(r)',20X,1PE17.10)
1300  FORMAT(/,' VALUES ',/,' K(r)',20X,1PE17.10,/,
     $' |GRAD(K(r))|',12X,1PE17.10,/,' GRAD(K(r))x',
     $13X,1PE17.10,/,' GRAD(K(r))y',13X,1PE17.10,/,
     $' GRAD(K(r))z',13X,1PE17.10,/,' Del**2(G(r))',12X,1PE17.10,/,
     $' Rho(r)',18X,1PE17.10,/,' DEL**2(Rho(r))',10X,1PE17.10,/,
     $' G(r)',20X,1PE17.10,/,' L(r)',20X,1PE17.10,/,' Vnuc(r)',
     $17X,1PE17.10,/,' V(r)',20X,1PE17.10)
1310  FORMAT(/,' VALUES ',/,' Vnuc(r)',17X,1PE17.10,/,
     $' |GRAD(Vnuc(r))|',9X,1PE17.10,/,' GRAD(Vnuc(r))x',
     $10X,1PE17.10,/,' GRAD(Vnuc(r))y',10X,1PE17.10,/,
     $' GRAD(Vnuc(r))z',10X,1PE17.10,/,' Del**2(Vnuc(r))',
     $9X,1PE17.10,/,' Rho(r)',18X,1PE17.10,/,' DEL**2(Rho(r))',
     $10X,1PE17.10,/,' G(r)',20X,1PE17.10,/,' K(r)',20X,
     $1PE17.10,/,' L(r)',20X,1PE17.10,/,' V(r)',20X,1PE17.10)
1320  FORMAT(/,' VALUES ',/,' V(r)',20X,1PE17.10,/,
     $' |GRAD(V(r))|',12X,1PE17.10,/,' GRAD(V(r))x',
     $13X,1PE17.10,/,' GRAD(V(r))y',13X,1PE17.10,/,
     $' GRAD(V(r))z',13X,1PE17.10,/,' Del**2(V(r))',
     $12X,1PE17.10,/,' Rho(r)',18X,1PE17.10,/,' DEL**2(Rho(r))',
     $10X,1PE17.10,/,' G(r)',20X,1PE17.10,/,' K(r)',20X,
     $1PE17.10,/,' L(r)',20X,1PE17.10,/,' Vnuc(r)',17X,1PE17.10)
C
      R = DSQRT(XYZ(1)**2+XYZ(2)**2+XYZ(3)**2)
      If(Icrit.eq.1)Then
      If(Iprint.eq.1)WRITE (IOUT,1000)
      WRITE (IWLP,1000)
      Else
      If(Iprint.eq.1)WRITE (IOUT,1010)
      WRITE (IWLP,1010)
      Endif
      If(Iprint.eq.1)WRITE (IOUT,1020) XYZ(1),XYZ(2),XYZ(3),R
      WRITE (IWLP,1020) XYZ(1),XYZ(2),XYZ(3),R
C
      If(Icrit.eq.1)Then
      If(Iprint.eq.1)WRITE (IOUT,1030)
      WRITE (IWLP,1030)
      Else
      If(Iprint.eq.1)WRITE (IOUT,1040)
      WRITE (IWLP,1040)
      Endif
C
      If(Iprint.eq.1)WRITE (IOUT,1050)
      WRITE (IWLP,1050)
      DO 10 I = 1,NCENT
        XREL=XYZ(1)-CO(IXC+I)
        YREL=XYZ(2)-CO(IYC+I)
        ZREL=XYZ(3)-CO(IZC+I)
        Dist=dsqrt(xrel**2+yrel**2+zrel**2)
        If(Iprint.eq.1)WRITE(IOUT,1060)ATnam(i),Xrel,Yrel,Zrel,Dist
        WRITE(IWLP,1060)ATnam(i),Xrel,Yrel,Zrel,Dist
10    Continue
C
      If(Iprint.eq.1)WRITE (IOUT,1070)
      WRITE (IWLP,1070)
      If(Iprint.eq.1)WRITE (IOUT,1080)
      WRITE (IWLP,1080)
      DO 20 I = 1,NCENT
        CALL GEOM (I,XYZ,RN,AYZ,AXZ,AXY)
        If(Iprint.eq.1)WRITE(IOUT,1090) ATNAM(I),AYZ,AXZ,AXY
        WRITE(IWLP,1090) ATNAM(I),AYZ,AXZ,AXY
20    CONTINUE
C
      If(IFUNC.eq.1)Then
      IWhat=2
      CALL GRDRHO(1,XYZ,RHO,WRHO,GRHO,HRHO,SG)
      CALL GRDD2R(1,XYZ,D2Rho,WD2Rho,GD2Rho,HD2Rho)
      CALL GRDKEG(1,XYZ,GEE,WGEE,GGEE,HGEE)
      CALL GRDKEK(1,XYZ,QUAY,WQUAY,GQUAY,HQUAY)
      CALL GRDVNE(1,XYZ,VNE,WVNE,GVNE,HVNE)
      CALL GRDV(1,XYZ,VEE,WVEE,GVEE,HVEE)
      CALL TRACE(HRHO,EV,WORK,3,IFAIL)
      If(ICRIT.EQ.1)Then
      EVSAVE(INP,1)=EV(1)
      EVSAVE(INP,2)=EV(2)
      EVSAVE(INP,3)=EV(3)
      Endif
      CALL TRACE(SG,EU,WORK,3,IFAIL)
      SV(1) = HRHO(1,3)
      SV(2) = HRHO(2,3)
      SV(3) = HRHO(3,3)
      XLag = -D2Rho/Four
      If(Iprint.eq.1)WRITE (IOUT,1100)
      WRITE (IWLP,1100)
      If(Iprint.eq.1)WRITE (IOUT,1110) (EV(I),I=1,3)
      WRITE (IWLP,1110) (EV(I),I=1,3)
      IF ((EV(1).LT.Zero).AND.(EV(2).LT.Zero)) THEN
        ELLIPT = (EV(1)/EV(2))-One
        If(Iprint.eq.1)WRITE (IOUT,1120) ELLIPT
        WRITE (IWLP,1120) ELLIPT
      ENDIF
      If(Iprint.eq.1)WRITE (IOUT,1130)
      WRITE (IWLP,1130)
      DO 150 I = 1,3
        If(Iprint.eq.1)WRITE(IOUT,1110)(HRHO(I,J),J=1,3)
        WRITE (IWLP,1110) (HRHO(I,J),J=1,3)
150   CONTINUE
C
      If(Iprint.eq.1)WRITE (IOUT,1140)
      WRITE (IWLP,1140)
      If(Iprint.eq.1)WRITE (IOUT,1110) (EU(I),I=1,3)
      WRITE (IWLP,1110) (EU(I),I=1,3)
      TRSG = EU(1)+EU(2)+EU(3)
      If(Iprint.eq.1)WRITE (IOUT,1150) TRSG
      WRITE (IWLP,1150) TRSG
      If(Iprint.eq.1)WRITE (IOUT,1160)
      WRITE (IWLP,1160)
      DO 160 I = 1,3
        If(Iprint.eq.1)WRITE (IOUT,1110) (SG(I,J),J=1,3)
        WRITE (IWLP,1110) (SG(I,J),J=1,3)
160   CONTINUE
      CALL DIVSTR(DSIG,DSIGM)
      If(Iprint.eq.1)WRITE (IOUT,1170)
      WRITE (IWLP,1170)
      If(Iprint.eq.1)WRITE (IOUT,1110) (DSIG(I),I=1,3)
      WRITE (IWLP,1110) (DSIG(I),I=1,3)
      If(Iprint.eq.1)WRITE (IOUT,1180) DSIGM
      WRITE (IWLP,1180) DSIGM
      If(Iprint.eq.1)Then
      WRITE (IOUT,1190) RHO,GRHO,WRHO(1),WRHO(2),WRHO(3),D2RHO,
     $  GEE,QUAY,XLag,VNE,VEE
      Endif
      WRITE (IWLP,1190) RHO,GRHO,WRHO(1),WRHO(2),WRHO(3),D2RHO,
     $  GEE,QUAY,XLag,VNE,VEE
C
      IF ((EV(1).LT.Zero) .AND. (EV(2).LT.Zero) 
     +.AND. (EV(3).GT.Zero).and.iwhole.eq.0.and.icrit.eq.1)THEN
        WRITE (IOUT,1200)
        READ (INPT,*) IWK
        IF (IWK .EQ. 1) THEN
          X1 = XYZ(1) + DXYZ*SV(1)
          Y1 = XYZ(2) + DXYZ*SV(2)
          Z1 = XYZ(3) + DXYZ*SV(3)
          CALL BOND(X1,Y1,Z1,BL1,NCAL,PN,IFunc,IWhole)
          If(Iprint.eq.1)WRITE (IOUT,1210) ATNAM(MINR),NCAL,EPSD,BL1
          WRITE (IWLP,1210) ATNAM(MINR),NCAL,EPSD,BL1
          If(Iprint.eq.1)WRITE (IOUT,1220) (PN(I),I=1,4)
          WRITE (IWLP,1220) (PN(I),I=1,4)
          If(Iprint.eq.1)WRITE (IOUT,1230)
          WRITE (IWLP,1230)
          If(Iprint.eq.1)WRITE (IOUT,1240) (BANGLE(I),BANGLE(I+3),I=1,3)
          WRITE (IWLP,1240) (BANGLE(I),BANGLE(I+3),I=1,3)
          IF (J1 .EQ. MINR) THEN
            JA = J1
            JB = J2
          ELSE
            JA = J2
            JB = J1
          END IF
          ANGLE(1,JA,JB) = BANGLE(1)
          ANGLE(2,JA,JB) = BANGLE(2)
          ANGLE(3,JA,JB) = BANGLE(3)
          X1 = XYZ(1) - DXYZ*SV(1)
          Y1 = XYZ(2) - DXYZ*SV(2)
          Z1 = XYZ(3) - DXYZ*SV(3)
          CALL BOND(X1,Y1,Z1,BL2,NCAL,PN)
          If(Iprint.eq.1)WRITE (IOUT,1210) ATNAM(MINR),NCAL,EPSD,BL2
          WRITE (IWLP,1210) ATNAM(MINR),NCAL,EPSD,BL2
          If(Iprint.eq.1)WRITE (IOUT,1220) (PN(I),I=1,4)
          WRITE (IWLP,1220) (PN(I),I=1,4)
          If(Iprint.eq.1)WRITE (IOUT,1230)
          WRITE (IWLP,1230)
          If(Iprint.eq.1)WRITE (IOUT,1240) (BANGLE(I),BANGLE(I+3),I=1,3)
          WRITE (IWLP,1240) (BANGLE(I),BANGLE(I+1),I=1,3)
          IF (J1 .EQ. MINR) THEN
            JA = J1
            JB = J2
          ELSE
            JA = J2
            JB = J1
          END IF
          ANGLE(1,JA,JB) = BANGLE(1)
          ANGLE(2,JA,JB) = BANGLE(2)
          ANGLE(3,JA,JB) = BANGLE(3)
          BL = BL1+BL2
          GBL = DSQRT((CO(IXC+J1)-CO(IXC+J2))**2 +
     +                (CO(IYC+J1)-CO(IYC+J2))**2 +
     +                (CO(IZC+J1)-CO(IZC+J2))**2)
          DBL = BL - GBL
          If(Iprint.eq.1)WRITE (IOUT,1250) BL
          WRITE (IWLP,1250) BL
          If(Iprint.eq.1)WRITE (IOUT,1260) GBL
          WRITE (IWLP,1260) GBL
          If(Iprint.eq.1)WRITE (IOUT,1270) DBL
          WRITE (IWLP,1270) DBL
        END IF
      END IF
C
      ElseIf(IFunc.eq.2)Then
      IWhat=3
      Call GrdD2r(0,XYZ,D2Rho,WD2Rho,GD2Rho,HD2Rho)
      CALL GRDRHO(2,XYZ,RHO,WRHO,GRHO,HRHO,SG)
      CALL GRDKEG(1,XYZ,GEE,WGEE,GGEE,HGEE)
      CALL GRDKEK(1,XYZ,QUAY,WQUAY,GQUAY,HQUAY)
      CALL GRDVNE(1,XYZ,VNE,WVNE,GVNE,HVNE)
      CALL GRDV(1,XYZ,VEE,WVEE,GVEE,HVEE)
      XLag=-D2RHo/Four
      HHD2Rh=HD2Rho(1,1)+HD2RHO(2,2)+HD2RHO(3,3)
      CALL TRACE(HD2RHO,EV,WORK,3,IFAIL)
      If(ICRIT.EQ.1)Then
      EVSAVE(INP,1)=EV(1)
      EVSAVE(INP,2)=EV(2)
      EVSAVE(INP,3)=EV(3)
      Endif
      If(Iprint.eq.1)WRITE (IOUT,1100)
      WRITE (IWLP,1100)
      If(Iprint.eq.1)WRITE (IOUT,1110) (EV(I),I=1,3)
      WRITE (IWLP,1110) (EV(I),I=1,3)
      If(Iprint.eq.1)WRITE (IOUT,1130)
      WRITE (IWLP,1130)
      DO 151 I = 1,3
        If(Iprint.eq.1)WRITE (IOUT,1110) (HD2RHO(I,J),J=1,3)
        WRITE (IWLP,1110) (HD2RHO(I,J),J=1,3)
151   CONTINUE
C
      If(Iprint.eq.1)Then
      WRITE (IOUT,1280) D2RHO,GD2RHO,WD2RHO(1),WD2RHO(2),WD2RHO(3),
     $  HHD2Rh,RHO,GEE,QUAY,XLag,VNE,VEE
      Endif
      WRITE (IWLP,1280) D2RHO,GD2RHO,WD2RHO(1),WD2RHO(2),WD2RHO(3),
     $  HHD2Rh,RHO,GEE,QUAY,XLag,VNE,VEE
C
      ElseIf(IFunc.eq.3)Then
      IWhat=2
      CALL GRDKEG(0,XYZ,GEE,WGEE,GGEE,HGEE)
      Call GrdD2r(1,XYZ,D2Rho,WD2Rho,GD2Rho,HD2Rho)
      CALL GRDRHO(2,XYZ,RHO,WRHO,GRHO,HRHO,SG)
      CALL GRDKEK(1,XYZ,QUAY,WQUAY,GQUAY,HQUAY)
      CALL GRDVNE(1,XYZ,VNE,WVNE,GVNE,HVNE)
      CALL GRDV(1,XYZ,VEE,WVEE,GVEE,HVEE)
      XLag=-D2RHo/Four
      HHGEE=HGEE(1,1)+HGEE(2,2)+HGEE(3,3)
      CALL TRACE(HGEE,EV,WORK,3,IFAIL)
      If(ICRIT.EQ.1)Then
      EVSAVE(INP,1)=EV(1)
      EVSAVE(INP,2)=EV(2)
      EVSAVE(INP,3)=EV(3)
      Endif
      If(Iprint.eq.1)WRITE (IOUT,1100)
      WRITE (IWLP,1100)
      If(Iprint.eq.1)WRITE (IOUT,1110) (EV(I),I=1,3)
      WRITE (IWLP,1110) (EV(I),I=1,3)
      If(Iprint.eq.1)WRITE (IOUT,1130)
      WRITE (IWLP,1130)
      DO 152 I = 1,3
        If(Iprint.eq.1)WRITE (IOUT,1110) (HGEE(I,J),J=1,3)
        WRITE (IWLP,1110) (HGEE(I,J),J=1,3)
152   CONTINUE
C
      If(Iprint.eq.1)Then
      WRITE (IOUT,1290) GEE,GGEE,WGEE(1),WGEE(2),WGEE(3),
     $  HHGEE,RHO,D2RHO,QUAY,XLag,VNE
      Endif
      WRITE (IWLP,1290) GEE,GGEE,WGEE(1),WGEE(2),WGEE(3),
     $  HHGEE,RHO,D2RHO,QUAY,XLag,VNE
C
      ElseIf(IFunc.eq.4)Then
      IWhat=3
      CALL GRDKEK(0,XYZ,QUAY,WQUAY,GQUAY,HQUAY)
      CALL GRDKEG(1,XYZ,GEE,WGEE,GGEE,HGEE)
      Call GrdD2r(1,XYZ,D2Rho,WD2Rho,GD2Rho,HD2Rho)
      CALL GRDRHO(2,XYZ,RHO,WRHO,GRHO,HRHO,SG)
      CALL GRDVNE(1,XYZ,VNE,WVNE,GVNE,HVNE)
      CALL GRDV(1,XYZ,VEE,WVEE,GVEE,HVEE)
      XLag=-D2RHo/Four
      HHQUAY=HQUAY(1,1)+HQUAY(2,2)+HQUAY(3,3)
      CALL TRACE(HQUAY,EV,WORK,3,IFAIL)
      If(ICRIT.EQ.1)Then
      EVSAVE(INP,1)=EV(1)
      EVSAVE(INP,2)=EV(2)
      EVSAVE(INP,3)=EV(3)
      Endif
      If(Iprint.eq.1)WRITE (IOUT,1100)
      WRITE (IWLP,1100)
      If(Iprint.eq.1)WRITE (IOUT,1110) (EV(I),I=1,3)
      WRITE (IWLP,1110) (EV(I),I=1,3)
      If(Iprint.eq.1)WRITE (IOUT,1130)
      WRITE (IWLP,1130)
      DO 153 I = 1,3
        If(Iprint.eq.1)WRITE (IOUT,1110) (HQUAY(I,J),J=1,3)
        WRITE (IWLP,1110) (HQUAY(I,J),J=1,3)
153   CONTINUE
C
      If(Iprint.eq.1)Then
      WRITE (IOUT,1300) QUAY,GQUAY,WQUAY(1),WQUAY(2),WQUAY(3),
     $  HHQUAY,RHO,D2RHO,GEE,XLag,VNE,VEE
      Endif
      WRITE (IWLP,1300) QUAY,GQUAY,WQUAY(1),WQUAY(2),WQUAY(3),
     $  HHQUAY,RHO,D2RHO,GEE,XLag,VNE,VEE
C
      ElseIf(IFunc.eq.5)Then
C
      IWhat=1
      CALL GRDRHO(1,XYZ,RHO,WRHO,GRHO,HRHO,SG)
      CALL GRDVNE(0,XYZ,VNE,WVNE,GVNE,HVNE)
      Call GrdD2r(1,XYZ,D2Rho,WD2Rho,GD2Rho,HD2Rho)
      CALL GRDKEK(1,XYZ,QUAY,WQUAY,GQUAY,HQUAY)
      CALL GRDKEG(1,XYZ,GEE,WGEE,GGEE,HGEE)
      CALL GRDV(1,XYZ,VEE,WVEE,GVEE,HVEE)
      XLag=-D2RHo/Four
      HHVNE=HVNE(1,1)+HVNE(2,2)+HVNE(3,3)
      CALL TRACE(HVNE,EV,WORK,3,IFAIL)
      If(ICRIT.EQ.1)Then
      EVSAVE(INP,1)=EV(1)
      EVSAVE(INP,2)=EV(2)
      EVSAVE(INP,3)=EV(3)
      Endif
      IF(Iprint.eq.1)WRITE (IOUT,1100)
      WRITE (IWLP,1100)
      IF(Iprint.eq.1)WRITE (IOUT,1110) (EV(I),I=1,3)
      WRITE (IWLP,1110) (EV(I),I=1,3)
      IF(Iprint.eq.1)WRITE (IOUT,1130)
      WRITE (IWLP,1130)
      DO 154 I = 1,3
        IF(Iprint.eq.1)WRITE (IOUT,1110) (HVNE(I,J),J=1,3)
        WRITE (IWLP,1110) (HVNE(I,J),J=1,3)
154   CONTINUE
C
      If(Iprint.eq.1)Then
      WRITE (IOUT,1310) VNE,GVNE,WVNE(1),WVNE(2),WVNE(3),
     $  HHVNE,RHO,D2RHO,GEE,QUAY,XLag,VEE
      Endif
      WRITE (IWLP,1310) VNE,GVNE,WVNE(1),WVNE(2),WVNE(3),
     $  HHVNE,RHO,D2RHO,GEE,QUAY,XLag,VEE
C
      ElseIf(IFunc.eq.6)Then
C
      IWhat=3
      CALL GRDV(0,XYZ,VEE,WVEE,GVEE,HVEE)
      CALL GRDRHO(2,XYZ,RHO,WRHO,GRHO,HRHO,SG)
      Call GrdD2r(1,XYZ,D2Rho,WD2Rho,GD2Rho,HD2Rho)
      CALL GRDKEK(1,XYZ,QUAY,WQUAY,GQUAY,HQUAY)
      CALL GRDKEG(1,XYZ,GEE,WGEE,GGEE,HGEE)
      CALL GRDVNE(1,XYZ,VNE,WVNE,GVNE,HVNE)
      XLag=-D2RHo/Four
      HHVEE=HVEE(1,1)+HVEE(2,2)+HVEE(3,3)
      CALL TRACE(HVEE,EV,WORK,3,IFAIL)
      If(ICRIT.EQ.1)Then
      EVSAVE(INP,1)=EV(1)
      EVSAVE(INP,2)=EV(2)
      EVSAVE(INP,3)=EV(3)
      Endif
      If(Iprint.eq.1)WRITE (IOUT,1100)
      WRITE (IWLP,1100)
      If(Iprint.eq.1)WRITE (IOUT,1110) (EV(I),I=1,3)
      WRITE (IWLP,1110) (EV(I),I=1,3)
      If(Iprint.eq.1)WRITE (IOUT,1130)
      WRITE (IWLP,1130)
      DO 155 I = 1,3
        If(Iprint.eq.1)WRITE (IOUT,1110) (HVEE(I,J),J=1,3)
        WRITE (IWLP,1110) (HVEE(I,J),J=1,3)
155   CONTINUE
C
      If(Iprint.eq.1)Then
      WRITE (IOUT,1320) VEE,GVEE,WVEE(1),WVEE(2),WVEE(3),
     $  HHVEE,RHO,D2RHO,GEE,QUAY,XLag,VNE
      Endif
      WRITE (IWLP,1320) VEE,GVEE,WVEE(1),WVEE(2),WVEE(3),
     $  HHVEE,RHO,D2RHO,GEE,QUAY,XLag,VNE
C
      Endif
C
      Return
      End
      SUBROUTINE RDWFN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*80 WFNTTL,JOBTTL
      CHARACTER*8 ATNAM
      CHARACTER*4 MODE
      PARAMETER (NTYPE=20,MaxOff=200000,MaxAtm=100)
      COMMON CO(MaxOff),IC(MaxOff),NCENT,NMO,NPRIMS
      COMMON /OFFSET/ ITYPE,ICENT,IMO,IEORB,IE,ICHARG,IXC,IYC,IZC,
     $  IXX, IYY, IZZ,IRR,IR2,IP,IPSI,IGX,IGY,IGZ,IGXX,IGXY,IGXZ,IGYY,
     $  IGYZ,IGZZ,IGXXX,IGXXY,IGXXZ,IGXYY,IGXZZ,IGXYZ,IGYYY,IGYYZ,
     $  IGYZZ,IGZZZ,IGXXXX,IGXXXY,IGXXXZ,IGXXYY,IGXXZZ,IGXYYY,IGXZZZ,
     $  IGXXYZ,IGXYYZ,IGXYZZ,IGYYYY,IGYYYZ,IGYYZZ,IGYZZZ,IGZZZZ,ICOFMX
      COMMON /STRING/ WFNTTL,JOBTTL,ATNAM(MaxAtm)
      COMMON /UNITS/  INPT,IOUT,IWFN,IWLP
      DATA ENDATA /8HEND DATA/,ZERO/0.d0/
C
      READ(IWFN,101) WFNTTL
C
      READ(IWFN,102) MODE,NMO,NPRIMS,NCENT
C
C     Set OffSets For Storage in IC:
C
C                             FUNCTION TYPES
      ITYPE = 0
C                             FUNCTION CENTRES
      ICENT = ITYPE+NPRIMS
C
C     Set OffSets For Storage in CO - MO Coefficients, Then others: 
C
      IMO=0
C                             ORBITAL ENERGIES OF EACH M.O.
      IEORB = IMO+NMO*NPRIMS
C                             EXPONENTS OF THE FUNCTIONS
      IE = IEORB+NMO
C                             NUCLEAR CHARGE OF EACH CENTRE
      ICHARG = IE+NPRIMS
C
C                             XYZ COORDS OF THE CENTRES IN
C                             ORIGINAL CARTESIAN SYSTEM
      IXC = ICHARG+NCENT
      IYC = IXC+NCENT
      IZC = IYC+NCENT
C
C                             XYZ COORDS OF CENTRES RELATIVE TO A
C                             CURRENT TEST POINT IN THE INTEGRATION
      IXX = IZC+NCENT
      IYY = IXX+NCENT
      IZZ = IYY+NCENT
C
C                             SQUARE OF DISTANCE FROM CENTRES TO
C                             CURRENT POINT
      IRR = IZZ+NCENT
C                             DISTANCE FROM CENTRES TO CURRENT POINT
      IR2 = IRR+NCENT
C                             OCCUPATION NUMBER OF EACH M.O
      IP = IR2+NCENT
C                             PSI VALUES FOR EACH M.O.
      IPSI = IP+NMO
C
      IGX = IPSI+NMO
      IGY = IGX+NMO
      IGZ = IGY+NMO
C
      IGXX = IGZ+NMO
      IGXY = IGXX+NMO
      IGXZ = IGXY+NMO
      IGYY = IGXZ+NMO
      IGYZ = IGYY+NMO
      IGZZ = IGYZ+NMO
C
      IGXXX=IGZZ+NMO
      IGXXY=IGXXX+NMO
      IGXXZ=IGXXY+NMO
      IGXYY=IGXXZ+NMO
      IGXZZ=IGXYY+NMO
      IGXYZ=IGXZZ+NMO
      IGYYY=IGXYZ+NMO
      IGYYZ=IGYYY+NMO
      IGYZZ=IGYYZ+NMO
      IGZZZ=IGYZZ+NMO
C
      IGXXXX=IGZZZ+NMO
      IGXXXY=IGXXXX+NMO
      IGXXXZ=IGXXXY+NMO
      IGXXYY=IGXXXZ+NMO
      IGXXZZ=IGXXYY+NMO
      IGXYYY=IGXXZZ+NMO
      IGXZZZ=IGXYYY+NMO
      IGXXYZ=IGXZZZ+NMO
      IGXYYZ=IGXXYZ+NMO
      IGXYZZ=IGXYYZ+NMO
      IGYYYY=IGXYZZ+NMO
      IGYYYZ=IGYYYY+NMO
      IGYYZZ=IGYYYZ+NMO
      IGYZZZ=IGYYZZ+NMO
      IGZZZZ=IGYZZZ+NMO
      ICOFMX=IGZZZZ+NMO
C
      DO 100 I = 1,NCENT
        READ (IWFN,103) ATNAM(I),J,CO(IXC+J),CO(IYC+J),CO(IZC+J),
     +                  CO(ICHARG+J)
100   CONTINUE
      READ (IWFN,104) (IC(ICENT+I),I=1,NPRIMS)
      READ (IWFN,104) (IC(ITYPE+I),I=1,NPRIMS)
      READ (IWFN,105) (CO(IE+I),I=1,NPRIMS)
C
      DO 130 I = 1,NPRIMS
       IF(IC(ITYPE+I).GT.NTYPE) GOTO 999
 130  CONTINUE
C
      DO 120 I = 1,NMO
        READ (IWFN,106) CO(IP+I),CO(IEORB+I)
        K = NPRIMS*(I-1)+IMO
        READ (IWFN,107) (CO(K+J),J=1,NPRIMS)
120   Continue
C
      Do 125 I=1,NPrims
      check=zero
      Do 126 J=1,NMO
      temp=dabs(CO(IMO+NPRIMS*(J-1)+I))
      If(temp.gt.check)check=temp
126   Continue
      CO(ICOFMX+I)=check
125   Continue
      
      READ (IWFN,108) CHECK
      IF (CHECK .NE. ENDATA) STOP ' RDWFN : END CARD NOT FOUND '
C
C    READ IN TOTAL SCF ENERGY AND -V/T
C
      READ (IWFN,109) TOTE,GAMMA
      GOTO 9999
  999 STOP 'EXTREME CANNOT WORK WITH g-, h- OR HIGHER PRIMITIVES'
 9999 CONTINUE
      RETURN
C
101   FORMAT (A80)
102   FORMAT (4X,A4,12X,3(I3,17X))
103   FORMAT(A8,11X,I3,2X,3F12.8,10X,F5.1)
104   FORMAT (20X,20I3)
105   FORMAT (10X,5E14.7)
106   FORMAT(35X,F12.8,15X,F12.8)
107   FORMAT (5E16.8)
108   FORMAT(A8)
109   FORMAT(17X,F20.12,18X,F13.8)
      END
      Subroutine Store(XYZ,X,Y,Z,J1,J2,Icancl,INP,JI1,JI2,NITER)
C
      Implicit Double Precision (A-H,O-Z)
C
      Parameter(MaxCrt=500)
      COMMON /UNITS/  INPT,IOUT,IWFN,IWLP
      Common /options/ Icut,Iprint,Eps,Epsnuc,Dmp,DmpNuc
      Dimension XYZ(3),X(MaxCrt),Y(MaxCrt),Z(MaxCrt),JI1(MaxCrt),
     $  JI2(MaxCrt)
      Save Close
      Data Close/1.d-6/
1000  FORMAT(' NEW CRITICAL POINT FOUND:  NUMBER ',I4)
1010  FORMAT(' REDUNDANT CRITICAL POINT FOUND:  SAME AS NUMBER ',I4)
1020  FORMAT(' NUMBER OF NEWTON-RAPHSON ITERATIONS : ',I4)
C
      ICancl=0
      If(INP.gt.0)Then
      DO 10 I=1,INP
      Dist=dsqrt((xyz(1)-X(I))**2+(xyz(2)-Y(I))**2+(XYZ(3)-Z(I))**2)
      If(Dist.lt.close)Icancl=1
      IF(Dist.lt.Close)ISame=I
10    Continue
      Endif
C
      If(Icancl.eq.0)Then
      INP = INP + 1
      X(INP) = XYZ(1)
      Y(INP) = XYZ(2)
      Z(INP) = XYZ(3)
      JI1(INP) = J1
      JI2(INP) = J2
      Write(iout,*)
      Write(IOUT,1000)INP
      If(Iprint.eq.1)Write(iout,*)
      If(Iprint.eq.1)Write(IOUT,1020)NITER
      Write(iwlp,*)
      Write(IWLP,1000)INP
      Write(iwlp,*)
      Write(IWLP,1020)NITER
      Else
      If(Iprint.eq.1)Write(iout,*)
      If(Iprint.eq.1)Write(iwlp,*)
      If(Iprint.eq.1)Write(IOUT,1010)Isame
      If(Iprint.eq.1)Write(IWLP,1010)Isame
      Endif
C
      Return 
      End
        SUBROUTINE      TQLGRM	(N, D, E, Z, IERR)
        IMPLICIT        DOUBLE PRECISION (A-H, O-Z)
        DIMENSION       D(*), E(*), Z(N,N)
	PARAMETER (AMACH = 16.0D-13)
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
        DIMENSION       D(*), E(*), Z(N,N)
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
