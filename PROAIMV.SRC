      PROGRAM PROAIMV
C
C    Version 94 Revision B
C
C    THIS PROGRAM DETERMINES PROPERTIES OF ATOMS IN MOLECULES FROM
C    AB-INITIO MOLECULAR WAVEFUNCTIONS.  THIS IS DONE BY INTEGRATING
C    CORRESPONDING PROPERTY DENSITIES OVER THE ATOMIC BASINS,
C    THE BOUNDARIES OF WHICH ARE SURFACES HAVING LOCAL ZERO-FLUX OF THE 
C    GRADIENT OF THE ELECTRON DENSITY (GRADRHO).  THE ATOMIC SURFACES
C    ARE FOUND EITHER BY WALKING ALONG TRAJECTORIES OF GRADRHO FROM 
C    THE (3,-1) BOND CRITICAL POINTS OF RHO WHICH ARE IN THE ATOMIC 
C    SURFACE ("PROAIM") OR BY FINDING THE POINTS ALONG EACH INTEGRATION
C    RAY WHERE THE ATTRACTOR OF THE GRADRHO TRAJECTORIES INTERSECTING
C    THE RAY CHANGES FROM THE NUCLEUS OF THE INTEGRATED ATOM TO ANOTHER
C    NUCLEUS, OR VICE VERSA.  THE LATTER SURFACE ALGORITHM IS CALLED 
C    "PROMEGA."
C
C    For Information on the Original PROAIM program and a description
C    of the output please see:
C    "CALCULATION OF THE AVERAGE PROPERTIES OF ATOMS IN 
C    MOLECULES. II"; F.W. Biegler Konig, R.F.W. Bader, T. Tang; Journal
C    of Computational Chemistry;  Volume 13 (No. 2); 1982
C
C     QUESTIONS AND SUGGESTIONS SHOULD BE DIRECTED TO:
C     Richard Bader   McMASTER UNIVERSITY, DEPT. OF CHEMISTRY,
C     HAMILTON, ONTARIO CANADA
C     BITNET ADDRESS:  BADER@MCMAIL.CIS.MCMASTER.CA
C     or 
C     Todd A. Keith
C     keith@babbage.chemistry.mcmaster.ca
C
C     JANUARY-MARCH 1981:
C     PROAIM WRITTEN BY F.W. BIEGLER-KOENIG AND J.A. DUKE MCMASTER UNIV.
C     NOVEMBER 1988-FEBRUARY 1989:
C     PROAIM MODIFIED AND STREAMLINED KEITH E LAIDIG MCMASTER UNIV.
C
C     INTEGRATION STRUCTURE HEAVILY MODIFIED,VECTORIZED --> "PROAIMV"
C     TODD A. KEITH:  McMASTER UNIVERSITY, HAMILTON ONTARIO 1991
C
C     PRIMITIVE CUTOFF ALGORITHM INCORPORATED 
C     TAK  JUNE 1, 1992
C
C     SEPARATION OF BETA SPHERE INTEGRATION FROM OUTER INTEGRATION
C     INCORPORATED.  TAK JUNE 10 1992
C
C     ADDITION OF F FUNCTIONS RICHARD BONE 1993
C
C     MAJOR CLEANUP AND ALLOWANCE OF ARBITRARY EVEN-ORDER QUADRATURE
C     FOR THETA, PHI AND RADIAL INTEGRATIONS.  ADDITION OF PROMEGA 
C     SURFACE ALGORITHM AND CHANGE OF INPUT/OUTPUT FORMAT.  (PROMEGA 
C     SURFACE ALGORITHM DEVELOPED BY TAK and James R. Cheeseman.)
C     TAK 12/93
C
C     INCORPORATION OF MAGNETIC PROPERTIES CAPABILITY 
C     TAK 3/94
C
C     Generalization of Atomic Overlap Matrix (AOM) and AOM derived
C     properties to ROHF, UHF and natural orbital wavefunctions.  
C     JRC and TAK 3/94
C
C     Fixed bug for F-functions in "gauscheck".  TAK 3/94
C
C     PROAIMV CAN HANDLE S,P,D(6) and F(10) TYPE GAUSSIAN FUNCTIONS
C              
C     The maximum number of Theta, Phi and Radial points are determined
C     by the parameters MaxTht, MaxPhi and MaxRad.  Presently, they are
C     set at 200, which is larger than anyone will probably ever want ...
C
C     The maximum number of nuclei, molecular orbitals and primitives
C     are determined by the parameters Mcent, MMO and MPRIMS.  Presently
C     they are set at 50, 100 and 500 respectively.
C
C**********************************************************************
C     TO LOWER MEMORY REQUIREMENTS, LOWER THE PARAMETER MPTS IN GAUS3,
C     FUNC3, INTARC and TRUDGE3 - BUT THE LOWER MPTS, THE SLOWER THE JOB.
C     DO NOT CHANGE THE PARAMETERS MPTS0 AND MPTSX IN GAUSCHECK UNLESS
C     THE CUTOFF ALGORITHM IS MODIFIED.
C**********************************************************************
C
C     To Run PROAIMV, two files are needed:  the AIMPAC wavefunction
C     file (with the extension ".wfn") and the integration input file
C     (with the extension ".inp").  PROAIMV produces an output file
C     with the extension ".int".  For example, to run PROAIMV on a 
C     carbon atom of c4h4 (with the executable "proaimv"):
C
C     proaimv c4h4_c1 c4h4
C
C     where c4h4_c1 refers to the file "c4h4_c1.inp" and c4h4 refers
C     to the wavefunction file "c4h4.wfn".
C
C     The input for PROAIMV is dependent upon which Surface Algorithm
C     is to be used.  The first surface algorithm (PROAIM) is usually
C     faster but requires more user input and may sometimes seriously
C     fail.  The second surface algorithm (PROMEGA) is slower but
C     rarely fails. 
C
C     As an Example, an Input File for a PROAIM Job is given below 
C     within the Starred (*) box for a carbon atom of 
C     tetrahedrane (C4H4):
C
C        *************************************************
CARD1    *C4H4_C1 RHF/6-31G**                            *
CARD2    *  C    1                                       *
CARD3    *PROAIM                                         *
CARD4    *  4 3 1                                        *
CARD5    *9.71207959E-10  1.98184301E-09  1.15576886E+00 *  (C1-C2 bcp)
CARD6    *1.15576886E+00 -4.37939987E-09 -2.58361695E-09 *  (C1-C3 bcp)
CARD7    *-8.22575222E-10  1.15576886E+00  5.95489085E-09*  (C1-C4 bcp)
CARD8    *1.72948569E+00  1.72948569E+00  1.72948568E+00 *  (C1-H5 bcp)
CARD9    *3.99523990E-01 -3.99523987E-01  3.99523978E-01 *  (C1-C2-C3 rcp)
CARD10   *-3.99523982E-01  3.99523994E-01  3.99523989E-01*  (C1-C2-C4 rcp)
CARD11   *3.99523989E-01  3.99523987E-01 -3.99523978E-01 *  (C1-C3-C4 rcp)
CARD12   *1.14033776E-08  2.52500697E-09 -3.50325120E-10 *  (C1-C2-C3-C4 ccp)
CARD13   *1 2 8 0                                        *
CARD14   *1 3 8 0                                        *
CARD15   *2 3 8 0                                        *
CARD16   *64 48 96                                       *
CARD17   *OPTIONS                                        *
CARD18   *INTEGER 1                                      *
CARD19   *6 1                                            *
CARD20   *REAL 2                                         *
CARD21   *1 9.0                                          *
CARD22   *4 1.0D-8                                       *
C        ************************************************* 
C
Card1 is the job title (up to 80 characters)
Card2 is the Atom Name in A4I4 Format
Card3 specifies that the PROAIM surface algorithm is to be used
Card4 specifies the number of bond critical points (bcp), ring critical
C     points (rcp) and cage critical points lying within the surface
C     of the carbon atom
Cards5-8 specify the coordinates of the bond critical points (in au).
cards9-11 specify the coordinates of the ring critical points (in au).
card12 specifies the coordinates of the cage critical point (in au).
Cards13-15 specify how the ring critical points are linked to the bond
C          critical points and the cage critical points.  Thus, Card13
C          says that the first rcp is connected to the first two bcp's
C          and to the cage critical point.  THE INFORMATION ABOUT THE
C          CRITICAL POINTS MUST BE DETERMINED BY THE PROGRAM "EXTREME".
Card16 Specifies three numerical integration parameters:  the number
C      of Phi planes, the number of theta planes and the number of 
C      radial points to be used per integration ray within the Beta Sphere.
Card17 Specifies whether the user wishes to change the default 
C      parameters of the program.  If only defualt values are to be
C      used then the OPTIONS card should be absent.
Card18 Specifies that one integer parameter is to be specified in the 
C      following cards.
Card19 specifies that integer parameter 6 is to be set to the value 1
Card20 Specifies that two real parameters are to be specified in the
C      following cards.
Card21 Specifies that the real parameter 1 is set to the value 9.0
Card22 Specifies that the real parameter 4 is set to the value 1.0D-8
C
C     A description of the parameters and their default values is as
C     follows if PROAIM is specified:
C
C INTEGER OPTIONS:
C    (1) = whether pre-job primitive cutoffs are to be used: 0/1 = No/Yes
C          Default is 1.
C    (2) = Number of points per gradrho path in the surface tracing
C          Default is 140
C    (3) = Number of basic gradrho paths used in the surface tracing
C          Default is 80
C    (4) = Multiple of the default number of radial points to be
C          used for integration outside the BEta sphere and the default
C          number of theta and phi planes to be used inside the Beta
C          Sphere.  Default is 1.
C    (5) = Maximum number of gradrho paths which can be inserted between
C          adjacent basic paths.  Default is 6
C
C    (6) = Whether to calculate the atomic overlap matrix: 0/1 = No/Yes
C          Default is 0.  For correlated wavefunctions, where there
C          are a relatively large number of Molecular orbitals, the
C          computation of the AOM can become dominant in terms of CPU
C          time.
C    (7) = Whether to calculate second-order magnetic properties - namely
C          the atomic contribution to the shielding tensors of the nuclei
C          and the atomic magnetic susceptibility tensor and
C          the atomic "net current tensor".  Note that for magnetic 
C          properties, the first-order wavefunctions for the Lx, Ly, Lz
C          Px, Py and Pz perturbations are required in addition to the
C          unperturbed wavefunction.  Default is no (0).
C
C    (8) = Type of Gauge Transformations to perform to calculate the
C          Current Distribution Within the atom, and hence the atom's
C          other magnetic properties.  0 = Use IGAIM method - gauge origin 
C          coincident with the nucleus of the integrated atom.  
C          1 = use another single gauge origin - the gauge origin should be
C          specified (X,Y,Z) on the following card.  2=Becke-Igaim.
C          Default is 0 (IGAIM).
C          
C    (9) = If this is an unrestricted wfn, the number of the first beta MO.
C          Default = 0.
C
C REAL OPTIONS:
C    (1) = Maximum distance from the nucleus of the integrated atom
C          to integrate to.  Default is 9.0 au
C    (2) = Value of the first rho isosurface.  Default is 0.001 au.
C    (3) = Value of the second rho isosurface.  Default is 0.002 au.
C    (4) = Cutoff Value for the primitive cutoff algorithms.  Smaller
C          is more accurate but more time-consuming.  Default is 1.0D-9.
C    (5) = Maximum Allowable Distance between ends of adjacent 
C          gradrho paths.  Default is 0.6 au.
C    (6) = Length of Grad Rho Paths in Surface Tracing.  
C          Default is 8.0 au.
C
C     The input file for a PROMEGA job is simpler.  An example is given
C     below in the starred (*) box for the same carbon atom:
C
C        **********************************************
CARD1    *C4H4_C1 RHF/6-31G**                         *
CARD2    *  C    1                                    *
CARD3    *PROMEGA                                     *
CARD4    *64 48 96                                    *
CARD5    *OPTIONS                                     *
CARD6    *INTEGER 1                                   *
CARD7    *3 4                                         *
CARD8    *REAL 2                                      *
CARD9    *6 0.0015                                    *
CARD10   *7 0.03                                      *
C        **********************************************
C
Card1 specifies the job title A80
Card2 specifies the atom name A4I4
Card3 specifies that the PROMEGA surface algorithm is to be used
Card4 is the number of phi, theta and radial points as in Card 16 above.
Card5 specifies whether optional parameters are to be supplied in the
C     following cards.  If all default values are to be used, this
C     OPTIONS card should not be specified
Card6 Specifies that one integer parameter is to be specified in the 
C     following cards.
Card7 specifies that integer parameter 3 is to be set to the value 4.
Card8 Specifies that two real parameters are to be specified in the
C      following cards.
Card9 Specifies that the real parameter 6 is set to the value 0.0015
Card10 Specifies that the real parameter 7 is set to the value 0.03
C     
C     A description of the parameters and their default values is as
C     follows if PROMEGA is specified:
C
C INTEGER OPTIONS:
C     (1) = whether second and third intersections of the integration
C           rays with the atomic surface are searched for 1/0 = No/Yes.
C           Default is 1 (search only for first intersections).
C     (2) = Number of small steps per regular step in tracing the
C           gradrho trajectories.  Default is 5
C     (3) = Order of Adam's Bashforth-Moulton predictor corrector method
C           to be used in calculating the gradrho trajectries.
C           Defualt is 6.
C     (4) = Multiple of the default number of radial points to be
C           used for integration outside the BEta sphere and the default
C           number of Theta and Phi planes to be used inside the Beta 
C           Sphere.  Default is 1.
C
C    (6) = Whether to calculate the atomic overlap matrix: 0/1 = No/Yes
C          Default is 0.  For correlated wavefunctions, where there
C          are a relatively large number of Molecular orbitals, the
C          computation of the AOM can become dominant in terms of CPU
C          time so this option should be considered for such cases.
C
C    (7) = Whether to calculate second-order magnetic properties - namely
C          the atomic contribution to the shielding tensors of the nuclei
C          and the atomic magnetic susceptibility tensor and
C          the atomic "net current tensor".  Note that for magnetic 
C          properties, the first-order wavefunctions for the Lx, Ly, Lz
C          Px, Py and Pz perturbations are required in addition to the
C          unperturbed wavefunction.
C
C    (8) = Type of Gauge origin to use in order to calculate the
C          Current Distribution Within the atom, and hence the atom's
C          other magnetic properties.  0 = Use IGAIM method - gauge origin 
C          coincident with the nucleus of the integrated atom.  
C          1 = use another single gauge origin - the gauge origin should be
C          specified (X,Y,Z) on the following card.  2=Becke-Igaim.
C          Default is 0 (IGAIM).
C          
C    (9) = If this is an unrestricted wfn, the number of the first beta MO.
C          Default = 0.
C
C REAL OPTIONS:
C     (1) = Maximum Distance from nucleus of integrated atom to 
C           integrate to.  Default is 9.0 au
C     (2) = Value of the first rho isosurface.  Default is 0.001 au.
C     (3) = Value of the second rho isosurface.  Default is 0.002 au.
C     (4) = Value for the primitive cutoff algorithms.  Smaller is more
C          accurate but more time-consuming.  Default is 1.0D-9.
C     (5) = How far out along the integration rays to search for
C           an intersection.  Default is 7.5.
C     (6) = How close to the atomic surface to get in intersections
C           search.  Default is 0.001.  Smaller is better but more time
C           consuming.
C     (7) = Step Size for tracing the gradrho trajectories.
C           Default is 0.025 au
C     (8) = Initial Step Size along integration rays in search for 
C           First intersections.  Default is 0.25 au
C     (9) = Initial Step Size along integration rays in search for
C           Second and Third intersections.  Default is 0.025 au
C
C     Atomic Units only are used throughout ...
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 WFNTTL,JOBTTL,LINE
      CHARACTER*40 INP,INT,WFN
      CHARACTER*8 AT,ATNAM,ANUC
      CHARACTER*7 OptStr,IStr,StrIn,UpStr
      CHARACTER*6 PAIM,PMEGA,Gstrng,Gupstr
      CHARACTER*4 FINP, FWFN, FINT,RSTR
      LOGICAL RHF,ROHF,ABNAT,RNAT
      PARAMETER (MCENT=50, MaxCrt=20,MxBcrt=10,MaxPhi=200,MaxTht=200,
     $MaxRad=200,MaxPrp=200,MMO=300)
      COMMON/UNITS/ ISRF,INPT,IOUT,IWFN,IDBG
      COMMON /DIRK/ Toll,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,IMagM,NFBETA,RHF,ROHF,
     $ABNAT,RNAT,XGO,YGO,ZGO
      COMMON/C7/  CT(3,3),  X,  Y,  Z, NCRNT
      COMMON/PARA/ NINS(Mxbcrt),IBETPTS,INMAX,NNN,NPATH,DIST,
     $DMAX,Amax,ADIFF,Tinf
      COMMON/ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON/STRING/ WFNTTL,JOBTTL,ATNAM(MCENT),NAT
      COMMON/VALUES/ THRESH1,THRESH2,GAMMA,TOTE
      COMMON/NCUT/ CUTOFF,NDOCUT,NPR, NPRA, NPRB, NZEROA, NZEROB
      COMMON /LIMITS/ RMAX, NATMX, NFLIM, NPROPS
      COMMON /ORBTL/ EORB(MMO),PO(MMO),ROOTPO(MMO),NMO
      DIMENSION CRIT(3,MaxCrt),NSRC(4,Mxbcrt),IOPAIM(20),IOPEGA(20),
     $OPAIM(20),OPEGA(20)
      DATA FINP /'.inp'/, FWFN /'.wfn'/, FINT /'.int'/
      DATA ANUC /'NON NUCL'/,PAIM/'PROAIM'/,PMEGA/'PROMEG'/
      DATA OptStr/'OPTIONS'/,Istr/'INTEGER'/,Rstr/'REAL'/
      Data Pt1/0.1d0/,One/1.0d0/,Two/2.0d0/,Thresh3/1.1d0/,
     $thresh4/0.1d0/
496   Format('Optional Parameters Read From Input')
497   Format('All Default Parameters Will be Used')
498   Format('No Non-nuclear attractors Using Promega Yet',
     $' - Use Proaim Option')
500   FORMAT(A80)
501   FORMAT(A7,1X,1I1)
502   FORMAT(A4,1X,1I1)
510   FORMAT(A8)
511   FORMAT(/,'NORMAL TERMINATION OF PROAIMV')
513   Format(/,'JobTim = ',I2)
789   Format('Will calculate magnetic properties')
800   FORMAT(' PROAIMV - Version 94 - Revision B')
810   FORMAT(/,A80)
830   FORMAT(/,A80)
840   FORMAT(' -V/T FOR THIS WAVEFUNCTION = ',1F20.11)
850   FORMAT(' MOLECULAR SCF ENERGY (AU)  = ',1F20.11)
859   Format('Unrestricted Wavefunction But First Beta Orbital',/,
     $' Not Specified - Atomic Overlap Matrix Will Not be Calculated')
860   FORMAT(' INTEGRATION IS OVER ATOM ',A8)
861   FORMAT(' INTEGRATION IS OVER A Non-Nuclear Attractor')
1999  Format(A6)
2333  Format('Need to specify either proaim or promega as surface',
     $' method on Card 3')
2999  Format('PROAIM SURFACE ALGORITHM USED')
3010  Format('Requested Number of Phi Planes Exceeds Maximum of '
     $,1I6)
3011  Format('Requested Number of Theta Planes Exceeds Maximum of '
     $,1I6)
3012  Format('Requested Number of Radial Points Exceeds Maximum of '
     $,1I6)
3999  Format('PROMEGA SURFACE ALGORITHM USED')
4999  Format('Critical Points in Atomic Surface:')
5999  Format(1I3,' Bond ',1PE16.8,1X,1PE16.8,1X,1PE16.8)
6999  Format(1I3,' Ring ',1PE16.8,1X,1PE16.8,1X,1PE16.8)
7999  Format(11X,'Connected to bonds ',1I3,1X,1I3,' and cages ',
     $1I3,1X,1I3)
8999  Format(1I3,' Cage ',1PE16.8,1X,1PE16.8,1X,1PE16.8)
C
      CALL MAKNAME(1,INP,ILEN,FINP)
      IF (ILEN .EQ. 0) STOP ' usage: proaimv inpfile wfnfile ' 
      CALL MAKNAME(1,INT,ILEN,FINT)
      IF (ILEN .EQ. 0) STOP ' usage: proaimv inpfile wfnfile ' 
      CALL MAKNAME(2,WFN,ILEN,FWFN)
      IF (ILEN .EQ. 0) STOP ' usage: proaimv inpfile wfnfile ' 
C
      OPEN (INPT,FILE=INP,status='unknown')
      OPEN (IOUT,FILE=INT,status='unknown')
      OPEN (IWFN,FILE=WFN,status='unknown')
C
      WRITE(IOUT,800)
      READ(INPT,500) JOBTTL
      CALL RDPSI
      WRITE(IOUT,830) WFNTTL
      WRITE(IOUT,840) GAMMA
      WRITE(IOUT,850) TOTE
      WRITE(IOUT,810) JOBTTL
      READ(INPT,510) AT
      IF (AT .EQ. ANUC) READ(INPT,*) XNN,YNN,ZNN
      READ(INPT,1999)GSTRNG
      Call CaseUp(Gstrng,GUpStr,6)
      IF(GUpSTR(1:6).EQ.PAIM)Then
      IMEGA=0
      ElseIF(GUpSTR(1:6).EQ.PMEGA)THEN
      IMEGA=1
      Else
      Write(Iout,2333)
      Stop 'Unrecognized Surface Method Specified '
      Endif
      IF(IMEGA.EQ.0)Write(Iout,2999)
      IF(IMEGA.EQ.1)Write(Iout,3999)
      IF(AT.EQ.ANUC.AND.IMEGA.EQ.1)Then
      Write(iout,498)
      Stop 
      Endif
      IF(IMEGA.EQ.0)THEN
      READ(INPT,*) N,NRING,NCAGE
      NT = N + NRING + NCAGE
      DO 100 I = 1,NT
      READ(INPT,*) (CRIT(J,I),J=1,3)
100   CONTINUE
      IF (NRING .NE. 0) THEN
      DO 105 I = 1,NRING
      II = N + I
      READ(INPT,*) (NSRC(J,I),J=1,4)
105   CONTINUE
      ENDIF
      WRITE(IOUT,4999)
      DO 110 I=1,N
      WRITE(IOUT,5999)I,CRIT(1,I),CRIT(2,I),CRIT(3,I)
110   CONTINUE
      DO 115 I=1,NRING
      WRITE(IOUT,6999)N+I,CRIT(1,N+I),CRIT(2,N+I),CRIT(3,N+I)
      WRITE(IOUT,7999)NSRC(1,I),NSRC(2,I),NSRC(3,I),NSRC(4,I)
115   CONTINUE
      DO 120 I=1,NCAGE
      WRITE(IOUT,8999)NT-NCAGE+I,CRIT(1,NT-NCAGE+I),
     $CRIT(2,NT-NCAGE+I),CRIT(3,NT-NCAGE+I)
120   CONTINUE
      ENDIF
      READ(INPT,*) IPHIPL,IT,IBETPTS
      IDIV=2
      IF(IMEGA.EQ.1)IDIV=4
      IPhipl=Iphipl/IDIV
      Iphipl=Iphipl*IDIV
      It=It/IDIV
      It=It*IDIV
      ibetpts=ibetpts/2
      ibetpts=ibetpts*2
      If(IPHIPL.GT.MaxPhi)Write(Iout,3010)MaxPhi
      If(IPHIPL.GT.MaxPhi)STOP 'Too many Phi Planes Requested'
      If(IT.GT.MaxTHt)Write(Iout,3011)MaxTht
      If(IT.GT.MaxTht)STOP 'Too many Theta Planes Requested'
      If(IBETPTS.GT.MaxRad)Write(Iout,3012)MaxRad
      If(IBETPTS.GT.MaxRad)STOP 'Too many Radial Points Requested'
C
      IOPAIM(1)=NDOCUT
      IOPAIM(2)=NNN
      IOPAIM(3)=NPATH
      IOPAIM(4)=NACC
      IOPAIM(5)=INMAX
      IOPAIM(6)=IDOAOM
      IOPAIM(7)=IDOMAG
      IOPAIM(8)=IMagM
      IOPAIM(9)=NFBETA
      OPAIM(1)=TINF
      OPAIM(2)=THRESH1
      OPAIM(3)=THRESH2
      OPAIM(4)=CUTOFF
      OPAIM(5)=DIST
      OPAIM(6)=DMAX
C
      IOPEGA(1)=NOSEC
      IOPEGA(2)=ISTEP
      IOPEGA(3)=NABMO
      IOPEGA(4)=NACC
      IOPEGA(6)=IDOAOM
      IOPEGA(7)=IDOMAG
      IOPEGA(8)=IMAGM
      IOPEGA(9)=NFBETA
      OPEGA(1)=TINF
      OPEGA(2)=THRESH1
      OPEGA(3)=THRESH2
      OPEGA(4)=CUTOFF
      OPEGA(5)=SIZE
      OPEGA(6)=CTF
      OPEGA(7)=STP
      OPEGA(8)=FSTP
      OPEGA(9)=SSTP
C
      IF(IMEGA.EQ.0)THEN
      READ(INPT,500,End=145)LINE
      Call CaseUp(Line,UpStr,7)
      If(UpStr(1:7).eq.OPTSTR)Then
      Read(inpT,501) Line,NOPT
      Call CaseUp(Line,UpStr,7)
      If(UpStr(1:7).ne.IStr) STOP 
     $ 'INTEGER N  Card Required when OPTION card is specified'
      DO 125 I=1,NOPT
      READ(INPT,*) NOP,IOPAIM(NOP)
      If(NOP.eq.8.and.iopaim(8).eq.1)Read(Inpt,*)xgo,ygo,zgo
125   CONTINUE
      Read(inpT,502,End=126) Line,NOPT
126   Call CaseUp(Line,UpStr,4)
      If(UpStr(1:4).ne.RStr) STOP 
     $  'REAL N  Card Required when OPTION card is specified'
      DO 130 I=1,NOPT
      READ(INPT,*) NOP,OPAIM(NOP)
130   CONTINUE
      Write(Iout,496)
      Goto 146
      Endif
      ELSEIF(IMEGA.EQ.1)THEN
      READ(INPT,500,END=145)LINE
      Call CaseUp(Line,UpStr,7)
      If(UpStr(1:7).eq.OPTSTR)Then
      Read(inpT,501) Line,NOPT
      Call CaseUp(Line,UpStr,7)
      If(UpStr(1:7).ne.IStr) STOP 
     $ 'INTEGER N  Card Required when OPTION card is specified'
      DO 135 I=1,NOPT
      READ(INPT,*) NOP,IOPEGA(NOP)
      If(NOP.eq.8.and.iopega(8).eq.1)Read(Inpt,*)xgo,ygo,zgo
135   CONTINUE
      Read(inpT,502,End=136) Line,NOPT
136   Call CaseUp(Line,UpStr,4)
      If(UpStr(1:4).ne.RStr) STOP 
     $  'REAL N  Card Required when OPTION card is specified'
      DO 140 I=1,NOPT
      READ(INPT,*) NOP,OPEGA(NOP)
140   CONTINUE
      Write(Iout,496)
      Goto 146
      Endif
      Endif
C
145   Continue
      Write(Iout,497)
146   Continue
C
      IF(IMEGA.EQ.0)THEN
      NDOCUT=IOPAIM(1)
      NNN=IOPAIM(2)
      NPATH=IOPAIM(3)
      NACC=IOPAIM(4)
      INMAX=IOPAIM(5)
      IDOAOM=IOPAIM(6)
      IDOMAG=IOPAIM(7)
      IMAGM=IOPAIM(8)
      NFBETA=IOPAIM(9)
      TINF=OPAIM(1)
      THRESH1=OPAIM(2)
      THRESH2=OPAIM(3)
      CUTOFF=OPAIM(4)
      DIST=OPAIM(5)
      DMAX=OPAIM(6)
      ELSE
      NOSEC=IOPEGA(1)
      ISTEP=IOPEGA(2)
      NABMO=IOPEGA(3)
      NACC=IOPEGA(4)
      IDOAOM=IOPEGA(6)
      IDOMAG=IOPEGA(7)
      IMAGM=IOPEGA(8)
      NFBETA=IOPEGA(9)
      TINF=OPEGA(1)
      THRESH1=OPEGA(2)
      THRESH2=OPEGA(3)
      CUTOFF=OPEGA(4)
      SIZE=OPEGA(5)
      CTF=OPEGA(6)
      STP=OPEGA(7)
      FSTP=OPEGA(8)
      SSTP=OPEGA(9)
      If(Size.Ge.Tinf)Size=Tinf-Pt1
      ENDIF
C
      If(Imega.eq.1)Ndocut=0
      lmo=nmo
      If(Idomag.eq.1)Then
      lmo=nmo/7
      NDoCut=0
      Write(iout,789)
      NProps=NProps+18+9*Ncent
      If(Nprops.gt.MaxPrp)Stop 'Redimension MaxPrp'
      Endif
C
C     figure out what type of wfn this is
C
      If(IDOAOM.eq.1)Then
      occmax=-Two
      Occmin=Two
      Do 1 I=1,LMO
      If(po(i).gt.occmax)occmax=po(i)
      If(po(i).lt.occmin)occmin=po(i)
1     Continue
      If((OccMax.le.Thresh3) .and. Nfbeta.gt.0) ABNat = .True. 
      If(OccMax.ge.thresh3.and.occmax.lt.Two.and.occmin.lt.thresh4)
     $RNAT=.true.
      If(OccMin.eq.Two) RHF = .True. 
      If(OccMin.eq.One .and. OccMax.eq.Two) ROHF = .True. 
      If(OccMax.le.Thresh3.and.nfbeta.eq.0)Then
      Write(iout,859)
      Idoaom=0
      Endif
      Endif
C
      IF (AT .EQ. ANUC) THEN
      Write(Iout,861)
      X = XNN
      Y = YNN
      Z = ZNN
      NCRNT = -1
      ELSE 
      DO 150 II = 1,NCENT
      IF (ATNAM(II) .EQ. AT) THEN
      Write(Iout,860)ATNAM(II)
      X = XC(II) 
      Y = YC(II) 
      Z = ZC(II) 
      NCRNT = II
      NAT = II
      GOTO 155
      END IF
150   CONTINUE
      STOP ' ATOM NAME NOT FOUND '
      END IF
155   CALL INTEG(BETA,IPHIPL,IT,CRIT,N,NRING,NCAGE,NSRC,IMEGA)
      CALL RESULT
      WRITE(IOUT,511)
C
      END
      DOUBLE PRECISION FUNCTION DASUM (N,DX,INCX)                       SHO18080
C                                                                       SHO18090
C                                  SPECIFICATIONS FOR ARGUMENTS         SHO18100
      DOUBLE PRECISION   DX(1)                                          SHO18110
      INTEGER            N,INCX                                         SHO18120
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   SHO18130
      INTEGER            I,M,MP1,NS                                     SHO18140
C                                  FIRST EXECUTABLE STATEMENT           SHO18150
      DASUM = 0.D0                                                      SHO18160
      IF (N.LE.0) RETURN                                                SHO18170
      IF (INCX.EQ.1) GO TO 10                                           SHO18180
C                                  CODE FOR INCREMENTS NOT EQUAL TO 1.  SHO18190
      NS = N*INCX                                                       SHO18200
      DO 5 I=1,NS,INCX                                                  SHO18210
         DASUM = DASUM+DABS(DX(I))                                      SHO18220
    5 CONTINUE                                                          SHO18230
      RETURN                                                            SHO18240
C                                  CODE FOR INCREMENTS EQUAL TO 1.      SHO18250
C                                    CLEAN-UP LOOP SO REMAINING VECTOR  SHO18260
C                                    LENGTH IS A MULTIPLE OF 6.         SHO18270
   10 M = N-(N/6)*6                                                     SHO18280
      IF (M.EQ.0) GO TO 20                                              SHO18290
      DO 15 I=1,M                                                       SHO18300
         DASUM = DASUM+DABS(DX(I))                                      SHO18310
   15 CONTINUE                                                          SHO18320
      IF (N.LT.6) RETURN                                                SHO18330
   20 MP1 = M+1                                                         SHO18340
      DO 25 I=MP1,N,6                                                   SHO18350
         DASUM = DASUM+DABS(DX(I))+DABS(DX(I+1))+DABS(DX(I+2))          SHO18360
     1   +DABS(DX(I+3))+DABS(DX(I+4))+DABS(DX(I+5))                     SHO18370
   25 CONTINUE                                                          SHO18380
      RETURN                                                            SHO18390
      END                                                               SHO18400
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
      SUBROUTINE  DCOPY(N,DX,INCX,DY,INCY)
      DOUBLE PRECISION DX(1),DY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
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
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(1),INFO
      DOUBLE PRECISION A(LDA,1)
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
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
      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
      INTEGER LDA,N,IPVT(1),JOB
      DOUBLE PRECISION A(LDA,1),B(1)
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DNRM2 (N,DX,INCX)                       SHO18810
C                                                                       SHO18820
C                                  SPECIFICATIONS FOR ARGUMENTS         SHO18830
      INTEGER            N,INCX                                         SHO18840
      DOUBLE PRECISION   DX(1)                                          SHO18850
C                                  SPECIFICATIONS FOR LOCAL VARIABLES   SHO18860
      INTEGER            I,J,NEXT,NN                                    SHO18870
      DOUBLE PRECISION   CUTLO,CUTHI,SUM,XMAX,ZERO,ONE,HITEST           SHO18880
      DATA               ZERO, ONE /0.0D0, 1.0D0/                       SHO18890
      DATA               CUTLO, CUTHI / 8.232D-11,  1.304D19 /          SHO18900
C                                  FIRST EXECUTABLE STATEMENT           SHO18910
      IF (N.GT.0) GO TO 5                                               SHO18920
      DNRM2 = ZERO                                                      SHO18930
      GO TO 70                                                          SHO18940
C                                                                       SHO18950
    5 ASSIGN 15 TO NEXT                                                 SHO18960
      SUM = ZERO                                                        SHO18970
      NN = N*INCX                                                       SHO18980
C                                  BEGIN MAIN LOOP                      SHO18990
      I = 1                                                             SHO19000
   10 GO TO NEXT, (15,20,35,40)                                         SHO19010
   15 IF (DABS(DX(I)).GT.CUTLO) GO TO 55                                SHO19020
      ASSIGN 20 TO NEXT                                                 SHO19030
      XMAX = ZERO                                                       SHO19040
C                                  PHASE 1. SUM IS ZERO                 SHO19050
   20 IF (DX(I).EQ.ZERO) GO TO 65                                       SHO19060
      IF (DABS(DX(I)).GT.CUTLO) GO TO 55                                SHO19070
C                                  PREPARE FOR PHASE 2.                 SHO19080
      ASSIGN 35 TO NEXT                                                 SHO19090
      GO TO 30                                                          SHO19100
C                                  PREPARE FOR PHASE 4.                 SHO19110
   25 I = J                                                             SHO19120
      ASSIGN 40 TO NEXT                                                 SHO19130
      SUM = (SUM/DX(I))/DX(I)                                           SHO19140
   30 XMAX = DABS(DX(I))                                                SHO19150
      GO TO 45                                                          SHO19160
C                                  PHASE 2. SUM IS SMALL. SCALE TO      SHO19170
C                                    AVOID DESTRUCTIVE UNDERFLOW.       SHO19180
   35 IF (DABS(DX(I)).GT.CUTLO) GO TO 50                                SHO19190
C                                  COMMON CODE FOR PHASES 2 AND 4. IN   SHO19200
C                                    PHASE 4 SUM IS LARGE. SCALE TO     SHO19210
C                                    AVOID OVERFLOW.                    SHO19220
   40 IF (DABS(DX(I)).LE.XMAX) GO TO 45                                 SHO19230
      SUM = ONE+SUM*(XMAX/DX(I))**2                                     SHO19240
      XMAX = DABS(DX(I))                                                SHO19250
      GO TO 65                                                          SHO19260
C                                                                       SHO19270
   45 SUM = SUM+(DX(I)/XMAX)**2                                         SHO19280
      GO TO 65                                                          SHO19290
C                                  PREPARE FOR PHASE 3.                 SHO19300
   50 SUM = (SUM*XMAX)*XMAX                                             SHO19310
C                                  FOR REAL OR D.P. SET HITEST =        SHO19320
C                                    CUTHI/N FOR COMPLEX SET HITEST =   SHO19330
C                                    CUTHI/(2*N)                        SHO19340
   55 HITEST = CUTHI/FLOAT(N)                                           SHO19350
C                                  PHASE 3. SUM IS MID-RANGE. NO        SHO19360
C                                    SCALING.                           SHO19370
      DO 60 J=I,NN,INCX                                                 SHO19380
         IF (DABS(DX(J)).GE.HITEST) GO TO 25                            SHO19390
   60 SUM = SUM+DX(J)**2                                                SHO19400
      DNRM2 = DSQRT(SUM)                                                SHO19410
      GO TO 70                                                          SHO19420
C                                                                       SHO19430
   65 CONTINUE                                                          SHO19440
      I = I+INCX                                                        SHO19450
      IF (I.LE.NN) GO TO 10                                             SHO19460
C                                  END OF MAIN LOOP. COMPUTE SQUARE     SHO19470
C                                    ROOT AND ADJUST FOR SCALING.       SHO19480
      DNRM2 = XMAX*DSQRT(SUM)                                           SHO19490
   70 CONTINUE                                                          SHO19500
      RETURN                                                            SHO19510
      END                                                               SHO19520
      SUBROUTINE DQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,JOB)
      INTEGER LDX,N,P,JOB
      INTEGER JPVT(P)
      DOUBLE PRECISION X(N,P),QRAUX(P),WORK(P)
C
C     DQRDC USES HOUSEHOLDER TRANSFORMATIONS TO COMPUTE THE QR
C     FACTORIZATION OF AN N BY P MATRIX X.  COLUMN PIVOTING
C     BASED ON THE 2-NORMS OF THE REDUCED COLUMNS MAY BE
C     PERFORMED AT THE USERS OPTION.
C
C     ON ENTRY
C
C        X       DOUBLE PRECISION(LDX,P), WHERE LDX .GE. N.
C                X CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO BE
C                COMPUTED.
C
C        LDX     INTEGER.
C                LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C        N       INTEGER.
C                N IS THE NUMBER OF ROWS OF THE MATRIX X.
C
C        P       INTEGER.
C                P IS THE NUMBER OF COLUMNS OF THE MATRIX X.
C
C        JPVT    INTEGER(P).
C                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION
C                OF THE PIVOT COLUMNS.  THE K-TH COLUMN X(K) OF X
C                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE
C                VALUE OF JPVT(K).
C
C                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL
C                                      COLUMN.
C
C                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE COLUMN.
C
C                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL COLUMN.
C
C                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL COLUMNS
C                ARE MOVED TO THE BEGINNING OF THE ARRAY X AND FINAL
C                COLUMNS TO THE END.  BOTH INITIAL AND FINAL COLUMNS
C                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY
C                FREE COLUMNS ARE MOVED.  AT THE K-TH STAGE OF THE
C                REDUCTION, IF X(K) IS OCCUPIED BY A FREE COLUMN
C                IT IS INTERCHANGED WITH THE FREE COLUMN OF LARGEST
C                REDUCED NORM.  JPVT IS NOT REFERENCED IF
C                JOB .EQ. 0.
C
C        WORK    DOUBLE PRECISION(P).
C                WORK IS A WORK ARRAY.  WORK IS NOT REFERENCED IF
C                JOB .EQ. 0.
C
C        JOB     INTEGER.
C                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.
C                IF JOB .EQ. 0, NO PIVOTING IS DONE.
C                IF JOB .NE. 0, PIVOTING IS DONE.
C
C     ON RETURN
C
C        X       X CONTAINS IN ITS UPPER TRIANGLE THE UPPER
C                TRIANGULAR MATRIX R OF THE QR FACTORIZATION.
C                BELOW ITS DIAGONAL X CONTAINS INFORMATION FROM
C                WHICH THE ORTHOGONAL PART OF THE DECOMPOSITION
C                CAN BE RECOVERED.  NOTE THAT IF PIVOTING HAS
C                BEEN REQUESTED, THE DECOMPOSITION IS NOT THAT
C                OF THE ORIGINAL MATRIX X BUT THAT OF X
C                WITH ITS COLUMNS PERMUTED AS DESCRIBED BY JPVT.
C
C        QRAUX   DOUBLE PRECISION(P).
C                QRAUX CONTAINS FURTHER INFORMATION REQUIRED TO RECOVER
C                THE ORTHOGONAL PART OF THE DECOMPOSITION.
C
C        JPVT    JPVT(K) CONTAINS THE INDEX OF THE COLUMN OF THE
C                ORIGINAL MATRIX THAT HAS BEEN INTERCHANGED INTO
C                THE K-TH COLUMN, IF PIVOTING WAS REQUESTED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     DQRDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2
C     FORTRAN DABS,DMAX1,MIN0,DSQRT
C
C     INTERNAL VARIABLES
C
      INTEGER J,JP,L,LP1,LUP,MAXJ,PL,PU
      DOUBLE PRECISION MAXNRM,DNRM2,TT
      DOUBLE PRECISION DDOT,NRMXL,T
      LOGICAL NEGJ,SWAPJ
C
C
      PL = 1
      PU = 0
      IF (JOB .EQ. 0) GO TO 60
C
C        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS
C        ACCORDING TO JPVT.
C
         DO 20 J = 1, P
            SWAPJ = JPVT(J) .GT. 0
            NEGJ = JPVT(J) .LT. 0
            JPVT(J) = J
            IF (NEGJ) JPVT(J) = -J
            IF (.NOT.SWAPJ) GO TO 10
               IF (J .NE. PL) CALL DSWAP(N,X(1,PL),1,X(1,J),1)
               JPVT(J) = JPVT(PL)
               JPVT(PL) = J
               PL = PL + 1
   10       CONTINUE
   20    CONTINUE
         PU = P
         DO 50 JJ = 1, P
            J = P - JJ + 1
            IF (JPVT(J) .GE. 0) GO TO 40
               JPVT(J) = -JPVT(J)
               IF (J .EQ. PU) GO TO 30
                  CALL DSWAP(N,X(1,PU),1,X(1,J),1)
                  JP = JPVT(PU)
                  JPVT(PU) = JPVT(J)
                  JPVT(J) = JP
   30          CONTINUE
               PU = PU - 1
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
C
C     COMPUTE THE NORMS OF THE FREE COLUMNS.
C
      IF (PU .LT. PL) GO TO 80
      DO 70 J = PL, PU
         QRAUX(J) = DNRM2(N,X(1,J),1)
         WORK(J) = QRAUX(J)
   70 CONTINUE
   80 CONTINUE
C
C     PERFORM THE HOUSEHOLDER REDUCTION OF X.
C
      LUP = MIN0(N,P)
      DO 200 L = 1, LUP
         IF (L .LT. PL .OR. L .GE. PU) GO TO 120
C
C           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT
C           INTO THE PIVOT POSITION.
C
            MAXNRM = 0.0D0
            MAXJ = L
            DO 100 J = L, PU
               IF (QRAUX(J) .LE. MAXNRM) GO TO 90
                  MAXNRM = QRAUX(J)
                  MAXJ = J
   90          CONTINUE
  100       CONTINUE
            IF (MAXJ .EQ. L) GO TO 110
               CALL DSWAP(N,X(1,L),1,X(1,MAXJ),1)
               QRAUX(MAXJ) = QRAUX(L)
               WORK(MAXJ) = WORK(L)
               JP = JPVT(MAXJ)
               JPVT(MAXJ) = JPVT(L)
               JPVT(L) = JP
  110       CONTINUE
  120    CONTINUE
         QRAUX(L) = 0.0D0
         IF (L .EQ. N) GO TO 190
C
C           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L.
C
            NRMXL = DNRM2(N-L+1,X(L,L),1)
            IF (NRMXL .EQ. 0.0D0) GO TO 180
               IF (X(L,L) .NE. 0.0D0) NRMXL = DSIGN(NRMXL,X(L,L))
               CALL DSCAL(N-L+1,1.0D0/NRMXL,X(L,L),1)
               X(L,L) = 1.0D0 + X(L,L)
C
C              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
C              UPDATING THE NORMS.
C
               LP1 = L + 1
               IF (P .LT. LP1) GO TO 170
               DO 160 J = LP1, P
                  T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
                  CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
                  IF (J .LT. PL .OR. J .GT. PU) GO TO 150
                  IF (QRAUX(J) .EQ. 0.0D0) GO TO 150
                     TT = 1.0D0 - (DABS(X(L,J))/QRAUX(J))**2
                     TT = DMAX1(TT,0.0D0)
                     T = TT
                     TT = 1.0D0 + 0.05D0*TT*(QRAUX(J)/WORK(J))**2
                     IF (TT .EQ. 1.0D0) GO TO 130
                        QRAUX(J) = QRAUX(J)*DSQRT(T)
                     GO TO 140
  130                CONTINUE
                        QRAUX(J) = DNRM2(N-L,X(L+1,J),1)
                        WORK(J) = QRAUX(J)
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
  170          CONTINUE
C
C              SAVE THE TRANSFORMATION.
C
               QRAUX(L) = X(L,L)
               X(L,L) = -NRMXL
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
      RETURN
      END

      SUBROUTINE DQRSL(X,LDX,N,K,QRAUX,Y,QY,QTY,B,RSD,XB,JOB,INFO)
      INTEGER LDX,N,K,JOB,INFO
      DOUBLE PRECISION X(N,K),QRAUX(K),Y(N),QY(N),QTY(N),B(K),RSD(N),
     *                 XB(N)
C
C     DQRSL APPLIES THE OUTPUT OF DQRDC TO COMPUTE COORDINATE
C     TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS.
C     FOR K .LE. MIN(N,P), LET XK BE THE MATRIX
C
C            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))
C
C     FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL
C     N X P MATRIX X THAT WAS INPUT TO DQRDC (IF NO PIVOTING WAS
C     DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR
C     ORIGINAL ORDER).  DQRDC PRODUCES A FACTORED ORTHOGONAL MATRIX Q
C     AND AN UPPER TRIANGULAR MATRIX R SUCH THAT
C
C              XK = Q * (R)
C                       (0)
C
C     THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS
C     X AND QRAUX.
C
C     ON ENTRY
C
C        X      DOUBLE PRECISION(LDX,P).
C               X CONTAINS THE OUTPUT OF DQRDC.
C
C        LDX    INTEGER.
C               LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C        N      INTEGER.
C               N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST
C               HAVE THE SAME VALUE AS N IN DQRDC.
C
C        K      INTEGER.
C               K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K
C               MUST NNOT BE GREATER THAN MIN(N,P), WHERE P IS THE
C               SAME AS IN THE CALLING SEQUENCE TO DQRDC.
C
C        QRAUX  DOUBLE PRECISION(P).
C               QRAUX CONTAINS THE AUXILIARY OUTPUT FROM DQRDC.
C
C        Y      DOUBLE PRECISION(N)
C               Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED
C               BY DQRSL.
C
C        JOB    INTEGER.
C               JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS
C               THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING
C               MEANING.
C
C                    IF A.NE.0, COMPUTE QY.
C                    IF B,C,D, OR E .NE. 0, COMPUTE QTY.
C                    IF C.NE.0, COMPUTE B.
C                    IF D.NE.0, COMPUTE RSD.
C                    IF E.NE.0, COMPUTE XB.
C
C               NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB
C               AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR
C               WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING
C               SEQUENCE.
C
C     ON RETURN
C
C        QY     DOUBLE PRECISION(N).
C               QY CONNTAINS Q*Y, IF ITS COMPUTATION HAS BEEN
C               REQUESTED.
C
C        QTY    DOUBLE PRECISION(N).
C               QTY CONTAINS TRANS(Q)*Y, IF ITS COMPUTATION HAS
C               BEEN REQUESTED.  HERE TRANS(Q) IS THE
C               TRANSPOSE OF THE MATRIX Q.
C
C        B      DOUBLE PRECISION(K)
C               B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM
C
C                    MINIMIZE NORM2(Y - XK*B),
C
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT
C               IF PIVOTING WAS REQUESTED IN DQRDC, THE J-TH
C               COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J)
C               OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO DQRDC.)
C
C        RSD    DOUBLE PRECISION(N).
C               RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B,
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS
C               ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE
C               ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK.
C
C        XB     DOUBLE PRECISION(N).
C               XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B,
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO
C               THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE
C               OF X.
C
C        INFO   INTEGER.
C               INFO IS ZERO UNLESS THE COMPUTATION OF B HAS
C               BEEN REQUESTED AND R IS EXACTLY SINGULAR.  IN
C               THIS CASE, INFO IS THE INDEX OF THE FIRST ZERO
C               DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED.
C
C     THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED
C     IF THEIR COMPUTATION IS NOT REQUESTED AND IN THIS CASE
C     CAN BE REPLACED BY DUMMY VARIABLES IN THE CALLING PROGRAM.
C     TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME
C     ARRAY FOR DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A
C     FREQUENTLY OCCURING EXAMPLE IS WHEN ONE WISHES TO COMPUTE
C     ANY OF B, RSD, OR XB AND DOES NOT NEED Y OR QTY.  IN THIS
C     CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD, OR XB, WHILE
C     PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE
C     COMPUTED.  THUS THE CALLING SEQUENCE
C
C          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
C
C     WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD
C     OVERWRITING Y.  MORE GENERALLY, EACH ITEM IN THE FOLLOWING
C     LIST CONTAINS GROUPS OF PERMISSIBLE IDENTIFICATIONS FOR
C     A SINGLE CALLINNG SEQUENCE.
C
C          1. (Y,QTY,B) (RSD) (XB) (QY)
C
C          2. (Y,QTY,RSD) (B) (XB) (QY)
C
C          3. (Y,QTY,XB) (B) (RSD) (QY)
C
C          4. (Y,QY) (QTY,B) (RSD) (XB)
C
C          5. (Y,QY) (QTY,RSD) (B) (XB)
C
C          6. (Y,QY) (QTY,XB) (B) (RSD)
C
C     IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO
C     THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     DQRSL USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     BLAS DAXPY,DCOPY,DDOT
C     FORTRAN DABS,MIN0,MOD
C
C     INTERNAL VARIABLES
C
      INTEGER I,J,JJ,JU,KP1
      DOUBLE PRECISION DDOT,T,TEMP
      LOGICAL CB,CQY,CQTY,CR,CXB
C
C
C     SET INFO FLAG.
C
      INFO = 0
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      CQY = JOB/10000 .NE. 0
      CQTY = MOD(JOB,10000) .NE. 0
      CB = MOD(JOB,1000)/100 .NE. 0
      CR = MOD(JOB,100)/10 .NE. 0
      CXB = MOD(JOB,10) .NE. 0
      JU = MIN0(K,N-1)
C
C     SPECIAL ACTION WHEN N=1.
C
      IF (JU .NE. 0) GO TO 40
         IF (CQY) QY(1) = Y(1)
         IF (CQTY) QTY(1) = Y(1)
         IF (CXB) XB(1) = Y(1)
         IF (.NOT.CB) GO TO 30
            IF (X(1,1) .NE. 0.0D0) GO TO 10
               INFO = 1
            GO TO 20
   10       CONTINUE
               B(1) = Y(1)/X(1,1)
   20       CONTINUE
   30    CONTINUE
         IF (CR) RSD(1) = 0.0D0
      GO TO 250
   40 CONTINUE
C
C        SET UP TO COMPUTE QY OR QTY.
C
         IF (CQY) CALL DCOPY(N,Y,1,QY,1)
         IF (CQTY) CALL DCOPY(N,Y,1,QTY,1)
         IF (.NOT.CQY) GO TO 70
C
C           COMPUTE QY.
C
            DO 60 JJ = 1, JU
               J = JU - JJ + 1
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 50
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -DDOT(N-J+1,X(J,J),1,QY(J),1)/X(J,J)
                  CALL DAXPY(N-J+1,T,X(J,J),1,QY(J),1)
                  X(J,J) = TEMP
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
         IF (.NOT.CQTY) GO TO 100
C
C           COMPUTE TRANS(Q)*Y.
C
            DO 90 J = 1, JU
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 80
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -DDOT(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)
                  CALL DAXPY(N-J+1,T,X(J,J),1,QTY(J),1)
                  X(J,J) = TEMP
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        SET UP TO COMPUTE B, RSD, OR XB.
C
         IF (CB) CALL DCOPY(K,QTY,1,B,1)
         KP1 = K + 1
         IF (CXB) CALL DCOPY(K,QTY,1,XB,1)
         IF (CR .AND. K .LT. N) CALL DCOPY(N-K,QTY(KP1),1,RSD(KP1),1)
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 120
            DO 110 I = KP1, N
               XB(I) = 0.0D0
  110       CONTINUE
  120    CONTINUE
         IF (.NOT.CR) GO TO 140
            DO 130 I = 1, K
               RSD(I) = 0.0D0
  130       CONTINUE
  140    CONTINUE
         IF (.NOT.CB) GO TO 190
C
C           COMPUTE B.
C
            DO 170 JJ = 1, K
               J = K - JJ + 1
               IF (X(J,J) .NE. 0.0D0) GO TO 150
                  INFO = J
C           ......EXIT
                  GO TO 180
  150          CONTINUE
               B(J) = B(J)/X(J,J)
               IF (J .EQ. 1) GO TO 160
                  T = -B(J)
                  CALL DAXPY(J-1,T,X(1,J),1,B,1)
  160          CONTINUE
  170       CONTINUE
  180       CONTINUE
  190    CONTINUE
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 240
C
C           COMPUTE RSD OR XB AS REQUIRED.
C
            DO 230 JJ = 1, JU
               J = JU - JJ + 1
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 220
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  IF (.NOT.CR) GO TO 200
                     T = -DDOT(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)
                     CALL DAXPY(N-J+1,T,X(J,J),1,RSD(J),1)
  200             CONTINUE
                  IF (.NOT.CXB) GO TO 210
                     T = -DDOT(N-J+1,X(J,J),1,XB(J),1)/X(J,J)
                     CALL DAXPY(N-J+1,T,X(J,J),1,XB(J),1)
  210             CONTINUE
                  X(J,J) = TEMP
  220          CONTINUE
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
      RETURN
      END

      SUBROUTINE  DSCAL(N,DA,DX,INCX)
      DOUBLE PRECISION DA,DX(1)
      INTEGER I,INCX,M,MP1,N,NINCX
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
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
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
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
      DOUBLE PRECISION FUNCTION EPSLON (X)
      DOUBLE PRECISION X
      DOUBLE PRECISION A,B,C,EPS
      A = 4.0D0/3.0D0
   10 B = A - 1.0D0
      C = B + B + B
      EPS = DABS(C-1.0D0)
      IF (EPS .EQ. 0.0D0) GO TO 10
      EPSLON = EPS*DABS(X)
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
      SUBROUTINE LLSQB (A,IA,M,N,B,TOL,KBASIS,X,H,IP,IER)               
C
C    LLSQB IS HEAVILY BASED UPON THE IMSL ROUTINE LLSQF
C    REQUIRED TO SOLVE LINEAR LEAST SQUARES IN RING
C
      DOUBLE PRECISION   A(IA,N),B(M),TOL,X(N),H(N)                     
      DOUBLE PRECISION   BB,DLOSS,DLOSSJ,RCOND,RCONDJ,RNORM,TMP,XNORM   
      DOUBLE PRECISION   DASUM,DDOT,DNRM2  
      INTEGER            I,IER,J,JCOL,JJ,JSTART,K,KP1,L,LDIAG,LMAX      
      INTEGER            IA,M,N,KBASIS,IP(N)                         
C
      LDIAG = MIN0(M,N)                                                 
      IER = 129                                                         
      IF (LDIAG.LE.0) GO TO 9000                                        
      IER = 130                                                         
      IF (TOL.GT.1.0D0) GO TO 9000                                      
      IER = 0                                                           
      JSTART = MAX0(KBASIS+1,1)                                         
      DO 35 J=1,LDIAG                                                   
         IP(J) = J                                                      
         IF (J.LE.KBASIS) GO TO 30                                      
         LMAX = J                                                       
         IF (J.EQ.JSTART) GO TO 10                                      
C                                  UPDATE COLUMN LENGTHS AND FIND LMAX  
         DLOSSJ = 1.0D0                                                 
         IF (BB.EQ.0.0D0) GO TO 30                                      
         TMP = BB                                                       
         BB = BB*DSQRT(DMAX1(1.0D0-(B(J-1)/BB)**2,0.0D0))               
         IF (BB.EQ.0.0D0) GO TO 30                                      
         DLOSSJ = BB/TMP                                                
         DO 5 L=J,N                                                     
            IF (H(L).EQ.0.0D0) GO TO 5                                  
            TMP = H(L)                                                  
            H(L) = H(L)*DSQRT(DMAX1(1.0D0-(A(J-1,L)/H(L))**2,0.0D0))    
            DLOSSJ = DMIN1(DLOSSJ,H(L)/TMP)                             
            TMP = X(L)                                                  
            X(L) = 0.0D0                                                
            IF (H(L).EQ.0.0D0) GO TO 5                                  
            X(L) = TMP-A(J-1,L)*B(J-1)                                  
            IF (H(LMAX).EQ.0.0D0) LMAX = L                              
            IF (DABS(X(L))/H(L).GT.DABS(X(LMAX))/H(LMAX)) LMAX = L      
    5    CONTINUE                                                       
         DLOSS = DLOSS*DLOSSJ                                           
         TMP = 10.0D0+DLOSS                                             
         IF (TMP.GT.10.0D0) GO TO 20                                    
C                                  COMPUTE COLUMN LENGTHS AND FIND LMAX 
   10    BB = DNRM2(M-J+1,B(J),1)                                       
         IF (BB.EQ.0.0D0) GO TO 30                                      
         DO 15 L=J,N                                                    
            H(L) = DNRM2(M-J+1,A(J,L),1)                                
            X(L) = 0.0D0                                                
            IF (H(L).EQ.0.0D0) GO TO 15                                 
            X(L) = DDOT(M-J+1,A(J,L),1,B(J),1)                          
            IF (H(LMAX).EQ.0.0D0) LMAX = L                              
            IF (DABS(X(L))/H(L).GT.DABS(X(LMAX))/H(LMAX)) LMAX = L      
   15    CONTINUE                                                       
         DLOSS = 1.0D0                                                  
C                                  LMAX HAS BEEN DETERMINED DO COLUMN   
C                                    INTERCHANGES IF NEEDED.            
   20    CONTINUE                                                       
         IP(J) = LMAX                                                   
         IF (LMAX.EQ.J) GO TO 30                                        
         DO 25 I=1,M                                                    
            TMP = A(I,J)                                                
            A(I,J) = A(I,LMAX)                                          
            A(I,LMAX) = TMP                                             
   25    CONTINUE                                                       
         H(LMAX) = H(J)                                                 
         X(LMAX) = X(J)                                                 
C                                  COMPUTE THE J-TH TRANSFORMATION AND  
C                                    APPLY IT TO A AND B.               
C                                                                       
   30    JCOL = MIN0(J+1,N)                                             
         CALL VHS12(1,J,J+1,M,A(1,J),1,H(J),A(1,JCOL),1,IA,N-J)         
         CALL VHS12(2,J,J+1,M,A(1,J),1,H(J),B,1,M,1)                    
   35 CONTINUE                                                          
C                                  DETERMINE THE NUMBER OF COLUMNS OF A 
C                                  TO BE INCLUDED IN THE BASIS SO THAT  
C                                  COND(AK) .LT. 1/TOL                  
C                                  AK = FIRST K COLUMNS OF A AFTER      
C                                       PIVOTING                        
C                                  RK = FIRST K COLUMNS OF R (Q*R = A)  
C                                  COND(AK) = NORM1(RK)*NORM1(RK**(-1)) 
      RCOND = 0.0D0                                                     
      K = 0                                                             
C                                  RNORM = NORM1(RK)                    
      RNORM = 0.0D0                                                     
C                                  XNORM = NORM1(RK**(-1))              
      XNORM = 0.0D0                                                     
      DO 55 J=1,LDIAG                                                   
         IF (DABS(A(J,J)).EQ.0.0D0) GO TO 60                            
         IF (TOL.LT.0.0D0) GO TO 50                                     
         RNORM = DMAX1(RNORM,DASUM(J,A(1,J),1))                         
         X(J) = 1.0D0/A(J,J)                                            
         IF (J.LT.2) GO TO 45                                           
         I = J                                                          
         DO 40 L=2,J                                                    
            I = I-1                                                     
            X(I) = -DDOT(J-I,X(I+1),1,A(I,I+1),IA)/A(I,I)               
   40    CONTINUE                                                       
   45    CONTINUE                                                       
         XNORM = DMAX1(XNORM,DASUM(J,X,1))                              
         RCONDJ = 1.0D0/(RNORM*XNORM)                                   
         IF (TOL.GE.RCONDJ) GO TO 60                                    
         RCOND = RCONDJ                                                 
   50    K = J                                                          
   55 CONTINUE                                                          
   60 KP1 = K+1                                                         
CJRC     KBASIS = K                                                        
      DO 65 J=1,N                                                       
   65 X(J) = 0.0D0                                                      
C                                  SPECIAL FOR KBASIS = 0               
      IF (KBASIS.EQ.0) GO TO 90                                         
C                                  SOLVE THE K BY K TRIANGULAR SYSTEM.  
      X(K) = B(K)/A(K,K)                                                
      IF (K.LT.2) GO TO 75                                              
      I = K                                                             
      DO 70 L=2,K                                                       
         I = I-1                                                        
         X(I) = (B(I)-DDOT(K-I,X(I+1),1,A(I,I+1),IA))/A(I,I)            
   70 CONTINUE                                                          
C                                  RE-ORDER THE SOLUTION VECTOR         
   75 J = LDIAG+1                                                       
      DO 80 JJ=1,LDIAG                                                  
         J = J-1                                                        
         L = IP(J)                                                      
         IF (L.EQ.J) GO TO 80                                           
         TMP = X(L)                                                     
         X(L) = X(J)                                                    
         X(J) = TMP                                                     
   80 CONTINUE                                                          
C                                  COMPUTE B - A*X                      
      DO 85 I=1,K                                                       
   85 B(I) = 0.0D0                                                      
   90 J = LDIAG+1                                                       
      DO 95 JJ=1,LDIAG                                                  
         J = J-1                                                        
         CALL VHS12(2,J,J+1,M,A(1,J),1,H(J),B,1,M,1)                    
   95 CONTINUE                                                          
      IF (TOL.GE.0.0D0) TOL = RCOND                                     
      GO TO 9005                                                        
 9000 CONTINUE                                                          
      WRITE (6,9002) IER
 9002 FORMAT(' IER = ',I4,' IN LLSQB ')
 9005 RETURN                                                            
      END                                                               
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      DOUBLE PRECISION A,B
      DOUBLE PRECISION P,R,S,T,U
      P = DMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DMIN1(DABS(A),DABS(B))/P)**2
   10 CONTINUE
         T = 4.0D0 + R
         IF (T .EQ. 4.0D0) GO TO 20
         S = R/T
         U = 1.0D0 + 2.0D0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
      INTEGER N,NM,IERR,MATZ
      DOUBLE PRECISION A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
   10 IF (MATZ .NE. 0) GO TO 20
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      DOUBLE PRECISION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
  100 E(I-1) = E(I)
      F = 0.0D0
      TST1 = 0.0D0
      E(N) = 0.0D0
      DO 240 L = 1, N
         J = 0
         H = DABS(D(L)) + DABS(E(L))
         IF (TST1 .LT. H) TST1 = H
         DO 110 M = L, N
            TST2 = TST1 + DABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
  110    CONTINUE
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
         R = PYTHAG(P,1.0D0)
         D(L) = E(L) / (P + DSIGN(R,P))
         D(L1) = E(L) * (P + DSIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
         DO 140 I = L2, N
  140    D(I) = D(I) - H
  145    F = F + H
         P = D(M)
         C = 1.0D0
         C2 = C
         EL1 = E(L1)
         S = 0.0D0
         MML = M - L
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            R = PYTHAG(P,E(I))
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
  200    CONTINUE
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + DABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
  300 CONTINUE
      GO TO 1001
 1000 IERR = L
 1001 RETURN
      END
        SUBROUTINE      TQLGRM	(N, D, E, Z, IERR)
        IMPLICIT        DOUBLE PRECISION (A-H, O-Z)
        DIMENSION       D(N), E(N), Z(N,N)
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
        RETURN
        END
      SUBROUTINE TQLRAT(N,D,E2,IERR)
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      DOUBLE PRECISION D(N),E2(N)
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
      F = 0.0D0
      T = 0.0D0
      E2(N) = 0.0D0
      DO 290 L = 1, N
         J = 0
         H = DABS(D(L)) + DSQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = EPSLON(T)
         C = B * B
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
  110    CONTINUE
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
         L1 = L + 1
         S = DSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * S)
         R = PYTHAG(P,1.0D0)
         D(L) = S / (P + DSIGN(R,P))
         H = G - D(L)
         DO 140 I = L1, N
  140    D(I) = D(I) - H
         F = F + H
         G = D(M)
         IF (G .EQ. 0.0D0) G = B
         H = G
         S = 0.0D0
         MML = M - L
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0D0) G = B
            H = G * P / R
  200    CONTINUE
         E2(L) = S * G
         D(L) = H
         IF (H .EQ. 0.0D0) GO TO 210
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0D0) GO TO 130
  210    P = D(L) + F
         IF (L .EQ. 1) GO TO 250
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
  250    I = 1
  270    D(I) = P
  290 CONTINUE
      GO TO 1001
 1000 IERR = L
 1001 RETURN
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

        CALL TREDIG	(N, E, W, H)
        CALL TQLGRM	(N, E, W, H, IERR)

        RETURN
        END
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
      INTEGER I,J,K,L,N,II,NM,JP1
      DOUBLE PRECISION A(NM,N),D(N),E(N),E2(N)
      DOUBLE PRECISION F,G,H,SCALE
      DO 100 I = 1, N
         D(I) = A(N,I)
         A(N,I) = A(I,I)
  100 CONTINUE
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(D(K))
         IF (SCALE .NE. 0.0D0) GO TO 140
         DO 125 J = 1, L
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = 0.0D0
  125    CONTINUE
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 300
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         IF (L .EQ. 1) GO TO 285
         DO 170 J = 1, L
  170    E(J) = 0.0D0
         DO 240 J = 1, L
            F = D(J)
            G = E(J) + A(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
            DO 200 K = JP1, L
               G = G + A(K,J) * D(K)
               E(K) = E(K) + A(K,J) * F
  200       CONTINUE
  220       E(J) = G
  240    CONTINUE
         F = 0.0D0
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
         H = F / (H + H)
         DO 250 J = 1, L
  250    E(J) = E(J) - H * D(J)
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
            DO 260 K = J, L
  260       A(K,J) = A(K,J) - F * E(K) - G * D(K)
  280    CONTINUE
  285    DO 290 J = 1, L
            F = D(J)
            D(J) = A(L,J)
            A(L,J) = A(I,J)
            A(I,J) = F * SCALE
  290    CONTINUE
  300 CONTINUE
      RETURN
      END
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
      INTEGER I,J,K,L,N,II,NM,JP1
      DOUBLE PRECISION A(NM,N),D(N),E(N),Z(NM,N)
      DOUBLE PRECISION F,G,H,HH,SCALE
      DO 100 I = 1, N
         DO 80 J = I, N
   80    Z(J,I) = A(J,I)
         D(I) = A(N,I)
  100 CONTINUE
      IF (N .EQ. 1) GO TO 510
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 2) GO TO 130
         DO 120 K = 1, L
  120    SCALE = SCALE + DABS(D(K))
         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = D(L)
         DO 135 J = 1, L
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
            Z(J,I) = 0.0D0
  135    CONTINUE
         GO TO 290
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
         F = D(L)
         G = -DSIGN(DSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         DO 170 J = 1, L
  170    E(J) = 0.0D0
         DO 240 J = 1, L
            F = D(J)
            Z(J,I) = F
            G = E(J) + Z(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
            DO 200 K = JP1, L
               G = G + Z(K,J) * D(K)
               E(K) = E(K) + Z(K,J) * F
  200       CONTINUE
  220       E(J) = G
  240    CONTINUE
         F = 0.0D0
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
         HH = F / (H + H)
         DO 250 J = 1, L
  250    E(J) = E(J) - HH * D(J)
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
            DO 260 K = J, L
  260       Z(K,J) = Z(K,J) - F * E(K) - G * D(K)
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
  280    CONTINUE
  290    D(I) = H
  300 CONTINUE
      DO 500 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0D0
         H = D(I)
         IF (H .EQ. 0.0D0) GO TO 380
         DO 330 K = 1, L
  330    D(K) = Z(K,I) / H
         DO 360 J = 1, L
            G = 0.0D0
            DO 340 K = 1, L
  340       G = G + Z(K,I) * Z(K,J)
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * D(K)
  360    CONTINUE
  380    DO 400 K = 1, L
  400    Z(K,I) = 0.0D0
  500 CONTINUE
  510 DO 520 I = 1, N
         D(I) = Z(N,I)
         Z(N,I) = 0.0D0
  520 CONTINUE
      Z(N,N) = 1.0D0
      E(1) = 0.0D0
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
      SUBROUTINE VHS12  (MODE,LP,L1,M,U,INCU,UP,C,INCC,ICV,NCV)         
C
C    VHSB IS QUITE SIMILAR TO VHS12 OF THE IMSL LIBRARY
C
      DOUBLE PRECISION   U(1),UP,C(1)                                   
      DOUBLE PRECISION   SM,B                                           
      DOUBLE PRECISION   ONE,CL,CLINV,SM1                               
      INTEGER            MODE,LP,L1,M,INCU,INCC,ICV,NCV                 
      INTEGER            IJ,ILP,IL1,IM,INCR,I2,I3,I4,J                  
C
      ONE = 1.D0                                                        
C                                                                       
      IF (0.GE.LP.OR.LP.GE.L1.OR.L1.GT.M) GO TO 9005                    
      ILP = (LP-1)*INCU+1                                               
      IL1 = (L1-1)*INCU+1                                               
      IM = (M-1)*INCU+1                                                 
      CL = DABS(U(ILP))                                                 
      IF (MODE.EQ.2) GO TO 15                                           
C                                  CONSTRUCT THE TRANSFORMATION.        
      DO 5 IJ=IL1,IM,INCU                                               
    5 CL = DMAX1(DABS(U(IJ)),CL)                                        
      IF (CL.LE.0.0D0) GO TO 9005                                       
      CLINV = ONE/CL                                                    
      SM = (U(ILP)*CLINV)**2                                            
      DO 10 IJ=IL1,IM,INCU                                              
   10 SM = SM+(U(IJ)*CLINV)**2                                          
C
      SM1 = SM                                                          
      CL = CL*DSQRT(SM1)                                                
      IF (U(ILP).GT.0.0D0) CL = -CL                                     
      UP = U(ILP)-CL                                                    
      U(ILP) = CL                                                       
      GO TO 20                                                          
C                                  APPLY THE TRANSFORMATION             
C                                    I+U*(U**T)/B TO C.                 
   15 IF (CL.LE.0.0D0) GO TO 9005                                       
   20 IF (NCV.LE.0) GO TO 9005                                          
      B = UP*U(ILP)                                                     
C
      IF (B.GE.0.0D0) GO TO 9005                                        
      B = ONE/B                                                         
      I2 = 1-ICV+INCC*(LP-1)                                            
      INCR = INCC*(L1-LP)                                               
      DO 35 J=1,NCV                                                     
         I2 = I2+ICV                                                    
         I3 = I2+INCR                                                   
         I4 = I3                                                        
         SM = C(I2)*UP                                                  
         DO 25 IJ=IL1,IM,INCU                                           
            SM = SM+C(I3)*U(IJ)                                         
            I3 = I3+INCC                                                
   25    CONTINUE                                                       
         IF (SM.EQ.0.0D0) GO TO 35                                      
         SM = SM*B                                                      
         C(I2) = C(I2)+SM*UP                                            
         DO 30 IJ=IL1,IM,INCU                                           
            C(I4) = C(I4)+SM*U(IJ)                                      
            I4 = I4+INCC                                                
   30    CONTINUE                                                       
   35 CONTINUE                                                          
 9005 RETURN                                                            
      END                                                               
      BLOCK DATA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Logical RHF,ROHF,ABNAT,RNAT
      Parameter (Mxbcrt=10,Mcent=50)
      COMMON /C2/ AB(9,9),AM(9,10)
      COMMON /UNITS/ ISRF,INPT,IOUT,IWFN,IDBG
      COMMON /PARA/ NINS(Mxbcrt),IBETPTS,INMAX,NNN,NPATH,
     $DIST,DMAX,AMAX,ADIFF,Tinf
      COMMON /LIMITS/RMAX, NATMX,  NFLIM,  NPROPS
      COMMON /DIRK/ TOLL,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,IMagm,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGo,YGo,ZGo
      COMMON /VALUES/ THRESH1,THRESH2,GAMMA,TOTE
      COMMON /NCUT/ CUTOFF,NDOCUT,NPR,NPRA,NPRB,NZEROA,NZEROB
      DATA INMAX,NNN,NPATH,DIST,DMAX,AMAX,ADIFF,Tinf
     + /6,140,80,0.6D0,8.0D0,1.D-12,0.0001D0,9.0d0/
      DATA ISRF,INPT,IOUT,IWFN,IDBG /2,5,6,10,13/
      DATA NATMX, NPROPS, NFLIM / 14, 43, 6 /
      DATA Stp/0.025d0/,Size/7.5d0/,FStp/0.25d0/,Ctf/0.001d0/,
     +SStp/0.025d0/,Nosec/1/,IStep/5/,NabMo/6/,Toll/0.2d0/,Nacc/1/,
     +IDOAOM/0/,IDOMAG/0/,IMagM/0/,nfbeta/0/,THRESH1/0.001d0/,
     +THRESH2/0.002d0/,NDOCUT/1/,CUTOFF/1.0d-9/,RHF/.false./,
     +ROHF/.False./,ABNAT/.false./,RNAT/.false./
C
C     AB AND AM CONTAIN THE COEFFICIENTS OF THE ADAMS-BASHFORTH-
C     MOULTON METHOD UP TO ORDER 10.
C
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
      SUBROUTINE CAGE(R1,C1,N1,QL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter (Mxbcrt=10,Maxpp=300)
      COMMON /C3/ WK(3,maxpp)
      COMMON /PARA/NINS(Mxbcrt),IBETPTS,INMAX,NNN,NPATH,DIST,
     $DMAX,AMAX,ADIFF,Tinf
      DIMENSION R1(3),C1(3,Mxbcrt)
      Data Four/4.0d0/
C
      NNN1 = NNN - 1
      ITEST = 0
C
C     CALCULATE THE APPROXIMATE LENGTH OF A PATH FROM A GIVEN
C     POINT TO A (3,+3) CRITICAL POINT.
C
      DO 100 I = 1,3
        WK(I,1) = R1(I)
100   CONTINUE
      QL0 = DSQRT((R1(1)-C1(1,N1))**2 + 
     +            (R1(2)-C1(2,N1))**2 +
     +            (R1(3)-C1(3,N1))**2)
      DS = DMIN1(DMAX/DFLOAT(NNN1),QL0/DMAX)
C
      NM = INT(QL0/DS)
10    CALL TRUDGE2(WK,DS,NM)
      NM1 = NM + 1
      QL = DSQRT((WK(1,NM1)-C1(1,N1))**2 +
     +           (WK(2,NM1)-C1(2,N1))**2 +
     +           (WK(3,NM1)-C1(3,N1))**2)
C
      IF (QL .LT. QL0) THEN
        QL = QL + DFLOAT(NM)*DS
        RETURN
      END IF
C
      IF (ITEST .GT. 0) THEN
        PRINT*,'QL.GE.QL0 IN SUBROUTINE CAGE'
        QL = QL0
        RETURN
      END IF
C
      ITEST = ITEST + 1
      DS = DMAX1(DS/Four,QL0/DFLOAT(NNN1))
      NM = INT(QL0/DS)
      GOTO 10
C
      END
      Subroutine CaseUp(InStr,OutStr,LenS)
      Implicit Integer(A-Z)
C
C     Translate a character string to upper case.
C
      Character*(*) InStr, OutStr
C
      IUpA = IChar('A')
      ILowA = IChar('a')
      ILowZ = IChar('z')
      Do  1 I = 1, LenS
        IX = IChar(InStr(I:I))
        If(IX.ge.ILowA.and.IX.le.ILowZ) IX = IX + IUpA - ILowA
        OutStr(I:I) = Char(IX)
1     Continue
      Return
      End
      SUBROUTINE EIGEN(C1,N)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter(MaxCrt=20)
      COMMON/EV/ V(3,MaxCrt),W(3,MaxCrt)
      DIMENSION R(3),H(3,3),WORK(3,3),C1(3,MaxCrt),FV1(3),FV2(3)
C
      DO 100 I = 1,N
        CALL HESS(C1(1,I),R,H)
        DO 20 II = 1,3
          DO 10 JJ = 1,3
            WORK(II,JJ) = H(II,JJ)
10        CONTINUE
20      CONTINUE
        CALL RS(3,3,WORK,R,1,H,FV1,FV2,IER)
C
        DO 100 J = 1,3
          V(J,I) = H(J,1)
          W(J,I) = H(J,2)
100   CONTINUE
C
      RETURN
      END
      SUBROUTINE FUNC3(R1,FR,IA,IBETA)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL RHF,ROHF,ABNAT,RNAT
      PARAMETER (MCENT=50, MMO=300, MPTS=500,MaxPrp=200)
      COMMON /VALUES/ THRESH1,THRESH2, GAMMA, TOTE
      COMMON /LIMITS/ RMAX, NATMX,  NFLIM,  NPROPS
      COMMON/ATOMS/XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),ROOTPO(MMO),NMO
       COMMON /ZZZZ/ PSI(MPTS,MMO),GX(MPTS,MMO),GY(MPTS,MMO),
     + GZ(MPTS,MMO),D2(MPTS,MMO)
      COMMON /DIST/ R2(MPTS,MCENT),DX(MPTS,MCENT),DY(MPTS,MCENT),
     +              DZ(MPTS,MCENT)
      COMMON /C7/  CT(3,3), X, Y, Z, NCRNT
      COMMON /DIRK/ TOLL,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,IMagm,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGo,YGo,ZGo
      DIMENSION R1(MPTS,3),FR(MPTS,MaxPrp),R(MPTS,3),W(MPTS,3),
     +VE(MPTS),RHO(MPTS),RMAX2(MPTS),CPX(MPTS),CPY(MPTS),CPZ(MPTS),
     +CDX(MPTS),CDY(MPTS),CDZ(MPTS),CX(MPTS),CY(MPTS),CZ(MPTS),
     +DDX(MPTS),DDY(MPTS),DDZ(MPTS),P(MPTS,MCENT)
      Save Zero,one,Two,Three,Four,Pt5,onept5,fine,DMill
      Data zero/0.0d0/,one/1.0d0/,Two/2.0d0/,Three/3.0d0/,
     $  Four/4.0d0/,Pt5/0.5d0/,Onept5/1.5d0/,fine/137.036d0/,
     $  DMill/1.0d6/
C
      LMO=NMO
      If(IDOMAG.eq.1)LMO=NMO/7

      DO 4 I=1,IA
      R(I,1)=R1(I,1)+X
      R(I,2)=R1(I,2)+Y
      R(I,3)=R1(I,3)+Z
      RMAX2(I)=DSQRT(R1(I,1)**2+R1(I,2)**2+R1(I,3)**2)
      RMAX = DMAX1(RMAX,RMAX2(I))
4     CONTINUE     
C
      CALL GAUS3(R,IA,IBETA)
C
      DO 95 J=1,IA
      FR(J,2)=Zero
      FR(J,3)=Zero
      VE(J)=Zero
      RHO(J)=Zero
      FR(J,37)=Zero
      FR(J,38)=Zero
      FR(J,42)=Zero
      FR(J,43)=Zero
      R(J,1)=Zero
      R(J,2)=Zero
      R(J,3)=Zero
      FR(J,29)=Zero
      FR(J,30)=Zero
      FR(J,31)=Zero
      FR(J,14)=Zero
      W(I,1)=Zero
      W(I,2)=Zero
      W(I,3)=Zero
      FR(J,26)=Zero
      FR(J,27)=Zero
      FR(J,28)=Zero
95    CONTINUE
C
      EN=Zero
C
      DO 100 I=1,LMO
      EN=EN+PO(I)
      DO 105 J=1,IA
      P1=PSI(J,I)*PSI(J,I)
C
      P2=PO(I)*PSI(J,I)
      W(J,1)=W(J,1)+P2*GX(J,I)
      W(J,2)=W(J,2)+P2*GY(J,I)
      W(J,3)=W(J,3)+P2*GZ(J,I)
      RHO(J)=RHO(J)+P1*PO(I)
      VE(J)=VE(J)+EORB(I)*P1*PO(I)
      FR(J,3)= FR(J,3) - PO(I)*(D2(J,I))*PSI(J,I)
      FR(J,2)= FR(J,2) + PO(I)*( GX(J,I)*GX(J,I)
     1           + GY(J,I)*GY(J,I)
     2           + GZ(J,I)*GZ(J,I))
      R(J,1)=R(J,1)-P2*Two*GX(J,I)
      R(J,2)=R(J,2)-P2*Two*GY(J,I)
      R(J,3)=R(J,3)-P2*Two*GZ(J,I)
105   CONTINUE
100   CONTINUE
C
       DO 180 I=1,NCENT
       DO 181 J=1,IA
       temp=charg(i)/(r2(j,i)**onept5)
       FR(J,29)=FR(J,29)+temp*DX(J,I)
       FR(J,30)=FR(J,30)+temp*DY(J,I)
       FR(J,31)=FR(J,31)+temp*DZ(J,I)
       FR(J,14)=FR(J,14)-CHARG(I)/DSQRT(R2(J,I))
 181   CONTINUE
 180   CONTINUE
C
      DO 121 J=1,IA 
      FR(J,1)=RHO(J)
      FR(J,2)= FR(J,2)*Pt5
      FR(J,3)= FR(J,3)*Pt5
      FR(J,4)=FR(J,3)-FR(J,2)
      FR(J,32)=FR(J,1)*FR(J,4)
      FR(J,33)=Four*(W(J,1)*W(J,1)+W(J,2)*W(J,2)+W(J,3)*W(J,3))
      FR(J,35)=RHO(J)*RMAX2(J)*RMAX2(J)*RMAX2(J)*RMAX2(J)
      FR(J,29)=RHO(J)*FR(J,29)
      FR(J,30)=RHO(J)*FR(J,30)
      FR(J,31)=RHO(J)*FR(J,31)
      FR(J,14)=RHO(J)*FR(J,14)
      FR(J,6)=RHO(J)/RMAX2(J)
      FR(J,7)=RHO(J)*RMAX2(J)
      FR(J,8)=RHO(J)*RMAX2(J)*RMAX2(J)
      FR(J,10)=-R1(J,1)*R(J,1)-R1(J,2)*R(J,2)-R1(J,3)*R(J,3)
      FR(J,9)=FR(J,10)/RMAX2(J)
      FR(J,11)=FR(J,10)*RMAX2(J)
      FR(J,12)=FR(J,10)*RMAX2(J)*RMAX2(J)
      FR(J,15)=Pt5*(FR(J,3)+FR(J,14)+VE(J))
      FR(J,16)=VE(J)-FR(J,15)
      FR(J,17)=-RHO(J)*R1(J,1)
      FR(J,18)=-RHO(J)*R1(J,2)
      FR(J,19)=-RHO(J)*R1(J,3)
      FR(J,20)=-RHO(J)*(Three*R1(J,1)*R1(J,1)-RMAX2(J)*RMAX2(J))
      FR(J,21)=-RHO(J)*THree*R1(J,1)*R1(J,2)
      FR(J,22)=-RHO(J)*Three*R1(J,1)*R1(J,3)
      FR(J,23)=-RHO(J)*(Three*R1(J,2)*R1(J,2)-RMAX2(J)*RMAX2(J))
      FR(J,24)=-RHO(J)*Three*R1(J,2)*R1(J,3)
      FR(J,25)=-RHO(J)*(Three*R1(J,3)*R1(J,3)-RMAX2(J)*RMAX2(J))
121   CONTINUE
C
      IF (NCRNT .EQ. -1) THEN
      CHG=Zero
      ELSE
      CHG=CHARG(NCRNT)
      ENDIF
C
      If(NCrnt.gt.0)Then
      DO 123 J=1,IA 
      FR(J,13)=-CHG*RHO(J)/DSQRT(R2(J,NCRNT))
      FR(J,26)=RHO(J)*CHG*R1(J,1)/R2(J,NCRNT)**OnePt5
      FR(J,27)=RHO(J)*CHG*R1(J,2)/R2(J,NCRNT)**OnePt5
      FR(J,28)=RHO(J)*CHG*R1(J,3)/R2(J,NCRNT)**OnePt5
123   Continue
      Endif
C
      DO 122 J=1,IA 
      temp=rho(j)/en
      IF (RHO(J).GE.THRESH1) FR(J,37)=One
      IF (RHO(J).GE.THRESH1) FR(J,42)=RHO(J)
      IF (RHO(J).GE.THRESH2) FR(J,38)=One
      IF (RHO(J).GE.THRESH2) FR(J,43)=RHO(J) 
      IF(temp.GT.Zero)FR(J,5)=-(temp)*DLOG(temp)
      IF(temp.LE.Zero)FR(J,5)=Zero
      FR(J,29)=FR(J,29)-FR(J,26)
      FR(J,30)=FR(J,30)-FR(J,27)
      FR(J,31)=FR(J,31)-FR(J,28)
122   CONTINUE
C
      IF(IDOMAG.EQ.1)Then
C
      Do 1555 J=1,IA
      R(J,1)=R1(J,1)+X
      R(J,2)=R1(J,2)+Y
      R(J,3)=R1(J,3)+Z
1555  Continue
C
      Const=-Dmill/(Fine**2)
C
      If(IMagM.eq.0)Then
C
      DO 1556 J=1,IA
      DDX(J)=X
      DDY(J)=Y
      DDZ(J)=Z
1556  Continue
C
      ElseIf(IMagM.eq.1)Then
C
      DO 1557 J=1,IA
      DDX(J)=XGo
      DDY(J)=YGo
      DDZ(J)=ZGo
1557  Continue
C
      ElseIf(Imagm.eq.2)Then
C
      Do 1558 I=1,NCent
      Do 1668 K=1,IA
      P(K,I)=One
1668  Continue
      Do 1559 J=1,NCent
      If(J.ne.I)Then
      Distij=dsqrt((xc(i)-xc(j))**2+(yc(i)-yc(j))**2+
     $             (zc(i)-zc(j))**2)
      Do 1600 K=1,IA
      Disti=dsqrt((R(K,1)-XC(I))**2+(R(K,2)-YC(I))**2+
     $            (R(K,3)-ZC(I))**2)
      Distj=dsqrt((R(K,1)-XC(J))**2+(R(K,2)-YC(J))**2+
     $            (R(K,3)-ZC(J))**2)
      dmu=(disti-distj)/distij
      f1=onept5*dmu-pt5*(dmu**3)
      f2=onept5*f1-pt5*(f1**3)
      f3=onept5*f2-pt5*(f2**3)
      ess=Pt5*(One-f3)
      P(K,I)=P(K,I)*ess
1600  Continue
      Endif
1559  Continue
1558  Continue
C       
      DO 1601 K=1,IA
      SMOVEX=Zero
      SMOVEY=Zero
      SMOVEZ=Zero
      Do 1602 I=1,NCent
      SMOVEX=SMOVEX+(R(K,1)-XC(I))*P(K,I)
      SMOVEY=SMOVEY+(R(K,2)-YC(I))*P(K,I)
      SMOVEZ=SMOVEZ+(R(K,3)-ZC(I))*P(K,I)
1602  Continue
      DDX(K)=R(K,1)-SMOVEX
      DDY(K)=R(K,2)-SMOVEY
      DDZ(K)=R(K,3)-SMOVEZ
1601  Continue
C
      Endif
C
        DO 300 II=1,3
C 
        DO 310 J=1,IA
        CPX(J)=Zero
        CPY(J)=Zero
        CPZ(J)=Zero
        CDX(J)=Zero
        CDY(J)=Zero
        CDZ(J)=Zero
310     CONTINUE
C
        DO 320 I = 1,LMO
        DO 322 J = 1,IA
      CPX(J) = CPX(J)-Pt5*PO(I)*(PSI(J,I)*GX(J,I+II*LMO)
     +-PSI(J,I+II*LMO)*GX(J,I))
      CPY(J) = CPY(J)-Pt5*PO(I)*(PSI(J,I)*GY(J,I+II*LMO)
     +-PSI(J,I+II*LMO)*GY(J,I))
      CPZ(J) = CPZ(J)-Pt5*PO(I)*(PSI(J,I)*GZ(J,I+II*LMO)
     +-PSI(J,I+II*LMO)*GZ(J,I))
322   CONTINUE
320   CONTINUE
C
      IF (II .EQ. 1) THEN
C
        DO 420 I = 1,LMO
        DO 422 J = 1,IA
      CPX(J)=CPX(J)+PO(I)*Pt5*(PSI(J,I)*DDY(J)*GX(J,I+LMO*6)-
     +PSI(J,I+LMO*6)*
     +         DDY(J)*GX(J,I)-PSI(J,I)*DDZ(J)*GX(J,I+LMO*5)+
     +         PSI(J,I+LMO*5)*DDZ(J)*GX(J,I))
      CPY(J)=CPY(J)+PO(I)*Pt5*(PSI(J,I)*DDY(J)*GY(J,I+LMO*6)-
     +PSI(J,I+LMO*6)*
     +         DDY(J)*GY(J,I)-PSI(J,I)*DDZ(J)*GY(J,I+LMO*5)+
     +         PSI(J,I+LMO*5)*DDZ(J)*GY(J,I))
      CPZ(J)=CPZ(J)+PO(I)*Pt5*(PSI(J,I)*DDY(J)*GZ(J,I+LMO*6)-
     +PSI(J,I+LMO*6)*
     +         DDY(J)*GZ(J,I)-PSI(J,I)*DDZ(J)*GZ(J,I+LMO*5)+
     +         PSI(J,I+LMO*5)*DDZ(J)*GZ(J,I))
422   CONTINUE
420   CONTINUE
C
      ELSE IF (II .EQ. 2) THEN
C
        DO 520 I=1,LMO
        DO 522 J=1,IA
      CPX(J)=CPX(J)+PO(I)*Pt5*(PSI(J,I)*DDZ(J)*GX(J,I+LMO*4)-
     +PSI(J,I+LMO*4)*
     +         DDZ(J)*GX(J,I)-PSI(J,I)*DDX(J)*GX(J,I+LMO*6)+
     +         PSI(J,I+LMO*6)*DDX(J)*GX(J,I))
      CPY(J)=CPY(J)+PO(I)*Pt5*(PSI(J,I)*DDZ(J)*GY(J,I+LMO*4)-
     +PSI(J,I+LMO*4)*
     +         DDZ(J)*GY(J,I)-PSI(J,I)*DDX(J)*GY(J,I+LMO*6)+
     +         PSI(J,I+LMO*6)*DDX(J)*GY(J,I))
      CPZ(J)=CPZ(J)+PO(I)*Pt5*(PSI(J,I)*DDZ(J)*GZ(J,I+LMO*4)-
     +PSI(J,I+LMO*4)*
     +         DDZ(J)*GZ(J,I)-PSI(J,I)*DDX(J)*GZ(J,I+LMO*6)+
     +         PSI(J,I+LMO*6)*DDX(J)*GZ(J,I))
522   CONTINUE
520   CONTINUE
C
      ELSE IF (II .EQ. 3) THEN
C
        DO 620 I=1,LMO
        DO 622 J=1,IA
      CPX(J)=CPX(J)+PO(I)*Pt5*(PSI(J,I)*DDX(J)*GX(J,I+LMO*5)-
     +PSI(J,I+LMO*5)*
     +         DDX(J)*GX(J,I)-PSI(J,I)*DDY(J)*GX(J,I+LMO*4)+
     +         PSI(J,I+LMO*4)*DDY(J)*GX(J,I))
      CPY(J)=CPY(J)+PO(I)*Pt5*(PSI(J,I)*DDX(J)*GY(J,I+LMO*5)-
     +PSI(J,I+LMO*5)*
     +         DDX(J)*GY(J,I)-PSI(J,I)*DDY(J)*GY(J,I+LMO*4)+
     +         PSI(J,I+LMO*4)*DDY(J)*GY(J,I))
      CPZ(J)=CPZ(J)+PO(I)*Pt5*(PSI(J,I)*DDX(J)*GZ(J,I+LMO*5)-
     +PSI(J,I+LMO*5)*
     +         DDX(J)*GZ(J,I)-PSI(J,I)*DDY(J)*GZ(J,I+LMO*4)+
     +         PSI(J,I+LMO*4)*DDY(J)*GZ(J,I))
622   CONTINUE
620   CONTINUE
C
        ENDIF
C
      IF (II .EQ. 1) THEN
	DO 2420 I=1,LMO
        DO 2422 J = 1,IA
      CDX(J)=Zero
      CDY(J)=CDY(J)+PO(I)*Pt5*(PSI(J,I)**2*(R(J,3)-DDZ(J)))
      CDZ(J)=CDZ(J)-PO(I)*Pt5*(PSI(J,I)**2*(R(J,2)-DDY(J)))
2422   CONTINUE
2420   CONTINUE
C
      ELSEIF (II .EQ. 2) THEN
	DO 2520 I=1,LMO
        DO 2522 J=1,IA
      CDX(J)=CDX(J)-PO(I)*Pt5*(PSI(J,I)**2*(R(J,3)-DDZ(J)))
      CDY(J)=Zero
      CDZ(J)=CDZ(J)+PO(I)*Pt5*(PSI(J,I)**2*(R(J,1)-DDX(J)))
 2522   CONTINUE
 2520   CONTINUE
C
       ELSEIF (II .EQ. 3) THEN
	 DO 2620 I=1,LMO
         DO 2622 J=1,IA
       CDX(J)=CDX(J)+PO(I)*Pt5*(PSI(J,I)**2*(R(J,2)-DDY(J)))
       CDY(J)=CDY(J)-PO(I)*Pt5*(PSI(J,I)**2*(R(J,1)-DDX(J)))
       CDZ(J)=Zero
 2622   CONTINUE
 2620   CONTINUE
C
      ENDIF 
C
      DO 400 J=1,IA
      CX(J)=CPX(J)+CDX(J)
      CY(J)=CPY(J)+CDY(J)
      CZ(J)=CPZ(J)+CDZ(J)
400   CONTINUE
C
      DO 323 J=1,IA
      FR(J,44+(II-1)*(6+NCENT*3)) = Pt5*(R1(J,2)
     +*CZ(J)-R1(J,3)*CY(J))
      FR(J,45+(II-1)*(6+NCent*3)) = Pt5*(R1(J,3)
     +*CX(J)-R1(J,1)*CZ(J))
      FR(J,46+(II-1)*(6+NCent*3)) = Pt5*(R1(J,1)
     +*CY(J)-R1(J,2)*CX(J))
C
      FR(J,47+(II-1)*(6+NCent*3)) = Pt5*CX(J)
      FR(J,48+(II-1)*(6+NCent*3)) = Pt5*CY(J)
      FR(J,49+(II-1)*(6+NCent*3)) = Pt5*CZ(J)
C
      DO 510 N=1,NCent
       XDD = R(J,1) - XC(N)
       YDD = R(J,2) - YC(N)
       ZDD = R(J,3) - ZC(N)
       DSQN=One/((DSQRT(XDD**2+YDD**2+ZDD**2))**3)
      FR(J,50+3*(N-1)+(II-1)*(6+NCent*3))=CONST*DSQN*
     +(YDD*CZ(J)-ZDD*CY(J))
      FR(J,51+3*(N-1)+(II-1)*(6+NCent*3))=CONST*DSQN*
     +(ZDD*CX(J)-XDD*CZ(J))
      FR(J,52+3*(N-1)+(II-1)*(6+NCent*3))=CONST*DSQN*
     +(XDD*CY(J)-YDD*CX(J))
510   CONTINUE
323   CONTINUE
300   CONTINUE
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE GAUS(XYZ)
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL RHF,ROHF,ABNAT,RNAT
      DOUBLE PRECISION NINE
      PARAMETER (MCENT=50, MMO=300, MPRIMS=700, NTYPE=20)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),ROOTPO(MMO),NMO
      COMMON /PRIMSA/ COOAMX(MMO),SUMA(MPRIMS),DIVA(MPRIMS),
     +COOA(MPRIMS,MMO),EXXA(MPRIMS),ICTA(MPRIMS),ITPA(NTYPE),NPRIMSA
      COMMON /ZZZZZ/ PSI(MMO),GX(MMO),GY(MMO),GZ(MMO),GXX(MMO),
     +GXY(MMO),GXZ(MMO),GYY(MMO),GYZ(MMO),GZZ(MMO)
      COMMON /DIRK/ TOLL,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,Imagm,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGO,YGO,ZGO
      DIMENSION XYZ(3),DX(MCENT),DY(MCENT),DZ(MCENT),R2(MCENT),
     $CHI(MPRIMS),CHIX(MPRIMS),CHIY(MPRIMS),CHIXX(MPRIMS),
     $CHIXY(MPRIMS),CHIXZ(MPRIMS),CHIYY(MPRIMS),CHIYZ(MPRIMS),
     $CHIZZ(MPRIMS),CHIZ(MPRIMS)
      Save Zero,one,two,three,four,five,six,seven,nine
      DATA ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,FOUR/4.0D0/,FIVE/5.0D0/,
     $SEVEN/7.0D0/,THREE/3.0D0/,SIX/6.0D0/,NINE/9.0D0/
C
      lmo=nmo
      If(Idomag.eq.1) lmo=nmo/7
      DO 110 J = 1,NCENT
        DX(J) = XYZ(1) - XC(J)
        DY(J) = XYZ(2) - YC(J)
        DZ(J) = XYZ(3) - ZC(J)
        R2(J)= DX(J)*DX(J)+DY(J)*DY(J)+DZ(J)*DZ(J)
110   CONTINUE
C
C       FOR S-TYPE
C
        DO 120 J = 1,ITPA(1) 
C
        A=SUMA(J)*DEXP(-EXXA(J)*R2(ICTA(J)))
C
        CHI(J)=A*DIVA(J)
        CHIX(J)=DX(ICTA(J))*A
        CHIY(J)=DY(ICTA(J))*A
        CHIZ(J)=DZ(ICTA(J))*A
        CHIXX(J)=(ONE+SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*A
        CHIYY(J)=(ONE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*A
        CHIZZ(J)=(ONE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*A
        CHIXY(J)=SUMA(J)*DX(ICTA(J))*DY(ICTA(J))*A
        CHIXZ(J)=SUMA(J)*DX(ICTA(J))*DZ(ICTA(J))*A
        CHIYZ(J)=SUMA(J)*DY(ICTA(J))*DZ(ICTA(J))*A
120     CONTINUE
C
C       FOR Px-TYPE
C
        DO 140 J=ITPA(1)+1,ITPA(2)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DX(ICTA(J))*A*SUMA(J)
C
        CHI(J)=A*DX(ICTA(J))
        CHIX(J)=A+DX(ICTA(J))*B
        CHIY(J)=DY(ICTA(J))*B
        CHIZ(J)=DZ(ICTA(J))*B
        CHIXX(J)=(THREE+SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*B
        CHIYY(J)=(ONE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*B
        CHIZZ(J)=(ONE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*B
        CHIXY(J)=(A+DX(ICTA(J))*B)*DY(ICTA(J))*SUMA(J)
        CHIXZ(J)=(A+DX(ICTA(J))*B)*DZ(ICTA(J))*SUMA(J)
        CHIYZ(J)=SUMA(J)*DY(ICTA(J))*DZ(ICTA(J))*B
140     CONTINUE
C
C       FOR Py-TYPE
C
        DO 160 J=ITPA(2)+1,ITPA(3)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DY(ICTA(J))*A*SUMA(J)
C
        CHI(J)=A*DY(ICTA(J))
        CHIX(J)=DX(ICTA(J))*B
        CHIY(J)=A+DY(ICTA(J))*B
        CHIZ(J)=DZ(ICTA(J))*B
        CHIXX(J)=(ONE+SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*B
        CHIYY(J)=(THREE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*B
        CHIZZ(J)=(ONE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*B
        CHIXY(J)=(A+DY(ICTA(J))*B)*DX(ICTA(J))*SUMA(J)
        CHIXZ(J)=SUMA(J)*DX(ICTA(J))*DZ(ICTA(J))*B
        CHIYZ(J)=(A+DY(ICTA(J))*B)*DZ(ICTA(J))*SUMA(J)
160     CONTINUE
C
C       FOR Pz-TYPE
C
        DO 180 J=ITPA(3)+1,ITPA(4)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DZ(ICTA(J))*A*SUMA(J)
C
        CHI(J)=A*DZ(ICTA(J))
        CHIX(J)=DX(ICTA(J))*B
        CHIY(J)=DY(ICTA(J))*B
        CHIZ(J)=A+DZ(ICTA(J))*B
        CHIXX(J)=(ONE+SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*B
        CHIYY(J)=(ONE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*B
        CHIZZ(J)=(THREE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*B
        CHIXY(J)=SUMA(J)*DX(ICTA(J))*DY(ICTA(J))*B
        CHIXZ(J)=(A+DZ(ICTA(J))*B)*DX(ICTA(J))*SUMA(J)
        CHIYZ(J)=(A+DZ(ICTA(J))*B)*DY(ICTA(J))*SUMA(J)
180     CONTINUE
C
C       FOR Dxx-TYPE
C
        DO 220 J=ITPA(4)+1,ITPA(5)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DX(ICTA(J))*DX(ICTA(J))*A*SUMA(J)
C
        CHI(J)=B*DIVA(J)
        CHIX(J)=(TWO*A+B)*DX(ICTA(J))
        CHIY(J)=DY(ICTA(J))*B
        CHIZ(J)=DZ(ICTA(J))*B
        CHIXX(J)=TWO*A+(FIVE+SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*B
        CHIYY(J)=(ONE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*B
        CHIZZ(J)=(ONE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*B
        CHIXY(J)=(TWO*A+B)*SUMA(J)*DX(ICTA(J))*DY(ICTA(J))
        CHIXZ(J)=(TWO*A+B)*SUMA(J)*DX(ICTA(J))*DZ(ICTA(J))
        CHIYZ(J)=SUMA(J)*DY(ICTA(J))*DZ(ICTA(J))*B
220     CONTINUE
C
C       FOR Dyy-TYPE
C
        DO 240 J=ITPA(5)+1,ITPA(6)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DY(ICTA(J))*DY(ICTA(J))*A*SUMA(J)
C
        CHI(J)=B*DIVA(J)
        CHIX(J)=DX(ICTA(J))*B
        CHIY(J)=(TWO*A+B)*DY(ICTA(J))
        CHIZ(J)=DZ(ICTA(J))*B
        CHIXX(J)=(ONE+SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*B
        CHIYY(J)=TWO*A+(FIVE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*B
        CHIZZ(J)=(ONE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*B
        CHIXY(J)=(TWO*A+B)*SUMA(J)*DX(ICTA(J))*DY(ICTA(J))
        CHIXZ(J)=SUMA(J)*DX(ICTA(J))*DZ(ICTA(J))*B
        CHIYZ(J)=(TWO*A+B)*SUMA(J)*DZ(ICTA(J))*DY(ICTA(J))
240     CONTINUE
C
C       FOR Dzz-TYPE
C
        DO 260 J=ITPA(6)+1,ITPA(7)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DZ(ICTA(J))*DZ(ICTA(J))*A*SUMA(J)
C
        CHI(J)=B*DIVA(J)
        CHIX(J)=DX(ICTA(J))*B
        CHIY(J)=DY(ICTA(J))*B
        CHIZ(J)=(TWO*A+B)*DZ(ICTA(J))
        CHIXX(J)=(ONE+SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*B
        CHIYY(J)=(ONE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*B
        CHIZZ(J)=TWO*A+(FIVE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*B
        CHIXY(J)=SUMA(J)*DX(ICTA(J))*DY(ICTA(J))*B
        CHIXZ(J)=(TWO*A+B)*SUMA(J)*DX(ICTA(J))*DZ(ICTA(J))
        CHIYZ(J)=(TWO*A+B)*SUMA(J)*DZ(ICTA(J))*DY(ICTA(J))
260     CONTINUE
C
C       FOR Dxy-TYPE
C
        DO 280 J=ITPA(7)+1,ITPA(8)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=SUMA(J)*DX(ICTA(J))*DY(ICTA(J))*A
C
        CHI(J)=B*DIVA(J)
        CHIX(J)=DX(ICTA(J))*B+DY(ICTA(J))*A
        CHIY(J)=DY(ICTA(J))*B+DX(ICTA(J))*A
        CHIZ(J)=DZ(ICTA(J))*B
        CHIXX(J)=(THREE+SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*B
        CHIYY(J)=(THREE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*B
        CHIZZ(J)=(ONE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*B
        CHIXY(J)=(ONE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*(A+
     +           A*SUMA(J)*DX(ICTA(J))*DX(ICTA(J))) 
        CHIXZ(J)=(ONE+SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*A*
     +           DY(ICTA(J))*DZ(ICTA(J))*SUMA(J)
        CHIYZ(J)=(ONE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*A*
     +           DX(ICTA(J))*DZ(ICTA(J))*SUMA(J)
280     CONTINUE
C
C       FOR Dxz-TYPE
C
        DO 320 J=ITPA(8)+1,ITPA(9)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DX(ICTA(J))*DZ(ICTA(J))*A*SUMA(J)
C
        CHI(J)=B*DIVA(J)
        CHIX(J)=DX(ICTA(J))*B+DZ(ICTA(J))*A
        CHIY(J)=DY(ICTA(J))*B
        CHIZ(J)=DZ(ICTA(J))*B+DX(ICTA(J))*A
        CHIXX(J)=(THREE+SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*B
        CHIYY(J)=(ONE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*B
        CHIZZ(J)=(THREE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*B
        CHIXZ(J)=(ONE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*(A+
     +           A*SUMA(J)*DX(ICTA(J))*DX(ICTA(J))) 
        CHIXY(J)=(ONE+SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*A*
     +           DZ(ICTA(J))*DY(ICTA(J))*SUMA(J)
        CHIYZ(J)=(ONE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*A*
     +           DX(ICTA(J))*DY(ICTA(J))*SUMA(J)
320     CONTINUE
C
C       FOR Dyz-TYPE
C
        DO 340 J=ITPA(9)+1,ITPA(10)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DY(ICTA(J))*DZ(ICTA(J))*A*SUMA(J)
C
        CHI(J)=B*DIVA(J)
        CHIX(J)=DX(ICTA(J))*B
        CHIY(J)=DY(ICTA(J))*B+DZ(ICTA(J))*A
        CHIZ(J)=DZ(ICTA(J))*B+DY(ICTA(J))*A
        CHIXX(J)=(ONE+SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*B
        CHIYY(J)=(THREE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*B
        CHIZZ(J)=(THREE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*B
        CHIYZ(J)=(ONE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*(A+
     +           A*SUMA(J)*DY(ICTA(J))*DY(ICTA(J))) 
        CHIXY(J)=(ONE+SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*A*
     +           DZ(ICTA(J))*DX(ICTA(J))*SUMA(J)
        CHIXZ(J)=(ONE+SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*A*
     +           DX(ICTA(J))*DY(ICTA(J))*SUMA(J)
340     CONTINUE
C
C
C       FOR Fxxx-TYPE
C
        DO 501 J=ITPA(10)+1,ITPA(11)
 
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        XSQ=DX(ICTA(J))*DX(ICTA(J))
        YSQ=DY(ICTA(J))*DY(ICTA(J))
        ZSQ=DZ(ICTA(J))*DZ(ICTA(J))
        B=XSQ*A
C
        CHI(J)=B*DX(ICTA(J))
        CHIX(J)=(THREE + SUMA(J)*XSQ)*B
        CHIY(J)=SUMA(J)*DY(ICTA(J))*CHI(J)
        CHIZ(J)=SUMA(J)*DZ(ICTA(J))*CHI(J)
C         CHIXX(J)=DX(ICTA(J))*(SIX + SEVEN*SUMA(J)*XSQ + 
C     1           SUMA(J)*SUMA(J)*XSQ*XSQ)*A
        CHIXX(J)=SIX*DX(ICTA(J))*A + SUMA(J)*(SEVEN+SUMA(J)*XSQ)*CHI(J)
        CHIYY(J)=SUMA(J)*CHI(J)*(ONE + SUMA(J)*YSQ)
        CHIZZ(J)=SUMA(J)*CHI(J)*(ONE + SUMA(J)*ZSQ)
        CHIXY(J)=SUMA(J)*DY(ICTA(J))*(THREE + SUMA(J)*XSQ)*B
        CHIXZ(J)=SUMA(J)*DZ(ICTA(J))*(THREE + SUMA(J)*XSQ)*B
        CHIYZ(J)=SUMA(J)*DY(ICTA(J))*CHIZ(J)
C
 501    CONTINUE
C
C       FOR Fyyy-TYPE
C
        DO 511 J=ITPA(11)+1,ITPA(12)
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        XSQ=DX(ICTA(J))*DX(ICTA(J))
        YSQ=DY(ICTA(J))*DY(ICTA(J))
        ZSQ=DZ(ICTA(J))*DZ(ICTA(J))
        B=YSQ*A
C
        CHI(J)=B*DY(ICTA(J))
        CHIX(J)=SUMA(J)*DX(ICTA(J))*CHI(J)
        CHIY(J)=(THREE + SUMA(J)*YSQ)*B
        CHIZ(J)=SUMA(J)*DZ(ICTA(J))*CHI(J)
        CHIXX(J)=SUMA(J)*CHI(J)*(ONE + SUMA(J)*XSQ)
        CHIYY(J)=SIX*DY(ICTA(J))*A + SUMA(J)*(SEVEN+SUMA(J)*YSQ)*CHI(J)
        CHIZZ(J)=SUMA(J)*CHI(J)*(ONE + SUMA(J)*ZSQ)
        CHIXY(J)=SUMA(J)*DX(ICTA(J))*(THREE + SUMA(J)*YSQ)*B
        CHIXZ(J)=SUMA(J)*DX(ICTA(J))*CHIZ(J)
        CHIYZ(J)=SUMA(J)*DZ(ICTA(J))*(THREE + SUMA(J)*YSQ)*B
C
 511    CONTINUE
C
C       FOR Fzzz-TYPE
C
        DO 521 J=ITPA(12)+1,ITPA(13)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        XSQ=DX(ICTA(J))*DX(ICTA(J))
        YSQ=DY(ICTA(J))*DY(ICTA(J))
        ZSQ=DZ(ICTA(J))*DZ(ICTA(J))
        B=ZSQ*A
C
        CHI(J)=B*DZ(ICTA(J))
        CHIX(J)=SUMA(J)*DX(ICTA(J))*CHI(J)
        CHIY(J)=SUMA(J)*DY(ICTA(J))*CHI(J)
        CHIZ(J)=(THREE + SUMA(J)*ZSQ)*B
        CHIXX(J)=SUMA(J)*CHI(J)*(ONE + SUMA(J)*XSQ)
        CHIYY(J)=SUMA(J)*CHI(J)*(ONE + SUMA(J)*YSQ)
        CHIZZ(J)=SIX*DZ(ICTA(J))*A + SUMA(J)*(SEVEN+SUMA(J)*ZSQ)*CHI(J)
        CHIXY(J)=SUMA(J)*DX(ICTA(J))*CHIY(J)
        CHIXZ(J)=SUMA(J)*DX(ICTA(J))*(THREE + SUMA(J)*ZSQ)*B
        CHIYZ(J)=SUMA(J)*DY(ICTA(J))*(THREE + SUMA(J)*ZSQ)*B
C
 521    CONTINUE
C
C       FOR Fxxy-TYPE
C
        DO 531 J=ITPA(13)+1,ITPA(14)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        XY=DX(ICTA(J))*DY(ICTA(J))
        XZ=DX(ICTA(J))*DZ(ICTA(J))
        XSQ=DX(ICTA(J))*DX(ICTA(J))
        YSQ=DY(ICTA(J))*DY(ICTA(J))
        ZSQ=DZ(ICTA(J))*DZ(ICTA(J))
C
        CHI(J)=XY*DX(ICTA(J))*A
        CHIX(J)=(TWO + SUMA(J)*XSQ)*XY*A
        CHIY(J)=(ONE + SUMA(J)*YSQ)*XSQ*A
        CHIZ(J)=SUMA(J)*DZ(ICTA(J))*CHI(J)
        CHIXX(J)=DY(ICTA(J))*A*(TWO + SUMA(J)*XSQ*(FIVE+XSQ*SUMA(J)))
        CHIYY(J)=SUMA(J)*(THREE + SUMA(J)*YSQ)*CHI(J)
        CHIZZ(J)=SUMA(J)*(ONE + SUMA(J)*ZSQ)*CHI(J)
        CHIXY(J)=DX(ICTA(J))*(TWO+SUMA(J)*(XSQ+YSQ*(TWO+SUMA(J)*XSQ)))*A
        CHIXZ(J)=SUMA(J)*XY*DZ(ICTA(J))*(TWO + SUMA(J)*XSQ)*A
        CHIYZ(J)=SUMA(J)*DZ(ICTA(J))*XSQ*(ONE+ SUMA(J)*YSQ)*A
C
 531    CONTINUE
C
C       FOR Fxxz-TYPE
C
        DO 541 J=ITPA(14)+1,ITPA(15)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        XY=DX(ICTA(J))*DY(ICTA(J))
        XZ=DX(ICTA(J))*DZ(ICTA(J))
        XSQ=DX(ICTA(J))*DX(ICTA(J))
        YSQ=DY(ICTA(J))*DY(ICTA(J))
        ZSQ=DZ(ICTA(J))*DZ(ICTA(J))
C
        CHI(J)=XZ*DX(ICTA(J))*A
        CHIX(J)=(TWO + SUMA(J)*XSQ)*XZ*A
        CHIY(J)=SUMA(J)*DY(ICTA(J))*CHI(J)
        CHIZ(J)=(ONE + SUMA(J)*ZSQ)*XSQ*A
        CHIXX(J)=DZ(ICTA(J))*A*(TWO + SUMA(J)*XSQ*(FIVE+XSQ*SUMA(J)))
        CHIYY(J)=SUMA(J)*(ONE + SUMA(J)*YSQ)*CHI(J)
        CHIZZ(J)=SUMA(J)*(THREE + SUMA(J)*ZSQ)*CHI(J)
        CHIXY(J)=SUMA(J)*XZ*DY(ICTA(J))*(TWO + SUMA(J)*XSQ)*A
        CHIXZ(J)=DX(ICTA(J))*(TWO+SUMA(J)*(XSQ+ZSQ*(TWO+SUMA(J)*XSQ)))*A
        CHIYZ(J)=SUMA(J)*DY(ICTA(J))*XSQ*(ONE+ SUMA(J)*ZSQ)*A
C
 541    CONTINUE
C
C       FOR Fyyz-TYPE
C
        DO 561 J=ITPA(15)+1,ITPA(16)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        XY=DX(ICTA(J))*DY(ICTA(J))
        YZ=DZ(ICTA(J))*DY(ICTA(J))
        XSQ=DX(ICTA(J))*DX(ICTA(J))
        YSQ=DY(ICTA(J))*DY(ICTA(J))
        ZSQ=DZ(ICTA(J))*DZ(ICTA(J))
C
        CHI(J)=YZ*DY(ICTA(J))*A
        CHIX(J)=SUMA(J)*DX(ICTA(J))*CHI(J)
        CHIY(J)=(TWO + SUMA(J)*YSQ)*YZ*A
        CHIZ(J)=(ONE + SUMA(J)*ZSQ)*YSQ*A
        CHIXX(J)=SUMA(J)*(ONE + SUMA(J)*XSQ)*CHI(J)
        CHIYY(J)=DZ(ICTA(J))*A*(TWO + SUMA(J)*YSQ*(FIVE+YSQ*SUMA(J)))
        CHIZZ(J)=SUMA(J)*(THREE + SUMA(J)*ZSQ)*CHI(J)
        CHIXY(J)=SUMA(J)*YZ*DX(ICTA(J))*(TWO + SUMA(J)*YSQ)*A
        CHIXZ(J)=SUMA(J)*DX(ICTA(J))*YSQ*(ONE+ SUMA(J)*ZSQ)*A
        CHIYZ(J)=DY(ICTA(J))*(TWO+SUMA(J)*(YSQ+ZSQ*(TWO+SUMA(J)*YSQ)))*A
C
 561    CONTINUE
C
C       FOR Fxyy-TYPE
C
        DO 551 J=ITPA(16)+1,ITPA(17)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        XY=DX(ICTA(J))*DY(ICTA(J))
        YZ=DY(ICTA(J))*DZ(ICTA(J))
        XSQ=DX(ICTA(J))*DX(ICTA(J))
        YSQ=DY(ICTA(J))*DY(ICTA(J))
        ZSQ=DZ(ICTA(J))*DZ(ICTA(J))
C
        CHI(J)=XY*DY(ICTA(J))*A
        CHIX(J)=(ONE + SUMA(J)*XSQ)*YSQ*A
        CHIY(J)=(TWO + SUMA(J)*YSQ)*XY*A
        CHIZ(J)=SUMA(J)*DZ(ICTA(J))*CHI(J)
        CHIXX(J)=SUMA(J)*(THREE + SUMA(J)*XSQ)*CHI(J)
        CHIYY(J)=DX(ICTA(J))*A*(TWO + SUMA(J)*YSQ*(FIVE+YSQ*SUMA(J)))
        CHIZZ(J)=SUMA(J)*(ONE + SUMA(J)*ZSQ)*CHI(J)
        CHIXY(J)=DY(ICTA(J))*(TWO+SUMA(J)*(YSQ+XSQ*(TWO+SUMA(J)*YSQ)))*A
        CHIYZ(J)=SUMA(J)*XY*DZ(ICTA(J))*(TWO + SUMA(J)*YSQ)*A
        CHIXZ(J)=SUMA(J)*DZ(ICTA(J))*YSQ*(ONE+ SUMA(J)*XSQ)*A
C
 551    CONTINUE
C
C       FOR Fxzz-TYPE
C
        DO 571 J=ITPA(17)+1,ITPA(18)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
C       XY=DY(ICTA(J))*DX(ICTA(J))
        XZ=DZ(ICTA(J))*DX(ICTA(J))
        YZ=DZ(ICTA(J))*DY(ICTA(J))
        XSQ=DX(ICTA(J))*DX(ICTA(J))
        YSQ=DY(ICTA(J))*DY(ICTA(J))
        ZSQ=DZ(ICTA(J))*DZ(ICTA(J))
C
        CHI(J)=XZ*DZ(ICTA(J))*A
        CHIX(J)=(ONE + SUMA(J)*XSQ)*ZSQ*A
        CHIY(J)=SUMA(J)*DY(ICTA(J))*CHI(J)
        CHIZ(J)=(TWO + SUMA(J)*ZSQ)*XZ*A
        CHIXX(J)=SUMA(J)*(THREE + SUMA(J)*XSQ)*CHI(J)
        CHIYY(J)=SUMA(J)*(ONE + SUMA(J)*YSQ)*CHI(J)
        CHIZZ(J)=DX(ICTA(J))*A*(TWO + SUMA(J)*ZSQ*(FIVE+ZSQ*SUMA(J)))
        CHIXY(J)=SUMA(J)*DY(ICTA(J))*ZSQ*(ONE+ SUMA(J)*XSQ)*A
        CHIXZ(J)=DZ(ICTA(J))*(TWO+SUMA(J)*(ZSQ+XSQ*(TWO+SUMA(J)*ZSQ)))*A
        CHIYZ(J)=SUMA(J)*XZ*DY(ICTA(J))*(TWO + SUMA(J)*ZSQ)*A
C
 571    CONTINUE
C
C       FOR Fyzz-TYPE
C
        DO 581 J=ITPA(18)+1,ITPA(19)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        XZ=DX(ICTA(J))*DZ(ICTA(J))
        YZ=DZ(ICTA(J))*DY(ICTA(J))
        XSQ=DX(ICTA(J))*DX(ICTA(J))
        YSQ=DY(ICTA(J))*DY(ICTA(J))
        ZSQ=DZ(ICTA(J))*DZ(ICTA(J))
C
        CHI(J)=YZ*DZ(ICTA(J))*A
        CHIX(J)=SUMA(J)*DX(ICTA(J))*CHI(J)
        CHIY(J)=(ONE + SUMA(J)*YSQ)*ZSQ*A
        CHIZ(J)=(TWO + SUMA(J)*ZSQ)*YZ*A
        CHIXX(J)=SUMA(J)*(ONE + SUMA(J)*XSQ)*CHI(J)
        CHIYY(J)=SUMA(J)*(THREE + SUMA(J)*YSQ)*CHI(J)
        CHIZZ(J)=DY(ICTA(J))*A*(TWO + SUMA(J)*ZSQ*(FIVE+ZSQ*SUMA(J)))
        CHIXY(J)=SUMA(J)*DX(ICTA(J))*ZSQ*(ONE+ SUMA(J)*YSQ)*A
        CHIYZ(J)=DZ(ICTA(J))*(TWO+SUMA(J)*(ZSQ+YSQ*(TWO+SUMA(J)*ZSQ)))*A
        CHIXZ(J)=SUMA(J)*YZ*DX(ICTA(J))*(TWO + SUMA(J)*ZSQ)*A
C
 581    CONTINUE
C
C       FOR Fxyz-TYPE
C
        DO 591 J=ITPA(19)+1,ITPA(20)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        XY=DX(ICTA(J))*DY(ICTA(J))
        YZ=DZ(ICTA(J))*DY(ICTA(J))
        XZ=DX(ICTA(J))*DZ(ICTA(J))
        XSQ=DX(ICTA(J))*DX(ICTA(J))
        YSQ=DY(ICTA(J))*DY(ICTA(J))
        XSQ=DZ(ICTA(J))*DZ(ICTA(J))
C
        CHI(J)=DX(ICTA(J))*YZ*A
        CHIX(J)=(ONE + SUMA(J)*XSQ)*YZ*A
        CHIY(J)=(ONE + SUMA(J)*YSQ)*XZ*A
        CHIZ(J)=(ONE + SUMA(J)*ZSQ)*XY*A
        CHIXX(J)=(THREE + XSQ*SUMA(J))*SUMA(J)*CHI(J)
        CHIYY(J)=(THREE + YSQ*SUMA(J))*SUMA(J)*CHI(J)
        CHIZZ(J)=(THREE + ZSQ*SUMA(J))*SUMA(J)*CHI(J)
        CHIXY(J)=DZ(ICTA(J))*(ONE+SUMA(J)*(XSQ+YSQ+SUMA(J)*XY*XY))*A
        CHIXZ(J)=DY(ICTA(J))*(ONE+SUMA(J)*(XSQ+ZSQ+SUMA(J)*XZ*XZ))*A
        CHIYZ(J)=DX(ICTA(J))*(ONE+SUMA(J)*(YSQ+ZSQ+SUMA(J)*YZ*YZ))*A
C
 591    CONTINUE
C
        DO 124 L = 1,LMO
          PSI(L) =ZERO
          GX(L) =ZERO
          GY(L) = ZERO
          GZ(L) =ZERO
          GXX(L) = ZERO
          GXY(L) =ZERO
          GXZ(L) =ZERO
          GYY(L) =ZERO
          GYZ(L) =ZERO
          GZZ(L) =ZERO
        DO 125 J = 1,NPRIMSA
          PSI(L) = PSI(L) + COOA(J,L)*CHI(J)
          GX(L) = GX(L) + COOA(J,L)*CHIX(J)
          GY(L) = GY(L) + COOA(J,L)*CHIY(J)
          GZ(L) = GZ(L) + COOA(J,L)*CHIZ(J)
          GXX(L) = GXX(L) + COOA(J,L)*CHIXX(J)
          GXY(L) = GXY(L) + COOA(J,L)*CHIXY(J)
          GXZ(L) = GXZ(L) + COOA(J,L)*CHIXZ(J)
          GYY(L) = GYY(L) + COOA(J,L)*CHIYY(J)
          GYZ(L) = GYZ(L) + COOA(J,L)*CHIYZ(J)
          GZZ(L) = GZZ(L) + COOA(J,L)*CHIZZ(J)
125      CONTINUE
124      CONTINUE
C
        RETURN
        END
      SUBROUTINE GAUS2(XYZ,IW)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL RHF,ROHF,ABNAT,RNAT
      DOUBLE PRECISION NINE
      PARAMETER (MCENT=50, MMO=300, MPRIMS=700, NTYPE=20)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),ROOTPO(MMO),NMO
      COMMON /PRIMSA/ COOAMX(MMO),SUMA(MPRIMS),DIVA(MPRIMS),
     +COOA(MPRIMS,MMO),EXXA(MPRIMS),ICTA(MPRIMS),ITPA(NTYPE),NPRIMSA
      COMMON /C7/  CT(3,3), X, Y, Z, NCRNT
      COMMON /DIRK/ TOLL,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,Imagm,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGO,YGO,ZGO
      COMMON /NCUT/ CUTOFF,NDOCUT,NPR,NPRA,NPRB,NZEROA,NZEROB
      DIMENSION R2(MCENT),DX(MCENT),DY(MCENT),DZ(MCENT),CHI(MPRIMS),
     $CHIX(MPRIMS),CHIZ(MPRIMS),CHIY(MPRIMS),PSI(MMO),GX(MMO),
     $GY(MMO),GZ(MMO),XYZ(3)
      Save Zero,one,two,three,four,five,six,seven,nine
      DATA ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,FOUR/4.0D0/,FIVE/5.0D0/,
     $SEVEN/7.0D0/,THREE/3.0D0/,SIX/6.0D0/,NINE/9.0D0/
C
      LMO=NMO
      If(IdoMag.eq.1)LMO=NMO/7
C
      If(IW.EQ.0)Then
      XYZ(1) = XYZ(1) + X
      XYZ(2) = XYZ(2) + Y
      XYZ(3) = XYZ(3) + Z
      Endif
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
        DO 120 J = 1,ITPA(1) 
C
        A=SUMA(J)*DEXP(-EXXA(J)*R2(ICTA(J)))
C
        CHI(J)=A*DIVA(J)
        CHIX(J)=DX(ICTA(J))*A
        CHIY(J)=DY(ICTA(J))*A
        CHIZ(J)=DZ(ICTA(J))*A
120     CONTINUE
C
C       FOR Px-TYPE
C
        DO 140 J=ITPA(1)+1,ITPA(2)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DX(ICTA(J))*A*SUMA(J)
C
        CHI(J)=A*DX(ICTA(J))
        CHIX(J)=A+DX(ICTA(J))*B
        CHIY(J)=DY(ICTA(J))*B
        CHIZ(J)=DZ(ICTA(J))*B
140     CONTINUE
C
C       FOR Py-TYPE
C
        DO 160 J=ITPA(2)+1,ITPA(3)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DY(ICTA(J))*A*SUMA(J)
C
        CHI(J)=A*DY(ICTA(J))
        CHIX(J)=DX(ICTA(J))*B
        CHIY(J)=A+DY(ICTA(J))*B
        CHIZ(J)=DZ(ICTA(J))*B
160     CONTINUE
C
C       FOR Pz-TYPE
C
        DO 180 J=ITPA(3)+1,ITPA(4)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DZ(ICTA(J))*A*SUMA(J)
C
        CHI(J)=A*DZ(ICTA(J))
        CHIX(J)=DX(ICTA(J))*B
        CHIY(J)=DY(ICTA(J))*B
        CHIZ(J)=A+DZ(ICTA(J))*B
180     CONTINUE
C
C       FOR Dxx-TYPE
C
        DO 220 J=ITPA(4)+1,ITPA(5)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DX(ICTA(J))*DX(ICTA(J))*A*SUMA(J)
C
        CHI(J)=B*DIVA(J)
        CHIX(J)=(TWO*A+B)*DX(ICTA(J))
        CHIY(J)=DY(ICTA(J))*B
        CHIZ(J)=DZ(ICTA(J))*B
220     CONTINUE
C
C       FOR Dyy-TYPE
C
        DO 240 J=ITPA(5)+1,ITPA(6)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DY(ICTA(J))*DY(ICTA(J))*A*SUMA(J)
C
        CHI(J)=B*DIVA(J)
        CHIX(J)=DX(ICTA(J))*B
        CHIY(J)=(TWO*A+B)*DY(ICTA(J))
        CHIZ(J)=DZ(ICTA(J))*B
240     CONTINUE
C
C       FOR Dzz-TYPE
C
        DO 260 J=ITPA(6)+1,ITPA(7)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DZ(ICTA(J))*DZ(ICTA(J))*A*SUMA(J)
C
        CHI(J)=B*DIVA(J)
        CHIX(J)=DX(ICTA(J))*B
        CHIY(J)=DY(ICTA(J))*B
        CHIZ(J)=(TWO*A+B)*DZ(ICTA(J))
260     CONTINUE
C
C       FOR Dxy-TYPE
C
        DO 280 J=ITPA(7)+1,ITPA(8)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DX(ICTA(J))*DY(ICTA(J))*A*SUMA(J)
C
        CHI(J)=B*DIVA(J)
        CHIX(J)=DX(ICTA(J))*B+DY(ICTA(J))*A
        CHIY(J)=DY(ICTA(J))*B+DX(ICTA(J))*A
        CHIZ(J)=DZ(ICTA(J))*B
280     CONTINUE
C
C       FOR Dxz-TYPE
C
        DO 320 J=ITPA(8)+1,ITPA(9)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DX(ICTA(J))*DZ(ICTA(J))*A*SUMA(J)
C
        CHI(J)=B*DIVA(J)
        CHIX(J)=DX(ICTA(J))*B+DZ(ICTA(J))*A
        CHIY(J)=DY(ICTA(J))*B
        CHIZ(J)=DZ(ICTA(J))*B+DX(ICTA(J))*A
320     CONTINUE
C
C       FOR Dyz-TYPE
C
        DO 340 J=ITPA(9)+1,ITPA(10)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DY(ICTA(J))*DZ(ICTA(J))*A*SUMA(J)
C
        CHI(J)=B*DIVA(J)
        CHIX(J)=DX(ICTA(J))*B
        CHIY(J)=DY(ICTA(J))*B+DZ(ICTA(J))*A
        CHIZ(J)=DZ(ICTA(J))*B+DY(ICTA(J))*A
340     CONTINUE
C
C
C       FOR Fxxx-TYPE
C
        DO 501 J=ITPA(10)+1,ITPA(11)
 
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DX(ICTA(J))*DX(ICTA(J))*A
C
        CHI(J)=B*DX(ICTA(J))
        CHIX(J)=(THREE + SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*B
        CHIY(J)=SUMA(J)*DY(ICTA(J))*CHI(J)
        CHIZ(J)=SUMA(J)*DZ(ICTA(J))*CHI(J)
 501    CONTINUE
C
C       FOR Fyyy-TYPE
C
        DO 511 J=ITPA(11)+1,ITPA(12)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DY(ICTA(J))*DY(ICTA(J))*A
C
        CHI(J)=B*DY(ICTA(J))
        CHIX(J)=SUMA(J)*DX(ICTA(J))*CHI(J)
        CHIY(J)=(THREE + SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*B
        CHIZ(J)=SUMA(J)*DZ(ICTA(J))*CHI(J)
C
 511    CONTINUE
C
C       FOR Fzzz-TYPE
C
        DO 521 J=ITPA(12)+1,ITPA(13)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        B=DZ(ICTA(J))*DZ(ICTA(J))*A
C
        CHI(J)=B*DZ(ICTA(J))
        CHIX(J)=SUMA(J)*DX(ICTA(J))*CHI(J)
        CHIY(J)=SUMA(J)*DY(ICTA(J))*CHI(J)
        CHIZ(J)=(THREE + SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*B
C
 521    CONTINUE
C
C       FOR Fxxy-TYPE
C
        DO 531 J=ITPA(13)+1,ITPA(14)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        BXY=DX(ICTA(J))*DY(ICTA(J))*A
        BXX=DX(ICTA(J))*DX(ICTA(J))*A
C
        CHI(J)=BXY*DX(ICTA(J))
        CHIX(J)=(TWO + SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*BXY
        CHIY(J)=(ONE + SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*BXX
        CHIZ(J)=SUMA(J)*DZ(ICTA(J))*CHI(J)
C
 531    CONTINUE
C
C       FOR Fxxz-TYPE
C
        DO 541 J=ITPA(14)+1,ITPA(15)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        BXZ=DX(ICTA(J))*DZ(ICTA(J))*A
        BXX=DX(ICTA(J))*DX(ICTA(J))*A
C
        CHI(J)=BXZ*DX(ICTA(J))
        CHIX(J)=(TWO + SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*BXZ
        CHIY(J)=SUMA(J)*DY(ICTA(J))*CHI(J)
        CHIZ(J)=(ONE + SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*BXX
C
 541    CONTINUE
C
C       FOR Fyyz-TYPE
C
        DO 561 J=ITPA(15)+1,ITPA(16)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        BYZ=DZ(ICTA(J))*DY(ICTA(J))*A
        BYY=DY(ICTA(J))*DY(ICTA(J))*A
C
        CHI(J)=BYZ*DY(ICTA(J))
        CHIX(J)=SUMA(J)*DX(ICTA(J))*CHI(J)
        CHIY(J)=(TWO + SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*BYZ
        CHIZ(J)=(ONE + SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*BYY
C
 561    CONTINUE
C
C       FOR Fxyy-TYPE
C
        DO 551 J=ITPA(16)+1,ITPA(17)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        BXY=DX(ICTA(J))*DY(ICTA(J))*A
        BYY=DY(ICTA(J))*DY(ICTA(J))*A
C
        CHI(J)=BXY*DY(ICTA(J))
        CHIX(J)=(ONE + SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*BYY
        CHIY(J)=(TWO + SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*BXY
        CHIZ(J)=SUMA(J)*DZ(ICTA(J))*CHI(J)
C
 551    CONTINUE
C
C       FOR Fxzz-TYPE
C
        DO 571 J=ITPA(17)+1,ITPA(18)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        BXZ=DZ(ICTA(J))*DX(ICTA(J))*A
        BZZ=DZ(ICTA(J))*DZ(ICTA(J))*A
C
        CHI(J)=BXZ*DZ(ICTA(J))
        CHIX(J)=(ONE + SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*BZZ
        CHIY(J)=SUMA(J)*DY(ICTA(J))*CHI(J)
        CHIZ(J)=(TWO + SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*BXZ
C
 571    CONTINUE
C
C       FOR Fyzz-TYPE
C
        DO 581 J=ITPA(18)+1,ITPA(19)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        BYZ=DZ(ICTA(J))*DY(ICTA(J))*A
        BZZ=DZ(ICTA(J))*DZ(ICTA(J))*A
C
        CHI(J)=BYZ*DZ(ICTA(J))
        CHIX(J)=SUMA(J)*DX(ICTA(J))*CHI(J)
        CHIY(J)=(ONE + SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*BZZ
        CHIZ(J)=(TWO + SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*BYZ
C
 581    CONTINUE
C
C       FOR Fxyz-TYPE
C
        DO 591 J=ITPA(19)+1,ITPA(20)
C
        A=DEXP(-EXXA(J)*R2(ICTA(J)))
        BXY=DX(ICTA(J))*DY(ICTA(J))*A
        BYZ=DZ(ICTA(J))*DY(ICTA(J))*A
        BXZ=DX(ICTA(J))*DZ(ICTA(J))*A
C
        CHI(J)=DX(ICTA(J))*BYZ
        CHIX(J)=(ONE + SUMA(J)*DX(ICTA(J))*DX(ICTA(J)))*BYZ
        CHIY(J)=(ONE + SUMA(J)*DY(ICTA(J))*DY(ICTA(J)))*BXZ
        CHIZ(J)=(ONE + SUMA(J)*DZ(ICTA(J))*DZ(ICTA(J)))*BXY
C
 591    CONTINUE
C
       chimax=zero
       Do 592 J=1,NPRIMSA
       check=dmax1(dabs(chi(j)),dabs(chix(j)),dabs(chiy(j)),
     $       dabs(chiz(j)),chimax)
       If(check.gt.chimax)chimax=check
592    Continue
       DO 345 L=1,LMO
       PSI(L)=ZERO
       GX(L)=ZERO
       GY(L)=ZERO
       GZ(L)=ZERO
       If(dabs(cooamx(l)*rootpo(l)*chimax).gt.cutoff)Then
       DO 346 J=1,NPRIMSA
       PSI(L)=PSI(L)+COOA(J,L)*CHI(J)
       GX(L)=GX(L)+COOA(J,L)*CHIX(J)
       GY(L)=GY(L)+COOA(J,L)*CHIY(J)
       GZ(L)=GZ(L)+COOA(J,L)*CHIZ(J)
346    CONTINUE
       Endif
345    CONTINUE
C
       DO 243 I = 1,3
       XYZ(I) = Zero
243    CONTINUE
C
C      CALCULATE THE (+/-) GRADIENT VECTOR AT THE POINT (/2).
C
       ISIGN=-1
       IF(IW.EQ.1)ISIGN=1
       DO 250 I =1,LMO
       XYZ(1) = XYZ(1) + ISIGN*PO(I)*PSI(I)*GX(I)
       XYZ(2) = XYZ(2) + ISIGN*PO(I)*PSI(I)*GY(I)
       XYZ(3) = XYZ(3) + ISIGN*PO(I)*PSI(I)*GZ(I)
250    CONTINUE
C
        RETURN
        END
      SUBROUTINE GAUS3(XYZ,PTS,IBETA)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL RHF,ROHF,ABNAT,RNAT
      DOUBLE PRECISION NINE
      INTEGER PTS
      PARAMETER (MCENT=50, MMO=300, MPRIMS=700, MPTS=500, NTYPE=20)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),ROOTPO(MMO),NMO
      COMMON /PRIMSA/ COOAMX(MMO),SUMA(MPRIMS),DIVA(MPRIMS),
     +COOA(MPRIMS,MMO),EXXA(MPRIMS),ICTA(MPRIMS),ITPA(NTYPE),NPRIMSA
      COMMON /PRIMSB/ COOBMX(MMO),SUMB(MPRIMS),DIVB(MPRIMS),
     +COOB(MPRIMS,MMO),EXXB(MPRIMS),ICTB(MPRIMS),ITPB(NTYPE),NPRIMSB
      COMMON /ZZZZ/ PSI(MPTS,MMO),GX(MPTS,MMO),GY(MPTS,MMO),
     + GZ(MPTS,MMO),D2(MPTS,MMO)
      COMMON /DIST/ R2(MPTS,MCENT),DX(MPTS,MCENT),DY(MPTS,MCENT),
     +DZ(MPTS,MCENT)
      COMMON /DER/ DXR(MPTS),DYR(MPTS),DZR(MPTS)
      COMMON /NCUT/ CUTOFF,NDOCUT,NPR,NPRA,NPRB,NZEROA,NZEROB
      COMMON /DIRK/ TOLL,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,Imagm,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGO,YGO,ZGO
      DIMENSION CHI(MPTS,MPRIMS),CHIX(MPTS,MPRIMS),CHIY(MPTS,MPRIMS),
     $CHIZ(MPTS,MPRIMS),CHID2(MPTS,MPRIMS),XYZ(MPTS,3),CHIMAX(MPRIMS)
      Save Zero,one,two,three,four,five,six,seven,nine
      DATA ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,FOUR/4.0D0/,FIVE/5.0D0/,
     $SEVEN/7.0D0/,THREE/3.0D0/,SIX/6.0D0/,NINE/9.0D0/,Hund/100.0d0/
C
      Cutoft=cutoff*Hund
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
      IF(IBETA.EQ.1)GOTO 555
      IF(IBETA.EQ.2)GOTO 5555
C
C       FOR S-TYPE
C
        DO 1120 J = 1,ITPA(1) 
        is=icta(j)
        DO 1122 I=1,PTS 
C
        A=SUMA(J)*DEXP(-EXXA(J)*R2(I,is))
C
        CHI(I,J)=A*DIVA(J)
        CHIX(I,J)=DX(I,is)*A
        CHIY(I,J)=DY(I,is)*A
        CHIZ(I,J)=DZ(I,is)*A
        CHID2(I,J)=(THREE+SUMA(J)*R2(I,is))*A
1122     CONTINUE
1120     CONTINUE
C
C       FOR Px-TYPE
C
        DO 1140 J=ITPA(1)+1,ITPA(2)
        is=icta(j)
        DO 1142 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        B=DX(I,is)*A*SUMA(J)
C
        CHI(I,J)=A*DX(I,is)
        CHIX(I,J)=A+DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=(FIVE+SUMA(J)*R2(I,is))*B
1142     CONTINUE
1140     CONTINUE
C
C       FOR Py-TYPE
C
        DO 1160 J=ITPA(2)+1,ITPA(3)
        is=icta(j)
        DO 1162 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        B=DY(I,is)*A*SUMA(J)
C
        CHI(I,J)=A*DY(I,is)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=A+DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=(FIVE+SUMA(J)*R2(I,is))*B
1162     CONTINUE
1160     CONTINUE
C
C       FOR Pz-TYPE
C
        DO 1180 J=ITPA(3)+1,ITPA(4)
        is=icta(j)
        DO 1182 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        B=DZ(I,is)*A*SUMA(J)
C
        CHI(I,J)=A*DZ(I,is)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=A+DZ(I,is)*B
        CHID2(I,J)=(FIVE+SUMA(J)*R2(I,is))*B
1182     CONTINUE
1180     CONTINUE
C
C       FOR Dxx-TYPE
C
        DO 2220 J=ITPA(4)+1,ITPA(5)
        is=icta(j)
        DO 2222 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        B=DX(I,is)*DX(I,is)*A*SUMA(J)
C
        CHI(I,J)=B*DIVA(J)
        CHIX(I,J)=(TWO*A+B)*DX(I,is)
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=TWO*A+(SEVEN+SUMA(J)*R2(I,is))*B
2222     CONTINUE
2220     CONTINUE
C
C       FOR Dyy-TYPE
C
        DO 2240 J=ITPA(5)+1,ITPA(6)
        is=icta(j)
        DO 2242 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        B=DY(I,is)*DY(I,is)*A*SUMA(J)
C
        CHI(I,J)=B*DIVA(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=(TWO*A+B)*DY(I,is)
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=TWO*A+(SEVEN+SUMA(J)*R2(I,is))*B
2242     CONTINUE
2240     CONTINUE
C
C       FOR Dzz-TYPE
C
        DO 2260 J=ITPA(6)+1,ITPA(7)
        is=icta(j)
        DO 2262 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        B=DZ(I,is)*DZ(I,is)*A*SUMA(J)
C
        CHI(I,J)=B*DIVA(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=(TWO*A+B)*DZ(I,is)
        CHID2(I,J)=TWO*A+(SEVEN+SUMA(J)*R2(I,is))*B
2262     CONTINUE
2260     CONTINUE
C
C       FOR Dxy-TYPE
C
        DO 2280 J=ITPA(7)+1,ITPA(8)
        is=icta(j)
        DO 2282 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        B=DX(I,is)*DY(I,is)*A*SUMA(J)
C
        CHI(I,J)=B*DIVA(J)
        CHIX(I,J)=DX(I,is)*B+DY(I,is)*A
        CHIY(I,J)=DY(I,is)*B+DX(I,is)*A
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=(SEVEN+SUMA(J)*R2(I,is))*B
2282     CONTINUE
2280     CONTINUE
C
C       FOR Dxz-TYPE
C
        DO 3320 J=ITPA(8)+1,ITPA(9)
        is=icta(j)
        DO 3322 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        B=DX(I,is)*DZ(I,is)*A*SUMA(J)
C
        CHI(I,J)=B*DIVA(J)
        CHIX(I,J)=DX(I,is)*B+DZ(I,is)*A
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B+DX(I,is)*A
        CHID2(I,J)=(SEVEN+SUMA(J)*R2(I,is))*B
3322     CONTINUE
3320     CONTINUE
C
C       FOR Dyz-TYPE
C
        DO 3340 J=ITPA(9)+1,ITPA(10)
        is=icta(j)
        DO 3342 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        B=DY(I,is)*DZ(I,is)*A*SUMA(J)
C
        CHI(I,J)=B*DIVA(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B+DZ(I,is)*A
        CHIZ(I,J)=DZ(I,is)*B+DY(I,is)*A
        CHID2(I,J)=(SEVEN+SUMA(J)*R2(I,is))*B
3342     CONTINUE
3340     CONTINUE
C
C
C       FOR Fxxx-TYPE
C
        DO 501 J=ITPA(10)+1,ITPA(11)
        is=icta(j)
        DO 502 I=1,PTS
 
        A=DEXP(-EXXA(J)*R2(I,is))
        B=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=B*DX(I,is)
        CHIX(I,J)=(THREE + SUMA(J)*DX(I,is)*DX(I,is))*B
        CHIY(I,J)=SUMA(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=SUMA(J)*DZ(I,is)*CHI(I,J)
        CHID2(I,J)=SIX*A*DX(I,is)+
     1             (NINE+SUMA(J)*R2(I,is))*CHI(I,J)*SUMA(J)
 502    CONTINUE
 501    CONTINUE
C
C       FOR Fyyy-TYPE
C
        DO 511 J=ITPA(11)+1,ITPA(12)
        is=icta(j)
        DO 512 I=1,PTS
        A=DEXP(-EXXA(J)*R2(I,is))
        B=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=B*DY(I,is)
        CHIX(I,J)=SUMA(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(THREE + SUMA(J)*DY(I,is)*DY(I,is))*B
        CHIZ(I,J)=SUMA(J)*DZ(I,is)*CHI(I,J)
        CHID2(I,J)=SIX*A*DY(I,is)+
     1             (NINE+SUMA(J)*R2(I,is))*CHI(I,J)*SUMA(J)
C
 512    CONTINUE
 511    CONTINUE
C
C       FOR Fzzz-TYPE
C
        DO 521 J=ITPA(12)+1,ITPA(13)
        is=icta(j)
        DO 523 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        B=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=B*DZ(I,is)
        CHIX(I,J)=SUMA(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=SUMA(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(THREE + SUMA(J)*DZ(I,is)*DZ(I,is))*B
        CHID2(I,J)=SIX*A*DZ(I,is)+
     1             (NINE+SUMA(J)*R2(I,is))*CHI(I,J)*SUMA(J)
C
 523    CONTINUE
 521    CONTINUE
C
C       FOR Fxxy-TYPE
C
        DO 531 J=ITPA(13)+1,ITPA(14)
        is=icta(j)
        DO 532 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BXX=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=BXY*DX(I,is)
        CHIX(I,J)=(TWO + SUMA(J)*DX(I,is)*DX(I,is))*BXY
        CHIY(I,J)=(ONE + SUMA(J)*DY(I,is)*DY(I,is))*BXX
        CHIZ(I,J)=SUMA(J)*DZ(I,is)*CHI(I,J)
        CHID2(I,J)=TWO*A*DY(I,is)+
     1             (NINE+SUMA(J)*R2(I,is))*CHI(I,J)*SUMA(J)
C
 532    CONTINUE
 531    CONTINUE
C
C       FOR Fxxz-TYPE
C
        DO 541 J=ITPA(14)+1,ITPA(15)
        is=icta(j)
        DO 543 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        BXZ=DX(I,is)*DZ(I,is)*A
        BXX=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=BXZ*DX(I,is)
        CHIX(I,J)=(TWO + SUMA(J)*DX(I,is)*DX(I,is))*BXZ
        CHIY(I,J)=SUMA(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(ONE + SUMA(J)*DZ(I,is)*DZ(I,is))*BXX
        CHID2(I,J)=TWO*A*DZ(I,is)+
     1             (NINE+SUMA(J)*R2(I,is))*CHI(I,J)*SUMA(J)
C
 543    CONTINUE
 541    CONTINUE
C
C       FOR Fyyz-TYPE
C
        DO 561 J=ITPA(15)+1,ITPA(16)
        is=icta(j)
        DO 563 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        BYZ=DZ(I,is)*DY(I,is)*A
        BYY=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=BYZ*DY(I,is)
        CHIX(I,J)=SUMA(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(TWO + SUMA(J)*DY(I,is)*DY(I,is))*BYZ
        CHIZ(I,J)=(ONE + SUMA(J)*DZ(I,is)*DZ(I,is))*BYY
        CHID2(I,J)=TWO*A*DZ(I,is)+
     1             (NINE+SUMA(J)*R2(I,is))*CHI(I,J)*SUMA(J)
C
 563    CONTINUE
 561    CONTINUE
C
C       FOR Fxyy-TYPE
C
        DO 551 J=ITPA(16)+1,ITPA(17)
        is=icta(j)
        DO 552 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BYY=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=BXY*DY(I,is)
        CHIX(I,J)=(ONE + SUMA(J)*DX(I,is)*DX(I,is))*BYY
        CHIY(I,J)=(TWO + SUMA(J)*DY(I,is)*DY(I,is))*BXY
        CHIZ(I,J)=SUMA(J)*DZ(I,is)*CHI(I,J)
        CHID2(I,J)=TWO*A*DX(I,is)+
     1             (NINE+SUMA(J)*R2(I,is))*CHI(I,J)*SUMA(J)
C
 552    CONTINUE
 551    CONTINUE
C
C       FOR Fxzz-TYPE
C
        DO 571 J=ITPA(17)+1,ITPA(18)
        is=icta(j)
        DO 572 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        BXZ=DZ(I,is)*DX(I,is)*A
        BZZ=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=BXZ*DZ(I,is)
        CHIX(I,J)=(ONE + SUMA(J)*DX(I,is)*DX(I,is))*BZZ
        CHIY(I,J)=SUMA(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(TWO + SUMA(J)*DZ(I,is)*DZ(I,is))*BXZ
        CHID2(I,J)=TWO*A*DX(I,is)+
     1             (NINE+SUMA(J)*R2(I,is))*CHI(I,J)*SUMA(J)
C
 572    CONTINUE
 571    CONTINUE
C
C       FOR Fyzz-TYPE
C
        DO 581 J=ITPA(18)+1,ITPA(19)
        is=icta(j)
        DO 583 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        BYZ=DZ(I,is)*DY(I,is)*A
        BZZ=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=BYZ*DZ(I,is)
        CHIX(I,J)=SUMA(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(ONE + SUMA(J)*DY(I,is)*DY(I,is))*BZZ
        CHIZ(I,J)=(TWO + SUMA(J)*DZ(I,is)*DZ(I,is))*BYZ
        CHID2(I,J)=TWO*A*DY(I,is)+
     1             (NINE+SUMA(J)*R2(I,is))*CHI(I,J)*SUMA(J)
C
 583    CONTINUE
 581    CONTINUE
C
C       FOR Fxyz-TYPE
C
        DO 591 J=ITPA(19)+1,ITPA(20)
        is=icta(j)
        DO 592 I=1,PTS
C
        A=DEXP(-EXXA(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BYZ=DZ(I,is)*DY(I,is)*A
        BXZ=DX(I,is)*DZ(I,is)*A
C
        CHI(I,J)=DX(I,is)*BYZ
        CHIX(I,J)=(ONE + SUMA(J)*DX(I,is)*DX(I,is))*BYZ
        CHIY(I,J)=(ONE + SUMA(J)*DY(I,is)*DY(I,is))*BXZ
        CHIZ(I,J)=(ONE + SUMA(J)*DZ(I,is)*DZ(I,is))*BXY
        CHID2(I,J)=(NINE+SUMA(J)*R2(I,is))*CHI(I,J)*SUMA(J)
C
 592    CONTINUE
 591    CONTINUE
C
      DO 100 J=1,NPRIMSA
      Temp=Zero
      DO 101 I=1,PTS
      Check=Dmax1(Dabs(Chi(I,J)),Dabs(Chix(I,J)),Dabs(Chiy(I,J)),
     $Dabs(ChiZ(I,J)),Dabs(Chid2(I,J)),temp)
      temp=check
101   Continue
      Chimax(j)=Temp
100   Continue
C
      DO 3346 L=1,NMO
      occ=two
      if(idomag.eq.0)occ=rootpo(l)
      DO 3347 I=1,PTS
C
      PSI(I,L)=ZERO
      GX(I,L)=ZERO
      GY(I,L)=ZERO
      GZ(I,L)=ZERO
      D2(I,L)=ZERO
C 
3347  CONTINUE
C
      DO 3348 J=1,NPRIMSA
C
      Check=Dabs(COOA(J,L)*CHIMAX(J)*occ)
      IF(Check.GT.CUTOFF)THEN
      TEMP=COOA(J,L)
      DO 3349 I=1,PTS
C
      PSI(I,L)=PSI(I,L)+TEMP*CHI(I,J)
      GX(I,L)=GX(I,L)+TEMP*CHIX(I,J)
      GY(I,L)=GY(I,L)+TEMP*CHIY(I,J)
      GZ(I,L)=GZ(I,L)+TEMP*CHIZ(I,J)
      D2(I,L)=D2(I,L)+TEMP*CHID2(I,J)
C
3349  CONTINUE
      ENDIF
3348  CONTINUE
3346  CONTINUE
C
      GOTO 9999
C
555   CONTINUE
C
C       FOR S-TYPE
C
        DO 120 J = 1,ITPB(1) 
        is=ictb(j)
        DO 122 I=1,PTS 
C
        A=SUMB(J)*DEXP(-EXXB(J)*R2(I,is))
C
        CHI(I,J)=A*DIVB(J)
        CHIX(I,J)=DX(I,is)*A
        CHIY(I,J)=DY(I,is)*A
        CHIZ(I,J)=DZ(I,is)*A
        CHID2(I,J)=(THREE+SUMB(J)*R2(I,is))*A
122     CONTINUE
120     CONTINUE
C
C       FOR Px-TYPE
C
        DO 140 J=ITPB(1)+1,ITPB(2)
        is=ictb(j)
        DO 142 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DX(I,is)*A*SUMB(J)
C
        CHI(I,J)=A*DX(I,is)
        CHIX(I,J)=A+DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=(FIVE+SUMB(J)*R2(I,is))*B
142     CONTINUE
140     CONTINUE
C
C       FOR Py-TYPE
C
        DO 160 J=ITPB(2)+1,ITPB(3)
        is=ictb(j)
        DO 162 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DY(I,is)*A*SUMB(J)
C
        CHI(I,J)=A*DY(I,is)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=A+DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=(FIVE+SUMB(J)*R2(I,is))*B
162     CONTINUE
160     CONTINUE
C
C       FOR Pz-TYPE
C
        DO 180 J=ITPB(3)+1,ITPB(4)
        is=ictb(j)
        DO 182 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DZ(I,is)*A*SUMB(J)
C
        CHI(I,J)=A*DZ(I,is)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=A+DZ(I,is)*B
        CHID2(I,J)=(FIVE+SUMB(J)*R2(I,is))*B
182     CONTINUE
180     CONTINUE
C
C       FOR Dxx-TYPE
C
        DO 220 J=ITPB(4)+1,ITPB(5)
        is=ictb(j)
        DO 222 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DX(I,is)*DX(I,is)*A*SUMB(J)
C
        CHI(I,J)=B*DIVB(J)
        CHIX(I,J)=(TWO*A+B)*DX(I,is)
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=TWO*A+(SEVEN+SUMB(J)*R2(I,is))*B
222     CONTINUE
220     CONTINUE
C
C       FOR Dyy-TYPE
C
        DO 240 J=ITPB(5)+1,ITPB(6)
        is=ictb(j)
        DO 242 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DY(I,is)*DY(I,is)*A*SUMB(J)
C
        CHI(I,J)=B*DIVB(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=(TWO*A+B)*DY(I,is)
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=TWO*A+(SEVEN+SUMB(J)*R2(I,is))*B
242     CONTINUE
240     CONTINUE
C
C       FOR Dzz-TYPE
C
        DO 260 J=ITPB(6)+1,ITPB(7)
        is=ictb(j)
        DO 262 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DZ(I,is)*DZ(I,is)*A*SUMB(J)
C
        CHI(I,J)=B*DIVB(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=(TWO*A+B)*DZ(I,is)
        CHID2(I,J)=TWO*A+(SEVEN+SUMB(J)*R2(I,is))*B
262     CONTINUE
260     CONTINUE
C
C       FOR Dxy-TYPE
C
        DO 280 J=ITPB(7)+1,ITPB(8)
        is=ictb(j)
        DO 282 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DX(I,is)*DY(I,is)*A*SUMB(J)
C
        CHI(I,J)=B*DIVB(J)
        CHIX(I,J)=DX(I,is)*B+DY(I,is)*A
        CHIY(I,J)=DY(I,is)*B+DX(I,is)*A
        CHIZ(I,J)=DZ(I,is)*B
        CHID2(I,J)=(SEVEN+SUMB(J)*R2(I,is))*B
282     CONTINUE
280     CONTINUE
C
C       FOR Dxz-TYPE
C
        DO 320 J=ITPB(8)+1,ITPB(9)
        is=ictb(j)
        DO 322 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DX(I,is)*DZ(I,is)*A*SUMB(J)
C
        CHI(I,J)=B*DIVB(J)
        CHIX(I,J)=DX(I,is)*B+DZ(I,is)*A
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B+DX(I,is)*A
        CHID2(I,J)=(SEVEN+SUMB(J)*R2(I,is))*B
322     CONTINUE
320     CONTINUE
C
C       FOR Dyz-TYPE
C
        DO 340 J=ITPB(9)+1,ITPB(10)
        is=ictb(j)
        DO 342 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DY(I,is)*DZ(I,is)*A*SUMB(J)
C
        CHI(I,J)=B*DIVB(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B+DZ(I,is)*A
        CHIZ(I,J)=DZ(I,is)*B+DY(I,is)*A
        CHID2(I,J)=(SEVEN+SUMB(J)*R2(I,is))*B
342     CONTINUE
340     CONTINUE
C
C
C       FOR Fxxx-TYPE
C
        DO 901 J=ITPB(10)+1,ITPB(11)
        is=ictb(j)
        DO 902 I=1,PTS
 
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=B*DX(I,is)
        CHIX(I,J)=(THREE + SUMB(J)*DX(I,is)*DX(I,is))*B
        CHIY(I,J)=SUMB(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=SUMB(J)*DZ(I,is)*CHI(I,J)
        CHID2(I,J)=SIX*A*DX(I,is)+
     1             (NINE+SUMB(J)*R2(I,is))*CHI(I,J)*SUMB(J)
 902    CONTINUE
 901    CONTINUE
C
C       FOR Fyyy-TYPE
C
        DO 911 J=ITPB(11)+1,ITPB(12)
        is=ictb(j)
        DO 912 I=1,PTS
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=B*DY(I,is)
        CHIX(I,J)=SUMB(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(THREE + SUMB(J)*DY(I,is)*DY(I,is))*B
        CHIZ(I,J)=SUMB(J)*DZ(I,is)*CHI(I,J)
        CHID2(I,J)=SIX*A*DY(I,is)+
     1             (NINE+SUMB(J)*R2(I,is))*CHI(I,J)*SUMB(J)
C
 912    CONTINUE
 911    CONTINUE
C
C       FOR Fzzz-TYPE
C
        DO 921 J=ITPB(12)+1,ITPB(13)
        is=ictb(j)
        DO 923 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=B*DZ(I,is)
        CHIX(I,J)=SUMB(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=SUMB(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(THREE + SUMB(J)*DZ(I,is)*DZ(I,is))*B
        CHID2(I,J)=SIX*A*DZ(I,is)+
     1             (NINE+SUMB(J)*R2(I,is))*CHI(I,J)*SUMB(J)
C
 923    CONTINUE
 921    CONTINUE
C
C       FOR Fxxy-TYPE
C
        DO 931 J=ITPB(13)+1,ITPB(14)
        is=ictb(j)
        DO 932 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BXX=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=BXY*DX(I,is)
        CHIX(I,J)=(TWO + SUMB(J)*DX(I,is)*DX(I,is))*BXY
        CHIY(I,J)=(ONE + SUMB(J)*DY(I,is)*DY(I,is))*BXX
        CHIZ(I,J)=SUMB(J)*DZ(I,is)*CHI(I,J)
        CHID2(I,J)=TWO*A*DY(I,is)+
     1             (NINE+SUMB(J)*R2(I,is))*CHI(I,J)*SUMB(J)
C
 932    CONTINUE
 931    CONTINUE
C
C       FOR Fxxz-TYPE
C
        DO 941 J=ITPB(14)+1,ITPB(15)
        is=ictb(j)
        DO 943 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BXZ=DX(I,is)*DZ(I,is)*A
        BXX=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=BXZ*DX(I,is)
        CHIX(I,J)=(TWO + SUMB(J)*DX(I,is)*DX(I,is))*BXZ
        CHIY(I,J)=SUMB(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(ONE + SUMB(J)*DZ(I,is)*DZ(I,is))*BXX
        CHID2(I,J)=TWO*A*DZ(I,is)+
     1             (NINE+SUMB(J)*R2(I,is))*CHI(I,J)*SUMB(J)
C
 943    CONTINUE
 941    CONTINUE
C
C       FOR Fyyz-TYPE
C
        DO 961 J=ITPB(15)+1,ITPB(16)
        is=ictb(j)
        DO 963 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BYZ=DZ(I,is)*DY(I,is)*A
        BYY=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=BYZ*DY(I,is)
        CHIX(I,J)=SUMB(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(TWO + SUMB(J)*DY(I,is)*DY(I,is))*BYZ
        CHIZ(I,J)=(ONE + SUMB(J)*DZ(I,is)*DZ(I,is))*BYY
        CHID2(I,J)=TWO*A*DZ(I,is)+
     1             (NINE+SUMB(J)*R2(I,is))*CHI(I,J)*SUMB(J)
C
 963    CONTINUE
 961    CONTINUE
C
C       FOR Fxyy-TYPE
C
        DO 951 J=ITPB(16)+1,ITPB(17)
        is=ictb(j)
        DO 952 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BYY=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=BXY*DY(I,is)
        CHIX(I,J)=(ONE + SUMB(J)*DX(I,is)*DX(I,is))*BYY
        CHIY(I,J)=(TWO + SUMB(J)*DY(I,is)*DY(I,is))*BXY
        CHIZ(I,J)=SUMB(J)*DZ(I,is)*CHI(I,J)
        CHID2(I,J)=TWO*A*DX(I,is)+
     1             (NINE+SUMB(J)*R2(I,is))*CHI(I,J)*SUMB(J)
C
 952    CONTINUE
 951    CONTINUE
C
C       FOR Fxzz-TYPE
C
        DO 971 J=ITPB(17)+1,ITPB(18)
        is=ictb(j)
        DO 972 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BXZ=DZ(I,is)*DX(I,is)*A
        BZZ=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=BXZ*DZ(I,is)
        CHIX(I,J)=(ONE + SUMB(J)*DX(I,is)*DX(I,is))*BZZ
        CHIY(I,J)=SUMB(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(TWO + SUMB(J)*DZ(I,is)*DZ(I,is))*BXZ
        CHID2(I,J)=TWO*A*DX(I,is)+
     1             (NINE+SUMB(J)*R2(I,is))*CHI(I,J)*SUMB(J)
C
 972    CONTINUE
 971    CONTINUE
C
C       FOR Fyzz-TYPE
C
        DO 981 J=ITPB(18)+1,ITPB(19)
        is=ictb(j)
        DO 983 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BYZ=DZ(I,is)*DY(I,is)*A
        BZZ=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=BYZ*DZ(I,is)
        CHIX(I,J)=SUMB(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(ONE + SUMB(J)*DY(I,is)*DY(I,is))*BZZ
        CHIZ(I,J)=(TWO + SUMB(J)*DZ(I,is)*DZ(I,is))*BYZ
        CHID2(I,J)=TWO*A*DY(I,is)+
     1             (NINE+SUMB(J)*R2(I,is))*CHI(I,J)*SUMB(J)
C
 983    CONTINUE
 981    CONTINUE
C
C       FOR Fxyz-TYPE
C
        DO 991 J=ITPB(19)+1,ITPB(20)
        is=ictb(j)
        DO 992 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BYZ=DZ(I,is)*DY(I,is)*A
        BXZ=DX(I,is)*DZ(I,is)*A
C
        CHI(I,J)=DX(I,is)*BYZ
        CHIX(I,J)=(ONE + SUMB(J)*DX(I,is)*DX(I,is))*BYZ
        CHIY(I,J)=(ONE + SUMB(J)*DY(I,is)*DY(I,is))*BXZ
        CHIZ(I,J)=(ONE + SUMB(J)*DZ(I,is)*DZ(I,is))*BXY
        CHID2(I,J)=(NINE+SUMB(J)*R2(I,is))*CHI(I,J)*SUMB(J)
C
 992    CONTINUE
 991    CONTINUE
C
      DO 1000 J=1,NPRIMSB
      Temp=Zero
      DO 1001 I=1,PTS
      Check=Dmax1(Dabs(Chi(I,J)),Dabs(Chix(I,J)),Dabs(Chiy(I,J)),
     $Dabs(ChiZ(I,J)),Dabs(Chid2(I,J)),temp)
      temp=check
1001  Continue
      Chimax(j)=Temp
1000  Continue
C
      DO 346 L=1,NMO
      occ=two
      if(idomag.eq.0)occ=rootpo(l)
      DO 347 I=1,PTS
C
      PSI(I,L)=ZERO
      GX(I,L)=ZERO
      GY(I,L)=ZERO
      GZ(I,L)=ZERO
      D2(I,L)=ZERO
C
347   CONTINUE
C
      DO 348 J=1,NPRIMSB
      Check=Dabs(COOB(J,L)*CHIMAX(J)*occ)
      IF(Check.Gt.Cutoff)THEN
      TEMP=COOB(J,L)
      DO 349 I=1,PTS
      PSI(I,L)=PSI(I,L)+TEMP*CHI(I,J)
      GX(I,L)=GX(I,L)+TEMP*CHIX(I,J)
      GY(I,L)=GY(I,L)+TEMP*CHIY(I,J)
      GZ(I,L)=GZ(I,L)+TEMP*CHIZ(I,J)
      D2(I,L)=D2(I,L)+TEMP*CHID2(I,J)
349   CONTINUE
      ENDIF
348   CONTINUE
346   CONTINUE
C
      GOTO 9999
C
5555  CONTINUE
C
        lmo=nmo
        If(Idomag.eq.1)lmo=nmo/7
C
C       FOR S-TYPE
C
        DO 8120 J = 1,ITPB(1)
        is=ictb(j)
        DO 8122 I=1,PTS 
C
        A=SUMB(J)*DEXP(-EXXB(J)*R2(I,is))
C
        CHI(I,J)=A*DIVB(J)
        CHIX(I,J)=DX(I,is)*A
        CHIY(I,J)=DY(I,is)*A
        CHIZ(I,J)=DZ(I,is)*A
8122     CONTINUE
8120     CONTINUE
C
C       FOR Px-TYPE
C
        DO 8140 J=ITPB(1)+1,ITPB(2)
        is=ictb(j)
        DO 8142 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DX(I,is)*A*SUMB(J)
C
        CHI(I,J)=A*DX(I,is)
        CHIX(I,J)=A+DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
8142     CONTINUE
8140     CONTINUE
C
C       FOR Py-TYPE
C
        DO 8160 J=ITPB(2)+1,ITPB(3)
        is=ictb(j)
        DO 8162 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DY(I,is)*A*SUMB(J)
C
        CHI(I,J)=A*DY(I,is)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=A+DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
8162     CONTINUE
8160     CONTINUE
C
C       FOR Pz-TYPE
C
        DO 8180 J=ITPB(3)+1,ITPB(4)
        is=ictb(j)
        DO 8182 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DZ(I,is)*A*SUMB(J)
C
        CHI(I,J)=A*DZ(I,is)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=A+DZ(I,is)*B
8182     CONTINUE
8180     CONTINUE
C
C       FOR Dxx-TYPE
C
        DO 8220 J=ITPB(4)+1,ITPB(5)
        is=ictb(j)
        DO 8222 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DX(I,is)*DX(I,is)*A*SUMB(J)
C
        CHI(I,J)=B*DIVB(J)
        CHIX(I,J)=(TWO*A+B)*DX(I,is)
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B
8222     CONTINUE
8220     CONTINUE
C
C       FOR Dyy-TYPE
C
        DO 8240 J=ITPB(5)+1,ITPB(6)
        is=ictb(j)
        DO 8242 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DY(I,is)*DY(I,is)*A*SUMB(J)
C
        CHI(I,J)=B*DIVB(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=(TWO*A+B)*DY(I,is)
        CHIZ(I,J)=DZ(I,is)*B
8242     CONTINUE
8240     CONTINUE
C
C       FOR Dzz-TYPE
C
        DO 8260 J=ITPB(6)+1,ITPB(7)
        is=ictb(j)
        DO 8262 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DZ(I,is)*DZ(I,is)*A*SUMB(J)
C
        CHI(I,J)=B*DIVB(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=(TWO*A+B)*DZ(I,is)
8262     CONTINUE
8260     CONTINUE
C
C       FOR Dxy-TYPE
C
        DO 8280 J=ITPB(7)+1,ITPB(8)
        is=ictb(j)
        DO 8282 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DX(I,is)*DY(I,is)*A*SUMB(J)
C
        CHI(I,J)=B*DIVB(J)
        CHIX(I,J)=DX(I,is)*B+DY(I,is)*A
        CHIY(I,J)=DY(I,is)*B+DX(I,is)*A
        CHIZ(I,J)=DZ(I,is)*B
8282     CONTINUE
8280     CONTINUE
C
C       FOR Dxz-TYPE
C
        DO 8320 J=ITPB(8)+1,ITPB(9)
        is=ictb(j)
        DO 8322 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DX(I,is)*DZ(I,is)*A*SUMB(J)
C
        CHI(I,J)=B*DIVB(J)
        CHIX(I,J)=DX(I,is)*B+DZ(I,is)*A
        CHIY(I,J)=DY(I,is)*B
        CHIZ(I,J)=DZ(I,is)*B+DX(I,is)*A
8322     CONTINUE
8320     CONTINUE
C
C       FOR Dyz-TYPE
C
        DO 8340 J=ITPB(9)+1,ITPB(10)
        is=ictb(j)
        DO 8342 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DY(I,is)*DZ(I,is)*A*SUMB(J)
C
        CHI(I,J)=B*DIVB(J)
        CHIX(I,J)=DX(I,is)*B
        CHIY(I,J)=DY(I,is)*B+DZ(I,is)*A
        CHIZ(I,J)=DZ(I,is)*B+DY(I,is)*A
8342     CONTINUE
8340     CONTINUE
C
C
C       FOR Fxxx-TYPE
C
        DO 8901 J=ITPB(10)+1,ITPB(11)
        is=ictb(j)
        DO 8902 I=1,PTS
 
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=B*DX(I,is)
        CHIX(I,J)=(THREE + SUMB(J)*DX(I,is)*DX(I,is))*B
        CHIY(I,J)=SUMB(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=SUMB(J)*DZ(I,is)*CHI(I,J)
 8902    CONTINUE
 8901    CONTINUE
C
C       FOR Fyyy-TYPE
C
        DO 8911 J=ITPB(11)+1,ITPB(12)
        is=ictb(j)
        DO 8912 I=1,PTS
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=B*DY(I,is)
        CHIX(I,J)=SUMB(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(THREE + SUMB(J)*DY(I,is)*DY(I,is))*B
        CHIZ(I,J)=SUMB(J)*DZ(I,is)*CHI(I,J)
8912    CONTINUE
8911    CONTINUE
C
C       FOR Fzzz-TYPE
C
        DO 8921 J=ITPB(12)+1,ITPB(13)
        is=ictb(j)
        DO 8923 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        B=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=B*DZ(I,is)
        CHIX(I,J)=SUMB(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=SUMB(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(THREE + SUMB(J)*DZ(I,is)*DZ(I,is))*B
 8923    CONTINUE
 8921    CONTINUE
C
C       FOR Fxxy-TYPE
C
        DO 8931 J=ITPB(13)+1,ITPB(14)
        is=ictb(j)
        DO 8932 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BXX=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=BXY*DX(I,is)
        CHIX(I,J)=(TWO + SUMB(J)*DX(I,is)*DX(I,is))*BXY
        CHIY(I,J)=(ONE + SUMB(J)*DY(I,is)*DY(I,is))*BXX
        CHIZ(I,J)=SUMB(J)*DZ(I,is)*CHI(I,J)
 8932    CONTINUE
 8931    CONTINUE
C
C       FOR Fxxz-TYPE
C
        DO 8941 J=ITPB(14)+1,ITPB(15)
        is=ictb(j)
        DO 8943 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BXZ=DX(I,is)*DZ(I,is)*A
        BXX=DX(I,is)*DX(I,is)*A
C
        CHI(I,J)=BXZ*DX(I,is)
        CHIX(I,J)=(TWO + SUMB(J)*DX(I,is)*DX(I,is))*BXZ
        CHIY(I,J)=SUMB(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(ONE + SUMB(J)*DZ(I,is)*DZ(I,is))*BXX
 8943    CONTINUE
 8941    CONTINUE
C
C       FOR Fyyz-TYPE
C
        DO 8961 J=ITPB(15)+1,ITPB(16)
        is=ictb(j)
        DO 8963 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BYZ=DZ(I,is)*DY(I,is)*A
        BYY=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=BYZ*DY(I,is)
        CHIX(I,J)=SUMB(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(TWO + SUMB(J)*DY(I,is)*DY(I,is))*BYZ
        CHIZ(I,J)=(ONE + SUMB(J)*DZ(I,is)*DZ(I,is))*BYY
 8963    CONTINUE
 8961    CONTINUE
C
C       FOR Fxyy-TYPE
C
        DO 8951 J=ITPB(16)+1,ITPB(17)
        is=ictb(j)
        DO 8952 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BYY=DY(I,is)*DY(I,is)*A
C
        CHI(I,J)=BXY*DY(I,is)
        CHIX(I,J)=(ONE + SUMB(J)*DX(I,is)*DX(I,is))*BYY
        CHIY(I,J)=(TWO + SUMB(J)*DY(I,is)*DY(I,is))*BXY
        CHIZ(I,J)=SUMB(J)*DZ(I,is)*CHI(I,J)
 8952    CONTINUE
 8951    CONTINUE
C
C       FOR Fxzz-TYPE
C
        DO 8971 J=ITPB(17)+1,ITPB(18)
        is=ictb(j)
        DO 8972 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BXZ=DZ(I,is)*DX(I,is)*A
        BZZ=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=BXZ*DZ(I,is)
        CHIX(I,J)=(ONE + SUMB(J)*DX(I,is)*DX(I,is))*BZZ
        CHIY(I,J)=SUMB(J)*DY(I,is)*CHI(I,J)
        CHIZ(I,J)=(TWO + SUMB(J)*DZ(I,is)*DZ(I,is))*BXZ
 8972    CONTINUE
 8971    CONTINUE
C
C       FOR Fyzz-TYPE
C
        DO 8981 J=ITPB(18)+1,ITPB(19)
        is=ictb(j)
        DO 8983 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BYZ=DZ(I,is)*DY(I,is)*A
        BZZ=DZ(I,is)*DZ(I,is)*A
C
        CHI(I,J)=BYZ*DZ(I,is)
        CHIX(I,J)=SUMB(J)*DX(I,is)*CHI(I,J)
        CHIY(I,J)=(ONE + SUMB(J)*DY(I,is)*DY(I,is))*BZZ
        CHIZ(I,J)=(TWO + SUMB(J)*DZ(I,is)*DZ(I,is))*BYZ
 8983    CONTINUE
 8981    CONTINUE
C
C       FOR Fxyz-TYPE
C
        DO 8991 J=ITPB(19)+1,ITPB(20)
        is=ictb(j)
        DO 8992 I=1,PTS
C
        A=DEXP(-EXXB(J)*R2(I,is))
        BXY=DX(I,is)*DY(I,is)*A
        BYZ=DZ(I,is)*DY(I,is)*A
        BXZ=DX(I,is)*DZ(I,is)*A
C
        CHI(I,J)=DX(I,is)*BYZ
        CHIX(I,J)=(ONE + SUMB(J)*DX(I,is)*DX(I,is))*BYZ
        CHIY(I,J)=(ONE + SUMB(J)*DY(I,is)*DY(I,is))*BXZ
        CHIZ(I,J)=(ONE + SUMB(J)*DZ(I,is)*DZ(I,is))*BXY
8992    CONTINUE
8991    CONTINUE
C
      DO 8000 J=1,NPRIMSB
      Temp=Zero
      DO 8001 I=1,PTS
      Check=Dmax1(Dabs(Chi(I,J)),Dabs(Chix(I,J)),Dabs(Chiy(I,J)),
     $Dabs(ChiZ(I,J)),temp)
      temp=check
8001  Continue
      Chimax(j)=Temp
8000  Continue
C
      DO 850 I=1,PTS
      DXR(I)=Zero
      DYR(I)=Zero
      DZR(I)=Zero
850   Continue
C
      DO 846 L=1,LMO
      occ=rootpo(L)
      DO 847 I=1,PTS
      PSI(I,1)=ZERO
      GX(I,1)=ZERO
      GY(I,1)=ZERO
      GZ(I,1)=ZERO
847   CONTINUE
C
      DO 848 J=1,NPRIMSB
      Check=Dabs(COOB(J,L)*CHIMAX(J)*occ)
      IF(Check.Gt.CutofT)THEN
      Temp=coob(j,l)
      DO 849 I=1,PTS
      PSI(I,1)=PSI(I,1)+temp*CHI(I,J)
      GX(I,1)=GX(I,1)+temp*CHIX(I,J)
      GY(I,1)=GY(I,1)+temp*CHIY(I,J)
      GZ(I,1)=GZ(I,1)+temp*CHIZ(I,J)
849   CONTINUE
      ENDIF
848   CONTINUE
      Do 859 I=1,PTS
      DXR(I)=DXR(I)+PO(L)*PSI(I,1)*GX(I,1)
      DYR(I)=DYR(I)+PO(L)*PSI(I,1)*GY(I,1)
      DZR(I)=DZR(I)+PO(L)*PSI(I,1)*GZ(I,1)
859   Continue
846   CONTINUE
C
9999  Continue
C
      RETURN
      END
      SUBROUTINE GAUSCHECK(BETA,NT,NBCP,CRIT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION NINE
      PARAMETER (MCENT=50, MMO=300, MPRIMS0=800, MPTS0=200, 
     $MPRIMS=700,MPTSX=30, NTYPE=20)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /UNITS/ ISRF,INPT,IOUT,IWFN,IDBG
      COMMON /ORBTL/ EORB(MMO),PO(MMO),ROOTPO(MMO),NMO
      COMMON /PRIMS0/ COOMX(MMO),SUM(MPRIMS0),DIV(MPRIMS0),
     +COO(MPRIMS0,MMO),EXX(MPRIMS0),ICT(MPRIMS0),ITP(NTYPE),NPRIMS,
     +ITYPE(MPRIMS0)
      COMMON /NCUT/ CUTOFF,NDOCUT,NPR,NPRA,NPRB,NZEROA,NZEROB
      COMMON /PRIMSA/ COOAMX(MMO),SUMA(MPRIMS),DIVA(MPRIMS),
     +COOA(MPRIMS,MMO),EXXA(MPRIMS),ICTA(MPRIMS),ITPA(NTYPE),NPRIMSA
      COMMON /PRIMSB/ COOBMX(MMO),SUMB(MPRIMS),DIVB(MPRIMS),
     +COOB(MPRIMS,MMO),EXXB(MPRIMS),ICTB(MPRIMS),ITPB(NTYPE),NPRIMSB
      COMMON /C7/ CT(3,3),X,Y,Z,NCRNT
      DIMENSION DX(MPTS0,MCENT),DY(MPTS0,MCENT),DZ(MPTS0,MCENT),
     $R2(MPTS0,MCENT),CHI(MPTS0,MPRIMS0),CHIX(MPTS0,MPRIMS0),
     $CHIY(MPTS0,MPRIMS0),CHIZ(MPTS0,MPRIMS0),CHID2(MPTS0,MPRIMS0),
     $CRIT(3,20),IDEL(MPRIMS0),XYZ(MPTS0,3),ITEMP(MPRIMS0,MMO),
     $ITYPEA(MPRIMS),ITYPEB(MPRIMS),STEP(MPTSX),ISTEP(MPTSX),
     $PSI(MPTS0,MMO),GX(MPTS0,MMO),GY(MPTS0,MMO),GZ(MPTS0,MMO),
     $GRX(MPTS0),GRY(MPTS0),GRZ(MPTS0),STPTX(21),STPTY(21),
     $STPTZ(21),BITA(21),chimax(mprims0)
      DATA ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,FOUR/4.0D0/,FIVE/5.0D0/,
     $SEVEN/7.0D0/,THREE/3.0D0/,SIX/6.0D0/,NINE/9.0D0/,PT1/0.1d0/,
     $Pt2/0.2d0/,Pt05/0.05d0/,IPHI/12/,ITHETA/6/,
     $DINC/0.20D0/,ITIMES/mptsx/
      DATA STEP /0.01D0,0.01D0,0.015D0,0.015D0,0.02D0,0.025D0,0.03D0,
     +     0.035D0,0.04D0,0.05D0,0.06D0,0.07D0,0.08D0,0.09D0,0.1D0,
     +     0.1D0,0.1D0,0.1D0,0.1D0,0.1D0,0.1D0,0.1D0,0.1D0,0.1D0,
     +     0.1D0,0.1D0,0.1D0,0.1D0,0.1D0,0.1D0/
      DATA ISTEP /7,7,7,7,6,6,6,6,5,5,4,4,4,3,3,3,3,3,3,3,3,3,3,3,
     +     3,3,3,3,3,3/
3300  FORMAT(1I3,' Beta sphere prims  Redimension MPRIMS')
3301  FORMAT(1I3,' prims outside Beta sphere  Redimension MPRIMS')
962   FORMAT(' PRE-INTEGRATION PRIMITIVE CUTOFF ALGORITHM NOT USED ')
969   FORMAT(' PRE-INTEGRATION PRIMITIVE CUTOFF ALGORITHM USED ')
970   FORMAT(' DYNAMIC CUTOFFS USED THROUGHOUT ')
1971  Format(' CUTOFF VALUE USED IS ',1PE8.2)
964   FORMAT(' TOTAL NUMBER OF PRIMITIVES = ',1I6)
965   FORMAT
     +(' NUMBER OF PRIMITIVES USED OUTSIDE BETA SPHERE= ',1I6)
966   FORMAT
     +(' NUMBER OF PRIMITIVES USED INSIDE BETA SPHERE= ',1I6)
967   FORMAT
     +(1I6,' OF THE ',1I6, ' PRIM COEFFS ZEROED OUTSIDE BETA SPHERE')
968   FORMAT
     +(1I6,' OF THE ', 1I6,' PRIM COEFFS ZEROED INSIDE BETA SPHERE')
C
      PI=DACOS(-one)
C
      IF(NDOCUT.EQ.1) GOTO 111
C
      DO 95 L=1,NMO
      COOAMX(L)=COOMX(L)
      COOBMX(L)=COOMX(L)
      DO 96 J=1,NPRIMS
      COOA(J,L)=COO(J,L)
      COOB(J,L)=COO(J,L)
96    CONTINUE
95    CONTINUE
C
      DO 97 J=1,NPRIMS
      SUMA(J)=SUM(J)
      SUMB(J)=SUM(J)
      ICTA(J)=ICT(J)
      ICTB(J)=ICT(J)
      DIVA(J)=DIV(J)
      DIVB(J)=DIV(J)
      EXXA(J)=EXX(J)
      EXXB(J)=EXX(J)
97    CONTINUE
C
      DO 98 I=1,NTYPE
      ITPA(I)=ITP(I)
      ITPB(I)=ITP(I)
98    CONTINUE
C
      NPRIMSA=NPRIMS
      NPRIMSB=NPRIMS
C
      GOTO 999
C
111   THINC=PI/ITHETA
      HTHINC=THINC/TWO
      PHINC=TWO*PI/IPHI
      HPHINC=PHINC/TWO
C     
      DO 185 L=1,NMO
      DO 186 J=1,NPRIMS
      IDEL(J)=0
      ITEMP(J,L)=0
186   CONTINUE
185   CONTINUE
C
      MARK=IDINT(BETA/DINC)
      NPTS=MARK+2
      DO 1 II=1,IPHI
      PH=II*PHINC-HPHINC
      DO 2 JJ=1,ITHETA
      TH=JJ*THINC-HTHINC
C
      XYZ(1,1)=X+PT05*DSIN(TH)*DCOS(PH)
      XYZ(1,2)=Y+PT05*DSIN(TH)*DSIN(PH)
      XYZ(1,3)=Z+PT05*DCOS(TH)
      XYZ(2,1)=X+BETA*DSIN(TH)*DCOS(PH)
      XYZ(2,2)=Y+BETA*DSIN(TH)*DSIN(PH)
      XYZ(2,3)=Z+BETA*DCOS(TH)
C
      DO 3 K=1,MARK
      OFF=K*DINC
      XYZ(K+2,1)=X+OFF*DSIN(TH)*DCOS(PH)
      XYZ(K+2,2)=Y+OFF*DSIN(TH)*DSIN(PH)
      XYZ(K+2,3)=Z+OFF*DCOS(TH)
3     CONTINUE
C
      DO 110 J = 1,NCENT
       DO 112 I=1,NPTS 
        DX(I,J) = XYZ(I,1) - XC(J)
        DY(I,J) = XYZ(I,2) - YC(J)
        DZ(I,J) = XYZ(I,3) - ZC(J)
        R2(I,J)= DX(I,J)*DX(I,J)+DY(I,J)*DY(I,J)+DZ(I,J)*DZ(I,J)
112   CONTINUE
110   CONTINUE
C
C    CALCULATE THE VALUE OF THE BASIS FUNCITONS AND
C    THEIR FIRST DERIVATIVES AT THE PTS POINTS. 
C
C
C       FOR S-TYPE
C
        DO 120 J = 1,ITP(1) 
        DO 122 I=1,NPTS 
C
        A=SUM(J)*DEXP(-EXX(J)*R2(I,ICT(J)))
C
        CHI(I,J)=A*DIV(J)
        CHIX(I,J)=DX(I,ICT(J))*A
        CHIY(I,J)=DY(I,ICT(J))*A
        CHIZ(I,J)=DZ(I,ICT(J))*A
        CHID2(I,J)=(THREE+SUM(J)*R2(I,ICT(J)))*A
122     CONTINUE
120     CONTINUE
C
C       FOR Px-TYPE
C
        DO 140 J=ITP(1)+1,ITP(2)
        DO 142 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DX(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=A*DX(I,ICT(J))
        CHIX(I,J)=A+DX(I,ICT(J))*B
        CHIY(I,J)=DY(I,ICT(J))*B
        CHIZ(I,J)=DZ(I,ICT(J))*B
        CHID2(I,J)=(FIVE+SUM(J)*R2(I,ICT(J)))*B
142     CONTINUE
140     CONTINUE
C
C       FOR Py-TYPE
C
        DO 160 J=ITP(2)+1,ITP(3)
        DO 162 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DY(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=A*DY(I,ICT(J))
        CHIX(I,J)=DX(I,ICT(J))*B
        CHIY(I,J)=A+DY(I,ICT(J))*B
        CHIZ(I,J)=DZ(I,ICT(J))*B
        CHID2(I,J)=(FIVE+SUM(J)*R2(I,ICT(J)))*B
162     CONTINUE
160     CONTINUE
C
C       FOR Pz-TYPE
C
        DO 180 J=ITP(3)+1,ITP(4)
        DO 182 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DZ(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=A*DZ(I,ICT(J))
        CHIX(I,J)=DX(I,ICT(J))*B
        CHIY(I,J)=DY(I,ICT(J))*B
        CHIZ(I,J)=A+DZ(I,ICT(J))*B
        CHID2(I,J)=(FIVE+SUM(J)*R2(I,ICT(J)))*B
182     CONTINUE
180     CONTINUE
C
C       FOR Dxx-TYPE
C
        DO 220 J=ITP(4)+1,ITP(5)
        DO 222 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DX(I,ICT(J))*DX(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=(TWO*A+B)*DX(I,ICT(J))
        CHIY(I,J)=DY(I,ICT(J))*B
        CHIZ(I,J)=DZ(I,ICT(J))*B
        CHID2(I,J)=TWO*A+(SEVEN+SUM(J)*R2(I,ICT(J)))*B
222     CONTINUE
220     CONTINUE
C
C       FOR Dyy-TYPE
C
        DO 240 J=ITP(5)+1,ITP(6)
        DO 242 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DY(I,ICT(J))*DY(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,ICT(J))*B
        CHIY(I,J)=(TWO*A+B)*DY(I,ICT(J))
        CHIZ(I,J)=DZ(I,ICT(J))*B
        CHID2(I,J)=TWO*A+(SEVEN+SUM(J)*R2(I,ICT(J)))*B
242     CONTINUE
240     CONTINUE
C
C       FOR Dzz-TYPE
C
        DO 260 J=ITP(6)+1,ITP(7)
        DO 262 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DZ(I,ICT(J))*DZ(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,ICT(J))*B
        CHIY(I,J)=DY(I,ICT(J))*B
        CHIZ(I,J)=(TWO*A+B)*DZ(I,ICT(J))
        CHID2(I,J)=TWO*A+(SEVEN+SUM(J)*R2(I,ICT(J)))*B
262     CONTINUE
260     CONTINUE
C
C       FOR Dxy-TYPE
C
        DO 280 J=ITP(7)+1,ITP(8)
        DO 282 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DX(I,ICT(J))*DY(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,ICT(J))*B+DY(I,ICT(J))*A
        CHIY(I,J)=DY(I,ICT(J))*B+DX(I,ICT(J))*A
        CHIZ(I,J)=DZ(I,ICT(J))*B
        CHID2(I,J)=(SEVEN+SUM(J)*R2(I,ICT(J)))*B
282     CONTINUE
280     CONTINUE
C
C       FOR Dxz-TYPE
C
        DO 320 J=ITP(8)+1,ITP(9)
        DO 322 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DX(I,ICT(J))*DZ(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,ICT(J))*B+DZ(I,ICT(J))*A
        CHIY(I,J)=DY(I,ICT(J))*B
        CHIZ(I,J)=DZ(I,ICT(J))*B+DX(I,ICT(J))*A
        CHID2(I,J)=(SEVEN+SUM(J)*R2(I,ICT(J)))*B
322     CONTINUE
320     CONTINUE
C
C       FOR Dyz-TYPE
C
        DO 340 J=ITP(9)+1,ITP(10)
        DO 342 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DY(I,ICT(J))*DZ(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,ICT(J))*B
        CHIY(I,J)=DY(I,ICT(J))*B+DZ(I,ICT(J))*A
        CHIZ(I,J)=DZ(I,ICT(J))*B+DY(I,ICT(J))*A
        CHID2(I,J)=(SEVEN+SUM(J)*R2(I,ICT(J)))*B
342     CONTINUE
340     CONTINUE
C
C       FOR Fxxx-TYPE
C
        DO 501 J=ITP(10)+1,ITP(11)
        DO 502 I=1,NPTS
 
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DX(I,ICT(J))*DX(I,ICT(J))*A
C
        CHI(I,J)=B*DX(I,ICT(J))
        CHIX(I,J)=(THREE + SUM(J)*DX(I,ICT(J))*DX(I,ICT(J)))*B
        CHIY(I,J)=SUM(J)*DY(I,ICT(J))*CHI(I,J)
        CHIZ(I,J)=SUM(J)*DZ(I,ICT(J))*CHI(I,J)
        CHID2(I,J)=SIX*A*DX(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
 502    CONTINUE
 501    CONTINUE
C
C       FOR Fyyy-TYPE
C
        DO 511 J=ITP(11)+1,ITP(12)
        DO 512 I=1,NPTS
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DY(I,ICT(J))*DY(I,ICT(J))*A
C
        CHI(I,J)=B*DY(I,ICT(J))
        CHIX(I,J)=SUM(J)*DX(I,ICT(J))*CHI(I,J)
        CHIY(I,J)=(THREE + SUM(J)*DY(I,ICT(J))*DY(I,ICT(J)))*B
        CHIZ(I,J)=SUM(J)*DZ(I,ICT(J))*CHI(I,J)
        CHID2(I,J)=SIX*A*DY(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 512    CONTINUE
 511    CONTINUE
C
C       FOR Fzzz-TYPE
C
        DO 521 J=ITP(12)+1,ITP(13)
        DO 523 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DZ(I,ICT(J))*DZ(I,ICT(J))*A
C
        CHI(I,J)=B*DZ(I,ICT(J))
        CHIX(I,J)=SUM(J)*DX(I,ICT(J))*CHI(I,J)
        CHIY(I,J)=SUM(J)*DY(I,ICT(J))*CHI(I,J)
        CHIZ(I,J)=(THREE + SUM(J)*DZ(I,ICT(J))*DZ(I,ICT(J)))*B
        CHID2(I,J)=SIX*A*DZ(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 523    CONTINUE
 521    CONTINUE
C
C       FOR Fxxy-TYPE
C
        DO 531 J=ITP(13)+1,ITP(14)
        DO 532 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BXY=DX(I,ICT(J))*DY(I,ICT(J))*A
        BXX=DX(I,ICT(J))*DX(I,ICT(J))*A
C
        CHI(I,J)=BXY*DX(I,ICT(J))
        CHIX(I,J)=(TWO + SUM(J)*DX(I,ICT(J))*DX(I,ICT(J)))*BXY
        CHIY(I,J)=(ONE + SUM(J)*DY(I,ICT(J))*DY(I,ICT(J)))*BXX
        CHIZ(I,J)=SUM(J)*DZ(I,ICT(J))*CHI(I,J)
        CHID2(I,J)=TWO*A*DY(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 532    CONTINUE
 531    CONTINUE
C
C       FOR Fxxz-TYPE
C
        DO 541 J=ITP(14)+1,ITP(15)
        DO 543 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BXZ=DX(I,ICT(J))*DZ(I,ICT(J))*A
        BXX=DX(I,ICT(J))*DX(I,ICT(J))*A
C
        CHI(I,J)=BXZ*DX(I,ICT(J))
        CHIX(I,J)=(TWO + SUM(J)*DX(I,ICT(J))*DX(I,ICT(J)))*BXZ
        CHIY(I,J)=SUM(J)*DY(I,ICT(J))*CHI(I,J)
        CHIZ(I,J)=(ONE + SUM(J)*DZ(I,ICT(J))*DZ(I,ICT(J)))*BXX
        CHID2(I,J)=TWO*A*DZ(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 543    CONTINUE
 541    CONTINUE
C
C       FOR Fyyz-TYPE
C
        DO 561 J=ITP(15)+1,ITP(16)
        DO 563 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BYZ=DZ(I,ICT(J))*DY(I,ICT(J))*A
        BYY=DY(I,ICT(J))*DY(I,ICT(J))*A
C
        CHI(I,J)=BYZ*DY(I,ICT(J))
        CHIX(I,J)=SUM(J)*DX(I,ICT(J))*CHI(I,J)
        CHIY(I,J)=(TWO + SUM(J)*DY(I,ICT(J))*DY(I,ICT(J)))*BYZ
        CHIZ(I,J)=(ONE + SUM(J)*DZ(I,ICT(J))*DZ(I,ICT(J)))*BYY
        CHID2(I,J)=TWO*A*DZ(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 563    CONTINUE
 561    CONTINUE
C
C       FOR Fxyy-TYPE
C
        DO 551 J=ITP(16)+1,ITP(17)
        DO 552 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BXY=DX(I,ICT(J))*DY(I,ICT(J))*A
        BYY=DY(I,ICT(J))*DY(I,ICT(J))*A
C
        CHI(I,J)=BXY*DY(I,ICT(J))
        CHIX(I,J)=(ONE + SUM(J)*DX(I,ICT(J))*DX(I,ICT(J)))*BYY
        CHIY(I,J)=(TWO + SUM(J)*DY(I,ICT(J))*DY(I,ICT(J)))*BXY
        CHIZ(I,J)=SUM(J)*DZ(I,ICT(J))*CHI(I,J)
        CHID2(I,J)=TWO*A*DX(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 552    CONTINUE
 551    CONTINUE
C
C       FOR Fxzz-TYPE
C
        DO 571 J=ITP(17)+1,ITP(18)
        DO 572 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BXZ=DZ(I,ICT(J))*DX(I,ICT(J))*A
        BZZ=DZ(I,ICT(J))*DZ(I,ICT(J))*A
C
        CHI(I,J)=BXZ*DZ(I,ICT(J))
        CHIX(I,J)=(ONE + SUM(J)*DX(I,ICT(J))*DX(I,ICT(J)))*BZZ
        CHIY(I,J)=SUM(J)*DY(I,ICT(J))*CHI(I,J)
        CHIZ(I,J)=(TWO + SUM(J)*DZ(I,ICT(J))*DZ(I,ICT(J)))*BXZ
        CHID2(I,J)=TWO*A*DX(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 572    CONTINUE
 571    CONTINUE
C
C       FOR Fyzz-TYPE
C
        DO 581 J=ITP(18)+1,ITP(19)
        DO 583 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BYZ=DZ(I,ICT(J))*DY(I,ICT(J))*A
        BZZ=DZ(I,ICT(J))*DZ(I,ICT(J))*A
C
        CHI(I,J)=BYZ*DZ(I,ICT(J))
        CHIX(I,J)=SUM(J)*DX(I,ICT(J))*CHI(I,J)
        CHIY(I,J)=(ONE + SUM(J)*DY(I,ICT(J))*DY(I,ICT(J)))*BZZ
        CHIZ(I,J)=(TWO + SUM(J)*DZ(I,ICT(J))*DZ(I,ICT(J)))*BYZ
        CHID2(I,J)=TWO*A*DY(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 583    CONTINUE
 581    CONTINUE
C
C       FOR Fxyz-TYPE
C
        DO 591 J=ITP(19)+1,ITP(20)
        DO 592 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BXY=DX(I,ICT(J))*DY(I,ICT(J))*A
        BYZ=DZ(I,ICT(J))*DY(I,ICT(J))*A
        BXZ=DX(I,ICT(J))*DZ(I,ICT(J))*A
C
        CHI(I,J)=DX(I,ICT(J))*BYZ
        CHIX(I,J)=(ONE + SUM(J)*DX(I,ICT(J))*DX(I,ICT(J)))*BYZ
        CHIY(I,J)=(ONE + SUM(J)*DY(I,ICT(J))*DY(I,ICT(J)))*BXZ
        CHIZ(I,J)=(ONE + SUM(J)*DZ(I,ICT(J))*DZ(I,ICT(J)))*BXY
        CHID2(I,J)=(NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 592    CONTINUE
 591    CONTINUE
C
        Do 593 J=1,NPRIMS
        temp=zero
        Do 594 I=1,NPTS
        check=Dmax1(Dabs(chi(i,j)),dabs(chix(i,j)),dabs(chiy(i,j)),
     $        dabs(chiz(i,j)),dabs(chid2(i,j)),temp)
        If(check.gt.temp)temp=check
594     Continue
        chimax(j)=temp
593     continue
C
        DO 124 L=1,NMO
        DO 125 J = 1,NPRIMS
        If(dabs(coo(j,l)*chimax(j)*rootpo(l)).lt.cutoff)then
          ITEMP(J,L)=ITEMP(J,L)+1
          IDEL(J)=IDEL(J)+1
          ENDIF
125      CONTINUE
124      CONTINUE
2        CONTINUE
1        CONTINUE
C
         ICHK0=IPHI*ITHETA
         ICHK=NMO*IPHI*ITHETA
         NPRIMSB=NPRIMS
         NCNT=0
         DO 128 J=1,NPRIMS
C
         IF(IDEL(J).EQ.ICHK.AND.ICT(J).NE.NCRNT) THEN
         NPRIMSB=NPRIMSB-1
         ELSE
         NCNT=NCNT+1
         DO 129 L=1,NMO
         IF(ITEMP(J,L).EQ.ICHK0.AND.ICT(J).NE.NCRNT) THEN
         COOB(NCNT,L)=ZERO
         NZEROB=NZEROB+1
         ELSE
         COOB(NCNT,L)=COO(J,L)
         ENDIF
129      CONTINUE
         ITYPEB(NCNT)=ITYPE(J)
         EXXB(NCNT)=EXX(J)
         SUMB(NCNT)=SUM(J)
         DIVB(NCNT)=DIV(J)
         ICTB(NCNT)=ICT(J)
         ENDIF
128      CONTINUE
C
          N=0
         DO 131 K=1,NTYPE
         DO 132 J=1,NPRIMSB
         IF(ITYPEB(J).EQ.K) THEN
         N=N+1
         ENDIF
132      CONTINUE
         ITPB(K)=N
131      CONTINUE
C
        DO 195 L=1,NMO
        DO 196 J=1,NPRIMS
        ITEMP(J,L)=0
        IDEL(J)=0
196     CONTINUE
195     CONTINUE
C
        STPTX(1)=X
        STPTY(1)=Y
        STPTZ(1)=Z
        BITA(1)=BETA-Pt1
        NCT=1
        DO 197 I=1,NBCP
        DST=DSQRT((CRIT(1,I)-X)**2+(CRIT(2,I)-Y)**2+(CRIT(3,I)-Z)**2)
        IF (DST.GT.(BETA+Pt2))THEN
        NCT=NCT+1
        STPTX(NCT)=((CRIT(1,I)-X)*(BETA+(DST-BETA)/TWO))/DST + X
        STPTY(NCT)=((CRIT(2,I)-Y)*(BETA+(DST-BETA)/TWO))/DST + Y
        STPTZ(NCT)=((CRIT(3,I)-Z)*(BETA+(DST-BETA)/TWO))/DST + Z
        BITA(NCT)=(DST-BETA)/THREE
        ENDIF
197     CONTINUE
C
        DO 198 JL=1,NCT
C
        N=0
        IF(JL.GT.1)THEN
        DO 700 II=1,IPHI/2
        PH=2*PHINC*II-PHINC
        DO 710 J=1,ITHETA/2
        TH=2*J*THINC-THINC
        N=N+1
      XYZ(N,1)=STPTX(JL)+BITA(JL)*DSIN(TH)*DCOS(PH)
      XYZ(N,2)=STPTY(JL)+BITA(JL)*DSIN(TH)*DSIN(PH)
      XYZ(N,3)=STPTZ(JL)+BITA(JL)*DCOS(TH)
710    CONTINUE
700    CONTINUE
       NPTS=ITHETA*IPHI/4
C
        ELSE
        DO 200 II=1,IPHI
        PH=PHINC*II-HPHINC
        DO 210 J=1,ITHETA
        TH=J*THINC-HTHINC
        N=N+1
      XYZ(N,1)=STPTX(JL)+BITA(JL)*DSIN(TH)*DCOS(PH)
      XYZ(N,2)=STPTY(JL)+BITA(JL)*DSIN(TH)*DSIN(PH)
      XYZ(N,3)=STPTZ(JL)+BITA(JL)*DCOS(TH)
210    CONTINUE
200    CONTINUE
        NPTS=ITHETA*IPHI
        ENDIF
      DO 303 II=1,ITIMES
      DO 304 JJ=1,ISTEP(II)
C
      IF(JL.EQ.NCT.AND.II.EQ.ITIMES.AND.JJ.EQ.ISTEP(ITIMES))THEN
      DO 307 I=1,NT
      XYZ(NPTS+I,1)=CRIT(1,I)
      XYZ(NPTS+I,2)=CRIT(2,I)
      XYZ(NPTS+I,3)=CRIT(3,I)
307   CONTINUE
      IF(NCT.GT.1)THEN
      DO 217 I=2,NCT
      XYZ(NPTS+NT+I,1)=STPTX(I)
      XYZ(NPTS+NT+I,2)=STPTY(I)
      XYZ(NPTS+NT+I,3)=STPTZ(I)
217   CONTINUE
      NPTS=NPTS+NT+NCT-1
      ELSE
      NPTS=NPTS+NT
      ENDIF
      ENDIF 
C
      DO 310 J = 1,NCENT
       DO 312 I=1,NPTS 
        DX(I,J) = XYZ(I,1) - XC(J)
        DY(I,J) = XYZ(I,2) - YC(J)
        DZ(I,J) = XYZ(I,3) - ZC(J)
        R2(I,J)= DX(I,J)*DX(I,J)+DY(I,J)*DY(I,J)+DZ(I,J)*DZ(I,J)
312   CONTINUE
310   CONTINUE
C
C    CALCULATE THE VALUE OF THE BASIS FUNCITONS AND
C    THEIR FIRST DERIVATIVES AT THE PTS POINTS. 
C
C
C       FOR S-TYPE
C
        DO 620 J = 1,ITP(1) 
        DO 622 I=1,NPTS 
C
        A=SUM(J)*DEXP(-EXX(J)*R2(I,ICT(J)))
C
        CHI(I,J)=A*DIV(J)
        CHIX(I,J)=DX(I,ICT(J))*A
        CHIY(I,J)=DY(I,ICT(J))*A
        CHIZ(I,J)=DZ(I,ICT(J))*A
        CHID2(I,J)=(THREE+SUM(J)*R2(I,ICT(J)))*A
622     CONTINUE
620     CONTINUE
C
C       FOR Px-TYPE
C
        DO 640 J=ITP(1)+1,ITP(2)
        DO 642 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DX(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=A*DX(I,ICT(J))
        CHIX(I,J)=A+DX(I,ICT(J))*B
        CHIY(I,J)=DY(I,ICT(J))*B
        CHIZ(I,J)=DZ(I,ICT(J))*B
        CHID2(I,J)=(FIVE+SUM(J)*R2(I,ICT(J)))*B
642     CONTINUE
640     CONTINUE
C
C       FOR Py-TYPE
C
        DO 360 J=ITP(2)+1,ITP(3)
        DO 362 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DY(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=A*DY(I,ICT(J))
        CHIX(I,J)=DX(I,ICT(J))*B
        CHIY(I,J)=A+DY(I,ICT(J))*B
        CHIZ(I,J)=DZ(I,ICT(J))*B
        CHID2(I,J)=(FIVE+SUM(J)*R2(I,ICT(J)))*B
362     CONTINUE
360     CONTINUE
C
C       FOR Pz-TYPE
C
        DO 380 J=ITP(3)+1,ITP(4)
        DO 382 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DZ(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=A*DZ(I,ICT(J))
        CHIX(I,J)=DX(I,ICT(J))*B
        CHIY(I,J)=DY(I,ICT(J))*B
        CHIZ(I,J)=A+DZ(I,ICT(J))*B
        CHID2(I,J)=(FIVE+SUM(J)*R2(I,ICT(J)))*B
382     CONTINUE
380     CONTINUE
C
C       FOR Dxx-TYPE
C
        DO 820 J=ITP(4)+1,ITP(5)
        DO 822 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DX(I,ICT(J))*DX(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=(TWO*A+B)*DX(I,ICT(J))
        CHIY(I,J)=DY(I,ICT(J))*B
        CHIZ(I,J)=DZ(I,ICT(J))*B
        CHID2(I,J)=TWO*A+(SEVEN+SUM(J)*R2(I,ICT(J)))*B
822     CONTINUE
820     CONTINUE
C
C       FOR Dyy-TYPE
C
        DO 840 J=ITP(5)+1,ITP(6)
        DO 842 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DY(I,ICT(J))*DY(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,ICT(J))*B
        CHIY(I,J)=(TWO*A+B)*DY(I,ICT(J))
        CHIZ(I,J)=DZ(I,ICT(J))*B
        CHID2(I,J)=TWO*A+(SEVEN+SUM(J)*R2(I,ICT(J)))*B
842     CONTINUE
840     CONTINUE
C
C       FOR Dzz-TYPE
C
        DO 660 J=ITP(6)+1,ITP(7)
        DO 662 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DZ(I,ICT(J))*DZ(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,ICT(J))*B
        CHIY(I,J)=DY(I,ICT(J))*B
        CHIZ(I,J)=(TWO*A+B)*DZ(I,ICT(J))
        CHID2(I,J)=TWO*A+(SEVEN+SUM(J)*R2(I,ICT(J)))*B
662     CONTINUE
660     CONTINUE
C
C       FOR Dxy-TYPE
C
        DO 680 J=ITP(7)+1,ITP(8)
        DO 682 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DX(I,ICT(J))*DY(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,ICT(J))*B+DY(I,ICT(J))*A
        CHIY(I,J)=DY(I,ICT(J))*B+DX(I,ICT(J))*A
        CHIZ(I,J)=DZ(I,ICT(J))*B
        CHID2(I,J)=(SEVEN+SUM(J)*R2(I,ICT(J)))*B
682     CONTINUE
680     CONTINUE
C
C       FOR Dxz-TYPE
C
        DO 420 J=ITP(8)+1,ITP(9)
        DO 422 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DX(I,ICT(J))*DZ(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,ICT(J))*B+DZ(I,ICT(J))*A
        CHIY(I,J)=DY(I,ICT(J))*B
        CHIZ(I,J)=DZ(I,ICT(J))*B+DX(I,ICT(J))*A
        CHID2(I,J)=(SEVEN+SUM(J)*R2(I,ICT(J)))*B
422     CONTINUE
420     CONTINUE
C
C       FOR Dyz-TYPE
C
        DO 440 J=ITP(9)+1,ITP(10)
        DO 442 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DY(I,ICT(J))*DZ(I,ICT(J))*A*SUM(J)
C
        CHI(I,J)=B*DIV(J)
        CHIX(I,J)=DX(I,ICT(J))*B
        CHIY(I,J)=DY(I,ICT(J))*B+DZ(I,ICT(J))*A
        CHIZ(I,J)=DZ(I,ICT(J))*B+DY(I,ICT(J))*A
        CHID2(I,J)=(SEVEN+SUM(J)*R2(I,ICT(J)))*B
442     CONTINUE
440     CONTINUE
C
C
C       FOR Fxxx-TYPE
C
        DO 901 J=ITP(10)+1,ITP(11)
        DO 902 I=1,NPTS
 
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DX(I,ICT(J))*DX(I,ICT(J))*A
C
        CHI(I,J)=B*DX(I,ICT(J))
        CHIX(I,J)=(THREE + SUM(J)*DX(I,ICT(J))*DX(I,ICT(J)))*B
        CHIY(I,J)=SUM(J)*DY(I,ICT(J))*CHI(I,J)
        CHIZ(I,J)=SUM(J)*DZ(I,ICT(J))*CHI(I,J)
        CHID2(I,J)=SIX*A*DX(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
 902    CONTINUE
 901    CONTINUE
C
C       FOR Fyyy-TYPE
C
        DO 911 J=ITP(11)+1,ITP(12)
        DO 912 I=1,NPTS
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DY(I,ICT(J))*DY(I,ICT(J))*A
C
        CHI(I,J)=B*DY(I,ICT(J))
        CHIX(I,J)=SUM(J)*DX(I,ICT(J))*CHI(I,J)
        CHIY(I,J)=(THREE + SUM(J)*DY(I,ICT(J))*DY(I,ICT(J)))*B
        CHIZ(I,J)=SUM(J)*DZ(I,ICT(J))*CHI(I,J)
        CHID2(I,J)=SIX*A*DY(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 912    CONTINUE
 911    CONTINUE
C
C       FOR Fzzz-TYPE
C
        DO 921 J=ITP(12)+1,ITP(13)
        DO 923 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        B=DZ(I,ICT(J))*DZ(I,ICT(J))*A
C
        CHI(I,J)=B*DZ(I,ICT(J))
        CHIX(I,J)=SUM(J)*DX(I,ICT(J))*CHI(I,J)
        CHIY(I,J)=SUM(J)*DY(I,ICT(J))*CHI(I,J)
        CHIZ(I,J)=(THREE + SUM(J)*DZ(I,ICT(J))*DZ(I,ICT(J)))*B
        CHID2(I,J)=SIX*A*DZ(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 923    CONTINUE
 921    CONTINUE
C
C       FOR Fxxy-TYPE
C
        DO 931 J=ITP(13)+1,ITP(14)
        DO 932 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BXY=DX(I,ICT(J))*DY(I,ICT(J))*A
        BXX=DX(I,ICT(J))*DX(I,ICT(J))*A
C
        CHI(I,J)=BXY*DX(I,ICT(J))
        CHIX(I,J)=(TWO + SUM(J)*DX(I,ICT(J))*DX(I,ICT(J)))*BXY
        CHIY(I,J)=(ONE + SUM(J)*DY(I,ICT(J))*DY(I,ICT(J)))*BXX
        CHIZ(I,J)=SUM(J)*DZ(I,ICT(J))*CHI(I,J)
        CHID2(I,J)=TWO*A*DY(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 932    CONTINUE
 931    CONTINUE
C
C       FOR Fxxz-TYPE
C
        DO 941 J=ITP(14)+1,ITP(15)
        DO 943 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BXZ=DX(I,ICT(J))*DZ(I,ICT(J))*A
        BXX=DX(I,ICT(J))*DX(I,ICT(J))*A
C
        CHI(I,J)=BXZ*DX(I,ICT(J))
        CHIX(I,J)=(TWO + SUM(J)*DX(I,ICT(J))*DX(I,ICT(J)))*BXZ
        CHIY(I,J)=SUM(J)*DY(I,ICT(J))*CHI(I,J)
        CHIZ(I,J)=(ONE + SUM(J)*DZ(I,ICT(J))*DZ(I,ICT(J)))*BXX
        CHID2(I,J)=TWO*A*DZ(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 943    CONTINUE
 941    CONTINUE
C
C       FOR Fyyz-TYPE
C
        DO 961 J=ITP(15)+1,ITP(16)
        DO 963 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BYZ=DZ(I,ICT(J))*DY(I,ICT(J))*A
        BYY=DY(I,ICT(J))*DY(I,ICT(J))*A
C
        CHI(I,J)=BYZ*DY(I,ICT(J))
        CHIX(I,J)=SUM(J)*DX(I,ICT(J))*CHI(I,J)
        CHIY(I,J)=(TWO + SUM(J)*DY(I,ICT(J))*DY(I,ICT(J)))*BYZ
        CHIZ(I,J)=(ONE + SUM(J)*DZ(I,ICT(J))*DZ(I,ICT(J)))*BYY
        CHID2(I,J)=TWO*A*DZ(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 963    CONTINUE
 961    CONTINUE
C
C       FOR Fxyy-TYPE
C
        DO 951 J=ITP(16)+1,ITP(17)
        DO 952 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BXY=DX(I,ICT(J))*DY(I,ICT(J))*A
        BYY=DY(I,ICT(J))*DY(I,ICT(J))*A
C
        CHI(I,J)=BXY*DY(I,ICT(J))
        CHIX(I,J)=(ONE + SUM(J)*DX(I,ICT(J))*DX(I,ICT(J)))*BYY
        CHIY(I,J)=(TWO + SUM(J)*DY(I,ICT(J))*DY(I,ICT(J)))*BXY
        CHIZ(I,J)=SUM(J)*DZ(I,ICT(J))*CHI(I,J)
        CHID2(I,J)=TWO*A*DX(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 952    CONTINUE
 951    CONTINUE
C
C       FOR Fxzz-TYPE
C
        DO 971 J=ITP(17)+1,ITP(18)
        DO 972 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BXZ=DZ(I,ICT(J))*DX(I,ICT(J))*A
        BZZ=DZ(I,ICT(J))*DZ(I,ICT(J))*A
C
        CHI(I,J)=BXZ*DZ(I,ICT(J))
        CHIX(I,J)=(ONE + SUM(J)*DX(I,ICT(J))*DX(I,ICT(J)))*BZZ
        CHIY(I,J)=SUM(J)*DY(I,ICT(J))*CHI(I,J)
        CHIZ(I,J)=(TWO + SUM(J)*DZ(I,ICT(J))*DZ(I,ICT(J)))*BXZ
        CHID2(I,J)=TWO*A*DX(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 972    CONTINUE
 971    CONTINUE
C
C       FOR Fyzz-TYPE
C
        DO 981 J=ITP(18)+1,ITP(19)
        DO 983 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BYZ=DZ(I,ICT(J))*DY(I,ICT(J))*A
        BZZ=DZ(I,ICT(J))*DZ(I,ICT(J))*A
C
        CHI(I,J)=BYZ*DZ(I,ICT(J))
        CHIX(I,J)=SUM(J)*DX(I,ICT(J))*CHI(I,J)
        CHIY(I,J)=(ONE + SUM(J)*DY(I,ICT(J))*DY(I,ICT(J)))*BZZ
        CHIZ(I,J)=(TWO + SUM(J)*DZ(I,ICT(J))*DZ(I,ICT(J)))*BYZ
        CHID2(I,J)=TWO*A*DY(I,ICT(J))+
     1             (NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 983    CONTINUE
 981    CONTINUE
C
C       FOR Fxyz-TYPE
C
        DO 991 J=ITP(19)+1,ITP(20)
        DO 992 I=1,NPTS
C
        A=DEXP(-EXX(J)*R2(I,ICT(J)))
        BXY=DX(I,ICT(J))*DY(I,ICT(J))*A
        BYZ=DZ(I,ICT(J))*DY(I,ICT(J))*A
        BXZ=DX(I,ICT(J))*DZ(I,ICT(J))*A
C
        CHI(I,J)=DX(I,ICT(J))*BYZ
        CHIX(I,J)=(ONE + SUM(J)*DX(I,ICT(J))*DX(I,ICT(J)))*BYZ
        CHIY(I,J)=(ONE + SUM(J)*DY(I,ICT(J))*DY(I,ICT(J)))*BXZ
        CHIZ(I,J)=(ONE + SUM(J)*DZ(I,ICT(J))*DZ(I,ICT(J)))*BXY
        CHID2(I,J)=(NINE+SUM(J)*R2(I,ICT(J)))*CHI(I,J)*SUM(J)
C
 992    CONTINUE
 991    CONTINUE
C
       DO 346 L=1,NMO
       DO 347 I=1,NPTS
C
       PSI(I,L)=ZERO
       GX(I,L)=ZERO
       GY(I,L)=ZERO
       GZ(I,L)=ZERO
C
347    CONTINUE
C
       DO 348 J=1,NPRIMS
       IF(COO(J,L).NE.ZERO)THEN
       DO 349 I=1,NPTS
C
       PSI(I,L)=PSI(I,L)+COO(J,L)*CHI(I,J)
       GX(I,L)=GX(I,L)+COO(J,L)*CHIX(I,J)
       GY(I,L)=GY(I,L)+COO(J,L)*CHIY(I,J)
       GZ(I,L)=GZ(I,L)+COO(J,L)*CHIZ(I,J)
C
349    CONTINUE
       ENDIF
348    CONTINUE
346    CONTINUE
C
        DO 301 I=1,NPTS
        GRX(I)=ZERO
        GRY(I)=ZERO
        GRZ(I)=ZERO
301     CONTINUE
C
        DO 305 L=1,NMO
        DO 306 I=1,NPTS
C
        GRX(I)=GRX(I)-PO(L)*GX(I,L)*PSI(I,L)
        GRY(I)=GRY(I)-PO(L)*GY(I,L)*PSI(I,L)
        GRZ(I)=GRZ(I)-PO(L)*GZ(I,L)*PSI(I,L)
C
306     CONTINUE
305     CONTINUE
C
        DO 309 I=1,NPTS
C
        DSQ=DSQRT(GRX(I)**2+GRY(I)**2+GRZ(I)**2)
C
        XYZ(I,1)=XYZ(I,1)+STEP(II)*GRX(I)/DSQ
        XYZ(I,2)=XYZ(I,2)+STEP(II)*GRY(I)/DSQ
        XYZ(I,3)=XYZ(I,3)+STEP(II)*GRZ(I)/DSQ
C
309     CONTINUE
C
304     CONTINUE
C
        Do 1593 J=1,NPRIMS
        temp=zero
        Do 1594 I=1,NPTS
        check=Dmax1(Dabs(chi(i,j)),dabs(chix(i,j)),dabs(chiy(i,j)),
     $        dabs(chiz(i,j)),dabs(chid2(i,j)),temp)
        If(check.gt.temp)temp=check
1594     Continue
         chimax(j)=temp
1593     continue
C
        DO 424 L=1,NMO
        DO 425 J = 1,NPRIMS
          IF(dabs(coo(j,l)*chimax(j)*rootpo(l)).lt.cutoff) Then
          ITEMP(J,L)=ITEMP(J,L)+1
          IDEL(J)=IDEL(J)+1
          ENDIF
425      CONTINUE
424      CONTINUE
303      CONTINUE
198      CONTINUE
C
         NPRIMSA=NPRIMS
         ICHK0=ITIMES*NCT
         ICHK=NMO*ITIMES*NCT
         NCNT=0
         DO 528 J=1,NPRIMS
C
         IF(IDEL(J).EQ.ICHK) THEN
         NPRIMSA=NPRIMSA-1
         ELSE
         NCNT=NCNT+1
         DO 529 L=1,NMO
         IF(ITEMP(J,L).EQ.ICHK0) THEN
         COOA(NCNT,L)=ZERO
         NZEROA=NZEROA+1
         ELSE
         COOA(NCNT,L)=COO(J,L)
         ENDIF
529      CONTINUE
         ITYPEA(NCNT)=ITYPE(J)
         EXXA(NCNT)=EXX(J)
         SUMA(NCNT)=SUM(J)
         DIVA(NCNT)=DIV(J)
         ICTA(NCNT)=ICT(J)
         ENDIF
528      CONTINUE
C
          N=0
         DO 631 K=1,NTYPE
         DO 632 J=1,NPRIMSA
         IF(ITYPEA(J).EQ.K) THEN
         N=N+1
         ENDIF
632      CONTINUE
         ITPA(K)=N
631      CONTINUE
C
       Do 833 L=1,NMO
       temp=zero
       Do 834 J=1,NPrimsA
       If(dabs(cooa(J,L)).gt.temp)temp=dabs(cooa(j,L))
834    continue
       cooamx(L)=temp
833    continue
C	
       Do 933 L=1,NMO
       temp=zero
       Do 934 J=1,NPrimsB
       If(dabs(coob(J,L)).gt.temp)temp=dabs(coob(j,L))
934    continue
       coobmx(L)=temp
933    continue
C	
999    CONTINUE
      NPR=NPRIMS
      NPRA=NPRIMSA
      NPRB=NPRIMSB
      IF(NPRA.GT.MPRIMS) THEN
      WRITE(IOUT,3301) NPRA
      STOP 'Too many primitives needed outside Beta sphere'
      ENDIF
      IF(NPRB.GT.MPRIMS) THEN
      WRITE(IOUT,3300) NPRB
      STOP 'Too many primitives needed inside Beta sphere'
      ENDIF
C
      WRITE(IOUT,970)
      Write(Iout,1971)Cutoff
      IF(NDOCUT.EQ.0)THEN
      WRITE(IOUT,962)
      ELSE
      WRITE(IOUT,969)
      WRITE(IOUT,964) NPR
      WRITE(IOUT,965) NPRA
      WRITE(IOUT,966) NPRB
      WRITE(IOUT,967) NZEROA, NMO*NPRA
      WRITE(IOUT,968) NZEROB, NMO*NPRB
      ENDIF
C
       RETURN
       END
      SUBROUTINE HESS(R1,R,H1)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL RHF,ROHF,ABNAT,RNAT
      PARAMETER (MMO=300,MCENT=50)
      COMMON /ORBTL/ EORB(MMO),PO(MMO),ROOTPO(MMO),NMO
      COMMON /ZZZZZ/ PSI(MMO),GX(MMO),GY(MMO),GZ(MMO),GXX(MMO),
     +GXY(MMO),GXZ(MMO),GYY(MMO),GYZ(MMO),GZZ(MMO)
      COMMON /C7/  CT(3,3), X, Y, Z, NCRNT
      COMMON /DIRK/ Toll,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,IMagM,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGO,YGO,ZGO
      DIMENSION W(3),R(3),H(3,3),R1(3),H1(3,3),XYZ(3)
      Save Zero,Two
      Data Zero/0.0d0/,Two/2.0d0/
C
      LMO=NMO
      If(idomag.eq.1)LMO=NMO/7
      R(1)=R1(1)
      R(2)=R1(2)
      R(3)=R1(3)
C
      R(1) = R(1) + X
      R(2) = R(2) + Y
      R(3) = R(3) + Z
C
      XYZ(1) = R(1)
      XYZ(2) = R(2)
      XYZ(3) = R(3)
C
      CALL GAUS(XYZ) 
      W(1)=Zero
      W(2)=Zero
      W(3)=Zero
      H(1,1)=Zero
      H(1,2)=Zero
      H(1,3)=Zero
      H(2,2)=Zero
      H(2,3)=Zero
      H(3,3)=Zero
C
      DO 20 I=1,LMO
      PROD=Two*PO(I)*PSI(I)
      W(1)=W(1)+GX(I)*PROD
      W(2)=W(2)+GY(I)*PROD
      W(3)=W(3)+GZ(I)*PROD
      H(1,1)=H(1,1)+PROD*GXX(I)+Two*PO(I)*GX(I)*GX(I)
      H(1,2)=H(1,2)+PROD*GXY(I)+Two*PO(I)*GX(I)*GY(I)
      H(1,3)=H(1,3)+PROD*GXZ(I)+Two*PO(I)*GX(I)*GZ(I)
      H(2,2)=H(2,2)+PROD*GYY(I)+Two*PO(I)*GY(I)*GY(I)
      H(2,3)=H(2,3)+PROD*GYZ(I)+Two*PO(I)*GY(I)*GZ(I)
      H(3,3)=H(3,3)+PROD*GZZ(I)+Two*PO(I)*GZ(I)*GZ(I)
20    CONTINUE
      DO 30 I=1,3
      DO 30 J=1,I
 30   H(I,J)=H(J,I)
      DO 1 I=1,3
      R(I)=W(I)
      DO 1 J=1,3
      H1(I,J)=H(I,J)
1     CONTINUE
      RETURN
      END
      SUBROUTINE INTARC(A,B,TH,PHI,FLINTS,IA,IT,AOM,ICT,IBETA)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL RHF,ROHF,ABNAT,RNAT
      PARAMETER (MCENT=50, MMO=300, MPTS=500,MaxTht=200,MaxPhi=200,
     $MaxPrp=200,MaxRad=200)
      COMMON /C1/  BASSIN(MaxPrp),TMP(MMO*(MMO+1)/2),naol
      COMMON /LIMITS/ RMAX, NATMX, NFLIM, NPROPS
      COMMON /ORBTL/ EORB(MMO),PO(MMO),ROOTPO(MMO),NMO
      COMMON /C5/ WP(MaxPhi),WTT(MaxTht)
      COMMON /ZZZZ/ PSI(MPTS,MMO),GX(MPTS,MMO),GY(MPTS,MMO),
     +GZ(MPTS,MMO),D2(MPTS,MMO)
      COMMON /DIRK/ Toll,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,IMagM,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGO,YGO,ZGO
      COMMON /NCUT/ CUTOFF,NDOCUT,NPR,NPRA,NPRB,NZEROA,NZEROB
      DIMENSION XA(MaxRad),WA(MaxRad),FU(MPTS,MaxPrp),
     $XXX(MaxTht*MaxRad),R(MaxTht*MaxRad,3),A(MaxTht),B(MaxTht),
     $TH(MaxTht),IA(MaxTht),ICT(MaxTht),FLINTS(MaxPrp),IVD(MPTS),
     $RR(MPTS,3),XX(MPTS),AOM(MMO*(MMO+1)/2),TEMP(MaxRad,3),
     $PSIMAX(MMO)
      SAVE ZERO,PT5,Hund
      DATA ZERO/0.0d0/,Pt5/0.5d0/,Hund/100.0d0/
C
      cutaom=cutoff*Hund
      LMO=NMO
      If(IDOMAg.eq.1)LMO=NMO/7
      KLA=0
      DO 2 I=1,IT
      S4=DSIN(TH(I))
      S1=S4*DCOS(PHI)
      S2=S4*DSIN(PHI)
      S3=DCOS(TH(I))
      H1=(B(I)+A(I))*PT5
      H2=(B(I)-A(I))*PT5
      CALL QUAD(TEMP,IA(I),XA,WA,H1,H2,MaxRad)
      DO 3 J=1,IA(I)
      KLA=KLA+1
      R(KLA,1)=XA(J)*S1
      R(KLA,2)=XA(J)*S2
      R(KLA,3)=XA(J)*S3
      XXX(KLA)=WA(J)*XA(J)*XA(J)*H2*WTT(ICT(I))
3     CONTINUE
2     CONTINUE
C
C     INTEGRATE THE IT RAYS IN BATCHES OF MPTS AND A REMAINDER (IMINUS)
C
      IDIV=KLA/MPTS
      IMINUS=KLA-MPTS*IDIV
      DO 5 J=1,IDIV
      IVD(J)=MPTS
5     CONTINUE 
      IVD(IDIV+1)=IMINUS
      ISEND=IDIV+1
      IF(IMINUS.EQ.0)ISEND=IDIV
C
      KK=0
      DO 8 JB=1,ISEND
      DO 9 J=1,IVD(JB)
      RR(J,1)=R(J+KK,1)
      RR(J,2)=R(J+KK,2)
      RR(J,3)=R(J+KK,3)
      XX(J)=XXX(J+KK) 
9     CONTINUE
      CALL FUNC3(RR,FU,IVD(JB),IBETA)
      DO 16 L=1,NProps
      DO 15 J=1,IVD(JB)
      FLINTS(L)=FLINTS(L)+FU(J,L)*XX(J)
15    CONTINUE
16    CONTINUE
      If(IDOAOM.Eq.1)Then
      Do 17 L=1,LMO
      temp1=zero
      Do 18 J=1,IVD(JB)
      temp2 = dmax1(dabs(psi(J,L)),temp1)
      If(temp2.gt.temp1)temp1=temp2
18    Continue
      psimax(L)=temp1
17    Continue
      wmax=zero
      Do 19 J=1,IVD(JB)
      If(dabs(xx(J)).gt.wmax)wmax=dabs(xx(j))
19    Continue
      NAOM=0
      DO 26 LA=1,LMO
      DO 27 LB=1,LA
      NAOM=NAOM+1
      AOMM=ZERO
      temp1=dabs(Psimax(LA)*Psimax(LB)*wmax)
      If(temp1.gt.cutaom)Then
      DO 25 J=1,IVD(JB)
      AOMM=AOMM+PSI(J,LA)*PSI(J,LB)*XX(J)
25    CONTINUE
      Endif
      AOM(NAOM)=AOM(NAOM)+AOMM
27    CONTINUE
26    CONTINUE
      Endif
      KK=KK+IVD(JB)
8     CONTINUE
C
      RETURN
      END
      SUBROUTINE INTEG(BETA,IPHIPL,IT,CRIT,N,NR,NC,NSRC,IMEGA)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL RHF,ROHF,ABNAT,RNAT
      INTEGER TOTINS
      PARAMETER (MMO=300,MaxTht=200,MaxPhi=200,MaxCrt=20,Mxbcrt=10,
     $MaxPrp=200,MCENT=50)
      COMMON /C1/  BASSIN(MaxPrp),AOM(MMO*(MMO+1)/2),NAOL
      COMMON /LIMITS/ RMAX, NATMX, NFLIM, NPROPS
      COMMON /BETINT/ IPHB,ITB
      COMMON/C5/WP(MaxPhi),WTT(MaxTht)
      COMMON /C7/  CT(3,3), X, Y, Z, NCRNT
      COMMON /S/ THETA(MaxTht),PHI(MaxPhi),SUR(3,MaxTht,MaxPhi)
      COMMON /UNITS/ ISRF,INPT,IOUT,IWFN,IDBG
      COMMON/PARA/NINS(Mxbcrt),IBETPTS,INMAX,NNN,NPATH,DIST,
     $DMAX,Amax,ADIFF,Tinf
      COMMON /VALUES/ THRESH1,THRESH2,GAMMA,TOTE
      COMMON/INSINF/ INSLIM(Mxbcrt)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /DIRK/ TOLL,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,NAcc,IDOAOM,IDOMAG,IMagm,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGO,YGO,ZGO
      COMMON /ORBTL/ EORB(MMO),PO(MMO),ROOTPO(MMO),NMO
      DIMENSION CRIT(3,MaxCrt),ANGLE(2,2,MxBcrt),TEMP(MaxTht,3),
     $NSRC(4,Mxbcrt)
      Save Zero,One,Two,Hund
      Data Zero/0.0d0/,One/1.0d0/,Two/2.0d0/,Hund/100.0d0/
865   FORMAT(1I6,' PHI AND ',1I6,' THETA PLANES IN BETA SPHERE')
870   FORMAT(1I6,' PHI AND ',1I6,' THETA PLANES OUTSIDE BETA SPHERE')
890   FORMAT(' RADIUS OF BETA SPHERE: ',1F8.4,' WITH ',1I6,' POINTS',
     +' PER RAY')
910   FORMAT(' VOL1 RHO CONTOUR THRESHOLD= ',1F8.4)
920   FORMAT(' VOL2 RHO CONTOUR THRESHOLD= ',1F8.4)
921   Format('Doing Beta Sphere Integration ...')
922   Format('Beta Sphere Integration is done ...')
927   Format(1I6,' Initial GradRho Trajectories Per Interatomic',
     $' Surface')
928   Format(1I6,' Points per GradRho Surface Trajectory')
929   Format('Max. Dist. Between Ends of Adjacent',
     $' GradRho Surface Trajectories = ',1PE8.2)
930   FORMAT(' INSERTION LIMIT USED          = ',1I4)
931   Format('Maximum Length of GradRho Surface Trajectories = ',
     $1PE8.2)
935   Format('Integrate Out to ',1PE8.2,' For Rays Intersecting',
     $' Surface at Infinity')
940   FORMAT(' INSERTION LIMIT REACHED ',1I6,' TIMES FOR ',
     +'SURFACE ',1I3)
950   FORMAT(' FOR SURFACE #',I2,' NUMBER OF INSERTED ',
     +'PATHS = ',1I6)
960   FORMAT(' TOTAL NUMBER OF INSERTED PATHS= ',1I5)
1999  Format('Doing Proaim Surface Routine ...')
2000  Format('Surface is done ...')
2999  Format('Doing Promega Surface Routine ...')
3001  Format('Order of ABM to use in Solving for GradRho',
     $' Basin Trajectories = ',1I1)
3002  Format(1I3,' Small Steps per regular Step in ABM')
3003  Format('Step Size used in in ABM = ',1PE8.2)
3004  Format('Maximum Distance from Nucleus to Search for',
     $' Intersections = ',1PE8.2)
3005  Format('Max. Dist. between Calculated and Actual Points',
     $' of Intersection = ',1PE8.2)
3006  Format('Initial Step Size Along Rays In First Intersection', 
     $' Search = ',1PE8.2)
3007  Format('Initial Step Size Along Rays In Second and Third', 
     $' Intersection Search = ',1PE8.2)
3999  Format('Radius of Capture Sphere for Atom ',I6,' is ',1PE8.2)
4000  Format('Search for first intersections only in intersec')
4100  Format('Search for first,second and third intersections 
     $in intersec')
4200  Format('Number of Second Intersections with Surface: ',1I6)
4300  Format('Number of Oscillations in Trudge: ',1I6)
5000  Format('Doing Integration Outside of Beta Sphere ...')
5001  Format('Default number of radial points used')
5002  Format(1I2,' x the default number of radial points used')
5003  Format('Default number of theta and phi planes used for Beta',
     $' Sphere')
5004  Format(1I2,' x default number of theta and phi planes used',
     $' for Beta Sphere')
5005  Format('Number of Theta Values Requested Exceeds Maximum.', 
     $ ' Will Use Maximum ')
5006  Format('Number of Phi Values Requested Exceeds Maximum.', 
     $ ' Will Use Maximum ')
6000  Format('Integration Outside of Beta Sphere is Done ...')
C
      lmo=nmo
      If(idomag.eq.1)lmo=nmo/7
      naol=lmo*(lmo+1)/2
C
      NT=N+NR+NC
      Beta=Zero
      IF(IMEGA.EQ.0)Then
      Beta=Hund
      Do 100 I=1,NT
      Tmp=Dsqrt((Crit(1,I)-X)**2+(Crit(2,I)-Y)**2+
     $(Crit(3,I)-Z)**2)
      If(Tmp.Lt.beta)Beta=Tmp
100   Continue
      If(Beta.eq.Hund.and.Ncent.eq.1)Beta=One
      Endif
      CALL GAUSCHECK(BETA,NT,N,CRIT)
      IF(IMEGA.EQ.1)Call SPhere(Beta)
C
      WRITE(IOUT,890) BETA, IBETPTS
C
      DO 105 I=1,NT
      CRIT(1,I)=CRIT(1,I)-X
      CRIT(2,I)=CRIT(2,I)-Y
105   CRIT(3,I)=CRIT(3,I)-Z
      DO 110 I=1,NPROPS
110   BASSIN(I)=Zero
      DO 115 LA=1,NAOL
115   AOM(LA)=Zero
C
      CALL EIGEN(CRIT,NT)
C
      CALL RING(CRIT,N,NR,ANGLE,NSRC)
C
      IF (BETA.LE.(0.2D0))THEN
      IPHB=4*NACC
      ITB=2*NACC
      ELSEIF(BETA.GT.(0.2D0).AND.BETA.LE.(0.4D0))THEN
      IPHB=6*NACC
      ITB=4*NACC
      ELSEIF(BETA.GT.(0.4D0).AND.BETA.LE.(0.6D0))THEN
      IPHB=10*NACC
      ITB=6*NACC
      ELSEIF(BETA.GT.(0.6D0).AND.BETA.LE.(0.8D0))THEN
      IPHB=16*NACC
      ITB=10*NACC
      ELSEIF(BETA.GT.(0.8D0).AND.BETA.LE.(1.0D0))THEN
      IPHB=22*NACC
      ITB= 14*NACC
      ELSEIF(BETA.GT.(1.0D0).AND.BETA.LE.(1.3D0))THEN
      IPHB=28*NACC
      ITB=18*NACC
      ELSEIF(BETA.GT.(1.3D0).AND.BETA.LE.(1.6D0))THEN
      IPHB=36*NACC
      ITB=24*NACC
      ELSEIF(BETA.GT.(1.6D0).AND.BETA.LE.(2.0D0))THEN
      IPHB=48*NACC
      ITB=32*NACC
      ELSEIF(BETA.GT.(2.0D0))THEN
      IPHB=64*NACC
      ITB=48*NACC
      ENDIF
C
      If(Nacc.eq.1)Write(Iout,5003)
      If(Nacc.gt.1)Write(Iout,5004)Nacc
C
      If(ITB.GT.MaxTht)Then
      Write(Iout,5005)MaxTht
      ITB=MaxTht
      Endif
C
      If(IPhB.GT.MaxPhi)Then
      Write(Iout,5006)MaxPhi
      IPhB=MaxPhi
      Endif
C
      WRITE(IOUT,865) IPHB,ITB
      WRITE(IOUT,870) IPHIPL,IT
      WRITE(IOUT,935) TINF
      WRITE(IOUT,910) THRESH1
      WRITE(IOUT,920) THRESH2
C
      PI=DACOS(-ONE)
      PHINC=TWO*PI/IPHB
      HPHINC=PHINC/TWO
      PHWGT=TWO*PI/IPHB
      DO 120 I=1,IPHB
      PHI(I)=I*PHINC-HPHINC
      WP(I)=PHWGT
120   CONTINUE
C
      DMULT=PI/TWO
      CALL QUAD(TEMP,ITB,THETA,WTT,DMULT,DMULT,MaxTht)
C
      DO 125 I=1,ITB
      WTT(I)=WTT(I)*DSIN(THETA(I))*DMULT
125   CONTINUE
C
      IBETA=1
      Write(iout,921)
      CALL PHINT(IBETA,BETA,BASSIN,IPHB,ITB,AOM)
      Write(iout,922)
C
      PI=DACOS(-ONE)
      PHINC=TWO*PI/IPHIPL
      HPHINC=PHINC/TWO
      PHWGT=TWO*PI/IPHIPL
      DO 130 I=1,IPHIPL
      PHI(I)=I*PHINC-HPHINC
      WP(I)=PHWGT
130   CONTINUE
C
      DMULT=PI/TWO
      CALL QUAD(TEMP,IT,THETA,WTT,DMULT,DMULT,MaxTht)
C
      DO 135 I=1,IT
      WTT(I)=WTT(I)*DSIN(THETA(I))*DMULT
135   CONTINUE
C
      If(IMEGA.EQ.0)THEN
      WRITE(IOUT,1999)
      WRITE(IOUT,927)NPATH
      WRITE(IOUT,928)NNN
      WRITE(IOUT,929)Dist
      WRITE(IOUT,931)DMAX
      WRITE(IOUT,930)INMAX
      CALL SURFACE(CRIT,N,NR,NC,ANGLE,NSRC,IPHIPL,IT)
      TOTINS = 0
      WRITE(IOUT,940) (INSLIM(II),II,II=1,N)
      WRITE(IOUT,950) (II,NINS(II),II=1,N)
      DO 140 II=1,N
      TOTINS = TOTINS+NINS(II)
140   CONTINUE
      WRITE(IOUT,960) TOTINS
      WRITE(IOUT,2000)
C
      ELSEIF(IMEGA.EQ.1)THEN
      WRITE(IOUT,2999)
      If(Nosec.eq.1)Write(Iout,4000)
      If(Nosec.eq.0)Write(Iout,4100)
      WRITE(IOUT,3001) NABMO
      WRITE(IOUT,3002) ISTEP
      Write(Iout,3003)Stp
      Write(Iout,3004)Size
      Write(Iout,3005)Ctf
      Write(Iout,3006)FStp
      If(Nosec.eq.0)Write(Iout,3007)SStp
      DO 145 IA=1,NCent
      Write(Iout,3999)IA,Dsqrt(toll2(ia))
145   Continue
      CALL INTERSEC(BETA,IPHIPL,IT)
      If(NOSEC.eq.0)Write(4200) ISECT
      Write(Iout,4300)ICP
      WRITE(IOUT,2000)
      ENDIF
C
      IBETA=0
      Write(Iout,5000)
      If(Nacc.eq.1)Write(Iout,5001)
      If(Nacc.gt.1)Write(Iout,5002)Nacc
      CALL PHINT(IBETA,BETA,BASSIN,IPHIPL,IT,AOM)
      Write(Iout,6000)
C
      RETURN
      END
      SUBROUTINE INTERSEC(BETA,IPHI,ITHETA)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL RHF,ROHF,ABNAT,RNAT
      PARAMETER (MaxTht=200,MaxPhi=200,Mxbcrt=10)
      PARAMETER (MCENT=50, MNPT=MaxTht*MaxPhi)
      COMMON /C7/ CT(3,3),X, Y, Z,NCRNT
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /S/ THETA(MaxTht),PHI(MaxPhi),SUR(3,MaxTht,MaxPhi)
      COMMON /DIRK/ TOLL,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,Imagm,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGO,YGO,ZGO
      COMMON /UNITS/ ISRF,INPT,IOUT,IWFN,IDBG
      COMMON/PARA/ NINS(Mxbcrt),IBETPTS,INMAX,NNN,NPATH,DSTP,
     $DMAX,AMax,ADIFF,Tinf
      DIMENSION ORIGIN(3), DELXYZ(MNPT,3),NGG(MNPT),XYZ(MNPT,3),
     $XYZT(MNPT,3),DISTF(MNPT),STEP(MNPT), RINC(MNPT),ITERMT(MNPT),
     $LFOR(MNPT),ITERM(MNPT), LBACK(MNPT), LDONE(MNPT),DIST(MNPT)
      Save Zero,One,Two,Four,IPhat,IThat,Scale,Dinc
      Data Zero/0.0d0/,One/1.0d0/,Two/2.0d0/,IPhat/16/,IThat/8/,
     $Scale/0.8d0/,Dinc/0.1d0/,FOur/4.0d0/
1999  Format('Found First Intersections of ',I6,
     $' Feeler Rays Out of ', I6, ' Total Rays ')
2999  Format(I6, ' First Intersections Remain to be Found ')
3999  Format('Searching For Second Intersections')
4999  Format('Searching For Third Intersections')
C
      Pi=DACOS(-ONE)
C
      ICOUNT=0
      ORIGIN(1) = X
      ORIGIN(2) = Y
      ORIGIN(3) = Z
C
      N=0
C
      IF(NOSEC.EQ.0)GOTO 777
C
      DO 10 I=1,(ITHETA-3),4
      DO 10 J=1,(IPHI-3),4
      N=N+1
      LDONE(N)=0
      LBACK(N)=0
      LFOR(N)=0
      DELXYZ(N,1) = DSIN(THETA(I))*DCOS(PHI(J))
      DELXYZ(N,2) = DSIN(THETA(I))*DSIN(PHI(J))
      DELXYZ(N,3) = DCOS(THETA(I))
10    CONTINUE
C
      NPOINTS = N
C
      FWR=Four
      DO 55 I =1,NPOINTS
      RINC(I) = FWR*FSTP
      STEP(I) = BETA+FWR*FSTP
55    CONTINUE
C
      DO 51 J=1,3
      DO 51 I=1,NPOINTS
      XYZ(I,J) = DELXYZ(I,J)*(STEP(I)) + ORIGIN(J)
51    CONTINUE
C
57    CONTINUE
C
      KK = 0
      DO 58 I =1,NPOINTS
      IF (LDONE(I).EQ.0) Then
      KK = KK + 1
      XYZT(KK,1) = XYZ(I,1)
      XYZT(KK,2) = XYZ(I,2)
      XYZT(KK,3) = XYZ(I,3)
      ENdif
58    CONTINUE
C
      CALL TRUDGE3(XYZT,KK,ITERMT)
C
      N = 1
      DO 59 I=1,NPOINTS
      IF (LDONE(I) .EQ. 0) Then
      ITERM(I) = ITERMT(N)
      N = N + 1
      ENdif
59    CONTINUE
C
      DO 60 I=1,NPOINTS
      IF (LDONE(I) .EQ. 1) Goto 60
C
      IF(((STEP(I)).GE.SIZE))THEN
      DISTF(I) = Size
      LDONE(I) = 1
      ICOUNT = ICOUNT + 1
      GOTO 60
      ENDIF
      IF(RINC(I).LE.CTF)THEN
      DISTF(I) = STEP(I)
      LDONE(I) = 1
      ICOUNT = ICOUNT + 1
      GOTO 60
      ENDIF
C
      IF ((ITERM(I) .EQ. 1) .AND. (LBACK(I) .EQ. 0)) THEN
      LFOR(I)=1 
      STEP(I) = STEP(I) + RINC(I)
      DO 110 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(STEP(I)) + ORIGIN(J)
110   CONTINUE
      GOTO 60
      ENDIF
C
      IF(ITERM(I) .EQ. 0)  THEN
      IF(LFOR(I).EQ.1) RINC(I) = RINC(I)/Two
      STEP(I) = STEP(I) - RINC(I)
      LBACK(I) = 1
      DO 120 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(STEP(I)) + ORIGIN(J)
120   CONTINUE
      GOTO 60
      ENDIF
C 
      IF ((ITERM(I) .EQ. 1) .AND. (LBACK(I) .EQ. 1)) THEN
      LFOR(I)=1 
      RINC(I) = RINC(I)/Two
      STEP(I) = STEP(I) + RINC(I)
      DO 130 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(STEP(I)) + ORIGIN(J)
130   CONTINUE
      GOTO 60
      ENDIF
C
60    CONTINUE
C
      IF (ICOUNT .EQ. NPOINTS) GOTO 88
C
      GOTO 57
C
88    CONTINUE
C
      Write(Iout,1999)NPOINTS,ITHETA*IPHI
777   CONTINUE
C
      NPOINTS = IPHI*ITHETA
C
      DO 550 I =1,NPOINTS
      LBACK(I) = 0
      LFOR(I) = 0
      RINC(I) = FSTP
      LDONE(I) = 0
      STEP(I)=BETA+FSTP
550   CONTINUE
C
      IF(NOSEC.EQ.1)THEN
C
      N=0
      DO 195 I=1,(ITHETA-3),4
      DO 196 J=1,(IPHI-3),4
      N=N+1
      MAR=(I-1)*IPHI+J
      LDONE(MAR)=1
      DIST(MAR)=DISTF(N)
196   CONTINUE
195   CONTINUE
C
      ENDIF
C
      N=1
      DO 100 I=1,ITHETA
      DO 100 J=1,IPHI
C
      DELXYZ(N,1) = DSIN(THETA(I))*DCOS(PHI(J))
      DELXYZ(N,2) = DSIN(THETA(I))*DSIN(PHI(J))
      DELXYZ(N,3) = DCOS(THETA(I))
      N=N+1
C
100   CONTINUE
C
      IF(NOSEC.EQ.1)THEN
      N=0
      DO 2123 I=1,ITHETA
      DO 2124 J=1,IPHI
      N=N+1
      IF(LDONE(N).NE.1)THEN
      STEP(N)=Zero
      PROJS=PI
      MM=0
      DO 2125 K=1,(ITHETA-3),4
      DO 2126 L=1,(IPHI-3),4
      M=(K-1)*IPHI+L
      PROJ=(DELXYZ(M,1)*DELXYZ(N,1)+DELXYZ(M,2)*DELXYZ(N,2)
     +    +DELXYZ(M,3)*DELXYZ(N,3))
      If(proj.ge.one)Then
      proj=Zero
      ElseIf(-proj.Ge.One)Then
      proj=pi
      Else
      proj=Dacos(proj)
      Endif
      IF(PROJ.LE.PROJS)THEN
      IF(PROJ.LT.PROJS)THEN
      STEP(N)=DIST(M)
      PROJS=PROJ
      MM=1
      ELSEIF(PROJ.EQ.PROJS) THEN
      STEP(N)=STEP(N)+DIST(M)
      MM=MM+1
      ENDIF
      ENDIF
2126  CONTINUE
2125  CONTINUE
      STEP(N)=(STEP(N)/MM)-FSTP
      ENDIF
2124  CONTINUE
2123  CONTINUE
C
      ENDIF
C
      DO 501 J=1,3
      DO 501 I=1,NPOINTS
      IF (LDONE(I).EQ.1)GOTO 501
      XYZ(I,J) = DELXYZ(I,J)*(STEP(I)) + ORIGIN(J)
501   CONTINUE
C
575   CONTINUE
      Write(Iout,2999)NPOINTS-ICount
      KK = 0
      DO 580 I =1,NPOINTS
      IF (LDONE(I) .EQ. 0) Then
      KK = KK + 1
      XYZT(KK,1) = XYZ(I,1)
      XYZT(KK,2) = XYZ(I,2)
      XYZT(KK,3) = XYZ(I,3)
      Endif
580   CONTINUE
C
      CALL TRUDGE3(XYZT,KK,ITERMT)
C
      N = 1
      DO 590 I=1,NPOINTS
      IF (LDONE(I) .EQ. 0) Then
      ITERM(I) = ITERMT(N)
      N = N + 1
      Endif
590   CONTINUE
C
      DO 600 I=1,NPOINTS
      IF (LDONE(I) .EQ. 1) GOTO 600
C
      IF(((STEP(I)).GE.SIZE))THEN
	DIST(I) = Size
        LDONE(I) = 1
	ICOUNT = ICOUNT + 1
        GOTO 600
	ENDIF
      IF(RINC(I).LE.CTF)THEN
	DIST(I) = STEP(I)
        LDONE(I) = 1
	ICOUNT = ICOUNT + 1
        GOTO 600
	ENDIF
C
      IF ((ITERM(I) .EQ. 1) .AND. (LBACK(I) .EQ. 0)) THEN
      LFOR(I)=1
      STEP(I) = STEP(I) + RINC(I)
      DO 1100 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(STEP(I)) + ORIGIN(J)
1100  CONTINUE
      GOTO 600
      ENDIF
C
      IF(ITERM(I) .EQ. 0)  THEN
      IF(LFOR(I).EQ.1) RINC(I) = RINC(I)/Two
      STEP(I) = STEP(I) - RINC(I)
      LBACK(I) = 1
      DO 1200 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(STEP(I)) + ORIGIN(J)
1200  CONTINUE
      GOTO 600
      ENDIF
C 
      IF ((ITERM(I) .EQ. 1) .AND. (LBACK(I) .EQ. 1)) THEN
      LFOR(I)=1
      RINC(I) = RINC(I)/Two
      STEP(I) = STEP(I) + RINC(I)
      DO 1300 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(STEP(I)) + ORIGIN(J)
1300  CONTINUE
      GOTO 600
      ENDIF
C
600   CONTINUE
C
      IF (ICOUNT .EQ. NPOINTS) GOTO 888
      GOTO 575
C
888   CONTINUE
C
      DO 890 I=1,ITHETA
      DO 890 J=1,IPHI
      SUR(1,I,J) = TINF
      SUR(2,I,J)=-One
      SUR(3,I,J)=TINF
890   CONTINUE
C
      N = 1
      DO 900 I=1,ITHETA
      DO 900 J=1,IPHI
      SUR(1,I,J) = DIST(N) 
      If(Sur(1,I,J).GE.Size)Sur(1,I,J)=Tinf
      N = N + 1
900   CONTINUE
C
      IF(NOSEC.EQ.1) GOTO 9999
C
      WRITE(IOUT,3999)
      ICOUNT=0
      DO 601 I=1,NPOINTS
      IF(DIST(I) .GE. Size) THEN
      LDONE(I) = 1
      NGG(I)=1
      ICOUNT=ICOUNT+1
      ELSE
      LDONE(I) = 0
      LBACK(I) = 0
      RINC(I) = SSTP 
      STEP(I) = SSTP 
      DO 602 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(DIST(I)+RINC(I)) + ORIGIN(J)
602   CONTINUE
      ENDIF
601   CONTINUE
C
      IF (ICOUNT .EQ. NPOINTS) GOTO 9999
C
675   CONTINUE
      KK = 0
      DO 680 I =1,NPOINTS
      IF (LDONE(I) .EQ. 0) Then
      KK = KK + 1
      XYZT(KK,1) = XYZ(I,1)
      XYZT(KK,2) = XYZ(I,2)
      XYZT(KK,3) = XYZ(I,3)
      Endif
680   CONTINUE
C
      CALL TRUDGE3(XYZT,KK,ITERMT)
C
      N = 1
      DO 690 I=1,NPOINTS
      IF (LDONE(I) .EQ. 0) Then
      ITERM(I) = ITERMT(N)
      N = N + 1
      Endif
690   CONTINUE
C
      DO 700 I=1,NPOINTS
      IF (LDONE(I) .EQ. 1) GOTO 700
C
      IF((DIST(I)+STEP(I)).GE.SIZE) THEN
      LDONE(I)=1
      ICOUNT=ICOUNT+1
      NGG(I)=1
      GOTO 700
      ENDIF
C
      IF (RINC(I) .LE. CTF) THEN 
      DIST(I) = DIST(I) + STEP(I)
      LDONE(I) = 1
      ICOUNT = ICOUNT + 1
      ISECT=ISECT+1
      GOTO 700
      ENDIF
C
      IF ((ITERM(I) .EQ. 0) .AND. (LBACK(I) .EQ. 0)) THEN
      RINC(I) = SSTP
      STEP(I) = STEP(I) + RINC(I)
      DO 2100 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(DIST(I)+STEP(I)) + ORIGIN(J)
2100  CONTINUE
      GOTO 700
      ENDIF
C
      IF ((ITERM(I) .EQ. 1)) THEN
      LBACK(I)=1
      RINC(I) = RINC(I)/Two
      STEP(I) = STEP(I) - RINC(I)
      DO 2300 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(DIST(I)+STEP(I)) + ORIGIN(J)
2300  CONTINUE
      GOTO 700
      ENDIF
C
      IF ((ITERM(I) .EQ. 0) .AND. (LBACK(I) .EQ. 1)) THEN
      RINC(I) = RINC(I)/Two
      STEP(I) = STEP(I) + RINC(I)
      DO 3100 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(DIST(I)+STEP(I)) + ORIGIN(J)
3100  CONTINUE
      GOTO 700
      ENDIF
C
700   CONTINUE
C
      IF (ICOUNT .NE. NPOINTS) GOTO 675
C
      N=1
      DO 9900 I=1,ITHETA
      DO 9910 J=1,IPHI
      IF(NGG(N).EQ.0) THEN
      SUR(2,I,J)=DIST(N)
      ENDIF
      N = N + 1
9910   CONTINUE
9900   CONTINUE
C
      ICOUNT=0
      DO 1601 I=1,NPOINTS
      IF(NGG(I) .EQ. 1) THEN
      LDONE(I) = 1
      ICOUNT=ICOUNT+1
      ELSE
      LDONE(I) = 0
      LBACK(I) = 0
      RINC(I) = SSTP 
      STEP(I) = SSTP 
      DO 1602 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(DIST(I)+RINC(I)) + ORIGIN(J)
1602  CONTINUE
      ENDIF
1601  CONTINUE
C
      IF (ICOUNT .EQ. NPOINTS) GOTO 9999
C
      WRITE(IOUT,4999)
1675   CONTINUE
      KK = 0
      DO 1680 I =1,NPOINTS
      IF (LDONE(I) .EQ. 0) Then
      KK = KK + 1
      XYZT(KK,1) = XYZ(I,1)
      XYZT(KK,2) = XYZ(I,2)
      XYZT(KK,3) = XYZ(I,3)
      Endif
1680   CONTINUE
C
      CALL TRUDGE3(XYZT,KK,ITERMT)
C
      N = 1
      DO 1690 I=1,NPOINTS
      IF (LDONE(I) .EQ. 0) Then
      ITERM(I) = ITERMT(N)
      N = N + 1
      Endif
1690   CONTINUE
C
      DO 1700 I=1,NPOINTS
      IF (LDONE(I) .EQ. 1) GOTO 1700
C
      IF((DIST(I)+STEP(I)).GE.SIZE) THEN
      LDONE(I)=1
      ICOUNT=ICOUNT+1
      DIST(I)=Size
      GOTO 1700
      ENDIF
C      
      IF (RINC(I) .LE. CTF) THEN 
      DIST(I) = DIST(I) + STEP(I)
      LDONE(I) = 1
      ICOUNT = ICOUNT + 1
      GOTO 1700
      ENDIF
C
      IF ((ITERM(I) .EQ. 1) .AND. (LBACK(I) .EQ. 0)) THEN
      RINC(I) = SSTP
      STEP(I) = STEP(I) + RINC(I)
      DO 12100 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(DIST(I)+STEP(I)) + ORIGIN(J)
12100  CONTINUE
      GOTO 1700
      ENDIF
C
      IF ((ITERM(I) .EQ. 0)) THEN
      LBACK(I)=1
      RINC(I) = RINC(I)/Two
      STEP(I) = STEP(I) - RINC(I)
      DO 12300 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(DIST(I)+STEP(I)) + ORIGIN(J)
12300  CONTINUE
      GOTO 1700
      ENDIF
C
      IF ((ITERM(I) .EQ. 1) .AND. (LBACK(I) .EQ. 1)) THEN
      RINC(I) = RINC(I)/Two
      STEP(I) = STEP(I) + RINC(I)
      DO 13100 J=1,3
      XYZ(I,J) = DELXYZ(I,J)*(DIST(I)+STEP(I)) + ORIGIN(J)
13100  CONTINUE
      GOTO 1700
      ENDIF
C
1700   CONTINUE
C
       IF (ICOUNT .NE. NPOINTS) GOTO 1675
C
      N=1
      DO 19900 I=1,ITHETA
      DO 19910 J=1,IPHI
      IF(NGG(N).EQ.0) THEN
      SUR(3,I,J)=DIST(N)
      If(Sur(3,I,J).GE.Size)Sur(3,I,J)=Tinf
      ENDIF
      N = N + 1
19910  CONTINUE
19900  CONTINUE
C
9999  Continue
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
      COMMON /C2/AB(9,9),AM(9,10)
      DIMENSION W1(3,I),W(3),R(3),R1(3)
      Save Zero
      Data Zero/0.0d0/
C
      DO 2 L=1,3
      H=Zero
      DO 106 J=1,I
106   H=H+AB(I,J)*W1(L,J)
      R(L)=W(L)+H*DS
2     R1(L)=R(L)
C
      CALL GAUS2(R1,0) 
      G=DSQRT(R1(1)**2+R1(2)**2+R1(3)**2)
      DO 5 L=1,3
      H=AM(I,1)*R1(L)/G
      DO 103 J=1,I
103   H=H+AM(I,J+1)*W1(L,J)
5     R(L)=W(L)+H*DS
      RETURN
      END
      SUBROUTINE PATH(IBETA,BEDA,JJ,FLINTS,IT,IPHIPL,AOMF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL RHF,ROHF,ABNAT,RNAT
      PARAMETER (MMO=300,MaxTht=200,MaxPhi=200,MaxPrp=200,
     $MxbCrt=10,Mcent=50,MaxRad=200)
      COMMON /C1/  BASSIN(MaxPrp),TMP(MMO*(MMO+1)/2),NAOL
      COMMON /LIMITS/ RMAX, NATMX, NFLIM, NPROPS
      COMMON /PARA/NINS(MxbCrt),IBETPTS,INMAX,NNN,NPATH,DIST,
     $DMAX,Amax,ADIFF,Tinf
      COMMON /DIRK/ Toll,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,IMagm,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGO,YGO,ZGO
      COMMON /S/ THETA(MaxTht),PHI(MaxPhi),SUR(3,MaxTht,MaxPhi)
      COMMON /UNITS/ ISRF,INPT,IOUT,IWFN,IDBG
      DIMENSION SURF(MaxTht),D1(MaxTht),D2(MaxTht),IA(MaxTht),
     $IA1(MaxTht),IA2(MaxTht),BETA(MAxTht),TZERO(MaxTht),
     $FLINTS(MaxPrp),THT(MaxTht),ICT(MaxTht),SU2(MaxTht),
     $SU3(MaxTht),IAA2(MaxTht),AOMF(MMO*(MMO+1)/2)
      Save Zero,One,Ptoone,Two
      Data Zero/0.0d0/,One/1.0d0/,Ptoone/0.01d0/,Two/2.0d0/
1999  Format('Requested Number of Radial Points Exceeds Maximum.',
     $'  Will Use Maximum.')
C
      DO 100 L=1,NProps
      FLINTS(L) = Zero
100   CONTINUE
      DO 105 LA=1,NAOL
      AOMF(LA)=Zero
105   CONTINUE
C
      IF(IBETA.EQ.1)GOTO 555
C
      DO 1 II=1,IT
      SURF(II)=SUR(1,II,JJ)
      BETA(II)=BEDA
      IF(SURF(II).LT.(Tinf-Ptoone))GOTO 6
      IF(II.LE.1.OR.II.GE.IT)GOTO 5
      IF(DABS(SURF(II)-SUR(1,II-1,JJ)).LE.(Tinf/2).OR.
     1 DABS(SURF(II)-SUR(1,II+1,JJ)).LE.(Tinf/2))GOTO 5
      SURF(II)=SUR(1,II-1,JJ)+(SUR(1,II+1,JJ)-SUR(1,II-1,JJ))*
     1 (THETA(II)-THETA(II-1))/(THETA(II+1)-THETA(II-1))
      GOTO 6
5     IF(JJ.LE.1.OR.JJ.GE.IPHIPL)GOTO 6
      IF(DABS(SURF(II)-SUR(1,II,JJ-1)).LE.(Tinf/2).OR.
     1 DABS(SURF(II)-SUR(1,II,JJ+1)).LE.(Tinf/2))GOTO 6
      SURF(II)=SUR(1,II,JJ-1)+(SUR(1,II,JJ+1)-SUR(1,II,JJ-1))*
     1 (PHI(JJ)-PHI(JJ-1))/(PHI(JJ+1)-PHI(JJ-1))
6     CONTINUE
      IF(SUR(2,II,JJ)-SURF(II).LT.Ptoone.AND.
     $SUR(3,II,JJ).GT.(Tinf-Ptoone))SUR(2,II,JJ)=-One
1     CONTINUE
C
C
      DO 17 II=1,IT
      D1(II)=SURF(II)-BEDA
      D2(II)=SUR(3,II,JJ)-SUR(2,II,JJ)
17     CONTINUE
C
      DO 9 II=1,IT
      IF(D1(II).LT.0.1)IA1(II)=NACC*8
      IF(D1(II).GE.0.1.AND.D1(II).LT.2.0)
     $IA1(II)=NACC*(8+2*IDINT(10*D1(II)))
      IF(D1(II).GE.2.0.AND.D1(II).LT.3.0)IA1(II)=NACC*54
      IF(D1(II).GE.3.0)IA1(II)=NACC*64
      IF(D2(II).LT.0.1)IA2(II)=NACC*8
      IF(D2(II).GE.0.1.AND.D2(II).LT.2.0)
     $IA2(II)=NACC*(8+2*IDINT(10*D2(II)))
      IF(D2(II).GE.2.0.AND.D2(II).LT.3.0)IA2(II)=NACC*54
      IF(D2(II).GE.3.0)IA2(II)=NACC*64
      ICT(II)=II
      If(IA1(II).gt.MaxRad)Then
      Write(Iout,1999)
      IA1(II)=MaxRad
      Endif
      If(IA2(II).gt.MaxRad)Then
      IA2(II)=MaxRad
      Write(Iout,1999)
      Endif
9     CONTINUE
C
555   CONTINUE
C
      IF(IBETA.EQ.1)THEN
      DO 18 II=1,IT
      IA(II) = IBETPTS
      TZERO(II) = Zero
      BETA(II)=BEDA
      ICT(II)=II
18     CONTINUE
      CALL INTARC(TZERO,BETA,THETA,PHI(JJ),FLINTS,IA,IT,AOMF,ICT,IBETA)
      GOTO 999
      ENDIF
C
      IF(IBETA.EQ.0)THEN  
      CALL INTARC(BETA,SURF,THETA,PHI(JJ),FLINTS,IA1,IT,AOMF,ICT,IBETA)
C
      K=0
      DO 10 II=1,IT
      IF(SUR(2,II,JJ).GT.SUR(1,II,JJ)) Then
      K=K+1
      ICT(K)=II
      SU2(K)=SUR(2,II,JJ)
      SU3(K)=SUR(3,II,JJ)
      THT(K)=THETA(II)
      IAA2(K)=IA2(II)
      Endif
10    CONTINUE
C
      IF (K .GT. 0) THEN
      CALL INTARC(SU2,SU3,THT,PHI(JJ),FLINTS,IAA2,K,AOMF,ICT,IBETA)
      ENDIF
      ENDIF 
C
999   RETURN
      END
      SUBROUTINE PHINT(IBETA,BETA,FINTPH,IPHIPL,IT,AOMT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MMO=300,MaxPhi=200,MaxTht=200,MaxPrp=200)
      COMMON /C5/ WP(MaxPhi),WTT(MaxTht)
      COMMON /LIMITS/ RMAX, NATMX, NFLIM, NPROPS
      COMMON /S/ THETA(MAxTht),PHI(MaxPhi),SUR(3,MaxTht,MaxPhi)
      COMMON /ORBTL/ EORB(MMO),PO(MMO),ROOTPO(MMO),NMO
      COMMON /C1/  BASSIN(Maxprp),TMP(MMO*(MMO+1)/2),NAOL
      DIMENSION FINTPH(MaxPrp),FU(MaxPrp),AOMF(MMO*(MMO+1)/2),
     $AOMT(MMO*(MMO+1)/2)
C
      DO 1 I=1,IPHIPL
      CALL PATH(IBETA,BETA,I,FU,IT,IPHIPL,AOMF)
      DO 2 J=1,NProps
2     FINTPH(J)=FINTPH(J)+WP(I)*FU(J)
      DO 3 LA=1,NAOL
3     AOMT(LA)=AOMT(LA)+WP(I)*AOMF(LA)
1     CONTINUE
      RETURN
      END
      Subroutine Quad(Temp,NTheta,theta,wtt,DMult1,Dmult2,MaxTht)
      Implicit DOUBLE PRECISION (A-H,O-Z)
C
      Dimension Temp(MaxTht,3),theta(MaxTht),wtt(MaxTht)
      Save Zero,One,Two,Eight,Tol
      Data   Zero,One,Two,Eight,Tol
     $     / 0.D0,1.D0,2.D0,8.D0,1.D-12 /
C
        Do 500 i = 2,NTheta
          wtt(i) = Two - One / i
  500     temp(i,3) = One - One / i
        Pi   = DACos(-One)
        Pi1=Pi/Two
        MaxM = NTheta / 2
        Do 510 m = 1,MaxM
          DTheta = Pi * (4*m-1) / (4*NTheta+2)
  510   temp(m,1)=DCOS(DTheta+One/(Eight*NTheta**2*DTAN(DTheta)))
        temp(MaxM+1,1) = Zero
        Do 520 m = 1,MaxM
          z   = temp(m,1)
  530       P   = z
            P1  = One
            dP  = One
            dP1 = Zero
            Do 540 i = 2,NTheta
              Q   = wtt(i)*z*P - temp(i,3)*P1
              dQ  = wtt(i)*(z*dP+P) - temp(i,3)*dP1
              P1  = P
              P   = Q
              dP1 = dP
  540         dP  = dQ
           z = z - P/dP
           If(DAbs(P).gt.Tol) goto 530
           temp(m,1) = z
  520      temp(m,2) = Two / ((One-z**2)*dP**2)
        Do 560 m = 1,MaxM
          theta(m)=DMult1-DMult2*temp(m,1)
          theta(Ntheta-m+1)=DMult1+DMult2*temp(m,1)
          wtt(m)=temp(m,2)
          wtt(Ntheta-m+1)=temp(m,2)
  560     Continue
      Return
      End
      SUBROUTINE RDPSI
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*80 WFNTTL,JOBTTL,LINE
      CHARACTER*8 ATNAM
      PARAMETER (MCENT=50, MMO=300, MPRIMS0=800, NTYPE=20)
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /ORBTL/ EORB(MMO),PO(MMO),ROOTPO(MMO),NMO
      COMMON /STRING/ WFNTTL,JOBTTL,ATNAM(MCENT),NAT
      COMMON /UNITS/ ISRF,INPT,IOUT,IWFN,IDBG
      COMMON /VALUES/ THRESH1,THRESH2,GAMMA,TOTE
      COMMON /PRIMS0/ COOMX(MMO),SUM(MPRIMS0),DIV(MPRIMS0),
     +COO(MPRIMS0,MMO),EXX(MPRIMS0),ICT(MPRIMS0),ITP(NTYPE),NPRIMS,
     +ITYPO(MPRIMS0)
      DIMENSION ITYPE(MPRIMS0),CO(MPRIMS0,MMO),ICENT(MPRIMS0),
     $EX(MPRIMS0)
      Save Zero,One,Two
      Data Zero/0.0d0/,One/1.0d0/,Two/2.0d0/
C
C    THE ITYPE ARRAY REPRESENTS THE FOLLOWING GAUSSIAN ORBITAL TYPES:
C     S
C     PX, PY, PZ 
C     DXX, DYY, DZZ, DXY, DXZ, DYZ
C     FXXX, FYYY, FZZZ, FXXY, FXXZ, FYYZ, FXYY, FXZZ, FYZZ, FXYZ
C
      READ(IWFN,101) WFNTTL
C
      READ(IWFN,102) MODE,NMO,NPRIMS,NCENT
C
      IF (NMO .GT. MMO) THEN
      WRITE(IOUT,3300) NMO
      STOP 'Too many molecular orbitals'
      ENDIF
      IF (NPRIMS .GT. MPRIMS0) THEN
      WRITE(IOUT,3301) NPRIMS
      STOP 'Too many total primitives'
      ENDIF
      IF (NCENT .GT. MCENT) THEN
      WRITE(IOUT,3302) NCENT
      STOP 'Too many atoms'
      ENDIF
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
      K=0
      DO 145 J=1,NPRIMS
      ICNT=0
      DO 146 L=1,NMO
      IF(CO(J,L).EQ.ZERO)ICNT=ICNT+1
146   CONTINUE
      IF(ICNT.NE.NMO)THEN
      K=K+1
      ITYPE(K)=ITYPE(J)
      ICENT(K)=ICENT(J)
      EX(K)=EX(J)
      DO 147 LL=1,NMO
      CO(K,LL)=CO(J,LL)
147   CONTINUE
      ENDIF
145   CONTINUE
      NPRIMS=K
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
      ITYPO(N)=K
180   CONTINUE
      ENDIF
170   CONTINUE
      ITP(K)=N
160   CONTINUE 
C
      Do 148 L=1,NMO
      temp=zero
      temp1=dabs(po(l))
      rootpo(l)=dsqrt(temp1)
      Do 149 J=1,Nprims
      if(dabs(coo(j,l)).gt.temp)temp=dabs(coo(j,l))
149   continue
      coomx(l)=temp
148   continue
C
101   FORMAT (1A80)
102   FORMAT (4X,A4,12X,3(I3,17X))
103   FORMAT (1A8,11X,1I3,2X,3F12.8,10X,F5.1)
104   FORMAT (20X,20I3)
105   FORMAT (10X,5E14.7)
106   FORMAT (35X,F12.8,15X,F12.8)
107   FORMAT (5E16.8)
109   FORMAT (17X,F20.12,18X,F13.8)
3300  FORMAT(1I3,' molecular orbitals  Redimension MMO')
3301  FORMAT(1I3,' total primitives  Redimension MPRIMS0')
3302  FORMAT(1I3,' atoms in molecule  Redimension MCENT')
      RETURN
      END
      SUBROUTINE RESULT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*80 WFNTTL,JOBTTL
      CHARACTER*8 ATNAM
      LOGICAL RHF,ROHF,ABNAT,RNAT
      PARAMETER (MCENT=50, MMO=300,MaxPrp=200)
      COMMON /LIMITS/ RMAX, NATMX, NFLIM, NPROPS
      COMMON /C1/  BAS(MaxPrp),AOM(MMO*(MMO+1)/2),naol
      COMMON/ATOMS/XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /UNITS/ ISRF,INPT,IOUT,IWFN,IDBG
      COMMON /STRING/ WFNTTL,JOBTTL,ATNAM(MCENT),NAT
      COMMON /VALUES/ THRESH1,THRESH2,GAMMA,TOTE
      COMMON /ORBTL/ EORB(MMO),PO(MMO),ROOTPO(MMO),NMO
      COMMON /DIRK/ Toll,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,IMAGM,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGO,YGO,ZGO
      COMMON /C7/  CT(3,3), X, Y, Z, NCRNT
      DIMENSION AOMT(MMO*(MMO+1)/2),QD(3,3),EV(3),WORK(3)
      Save Zero,One,Two
      Data Zero/0.0d0/,One/1.0d0/,Two/2.0d0/
C
      LMO=NMO
      If(idomag.eq.1)lmo=nmo/7
C
C    DIAGONALIZE THE ATOMIC QUADRAPOLE MOMENT
C
      QD(1,1) = BAS(20)
      QD(1,2) = BAS(21)
      QD(2,1) = BAS(21)
      QD(1,3) = BAS(22)
      QD(3,1) = BAS(22)
      QD(2,2) = BAS(23)
      QD(2,3) = BAS(24)
      QD(3,2) = BAS(24)
      QD(3,3) = BAS(25)
C
      CALL TRACE(QD,EV,WORK,3,IERR)
C
      WRITE(IOUT,970)
970   FORMAT(/' RESULTS OF THE INTEGRATION')
      CHARGE = CHARG(NAT) - BAS(1)
      WRITE(IOUT,1000) BAS(1),CHARGE
1000  FORMAT(5X,'         N ',1PE22.14,5X,'NET CHARGE ',1PE22.14)
      WRITE(IOUT,1020) BAS(2)
1020  FORMAT(5X,'         G ',1PE22.14)
      EOFATOM=BAS(3)*(ONE-GAMMA)
      WRITE(IOUT,1030) BAS(3),EOFATOM
1030  FORMAT(5X,'         K ',1PE22.14,5X,'   E(ATOM) ',1PE22.14)
      WRITE(IOUT,1040) BAS(4)
1040  FORMAT(5X,'         L ',1PE22.14)
      WRITE(IOUT,1050) BAS(5)
1050  FORMAT(5X,'         I ',1PE22.14)
      WRITE(IOUT,1060) BAS(6)
1060  FORMAT(5X,'     R(-1) ',1PE22.14)
      WRITE(IOUT,1070) BAS(7)
1070  FORMAT(5X,'        R1 ',1PE22.14)
      WRITE(IOUT,1080) BAS(8)
1080  FORMAT(5X,'        R2 ',1PE22.14)
      WRITE(IOUT,1370) BAS(35)
1370  FORMAT(5X,'        R4 ',1PE22.14)
      WRITE(IOUT,1090) BAS(9)
1090  FORMAT(5X,'    GR(-1) ',1PE22.14)
      WRITE(IOUT,1100) BAS(10)
1100  FORMAT(5X,'       GR0 ',1PE22.14)
      WRITE(IOUT,1110) BAS(11)
1110  FORMAT(5X,'       GR1 ',1PE22.14)
      WRITE(IOUT,1120) BAS(12)
1120  FORMAT(5X,'       GR2 ',1PE22.14)
      VNEO2 = (-TWO*(ONE-GAMMA)/GAMMA)*BAS(13)
      WRITE(IOUT,1130) BAS(13),VNEO2
1130  FORMAT(5X,'      VNEO ',1PE22.14,5X,' VNEO(COR) ',1PE22.14)
      VNET2 = (-TWO*(ONE-GAMMA)/GAMMA)*BAS(14)
      WRITE(IOUT,1140) BAS(14),VNET2
1140  FORMAT(5X,'      VNET ',1PE22.14,5X,' VNET(COR) ',1PE22.14)
      VEET2 = (-TWO*(ONE-GAMMA)/GAMMA)*BAS(16)
      WRITE(IOUT,1150) BAS(16),VEET2
1150  FORMAT(5X,'      VEET ',1PE22.14,5X,' VEET(COR) ',1PE22.14)
      WRITE(IOUT,1160) BAS(15)
1160  FORMAT(5X,'       EHF ',1PE22.14)
      VNN2 = EOFATOM - BAS(15)
C      VNN2 = TWO*EOFATOM - VNET2 - VEET2 
C      WRITE(IOUT,1170) VNN2
C1170  FORMAT(5X,'           ',22X,5X,'  VNN(COR) ',1PE22.14)
      VR2 = VEET2 + VNN2
      WRITE(IOUT,1180) VR2
1180  FORMAT(5X,'           ',22X,5X,' VREP(COR) ',1PE22.14)
      VTOT2 = VNET2 + VR2
      WRITE(IOUT,1190) VTOT2
1190  FORMAT(5X,'           ',22X,5X,'   V(ATOM) ',1PE22.14)
      WRITE(IOUT,1200) BAS(17)
1200  FORMAT(5X,'     EL DX ',1PE22.14)
      WRITE(IOUT,1210) BAS(18)
1210  FORMAT(5X,'     EL DY ',1PE22.14)
      WRITE(IOUT,1220) BAS(19)
1220  FORMAT(5X,'     EL DZ ',1PE22.14)
      AMU = DSQRT(BAS(17)*BAS(17) +
     +            BAS(18)*BAS(18) +
     +            BAS(19)*BAS(19))
      WRITE(IOUT,1230) AMU
1230  FORMAT(2X,'EL DIPOLE MAG ',1PE22.14)
C
      WRITE(IOUT,*)
      WRITE(IOUT,1235)
1235  FORMAT(5X,'ATOMIC QUADRUPOLE MOMENT TENSOR')
      WRITE(IOUT,1240) BAS(20)
1240  FORMAT(5X,'       QXX ',1PE22.14)
      WRITE(IOUT,1250) BAS(21)
1250  FORMAT(5X,'       QXY ',1PE22.14)
      WRITE(IOUT,1260) BAS(22)
1260  FORMAT(5X,'       QXZ ',1PE22.14)
      WRITE(IOUT,1270) BAS(23)
1270  FORMAT(5X,'       QYY ',1PE22.14)
      WRITE(IOUT,1280) BAS(24)
1280  FORMAT(5X,'       QYZ ',1PE22.14)
      WRITE(IOUT,1290) BAS(25)
1290  FORMAT(5X,'       QZZ ',1PE22.14)
C
      WRITE(IOUT,*)
      WRITE(IOUT,2089)
2089  FORMAT(5X,'EIGENVALUES OF QUADRUPOLE MOMENT TENSOR:')
      WRITE(IOUT,2090) EV(1), EV(2), EV(3)
2090  FORMAT(5X,1PE22.14,2X,1PE22.14,2X,1PE22.14)
C
      WRITE(IOUT,*)
      WRITE(IOUT,2091)
2091  FORMAT(5X,'EIGENVECTORS OF QUADRUPOLE MOMENT TENSOR:')
      WRITE(IOUT,2092) QD(1,1), QD(1,2), QD(1,3)
      WRITE(IOUT,2092) QD(2,1), QD(2,2), QD(2,3)
      WRITE(IOUT,2092) QD(3,1), QD(3,2), QD(3,3)
2092  FORMAT(5X,1PE22.14,2X,1PE22.14,2X,1PE22.14)
      WRITE(IOUT,*)
C
      WRITE(IOUT,1300) BAS(26)
1300  FORMAT(5X,'      FAXA ',1PE22.14)
      WRITE(IOUT,1310) BAS(27)
1310  FORMAT(5X,'      FAYA ',1PE22.14)
      WRITE(IOUT,1320) BAS(28)
1320  FORMAT(5X,'      FAZA ',1PE22.14)
      WRITE(IOUT,1330) BAS(29)
1330  FORMAT(5X,'      FBXA ',1PE22.14)
      WRITE(IOUT,1340) BAS(30)
1340  FORMAT(5X,'      FBYA ',1PE22.14)
      WRITE(IOUT,1350) BAS(31)
1350  FORMAT(5X,'      FBZA ',1PE22.14)
      WRITE(IOUT,1360) BAS(32)
1360  FORMAT(5X,'     RHO*L ',1PE22.14)
      WRITE(IOUT,1380) BAS(37)
1380  FORMAT(5X,'      VOL1 ',1PE22.14)
      WRITE(IOUT,1390) BAS(38)
1390  FORMAT(5X,'      VOL2 ',1PE22.14)
      WRITE(IOUT,1400) BAS(42)
1400  FORMAT(5X,'   N(VOL1) ',1PE22.14)
      WRITE(IOUT,1410) BAS(43)
1410  FORMAT(5X,'   N(VOL2) ',1PE22.14)
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      WRITE(IOUT,1580) RMAX
1580  FORMAT(/,' MAXIMUM DISTANCE REACHED FROM NUCLEUS =',1PE18.10)
C
      WRITE(IOUT,*)
      If(Idoaom.eq.0)Then
      Write(Iout,1599)
1599  Format('Atomic Overlap Matrix Not Calculated')
      Else
      WRITE(IOUT,1600)
1600  FORMAT(10X,'The Atomic Overlap Matrix')
      WRITE(IOUT,*)
C
      If(RHF.or.RNAT)Write(iout,1604)
1604  Format('Restricted Closed-Shell Wavefunction')
      If(ROHF)Write(iout,1602)
1602  Format('Restricted Open-Shell Wavefunction')
      If(ABNAT)Write(iout,1603)
1603  Format('Unrestricted Wavefunction')
C
      WRITE(IOUT,*)
C
C     If this is an unrestricted wfn (with alpha mos first and beta 
C     mos second and NfBeta is the first Beta MO), then zero the 
C     alpha-beta block of the atomic overlap matrix.
C
      If(ABNat) then
        Do 1705 I = 1, NfBeta-1 
          Do 1705 J = NfBeta, LMO 
            IJ = (J*(J-1))/2 + I
 1705       AOM(IJ) = Zero 
        endif
C
      K=0 
      DO 1710 LA=1,LMO
      DO 1720 LB=1,LA
      K=K+1
      AOMT(LB)=AOM(K)
1720  CONTINUE
      WRITE(IOUT,1700) (AOMT(J),J=1,LA)
1710  CONTINUE
1700  FORMAT(8F10.6)
C
      If(ABNat) then
        FOOA = Zero
        FOOB = Zero
        AN = Zero
        BN = Zero
        K=0
        Do I=1,NfBeta-1
          AN = AN+AOM(K+I)*PO(I)
          K=K+I
          enddo
        Do I=NfBeta,LMO
          BN = BN+AOM(K+I)*PO(I)
          K=K+I
          enddo
        K=0
        DO 2920 I = 1,NfBeta-1
          DO 2920 J = 1,I
            K=K+1
            HH = Two
            IF(I.EQ.J) HH = One
            ANMIN = dsqrt(PO(I)*PO(J))
            FOOA = FOOA-HH*AOM(K)**2*ANMIN
2920        CONTINUE 
C
        DO 2930 I = NfBeta,LMO
          DO 2930 J = NfBeta,I
            HH = Two
            IF(I.EQ.J) HH = One
            K = (I*(I-1))/2 + J
            BNMIN = dsqrt(PO(I)*PO(J))
            FOOB = FOOB-HH*AOM(K)**2*BNMIN
2930        CONTINUE 
        endif
C
      If(RHF.or.RNAT) then
        FOOA = Zero
        FOOB = Zero
        AN = Zero
        BN = Zero
        K=0
        DO 2940 I = 1,LMO
          AN = AN+AOM(K+I)*PO(I)
          BN = BN+AOM(K+I)*PO(I)
          K=K+I
2940      CONTINUE
         AN=AN/TWO
         BN=BN/TWO
C
        K=0
        DO 2950 I = 1,LMO
          DO 2950 J = 1,I
            K=K+1
            HH = Two
            IF(I.EQ.J)HH = One
            ANMIN=sqrt(po(i)*po(j))/Two
            FOOA = FOOA-HH*ANMIN*AOM(K)**2
            FOOB = FOOB-HH*ANMIN*AOM(K)**2
2950        CONTINUE
        endif
C
      If(ROHF) then
        Do I=1,LMO
          If(PO(I).eq.One) then
            Alpha1 = I
            goto 1959
            endif
          enddo
1959    Continue
        FOOA = Zero
        FOOB = Zero
        AN = Zero
        BN = Zero
        K=0
        DO I=1,Alpha1-1
          AN = AN+AOM(K+I)
          BN = BN+AOM(K+I)
          K=K+I
          enddo
        Do I=Alpha1,LMO
          AN = AN+AOM(K+I)
          K=K+I
          enddo
C
        K=0
        DO 2960 I = 1,Alpha1-1
          DO 2960 J = 1,I
            K=K+1
            HH = Two
            IF(I.EQ.J)HH = One
            FOOA1 = FOOA1-HH*AOM(K)**2
            FOOB1 = FOOB1-HH*AOM(K)**2
2960        CONTINUE
        DO 2965 I = 1,Alpha1-1
          DO 2965 J = Alpha1,LMO
            K=K+1
            FOOA2 = FOOA2-AOM(K)**2
2965        CONTINUE
        DO 2970 I = Alpha1,LMO
          DO 2970 J = Alpha1,I
            K=K+1
            HH = Two
            IF(I.EQ.J)HH = One
            FOOA3 = FOOA-HH*AOM(K)**2
2970        CONTINUE
           FOOA=FOOA1+FOOA2+FOOA3
           FOOB=FOOB1
        endif
C
      If(RHF.or.ROHF.or.ABNat.or.Rnat) then
        ALOC=DABS(FOOA)/AN
        BLOC=DABS(FOOB)/BN
        FLA=AN+FOOA
        FLB=BN+FOOB
        CN=AN+BN
        WRITE(IOUT,*)
        WRITE(IOUT,1460)AN
1460  FORMAT(5X,' ALPHA ELECTRONS (NA)             ',1PE22.14)
        WRITE(IOUT,1470) BN
1470  FORMAT(5X,' BETA ELECTRONS (NB)              ',1PE22.14)
        WRITE(IOUT,1480) CN
1480  FORMAT(5X,' TOTAL ELECTRONS (N)              ',1PE22.14)
        WRITE(IOUT,1490) FOOA
1490  FORMAT(5X,' ALPHA FERMI CORRELATION (FOOA)   ',1PE22.14)
        WRITE(IOUT,1500) FOOB
1500  FORMAT(5X,' BETA FERMI CORRELATION (FOOB)    ',1PE22.14)
        WRITE(IOUT,1510) ALOC
1510  FORMAT(5X,' ALPHA LOCALIZATION (ALOC)        ',1PE22.14)
        WRITE(IOUT,1520) BLOC
1520  FORMAT(5X,' BETA LOCALIZATION (BLOC)         ',1PE22.14)
        WRITE(IOUT,1530) FLA
1530  FORMAT(5X,' ALPHA FLUCTUATION (FLA)          ',1PE22.14)
        WRITE(IOUT,1540) FLB
1540  FORMAT(5X,' BETA FLUCTUATION (FLB)           ',1PE22.14)
        endif
C
      Endif
C
      IF (IDOMAG .EQ. 1) Then
C
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      WRITE(IOUT,1585)
1585  FORMAT(20X,'**********************************')
      WRITE(IOUT,1586) ATNAM(NAT)
1586  FORMAT(20X,'MAGNETIC PROPERTIES OF ',A8)
      WRITE(IOUT,1585)
      WRITE(IOUT,*)
      If(Imagm.eq.0)write(Iout,1541)
1541  Format('Used IGAIM Method')
      If(Imagm.eq.1)write(Iout,1542)
1542  Format('Used Single-Gauge-Origin Method, Gauge Origin is: ')
      If(Imagm.eq.1)write(Iout,*)Xgo,YGo,ZGo
      If(Imagm.eq.2)write(Iout,1543)
1543  Format('Used Becke-Igaim Method')
      WRITE(IOUT,*)
      WRITE(IOUT,1588)
1588  FORMAT('MAGNETIC SUSCEPTIBILITY IN ATOMIC UNITS')
      WRITE(IOUT,*)
      WRITE(IOUT,1589)
1589  FORMAT('The Internal Contribution')
      WRITE(IOUT,*)
      WRITE(IOUT,1590) BAS(44)
1590  FORMAT(10X,'CHIINTXX ',1PE22.14)
      WRITE(IOUT,1601) BAS(45)
1601  FORMAT(10X,'CHIINTYX ',1PE22.14)
      WRITE(IOUT,1611) BAS(46)
1611  FORMAT(10X,'CHIINTZX ',1PE22.14)
      WRITE(IOUT,1621) BAS(44+6+NCENT*3)
1621  FORMAT(10X,'CHIINTXY ',1PE22.14)
      WRITE(IOUT,1630) BAS(45+6+NCENT*3)
1630  FORMAT(10X,'CHIINTYY ',1PE22.14)
      WRITE(IOUT,1640) BAS(46+6+NCENT*3)
1640  FORMAT(10X,'CHIINTZY ',1PE22.14)
      WRITE(IOUT,1650) BAS(44+2*(6+NCENT*3))
1650  FORMAT(10X,'CHIINTXZ ',1PE22.14)
      WRITE(IOUT,1660) BAS(45+2*(6+ncent*3))
1660  FORMAT(10X,'CHIINTYZ ',1PE22.14)
      WRITE(IOUT,1670) BAS(46+2*(6+ncent*3))
1670  FORMAT(10X,'CHIINTZZ ',1PE22.14)
      WRITE(IOUT,*)
      A=(BAS(44)+BAS(45+6+ncent*3)+BAS(46+2*(6+ncent*3)))/3
      WRITE(IOUT,1672) A
1672  FORMAT(10X,'CHIINT(MEAN) ',1PE22.14)
      WRITE(IOUT,*)
      WRITE(IOUT,1675)
1675  FORMAT('The Nuclear-weighted Current Contribution.')
      WRITE(IOUT,*)
      A1=Y*BAS(49)-Z*BAS(48)
      WRITE(IOUT,1681) A1 
1681  FORMAT(10X,'CHIOTXX ',1PE22.14)
      A2=Z*BAS(47)-X*BAS(49)
      WRITE(IOUT,1691) A2 
1691  FORMAT(10X,'CHIOTYX ',1PE22.14)
      A3=X*BAS(48)-Y*BAS(47)
      WRITE(IOUT,1701) A3 
1701  FORMAT(10X,'CHIOTZX ',1PE22.14)
      A4=Y*BAS(49+6+ncent*3)-Z*BAS(48+6+ncent*3)
      WRITE(IOUT,1711) A4 
1711  FORMAT(10X,'CHIOTXY ',1PE22.14)
      A5=Z*BAS(47+6+ncent*3)-X*BAS(49+6+ncent*3)
      WRITE(IOUT,1721) A5 
1721  FORMAT(10X,'CHIOTYY ',1PE22.14)
      A6=X*BAS(48+6+ncent*3)-Y*BAS(47+6+ncent*3)
      WRITE(IOUT,1731) A6 
1731  FORMAT(10X,'CHIOTZY ',1PE22.14)
      A7=Y*BAS(49+2*(ncent*3+6))-Z*BAS(48+2*(6+ncent*3))
      WRITE(IOUT,1741) A7 
1741  FORMAT(10X,'CHIOTXZ ',1PE22.14)
      A8=Z*BAS(47+2*(6+ncent*3))-X*BAS(49+2*(6+ncent*3))
      WRITE(IOUT,1751) A8 
1751  FORMAT(10X,'CHIOTYZ ',1PE22.14)
      A9=X*BAS(48+2*(6+ncent*3))-Y*BAS(47+2*(6+ncent*3))
      WRITE(IOUT,1761) A9 
1761  FORMAT(10X,'CHIOTZZ ',1PE22.14)
      WRITE(IOUT,*)
      WRITE(IOUT,1762)(A1+A5+A9)/3
1762  FORMAT(10X,'CHIOT(MEAN) ',1PE22.14)
      WRITE(IOUT,*)
      WRITE(IOUT,1765)
1765  FORMAT('The Total Atomic Susceptibility.')
      WRITE(IOUT,*)
      WRITE(IOUT,1770) BAS(44)+A1
1770  FORMAT(10X,'CHITOTXX ',1PE22.14)
      WRITE(IOUT,1780) BAS(45)+A2
1780  FORMAT(10X,'CHITOTYX ',1PE22.14)
      WRITE(IOUT,1790) BAS(46)+A3
1790  FORMAT(10X,'CHITOTZX ',1PE22.14)
      WRITE(IOUT,1801) BAS(44+6+ncent*3)+A4
1801  FORMAT(10X,'CHITOTXY ',1PE22.14)
      WRITE(IOUT,1811) BAS(45+6+ncent*3)+A5
1811  FORMAT(10X,'CHITOTYY ',1PE22.14)
      WRITE(IOUT,1821) BAS(46+6+ncent*3)+A6
1821  FORMAT(10X,'CHITOTZY ',1PE22.14)
      WRITE(IOUT,1831) BAS(44+2*(6+ncent*3))+A7
1831  FORMAT(10X,'CHITOTXZ ',1PE22.14)
      WRITE(IOUT,1841) BAS(45+2*(6+ncent*3))+A8
1841  FORMAT(10X,'CHITOTYZ ',1PE22.14)
      WRITE(IOUT,1851) BAS(46+2*(6+ncent*3))+A9
1851  FORMAT(10X,'CHITOTZZ ',1PE22.14)
      WRITE(IOUT,*)
      WRITE(IOUT,1852)A+(A1+A5+A9)/3
1852  FORMAT(10X,'CHITOT(MEAN) ',1PE22.14)
      WRITE(IOUT,*)
      WRITE(IOUT,1855)
1855  FORMAT('Atomic Integral of the Current Density - Half Unit Field')
      WRITE(IOUT,*)
      WRITE(IOUT,2800) BAS(47)
2800  FORMAT(10X,'CURINTXX ',1PE22.14)
      WRITE(IOUT,2801) BAS(48)
2801  FORMAT(10X,'CURINTYX ',1PE22.14)
      WRITE(IOUT,2802) BAS(49)
2802  FORMAT(10X,'CURINTZX ',1PE22.14)
      WRITE(IOUT,2803) BAS(47+6+ncent*3)
2803  FORMAT(10X,'CURINTXY ',1PE22.14)
      WRITE(IOUT,2804) BAS(48+6+ncent*3)
2804  FORMAT(10X,'CURINTYY ',1PE22.14)
      WRITE(IOUT,2805) BAS(49+6+ncent*3)
2805  FORMAT(10X,'CURINTZY ',1PE22.14)
      WRITE(IOUT,2806) BAS(47+2*(6+ncent*3))
2806  FORMAT(10X,'CURINTXZ ',1PE22.14)
      WRITE(IOUT,2807) BAS(48+2*(6+ncent*3))
2807  FORMAT(10X,'CURINTYZ ',1PE22.14)
      WRITE(IOUT,2808) BAS(49+2*(6+ncent*3))
2808  FORMAT(10X,'CURINTZZ ',1PE22.14)
      WRITE(IOUT,*)
      WRITE(IOUT,1860)
1860  FORMAT(10X,'Atomic Contribution to Shielding Tensors in ppm')
      WRITE(IOUT,*)
      DO 3000 II=1,ncent
      WRITE(IOUT,*)
      WRITE(IOUT,1870) ATNAM(II)
1870  FORMAT(10X,A8)
      WRITE(IOUT,1880) BAS(50+(II-1)*3)
1880  FORMAT(10X,'SIGMAXX ',1PE22.14)
      WRITE(IOUT,1890) BAS(51+(II-1)*3)
1890  FORMAT(10X,'SIGMAYX ',1PE22.14)
      WRITE(IOUT,1900) BAS(52+(II-1)*3)
1900  FORMAT(10X,'SIGMAZX ',1PE22.14)
      WRITE(IOUT,1910) BAS(50+6+3*ncent+(II-1)*3)
1910  FORMAT(10X,'SIGMAXY ',1PE22.14)
      WRITE(IOUT,1921) BAS(51+6+3*ncent+(II-1)*3)
1921  FORMAT(10X,'SIGMAYY ',1PE22.14)
      WRITE(IOUT,1930) BAS(52+6+3*ncent+(II-1)*3)
1930  FORMAT(10X,'SIGMAZY ',1PE22.14)
      WRITE(IOUT,1940) BAS(50+2*(6+3*ncent)+(II-1)*3)
1940  FORMAT(10X,'SIGMAXZ ',1PE22.14)
      WRITE(IOUT,1950) BAS(51+2*(6+3*ncent)+(II-1)*3)
1950  FORMAT(10X,'SIGMAYZ ',1PE22.14)
      WRITE(IOUT,1960) BAS(52+2*(6+3*ncent)+(II-1)*3)
1960  FORMAT(10X,'SIGMAZZ ',1PE22.14)
      WRITE(IOUT,*)
      WRITE(IOUT,1970) (BAS(50+(II-1)*3)+BAS(51+6+3*ncent+
     +(II-1)*3)+BAS(52+2*(6+3*ncent)+(II-1)*3))/3
1970  FORMAT(10X,'SIGMA(MEAN) ',1PE22.14)
C
3000  CONTINUE
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE RING(CRIT,N,NR,ANGLE,NSRC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER(MaxCrt=20,Mxbcrt=10,MaxTht=200,MaxPhi=200,
     $Maxpp=300)
      COMMON /C3/ WK(3*Maxpp)
      COMMON /EV/ V(3,MaxCrt),W(3,MaxCrt)
      COMMON/PARA/NINS(Mxbcrt),IBETPTS,INMAX,NNN,NPATH,DIST,
     $DMAX,AMAX,ADIFF,tinf
      COMMON /S/ XI(MaxTht),YI(MaxPhi),ZI(3,MaxTht*MaxPhi)
      DIMENSION  CRIT(3,MaxCrt),XB(3,2),XX(3),AA(3,2),IP(2),
     $ANGLE(2,2,Mxbcrt),NSRC(4,Mxbcrt)
      Save Zero,One,Two,Pt8,Otmfiv,Otmfr,otmegt
      Data Zero/0.0d0/,One/1.0d0/,Two/2.0d0/,Pt8/0.8d0/,
     $otmfr/1.0d-4/,Otmfiv/1.0d-5/,otmegt/1.0d-8/
C
      PI=DACOS(-one)
      TPI=TWO*PI
C
      IF (NR .EQ. 0) RETURN
      DO 1000 I = 1,NR
C
C     FOR EACH RING CRIT. POINT COMPUTE THE PATHS TO THE
C     TWO NEIGHBOURING (3,-1) CRIT. POINTS.
C     IMSL ROUTINE IS USED.
C
      N1 = N + I
      DO 1010 J = 1,2
      FF = Zero
      IZZ = 1
      N2 = NSRC(J,I)
      IF (N2 .EQ. 0) GOTO 1010
      QL = DSQRT((CRIT(1,N2)-CRIT(1,N1))**2 +
     +           (CRIT(2,N2)-CRIT(2,N1))**2 + 
     +           (CRIT(3,N2)-CRIT(3,N1))**2) - Otmfr
      NM = INT(DFLOAT(NNN)*QL/DMAX)
      NM1 = NM + 1
      DS = QL/DFLOAT(NM)
C
C     STORE THE LENGTH OF THE STRAIGHT LINE FROM THE (3,-1) TO THE
C     (3,+1) CRITICAL POINT AND ITS ANGLE AT THE (3,-1) CRIT.
C     POINT AS INITIAL GUESSES.
C
      DO 100 K = 1,3
        AA(K,1) = V(K,N2)
        AA(K,2) = W(K,N2)
        WK(K) = CRIT(K,N1) - CRIT(K,N2)
100   CONTINUE
      G = DSQRT(WK(1)**2 + WK(2)**2 + WK(3)**2)
C
C  SCALE WK
C
      DO 110 K = 1,3
        WK(K) = WK(K)/G
110   CONTINUE
C
C    CALL TO LLSQB TO SOLVE LINEAR LEAST SQUARES AA*X = WK
C
      CALL LLSQB(AA,3,3,2,WK(1),-One,2,WK(4),WK(6),
     +           IP,IER)
C
      CO = DSQRT(WK(4)**2 + WK(5)**2)
      SI = WK(5)/CO
      CO = WK(4)/CO
      IF (CO .LT. -One) CO =-One
      IF (CO .GT.  One) CO = One
      ANGLE(1,J,I) = DACOS(CO)
      IF (SI .LT. Zero) ANGLE(1,J,I) = Tpi - ANGLE(1,J,I)
      ALP = ANGLE(1,J,I)
C
C     COMPUTE THE PATH WITH ANGLE ALP (INITIAL GUESS).
C
      DO 120 K = 1,3
        ZI(K,1) = CRIT(K,N2) + 
     +            Otmfr*(DCOS(ALP)*V(K,N2) + DSIN(ALP)*W(K,N2))
120   CONTINUE
C
      CALL TRUDGE2(ZI,DS,NM)
C
      D = (CRIT(1,N1)-ZI(1,NM1))**2 + 
     +    (CRIT(2,N1)-ZI(2,NM1))**2 +
     +    (CRIT(3,N1)-ZI(3,NM1))**2
      DO 130 K = 1,3
        XX(K) = ZI(K,NM1)
130   CONTINUE
      IF (D .LE. Otmegt) GOTO 4
C
C     THE SOLUTION PATH IS ASSUMED TO HAVE AN ANGLE ALP WITH
C     ALP1 < ALP < ALP2.
C
      ALP1=ALP-Pt8
      ALP2=ALP+Pt8
      DO 5 K=1,3
      ZI(K,201)=CRIT(K,N2)+Otmfr*(DCOS(ALP2)*V(K,N2)+
     1 DSIN(ALP2)*W(K,N2))
5     ZI(K,1)=CRIT(K,N2)+Otmfr*(DCOS(ALP1)*V(K,N2)+DSIN(ALP1)*W(K,N2))
      IZZ=IZZ+2
      CALL TRUDGE2(ZI,DS,NM)
      CALL TRUDGE2(ZI(1,201),DS,NM)
      H=Zero
      H1=Zero
      G=Zero
      G1=Zero
      DO 9 K=1,3
      XB(K,1)=CRIT(K,N1)-V(K,N1)
      XB(K,2)=CRIT(K,N1)+V(K,N1)
      H=H+(XB(K,1)-ZI(K,200+NM1))**2
      H1=H1+(XB(K,2)-ZI(K,200+NM1))**2
      G=G+(XB(K,1)-ZI(K,NM1))**2
      G1=G1+(XB(K,2)-ZI(K,NM1))**2
9     CONTINUE
      F1=-DABS(DSQRT(G)-DSQRT(G1))
      F2=DABS(DSQRT(H)-DSQRT(H1))
      IF(G1.GE.G)GOTO 10
      DO 11 K=1,3
      XB(K,1)=XB(K,2)
11    XB(K,2)=CRIT(K,N1)-V(K,N1)
C
C     USE THE PROPERTY THAT THE RING PATH IS PERPENDICULAR TO THE
C     EIGENVECTOR OF TH NEGATIVE EIGENVALUE OF THE HESSIAN AT THE
C     (3,+1) CRITICAL POINT (THIS PROPERTY IS EQUIVALENT TO FF=0.).
C     SET UP A REGULA FALSI SEARCH PROCESS FOR THE PATH WITH FF=0.
C
10    A1=(XX(1)-XB(1,1))**2+(XX(2)-XB(2,1))**2+(XX(3)-XB(3,1))**2
      A2=(XX(1)-XB(1,2))**2+(XX(2)-XB(2,2))**2+(XX(3)-XB(3,2))**2
      FF=DSQRT(A1)-DSQRT(A2)
      IF(DABS(FF).LE.Otmfr)GOTO 4
      IF(A1.LT.A2)ALP1=ALP
      IF(A1.LT.A2)F1=FF
      IF(A1.GT.A2)ALP2=ALP
      IF(A1.GT.A2)F2=FF
C
C     FORMULA OF THE REGULA FALSI FOR THE NEXT APPROXIMATION :
C
      ALP=(ALP1*F2-ALP2*F1)/(F2-F1)
      IF(DABS(ALP1-ALP2).LE.AMAX)GOTO 4
      DO 13 K=1,3
13    ZI(K,1)=CRIT(K,N2)+(DCOS(ALP)*V(K,N2)+DSIN(ALP)*W(K,N2))*Otmfr
      IZZ=IZZ+1
      CALL TRUDGE2(ZI,DS,NM)
      DO 14 K=1,3
14    XX(K)=ZI(K,NM1)
      D=(CRIT(1,N1)-XX(1))**2+(CRIT(2,N1)-XX(2))**2
     1 +(CRIT(3,N1)-XX(3))**2
      IF(D.LE.Otmegt)GOTO 4
      GOTO 10
4     ANGLE(1,J,I)=ALP
      ANGLE(2,J,I)=QL+DSQRT(D)-Otmfiv
C
1010  CONTINUE
1000  CONTINUE
      RETURN
      END
      SUBROUTINE RUNGE(R0,R1,GR0,DS)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION R0(3),R1(3),GR0(3),GR1(3)
      Save Pt5,Three,Six
      Data pt5/0.5d0/,Three/3.0d0/,Six/6.0d0/
C
      DO 20 I=1,3
      GR1(I)=DS*GR0(I)
      R1(I)=R0(I)+GR1(I)/Six
20    GR1(I)=R0(I)+Pt5*GR1(I)
C
      CALL GAUS2(GR1,0)
C
      G=DSQRT(GR1(1)**2+GR1(2)**2+GR1(3)**2)
      DO 30 I=1,3
      GR1(I)=DS*GR1(I)/G
      R1(I)=R1(I)+GR1(I)/Three
30     GR1(I)=R0(I)+Pt5*GR1(I)
      CALL GAUS2(GR1,0)
C
      G=DSQRT(GR1(1)**2+GR1(2)**2+GR1(3)**2)
      DO 40 I=1,3
      GR1(I)=DS*GR1(I)/G
      R1(I)=R1(I)+GR1(I)/Three
40     GR1(I)=R0(I)+GR1(I)
      CALL GAUS2(GR1,0)
C
      G=DSQRT(GR1(1)**2+GR1(2)**2+GR1(3)**2)
      DO 50 I=1,3
50    R1(I)=R1(I)+DS*GR1(I)/(Six*G)
      RETURN
      END
      Subroutine Sphere(Beta)
C
      Implicit Double Precision (A-H,O-Z)
C
      LOGICAL RHF,ROHF,ABNAT,RNAT
      Parameter(MaxTht=200,MaxPhi=200,Mcent=50,MPhiSm=16,MThtSm=8)
      COMMON /C7/ CT(3,3),X, Y, Z,NCRNT
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /DIRK/ Toll,TOLL2(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,Imagm,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGO,YGO,ZGO
      COMMON /UNITS/ ISRF,INPT,IOUT,IWFN,IDBG
      DIMENSION XYZ(MaxTht*Maxphi,3),ITermT(MaxTht*MaxPhi)
      Save One,Two,IPhat,IThat,Scale,Dinc
      Data One/1.0d0/,Two/2.0d0/,IPhat/MPhism/,IThat/MThtsm/,
     $Scale/0.8d0/,Dinc/0.1d0/
1999  Format('Radius of Capture Sphere for Atom ',I6,' is ',1PE8.2)
C
      Pi=Dacos(-one)
      NCRNTT=NCRNT
C
      DO 1 IA=1,NCent
      Toll2(IA)=Toll**2
1     Continue
C
      PHINC=TWO*PI/IPHAT
      HPHINC=PHINC/TWO
      THINC=PI/ITHAT
      HTHINC=THINC/TWO
      DO 2 IA=1,NCent
      NCRNT=IA
      N=0
      K=0
      Trad=Toll
      DO 3 I=1,IPHAT
      PH=PHINC*I-HPHINC
      DO 4 J=1,ITHAT
      N=N+1
      TH=THINC*J-HTHINC
      XYZ(N,1)=XC(IA)+Trad*DSin(TH)*DCos(Ph)
      XYZ(N,2)=YC(IA)+Trad*DSin(TH)*DSin(Ph)
      XYZ(N,3)=ZC(IA)+Trad*DCos(TH)
4     Continue
3     Continue
C
5     K=K+1
C
      CALL TRUDGE3(XYZ,N,ITERMT)
C
      IQuit=0
      Do 6 I=1,N
      If(ITermT(I).Eq.0)IQuit=1
6     Continue
C
      If(Iquit.eq.0)Then
      N=0
      DO 7 I=1,IPHAT
      PH=PHINC*I-HPHINC
      DO 8 J=1,ITHAT
      N=N+1
      TH=THINC*J-HTHINC
      XYZ(N,1)=XYZ(N,1)+DINC*DSin(TH)*DCos(Ph)
      XYZ(N,2)=XYZ(N,2)+DINC*DSin(TH)*DSin(Ph)
      XYZ(N,3)=XYZ(N,3)+DINC*DCos(TH)
8     Continue
7     Continue
      Trad=Trad+Dinc
      Toll2(ia)=((Trad-Stp)*Scale)**2
      Goto 5
      Endif
C
2     Continue
C
      NCRNT=NCRNTT
      BETA=Dsqrt(toll2(ncrnt))
C
      Return
      End
      SUBROUTINE SURFACE(C1,N,NR,NC,ANGLE,NSRC,IPHIPL,IT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter(Maxcrt=20,Mxbcrt=10,MaxTht=200,MaxPhi=200,Maxpp=300)
      COMMON/INSINF/ INSLIM(Mxbcrt)
      COMMON/PARA/ NINS(Mxbcrt),IBETPTS,INMAX,NNN,NPATH,DIST,
     $DMAX,AMax,ADIFF,Tinf
      COMMON /EV/ V(3,MaxCrt),W(3,MaxCrt)
      COMMON /S/ XI(MaxTht),YI(MaxPhi),ZI(3,MaxTht,MaxPhi)
      DIMENSION C1(3,Mxbcrt),XX0(3,Maxpp),XX1(3,Maxpp),
     $XS(3,MaxTht,MaxPhi),XX3(3,Maxpp),XS1(3,MaxTht*MaxPhi),
     $ANGLE(2,2,Mxbcrt),CC(3),IZ(MAxTht,MaxPhi),ICAGE(2,3),
     $NSRC(4,Mxbcrt),NSRC1(2,Mxbcrt),AA(3,3),XX2(3,Maxpp)
      Save Zero,one,Two,otmfiv,otmfr
      Data Zero/0.0d0/,One/1.0d0/,Two/2.0d0/,Otmfiv/1.0d-5/,
     $otmfr/1.0d-4/,Pt7/0.7d0/
C
      PI=DACOS(-ONE)
      TPI=Two*Pi
      TPiM=TPi-Otmfiv
      NPATH1=NPATH+1
      NNN1=NNN+1
      NNN2=NNN-1
C
C     INITIALIZE THE ARRAYS. XS CONTAINS THE POINTS ON THE RAYS WITH
C     DISTANCE 1 FROM THE NUCLEUS.
C
      DO 7 I=1,IT
      DO 7 J=1,IPHIPL
      IZ(I,J)=0
      ZI(1,I,J)=Tinf
      ZI(2,I,J)=-One
      ZI(3,I,J)=Tinf
      XS(1,I,J)=DSIN(XI(I))*DCOS(YI(J))
      XS(2,I,J)=DSIN(XI(I))*DSIN(YI(J))
      XS(3,I,J)=DCOS(XI(I))
7     CONTINUE
C
C     LOOP OVER THE NUMBER OF SURFACES
C
      DO 999 I=1,N
C---
      NINS(I) = 0
      INSLIM(I) = 0
C---
      DO 107 JJ=1,2
      DO 107 KK=1,5
107   NSRC1(JJ,KK)=0
C
C     SET UP THE INITIAL VALUES FOR EACH SURFACE.
C
      ALPO=Zero
      NRR=0
      NCC=0
      ITEST2=0
      ICAGE(1,1)=100
      ICAGE(1,2)=100
      ICAGE(1,3)=100
      ICAGE(2,1)=100
      ICAGE(2,2)=100
      ICAGE(2,3)=100
      DO 25 L=1,2
      IF(NR.EQ.0)GOTO 25
C
C     LOOK FOR RING CRIT. POINTS IN THE SURFACE.
C
      DO 26 K=1,NR
      IF(I.NE.NSRC(L,K))GOTO 26
C
C     LOOK FOR CAGE CRITICAL POINTS IN THE SURFACE.
C
      NRR=NRR+1
      NSRC1(1,NRR)=K
      NSRC1(2,NRR)=L
      IF(NC.EQ.0)GOTO 26
      DO 35 KK=1,NC
      K1=N+NR+KK
      IF(NSRC(3,K).NE.K1.AND.NSRC(4,K).NE.K1)GOTO 35
      NCC=NCC+1
C
C     STORE THE NUMBERS OF THE RING CRITICAL POINTS BORDERING
C     A CAGE.
C
      IF(ICAGE(1,KK).LT.99)ICAGE(2,KK)=NRR
      IF(ICAGE(1,KK).GE.99)ICAGE(1,KK)=NRR
35    CONTINUE
26    CONTINUE
25    CONTINUE
      IF(NCC.EQ.0)GOTO 37
      DO 50 J=1,NC
      IF(ICAGE(1,J).GE.99)GOTO 50
      DO 51 K=1,NRR
      IF(ICAGE(1,J).EQ.K)PH=ANGLE(1,NSRC1(2,K),NSRC1(1,K))
      IF(ICAGE(2,J).EQ.K)PH1=ANGLE(1,NSRC1(2,K),NSRC1(1,K))
51    CONTINUE
      IF(PH.LE.PH1)GOTO 52
      IF(PH1-PH+TPI.LT.PH-PH1)ITEST2=1
      IF(PH1-PH+TPI.LT.PH-PH1)GOTO 50
53    II=ICAGE(1,J)
      ICAGE(1,J)=ICAGE(2,J)
      ICAGE(2,J)=II
      GOTO 50
52    CONTINUE
C
C     ITEST2=1 INDICATES THAT WE ARE STARTING INSIDE A CAGE.
C
      IF(PH-PH1+TPI.LT.PH1-PH)ITEST2=1
      IF(PH-PH1+TPI.LT.PH1-PH)GOTO 53
50    IF(ITEST2.EQ.1)K3=J
      DO 90 L=1,3
90    IF(ITEST2.EQ.1)CC(L)=C1(L,N+NR+K3)
37    CONTINUE
      DO 1 J=1,NPATH1
C
C     NPATH IS THE NUMBER OF BASIC PATHS IN THE SURFACE WITH EQUIDISTANT
C     STARTING ANGLES.
C
31    ITEST1=0
      ITEST=0
      IF(J.NE.NPATH1)GOTO 14
C
C     THE LAST PATH IN THE SURFACE IS EQUAL TO THE FIRST.
C
      ALPN=TPI
      ALPB=TPI
      DO 15 K=1,NNN
      DO 15 L=1,3
15    XX1(L,K)=XX0(L,K)
      GOTO 24
C
C     COMPUTE THE INITIAL VALUES FOR THE NEXT PATH.
C
14    ALPN=DFLOAT(J-1)*TPI/DFLOAT(NPATH)
      ALPB=DFLOAT(J-1)*TPI/DFLOAT(NPATH)
      ALP=DFLOAT(J-1)*TPI/DFLOAT(NPATH)
24    IF(NRR.EQ.0.AND.J.EQ.NPATH1)GOTO 19
      IF(NRR.EQ.0)GOTO 17
C
C     DOES THE RING LIE BETWEEN THE TWO CONSIDERED PATHS ?
C
      DO 27 K=1,NRR
      IF(DABS(ALPN-ANGLE(1,NSRC1(2,K),NSRC1(1,K))).LT.Otmfiv)
     1 ANGLE(1,NSRC1(2,K),NSRC1(1,K))=ALPN
      PH=ANGLE(1,NSRC1(2,K),NSRC1(1,K))
      IF(PH.LT.Otmfiv.OR.PH.GT.TPIM)ANGLE(1,NSRC1(2,K),NSRC1(1,K))
     1 =TPI
      IF(PH.LT.OTMfiv.OR.PH.GT.TPIM)PH=Zero
      IF(ALPN.EQ.ZEro.AND.PH.EQ.Zero)GOTO 32
      IF(ALPN.GT.TPIM.AND.PH.EQ.Zero)GOTO 32
      IF(ALPN.LT.PH.OR.ALPO.GE.PH)GOTO 27
      ALPN=PH
32    ITEST1=3
      K4=K
C
C     THE HANDLING OF A RING CRIT. POINT .
C
      A=Otmfr*DCOS(PH)
      B=otmfr*DSIN(PH)
      DO 28 L=1,3
28    XX1(L,1)=A*V(L,I)+B*W(L,I)+C1(L,I)
      DS=DMAX/DFLOAT(NNN)
      IF(PH.EQ.Zero.AND.ALPN.EQ.Zero)GOTO 75
      IF(NCC.EQ.0)GOTO 75
      DO 74 L=1,NC
C
C     IS THIS RING CRITICAL POINT THE END OF A CAGE ?
C
      IF(ICAGE(2,L).NE.K)GOTO 74
      DO 76 LL=1,3
76    CC(LL)=C1(LL,N+NR+L)
      NN=NNN-4
      GOTO 63
74    CONTINUE
75    CONTINUE
C
C     COMPUTE POINTS ON THE PATH CONNECTING THE (3,-1) WITH
C     THE (3,+1) CRITICAL POINT.
C
      NN=INT(ANGLE(2,NSRC1(2,K),NSRC1(1,K))/DS)
      CALL TRUDGE2(XX1,DS,NN)
63    K2=NN+3
      K1=NSRC1(1,K)+N
      DO 42 L=1,3
42    IF(PH.EQ.Zero.AND.ALPN.EQ.Zero)XX2(L,K2)=C1(L,I)+
     $DCOS(Pt7)*V(L,I)+
     1 DSIN(Pt7)*W(L,I)
      A=(V(1,K1)-XX2(1,K2)+C1(1,K1))**2+(V(2,K1)-
     $XX2(2,K2)+C1(2,K1))**2
     1 +(V(3,K1)-XX2(3,K2)+C1(3,K1))**2
      B=(V(1,K1)+XX2(1,K2)-C1(1,K1))**2+(V(2,K1)+
     $XX2(2,K2)-C1(2,K1))**2
     1 +(V(3,K1)+XX2(3,K2)-C1(3,K1))**2
      VZ=One
      IF(B.LT.A)VZ=-One
      DO 29 L=1,3
29    XX1(L,NN+1)=Otmfr*VZ*V(L,K1)+C1(L,K1)
      NN1=NNN2-NN
      IF(PH.EQ.Zero.AND.ALPN.EQ.Zero)GOTO 55
      IF(PH.EQ.Zero)PH=Tpi
      IF(NCC.EQ.0)GOTO 55
      DO 54 L=1,NC
C
C     IF THE RING IS THE END OF A CAGE, COMPUTE POINTS ON THE
C     PATH CONNECTING THE (3,-1) WITH THE (3,+1) CRITICAL
C     POINT USING A NEW STEPSIZE.
C
      IF(ICAGE(2,L).NE.K)GOTO 54
      NNRL=N+NR+L
      CALL CAGE(XX1(1,NN+1),C1,NNRL,QL)
      DS=(QL+ANGLE(2,NSRC1(2,K),NSRC1(1,K)))/DFLOAT(NNN1)
      NN=INT(ANGLE(2,NSRC1(2,K),NSRC1(1,K))/DS)
      NN1=NNN2-NN
      NN2=NN-1
      CALL TRUDGE2(XX1,DS,NN2)
      DO 78 LL=1,3
78    XX1(LL,NN+1)=XX1(LL,NNN-3)
54    CONTINUE
55    CONTINUE
      IF(PH.NE.Zero)GOTO 62
      IF(NCC.EQ.0)GOTO 62
      DO 60 L=1,NC
C
C     IS THERE A RING WITH ANGLE 0. AND IS IT THE
C     BEGINNING OF A CAGE ?
C
      IF(ICAGE(1,L).NE.K)GOTO 60
      ITEST2=1
      DO 61 LL=1,3
61    CC(LL)=C1(LL,N+NR+L)
      K3=L
      NNRL=N+NR+L
      CALL CAGE(XX1(1,NN+1),C1,NNRL,QL)
      DS=(QL+ANGLE(2,NSRC1(2,K),NSRC1(1,K)))/DFLOAT(NNN1)
      NNN=NN
      NN=INT(ANGLE(2,NSRC1(2,K),NSRC1(1,K))/DS)
      DO 77 LL=1,3
77    XX1(LL,NN+1)=XX1(LL,NNN+1)
      NN2=NN-1
      CALL TRUDGE2(XX1,DS,NN2)
      GOTO 62
60    CONTINUE
      ITEST2=0
62    CONTINUE
C
C     COMPUTE POINTS ON THE PATH FROM THE RING CRITICAL POINT
C     TO INFINITY OR TO A CAGE CRITICAL POINT.
C
      CALL TRUDGE2(XX1(1,NN+1),DS,NN1)
      DO 30 L=1,3
      IF(ITEST2.EQ.1)XX1(L,NNN)=CC(L)
30    XX1(L,NN+1)=C1(L,K1)
      IF(PH.EQ.Zero)GOTO 43
      GOTO 19
27    CONTINUE
      IF(J.EQ.NPATH1)GOTO 19
C
C     SET UP INITIAL VALUES FOR A PATH.
C
17    A=Otmfr*DCOS(ALP)
      B=Otmfr*DSIN(ALP)
      DO 3 L=1,3
      XX1(L,1)=A*V(L,I)+B*W(L,I)+C1(L,I)
C----
 3    CONTINUE
C
C     IS THE CONSIDERED PATH INSIDE A CAGE ?
C
      DS=DMAX/DFLOAT(NNN)
      IF(ITEST2.EQ.0)GOTO 36
      NNRL=N+NR+K3
      CALL CAGE(XX1,C1,NNRL,QL)
      DS=QL/DFLOAT(NNN)
36    CALL TRUDGE2(XX1,DS,NNN2)
      DO 38 L=1,3
38    IF(ITEST2.EQ.1)XX1(L,NNN)=CC(L)
      IF(J.NE.1)GOTO 19
C
C     STORE THE FIRST (=LAST) PATH.
C
43    DO 5 K=1,NNN
      DO 5 L=1,3
5     XX0(L,K)=XX1(L,K)
      ITEST1=0
      GOTO 6
C
C     DO WE HAVE TO INSERT ANOTHER PATH BETWEEN THE TWO CONSIDERED ?
C
19    G=(XX2(1,NNN)-XX1(1,NNN))**2+(XX2(2,NNN)-XX1(2,NNN))**2+
     1 (XX2(3,NNN)-XX1(3,NNN))**2
      G=DSQRT(G)
      IF(G.GT.DIST.AND.ITEST.GT.0)GOTO 23
      IF(G.LE.DIST.AND.ITEST.EQ.0)GOTO 4
      IF(G.LE.DIST.AND.ITEST.GT.0)ITEST=-1
      IF(ITEST.LT.0)GOTO 4
C
C     STORE THE ALREADY COMPUTED PATH
C
      DO 16 K=1,NNN
      DO 16 L=1,3
16    XX3(L,K)=XX1(L,K)
23    ITEST=ITEST+1
C
C     LIMIT THE MAXIMAL NUMBER OF PATHS TO BE INSERTED BETWEEN
C     TWO BASIC PATHS.
C
      IF(ITEST.GT.INMAX) THEN
       ITEST = -1
       INSLIM(I) = INSLIM(I)+1
      ENDIF
      IF(DABS(ALPN-ALPO).LE.ADIFF)ITEST=-1
      IF(ITEST.LT.0)GOTO 4
C
C     SET UP THE INITIAL VALUES FOR THE INSERTED PATH.
C
      ALP=ALPO+(ALPN-ALPO)/Two**(DFLOAT(ITEST))
C
      NINS(I) =  NINS(I) + 1
C
      GOTO 17
C
C     PARTITION THE AREA BETWEEN TWO PATH INTO TRIANGLES.
C---
   4  CONTINUE
C
      DO 8 L=1,3
      AA(L,1)=C1(L,I)
      AA(L,2)=XX1(L,1)
      AA(L,3)=XX2(L,1)
8     CONTINUE
      CALL TRIA(AA,XS,XS1,IZ,IPHIPL,IT)
      DO 10 K=1,NNN2
      DO 11 L=1,3
      AA(L,1)=XX1(L,K)
      AA(L,2)=XX1(L,K+1)
11    AA(L,3)=XX2(L,K)
      CALL TRIA(AA,XS,XS1,IZ,IPHIPL,IT)
      IF(K+1.EQ.NNN)GOTO 6
      DO 12 L=1,3
      AA(L,1)=XX1(L,K+1)
      AA(L,2)=XX2(L,K)
12    AA(L,3)=XX2(L,K+1)
10    CALL TRIA(AA,XS,XS1,IZ,IPHIPL,IT)
C
C     FORGET THE FIRST PATH AND STORE THE SECOND.
C
6     DO 13 K=1,NNN
      DO 13 L=1,3
      XX2(L,K)=XX1(L,K)
13    IF(ITEST.EQ.-1)XX1(L,K)=XX3(L,K)
      IF(ITEST.EQ.0)GOTO 18
C
C     HANDLE THE CASE OF RING AND CAGE CRIT. POINTS.
C
      ITEST=0
      ALPO=ALP
      GOTO 19
18    ALPO=ALPN
      IF(ITEST1.EQ.0)GOTO 1
      IF(PH.GT.TPIM)GOTO 1
      DO 33 L=1,3
33    XX2(L,NN+1)=C1(L,K1)-OtmFr*VZ*V(L,K1)
      DS=DMAX/DFLOAT(NNN)
      IF(NCC.EQ.0)GOTO 59
      K=K4
      DO 57 L=1,NC
C
C     IS THE RING THE BEGINNING OF A CAGE ?
C
      IF(ICAGE(1,L).NE.K)GOTO 57
      ITEST2=1
      DO 58 LL=1,3
58    CC(LL)=C1(LL,N+NR+L)
      K3=L
      NNRL=N+NR+L
      CALL CAGE(XX2(1,NN+1),C1,NNRL,QL)
      DS=(QL+ANGLE(2,NSRC1(2,K),NSRC1(1,K)))/DFLOAT(NNN1)
      NN=INT(ANGLE(2,NSRC1(2,K),NSRC1(1,K))/DS)
      NN1=NNN2-NN
      NN2=NN-1
      DO 79 LL=1,3
79    XX2(LL,NN+1)=C1(LL,K1)-Otmfr*VZ*V(LL,K1)
      CALL TRUDGE2(XX2,DS,NN2)
      GOTO 59
57    CONTINUE
      IF(ITEST2.EQ.0)GOTO 59
C
C     IF THE RING IS THE END OF A CAGE, CHANGE THE STEPSIZE
C     BACK TO NORMAL.
C
      ITEST2=0
      NN=INT(ANGLE(2,NSRC1(2,K),NSRC1(1,K))/DS)
      NN1=NNN2-NN
      NN2=NN-1
      CALL TRUDGE2(XX2,DS,NN2)
      DO 80 LL=1,3
80    XX2(LL,NN+1)=C1(LL,K1)-Otmfr*VZ*V(LL,K1)
59    CONTINUE
      CALL TRUDGE2(XX2(1,NN+1),DS,NN1)
      DO 34 L=1,3
      IF(ITEST2.EQ.1)XX2(L,NNN)=CC(L)
34    XX2(L,NN+1)=C1(L,K1)
      IF(DABS(ALPB-ALPN).LT.Otmfiv)GOTO 1
      GOTO 31
1     CONTINUE
 999  CONTINUE
C
C     ORDER THE POINTS OF INTERSECTION ACCORDING To THEIR
C     MAGNITUDE.
C
      DO 2 I=1,IT
      DO 2 J=1,IPHIPL
      IF(IZ(I,J).LE.1)GOTO 2
      DO 41 L=1,2
      LL=3-L
      DO 40 K=1,LL
      IF(ZI(K,I,J).LE.ZI(K+1,I,J))GOTO 40
      A=ZI(K,I,J)
      ZI(K,I,J)=ZI(K+1,I,J)
      ZI(K+1,I,J)=A
40    CONTINUE
41    CONTINUE
2     CONTINUE
      RETURN
      END
      SUBROUTINE TRIA(A,XS,XS1,IZ,IP,IT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter(MaxTht=200,MaxPhi=200)
      COMMON /S/ XI(MaxTht),YI(MaxPhi),ZI(3,MaxTht,MaxPhi)
      COMMON /UNITS/ ISRF,INPT,IOUT,IWFN,IDBG
      DIMENSION A(3,3),XS(3,MaxTht,MaxPhi),XS1(3,MaxTht*MaxPhi),
     $AB(3,3),IZ(MaxTht,MaxPhi),LAA(MaxTht*MaxPhi),
     $MAA(MaxTht*MaxPhi),IPVT(3),Anew(3,3)
      Data alzero/1.0D-13/,three/3.0D0/,zero/0.D0/,zp1/0.1D0/,
     $d1min9/1.0D-9/,one/1.0D0/,two/2.0d0/,seven/7.0d0/,eight/8.0d0/
1999  Format('Degenerate Triangle in Tria ...')
c
      do 1234 j = 1,3
        do 1235 k = 1,3
          anew(k,j) = a(j,k)
1235    continue
1234  continue
C
      CALL DGEFA(A,3,3,IPVT,IER)
      IF(IER.GT.100)Then
      Write(Iout,1999)
      Return
      Endif
C
      PI = DACOS(-ONE)
      PI2 = TWO*PI
      PI3 = PI/TWO
C
      DO 4 L=1,3
      AB(L,3)=alzero+DSQRT(Anew(L,1)*anew(L,1)+
     $                     Anew(L,2)*anew(L,2)+Anew(L,3)*anew(L,3))
      AB(L,1)=DACOS(Anew(L,3)/AB(L,3))
      IF(DABS(Anew(L,1)).LT.d1min9) Then
      AB(L,2)= pi3
      IF(Anew(L,2).LT.zero)AB(L,2)=three*AB(L,2)
      Else
      AB(L,2)=ATAN(Anew(L,2)/Anew(L,1))
      IF(Anew(L,1).LT.zero) AB(L,2)=AB(L,2)+PI
      IF(Anew(L,1).GE.zero.AND.Anew(L,2).LT.zero) AB(L,2)=AB(L,2)+PI2
      Endif
4     CONTINUE
C
      DT=zero
      DP=zero
      DT1=PI
      DP1=PI2
      PI4=Seven*PI/Eight
C
      IF(AB(1,1).LE.PI4.AND.AB(2,1).LE.PI4.AND.AB(3,1).LE.PI4)
     1 DT1=DMAX1(AB(1,1),AB(2,1),AB(3,1))
      PI4=One*Pi/Eight
      IF(AB(1,1).GE.PI4.AND.AB(2,1).GE.PI4.AND.AB(3,1).GE.PI4)
     1 DT=DMIN1(AB(1,1),AB(2,1),AB(3,1))
      IF(AB(1,2).LE.PI.AND.AB(2,2).LE.PI.AND.AB(3,2).LE.PI)Then
      DP1=DMAX1(AB(1,2),AB(2,2),AB(3,2))
      DP=DMIN1(AB(1,2),AB(2,2),AB(3,2))
      Endif
      IF(AB(1,2).GE.PI.AND.AB(2,2).GE.PI.AND.AB(3,2).GE.PI)Then
      DP=DMIN1(AB(1,2),AB(2,2),AB(3,2))
      DP1=DMAX1(AB(1,2),AB(2,2),AB(3,2))
      Endif
C
      KK=0
      DO 7 J=1,IP
      IF(DP.LE.YI(J).AND.DP1.GE.YI(J))THEN
      DO 88 I=1,IT
      IF(DT.LE.XI(I).AND.DT1.GE.XI(I))THEN
      KK=KK+1
      LAA(KK)=I
      MAA(KK)=J 
      ENDIF
88    CONTINUE
      ENDIF
7     CONTINUE
C
      IF (KK .EQ. 0) RETURN
C
      DO 8 I=1,KK
      DO 9 L=1,3
9     XS1(L,I)=XS(L,LAA(I),MAA(I))
      CALL DGESL(A,3,3,IPVT,XS1(1,I),0)
      IF(XS1(1,I).GE.zero.and.XS1(2,I).Ge.Zero.and.
     $XS1(3,I).GE.Zero)Then
      H=XS1(1,I)+XS1(2,I)+XS1(3,I)
      recH = one/H
      IF(H.GT.zp1.AND.IZ(LAA(I),MAA(I)).EQ.0)
     +ZI(1,LAA(I),MAA(I))=recH
      IF(H.GT.zp1.AND.IZ(LAA(I),MAA(I)).EQ.1)
     +ZI(2,LAA(I),MAA(I))=recH
      IF(H.GT.zp1.AND.IZ(LAA(I),MAA(I)).EQ.2)
     +ZI(3,LAA(I),MAA(I))=recH
      IF(H.GT.zp1)IZ(LAA(I),MAA(I))=IZ(LAA(I),MAA(I))+1
      Endif
8     CONTINUE
C
      RETURN
      END
      SUBROUTINE TRUDGE2(XYZ,DS,NNN)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      Parameter (Maxpp=300)
      DIMENSION XYZ(3,Maxpp),W(3,Maxpp),W1(3,6),R(3)
C
      DO 100 M=1,5
      R(1)=XYZ(1,M)
      R(2)=XYZ(2,M)
      R(3)=XYZ(3,M)
      CALL GAUS2(R,0) 
      DSQ = DSQRT(R(1)**2 + R(2)**2 + R(3)**2)
      W(1,M)=R(1)/DSQ
      W(2,M)=R(2)/DSQ
      W(3,M)=R(3)/DSQ
      CALL RUNGE(XYZ(1,M),XYZ(1,M+1),W(1,M),DS)
100   CONTINUE
C
      DO 200 M=6,NNN
      R(1)=XYZ(1,M)
      R(2)=XYZ(2,M)
      R(3)=XYZ(3,M)
      CALL GAUS2(R,0) 
      DSQ= DSQRT(R(1)**2 + R(2)**2 + R(3)**2)
      W(1,M)=R(1)/DSQ
      W(2,M)=R(2)/DSQ
      W(3,M)=R(3)/DSQ
      DO 250 J=1,6
      W1(1,J)=W(1,M-J+1)
      W1(2,J)=W(2,M-J+1)
      W1(3,J)=W(3,M-J+1)
250   CONTINUE
      CALL MULTI1(XYZ(1,M),XYZ(1,M+1),W1,6,DS,M*DS)
200   CONTINUE
C
      RETURN
      END
      SUBROUTINE TRUDGE3(XYZA,NNPTH,ITERMT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL RHF,ROHF,ABNAT,RNAT
      Parameter (MaxTht=200,MaxPhi=200)
      PARAMETER (MCENT=50,MNPT=MaxTht*MaxPhi,MPTS=500)
      COMMON /DIRK/ TOLL,TOLLA(MCENT),STP,SIZE,FSTP,CTF,SSTP,NOSEC,
     $ISECT,ICP,ISTEP,NABMO,Nacc,IDOAOM,IDOMAG,IMagm,nfbeta,RHF,ROHF,
     $ABNAT,RNAT,XGO,YGO,ZGO
      COMMON /C2/ AB(9,9),AM(9,10)
      COMMON /C7/ CT(3,3),X, Y, Z,NCRNT
      COMMON /ATOMS/ XC(MCENT),YC(MCENT),ZC(MCENT),CHARG(MCENT),NCENT
      COMMON /DER/ DXR(MPTS),DYR(MPTS),DZR(MPTS)
      COMMON /UNITS/ ISRF,INPT,IOUT,IWFN,IDBG
      DIMENSION XYZA(MNPT,3),ITERM(MPTS),XYZG(MPTS,3),XYZ(MPTS,3),
     $DX(MPTS,9),DY(MPTS,9),DZ(MPTS,9),DX0(MPTS),DY0(MPTS),OSC(MPTS),
     $DZ0(MPTS),L(MPTS),MPS(MPTS),ITERMT(MNPT),MX(MPTS),XYZT(3),
     $Temp(3)
      Save Zero,One,MaxOsc
      Data Zero/0.0d0/,One/1.0d0/,MaxOsc/10/
1999  Format('Excessive Oscillations near: X= ',1PE8.2,' Y= ',
     $1PE8.2,' Z= ',1PE8.2)
2999  Format('Nonnuclear Attractor? ... Or too far out?')
C
      STPR=STP/ISTEP
      NabMo1=Nabmo-1
C
      DO 111 I=1,NNPTH
      ITERMT(I)=0
111   CONTINUE
C
      IDIV=NNPTH/MPTS
      IMINUS=NNPTH-IDIV*MPTS 
      DO 98 I=1,IDIV
      MPS(I)=MPTS
98    CONTINUE
      MPS(IDIV+1)=IMINUS
      ISEND=IDIV+1
      IF(IMINUS.EQ.0)ISEND=IDIV
C
      DO 1 JJ=1,ISEND
      JL=(JJ-1)*MPTS
      DO 2 KK=1,MPS(JJ)
      XYZ(KK,1)=XYZA(JL+KK,1)
      XYZ(KK,2)=XYZA(JL+KK,2)
      XYZ(KK,3)=XYZA(JL+KK,3)
      MX(KK)=0
      OSC(KK)=0
      DO 3 JR=1,NabMo
      DX(KK,JR)=Zero
      DY(KK,JR)=Zero
      DZ(KK,JR)=Zero
3     CONTINUE
2     CONTINUE
      NPTH=MPS(JJ)
      LC=0
      K=0
C
      DO 1100 IK=NabMo,1,-1
      DO 1105 JR=1,NPTH
      L(JR)=0
      ITERM(JR)=0
1105  Continue
      DO 1110 IL=1,ISTEP
C
      CALL GAUS3(XYZ,NPTH,2)
C
      DO 1128 I = 1,NPTH
      DSQ=STPR/(DSQRT(DXR(I)**2+DYR(I)**2+DZR(I)**2))
      XYZ(I,1) = XYZ(I,1) + DSQ*DXR(I) 
      XYZ(I,2) = XYZ(I,2) + DSQ*DYR(I) 
      XYZ(I,3) = XYZ(I,3) + DSQ*DZR(I) 
1128  CONTINUE
C
1110  CONTINUE
C
      DO 1127 I = 1,NPTH
      DO 1129  JJJ=1,NCENT
      DOTHER=(XYZ(I,1)-XC(JJJ))**2+
     +(XYZ(I,2)-YC(JJJ))**2+
     +(XYZ(I,3)-ZC(JJJ))**2
      IF(DOTHER.LE.TOLLA(JJJ))THEN
      L(I) = 1
      If(JJJ.Eq.NCrnt)ITERM(I) = 1
      If(JJJ.Ne.NCrnt)ITERM(I) = 0
      K = K +1
      LC=LC+1
      GOTO 1127
      ENDIF
1129  CONTINUE
1127  CONTINUE 
C
      NLM=0
      NLP=0
      DO 1500 J=1,NPTH
      IF(L(J).EQ.0) THEN
      NLP=NLP+1
      XYZ(NLP,1)=XYZ(J,1) 
      XYZ(NLP,2)=XYZ(J,2) 
      XYZ(NLP,3)=XYZ(J,3) 
      DO 1501 KB=1,NabMo
      DX(NLP,KB)=DX(J,KB)
      DY(NLP,KB)=DY(J,KB)
      DZ(NLP,KB)=DZ(J,KB)
1501  CONTINUE
      MX(NLP)=MX(J)+NLM
      ELSE
      NLM=NLM+1
      ITERMT(JL+J+MX(J))=ITERM(J)
      ENDIF
1500  CONTINUE
C
      NPTH=NPTH-LC
      LC=0
      IF(K.EQ.MPS(JJ))GOTO 1
C
      CALL GAUS3(XYZ,NPTH,2)
C
      DO 1200 I=1,NPTH
      Dsq=One/(Dsqrt(DXR(I)**2+DYR(I)**2+DZR(I)**2))
      DX(I,IK)=DXR(I)*Dsq
      DY(I,IK)=DYR(I)*Dsq
      DZ(I,IK)=DZR(I)*Dsq
1200  CONTINUE
C
1100  CONTINUE
C
      DO 26 I=1,NPTH
      DX0(I) = Zero
      DY0(I) = Zero
      DZ0(I) = Zero
      XYZG(I,1)=XYZ(I,1)
      XYZG(I,2)=XYZ(I,2)
      XYZG(I,3)=XYZ(I,3)
26    CONTINUE
C
30    CONTINUE
C
      CALL GAUS3(XYZ,NPTH,2)
C
      DO 116 J=1,NPTH
      DSQ=One/(DSQRT(DXR(J)**2+DYR(J)**2+DZR(J)**2))
      DX(J,1) = DXR(J)*DSQ
      DY(J,1) = DYR(J)*DSQ
      DZ(J,1) = DZR(J)*DSQ
      L(J) = 0
      ITERM(J)=0
116   CONTINUE
C
      DO 517 I=1,NPTH
      IF((DX0(I)*DX(I,1)+DY0(I)*DY(I,1)+DZ0(I)*DZ(I,1)).LT.Zero)
     +THEN
      ICP=ICP+1
      OSC(I)=OSC(I)+1
      If(OSC(I).GT.MaxOsc)Then
      Write(Iout,1999)Xyzg(I,1),Xyzg(I,2),Xyzg(I,3)
      Write(Iout,2999)
      Stop 'Aborted in trudge'
      Endif
      XYZT(1)=XYZG(I,1)
      XYZT(2)=XYZG(I,2)
      XYZT(3)=XYZG(I,3)
C
      DO 58 LL=1,ISTEP
      Temp(1)=Xyzt(1)
      Temp(2)=Xyzt(2)
      Temp(3)=Xyzt(3)
      CALL GAUS2(XYZT,1)
      DS = STPR/(DSQRT(XYZT(1)**2 + XYZT(2)**2 + XYZT(3)**2))
      XYZT(1) = Temp(1) + DS*XYZT(1)
      XYZT(2) = Temp(2) + DS*XYZT(2)
      XYZT(3) = Temp(3) + DS*XYZT(3)
58    CONTINUE
      XYZ(I,1)=XYZT(1)
      XYZ(I,2)=XYZT(2)
      XYZ(I,3)=XYZT(3)
      CALL GAUS2(XYZT,1)
      DS = One/(DSQRT(XYZT(1)**2 + XYZT(2)**2 + XYZT(3)**2))
      DX(I,1)=XYZT(1)*DS
      DY(I,1)=XYZT(2)*DS
      DZ(I,1)=XYZT(3)*DS
      ENDIF
517   CONTINUE
C
      DO 127 I = 1,NPTH
      XYZG(I,1) = XYZ(I,1)
      XYZG(I,2) = XYZ(I,2)
      XYZG(I,3) = XYZ(I,3)
      DO 129  JJJ=1,NCENT
      DOTHER=(XYZ(I,1)-XC(JJJ))**2+
     +(XYZ(I,2)-YC(JJJ))**2+
     +(XYZ(I,3)-ZC(JJJ))**2
      IF(DOTHER.LE.TOLLA(JJJ))THEN
      L(I) = 1
      IF(JJJ.EQ.Ncrnt)ITERM(I) = 1
      IF(JJJ.Ne.Ncrnt)ITERM(I) = 0
      K = K +1
      LC=LC+1
      GOTO 127
      ENDIF
129   CONTINUE
127   CONTINUE 
C
      NLM=0
      NLP=0
      DO 500 J=1,NPTH
      IF(L(J).EQ.0) THEN
      NLP=NLP+1
      XYZ(NLP,1)=XYZ(J,1) 
      XYZ(NLP,2)=XYZ(J,2) 
      XYZ(NLP,3)=XYZ(J,3) 
      XYZG(NLP,1)=XYZG(J,1) 
      XYZG(NLP,2)=XYZG(J,2) 
      XYZG(NLP,3)=XYZG(J,3) 
      DX0(NLP)=DX(J,1)
      DY0(NLP)=DY(J,1)
      DZ0(NLP)=DZ(J,1)
      OSC(NLP)=OSC(J)
      DO 501 KB=1,NabMo
      DX(NLP,KB)=DX(J,KB)
      DY(NLP,KB)=DY(J,KB)
      DZ(NLP,KB)=DZ(J,KB)
501   CONTINUE
      MX(NLP)=MX(J)+NLM
      ELSE
      NLM=NLM+1
      ITERMT(JL+J+MX(J))=ITERM(J)
      ENDIF
500   CONTINUE
C
      NPTH=NPTH-LC
      LC=0
      IF(K.EQ.MPS(JJ))GOTO 1
C
      DO 128 I = 1,NPTH
      HX=Zero
      HY=Zero
      HZ=Zero
      DO 228 KW=1,Nabmo
      HX=HX+DX(I,KW)*AB(NabMo,KW)
      HY=HY+DY(I,KW)*AB(NabMo,KW)
      HZ=HZ+DZ(I,KW)*AB(NabMo,KW)
228   CONTINUE
      XYZ(I,1) = XYZ(I,1) + STP*HX
      XYZ(I,2) = XYZ(I,2) + STP*HY
      XYZ(I,3) = XYZ(I,3) + STP*HZ
128   CONTINUE
C
      CALL GAUS3(XYZ,NPTH,2)
C
      DO 1428 I = 1,NPTH
      DSQ=AM(NABMO,1)/(DSQRT(DXR(I)**2+DYR(I)**2+DZR(I)**2))
      HX=DXR(I)*DSQ
      HY=DYR(I)*DSQ
      HZ=DZR(I)*DSQ
      DO 1228 KW=1,NabMo
      HX=HX+DX(I,KW)*AM(NabMo,KW+1)
      HY=HY+DY(I,KW)*AM(NabMo,KW+1)
      HZ=HZ+DZ(I,KW)*AM(NabMo,KW+1)
1228   CONTINUE
C
      XYZ(I,1) = XYZG(I,1) + STP*HX
      XYZ(I,2) = XYZG(I,2) + STP*HY
      XYZ(I,3) = XYZG(I,3) + STP*HZ
1428  CONTINUE
C
      DO 1528 I=1,NPTH
      DO 1628 KW=NabMo1,1,-1
      DX(I,KW+1)=DX(I,KW)
      DY(I,KW+1)=DY(I,KW)
      DZ(I,KW+1)=DZ(I,KW)
1628  CONTINUE
1528  CONTINUE

      GOTO 30
1     CONTINUE    
C
      RETURN
      END

