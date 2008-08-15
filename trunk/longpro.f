C    **************************PROGRAM LONGPRO*************************
C    Program LONGPRO calculates the dynamical evolution of a
C    stream's longitudinal profile
C**********************************************************************
      IMPLICIT NONE
      REAL*8 DELX,RLONG,DELT,XPRINT,TPRINT,WID(400),ELEV(400),MFEED
      REAL*8 STRDATA(6,400),VISC,RHO,VISKIN,DIMID,Q(400),TOX(400),G
      REAL*8 WSEBC,MANN,MFLUX(400),DECREM,TIME,FACTOR,INITELEV
      REAL*8 RHOM,RHOA,RHOC,VOLCON
      REAL*8 SIGMA,ALPHA,RATE,QMAX,OMEGA,REY,SEDCON
      INTEGER*4 NTIM,PRINT(5),NT,M,K,TURNON1,TURNOFF1,
     *          TURNON2,TURNOFF2,HINGE
      CHARACTER*24 FIN,FOUT1,FOUT
      CHARACTER*25 RUN
      CHARACTER*4 RUNNAME
C----------------------------------------------------------------------
C     COMMON BLOCK DECLARATIONS
C----------------------------------------------------------------------
      COMMON /COM1/M,NT,DELX,STRDATA,WID,ELEV,WSEBC,MANN
      COMMON /COM2/Q,DELT,XPRINT,TPRINT,TOX,RHO,G,RUN
      COMMON /COM3/VISKIN,DIMID,OMEGA,REY
      COMMON /COM4/SIGMA,VOLCON,RHOC,RHOM,RHOA,ALPHA
      COMMON /TECTO/RATE,TURNON1,TURNOFF1,TURNON2,TURNOFF2,HINGE
C----------------------------------------------------------------------
C     FILE MANAGEMENT PORTION : WHICH FILES ARE OPENED AS INPUT        
C     OR OUTPUT FILES                                                  
C----------------------------------------------------------------------
      OPEN(3, FILE='./files')
      READ (3, 109) FIN
      READ (3, 109) FOUT
      OPEN(4, FILE='./'//FIN)
      READ (4, 107) RUNNAME
      DO 1 K = 1, 24
         IF (FOUT(K:K) .EQ. ' ') GO TO 2
    1 CONTINUE
    2 CONTINUE
      K = K - 1
      FOUT1 = FOUT(1:K)//'A'
      OPEN(8, FILE='./'//FOUT1)
C----------------------------------------------------------------------
C     READ INPUT PARAMETERS                                            
C----------------------------------------------------------------------
      READ (4, 101) RLONG, DELX, DELT
      READ (4, 102) NTIM
      READ (4, 104) XPRINT, TPRINT
      READ (4, 100) DIMID
      READ (4, 108) RUN
      READ (4, 104) WSEBC, DECREM
      READ (4, 104) QMAX, MANN
      READ (4, 104) FACTOR, INITELEV
      READ (4, 106) TURNON1, TURNOFF1, TURNON2, TURNOFF2
      READ (4, 102) HINGE
      READ (4, 101) RATE, SEDCON, MFEED
  100   FORMAT(59x,f10.0)
  101   FORMAT(2(59X,F10.0,/),59X,F10.0)                    
  102   FORMAT(59X,I10)                                    
  103   FORMAT(/,5I10)                                    
  104   FORMAT(59X,F10.2,/,59X,F10.2)
  105   FORMAT(59X,E10.2)
  106   FORMAT(3(59X,I10,/),59X,I10)
  107   FORMAT(/,A4)
  108   FORMAT(/,A25)
  109   FORMAT(A24)
C----------------------------------------------------------------------
C     SET SOME INITIAL CONSTANTS
C----------------------------------------------------------------------
      DELT = DELT*3600.0
      DELX = DELX*1000.0
      XPRINT = XPRINT*1000.0
      TPRINT = TPRINT*3600.0
      RLONG = RLONG*1000.0
      M = INT(RLONG/DELX+0.1) + 1
      VISC = 0.0018
      RHO = 1000.0
      G = 9.80
      VISKIN = VISC/RHO
      SIGMA = 2650.0
      VOLCON = 0.70
      RHOC = 2650.0
      RHOA = 0.0
      RHOM = 3300.0
C----------------------------------------------------------------------
C    SET THE INITIAL BED ELEVATIONS, WIDTHS, AND SETTLING VEL
C----------------------------------------------------------------------
      CALL SETUP (INITELEV, DECREM, QMAX)
      CALL SETTLE
      REY = OMEGA*DIMID/VISKIN
C----------------------------------------------------------------------
C     WRITE OUT INITIAL PARAMETERS                                     
C----------------------------------------------------------------------
      IF (PRINT(2) .EQ. 1) THEN
         WRITE (7, 110) RUN
  110     FORMAT(' ',A25//)                                           
      ENDIF
C----------------------------------------------------------------------
C     BEGIN TIME LOOP
C----------------------------------------------------------------------
      DO 3 NT = 1, NTIM
         OPEN(19, FILE='./tst')
         WRITE (19, 111) NT
  111    FORMAT(1X,'TIMESTEP   ',I9)                                  
         CLOSE(19)
C----------------------------------------------------------------------
C     CALCULATE TECTONICS
C----------------------------------------------------------------------
         CALL TECTON
C----------------------------------------------------------------------
C     CALCULATE GRADUALLY VARIED FLOW
C----------------------------------------------------------------------
         CALL FLDTA
C----------------------------------------------------------------------
C     CALCULATE SED TRANSPORT AND AMOUNT OF EROSION OR DEPOSITION
C----------------------------------------------------------------------
         CALL ERODEP90 (MFLUX, FACTOR, SEDCON, MFEED)
C----------------------------------------------------------------------
C     WRITE OUT VARIABLES
C----------------------------------------------------------------------
         TIME = DELT*NT/TPRINT + 0.00001
         IF (TIME - FLOAT(INT(TIME)) .LE. 0.0001) CALL WRFLDTA
    3 CONTINUE
      STOP 
      END
 
 
***********************************************************************
*               SET INITIAL ELEVATION, DISCHARGES, AND WIDTHS         *
***********************************************************************
      SUBROUTINE SETUP(INITELEV,DECREM,QMAX)
      IMPLICIT NONE
      REAL*8 RHO,G,Q(400),DELT,XPRINT,TPRINT,TOX(400),DELX,
     *       STRDATA(6,400),WID(400),ELEV(400),WSEBC,MANN
      REAL*8 DECREM,QMAX,WMAX,INITELEV
      INTEGER*4 M,NT,JN
      CHARACTER*25 RUN
C----------------------------------------------------------------------
C     COMMON BLOCK DECLARATIONS
C----------------------------------------------------------------------
      COMMON /COM1/M,NT,DELX,STRDATA,WID,ELEV,WSEBC,MANN
      COMMON /COM2/Q,DELT,XPRINT,TPRINT,TOX,RHO,G,RUN
C***********************************************************************
C     ASSIGN INITIAL ELEVATIONS
C***********************************************************************
      DO JN = 1, M
         ELEV(JN) = INITELEV - DECREM*JN*DELX
      END DO
C***********************************************************************
C     DEFINE WIDTHS AND DISCHARGES AS FUNCTIONS OF X
C***********************************************************************
      WMAX = QMAX/4.0
      DO JN = 1, M
C----------------------------------------------------------------------
C     UNCOMMENT THE FOLLOWING LINE FOR A LINEARLY INCREASING 
C     DISCHARGE DOWNSTREAM
C----------------------------------------------------------------------
C       Q(JN)=QMAX*(FLOATJ(JN)/M)
        Q(JN) = QMAX
        WID(JN) = WMAX*(Q(JN)/QMAX)**0.5
      END DO
      RETURN 
      END
 
 
 
***********************************************************************
*                1-D GRADUALLY VARIED CHANNEL FLOW                    *
***********************************************************************
*  This subroutine calculates flow depth and cross-sectional average  *
*  flow velocity at n rectangular cross sections along a single       *
*  thread river using the standard step solution method for the       *
*  gradually varied flow equation.  See Henderson (1966), page 143,   *
*  for a discussion of the methodology. --R Sling; Fall 1989          *
***********************************************************************
*  The set up constants are:                                          *
*    M      ---Number of cross sections                               *
*    NT     ---Number of time steps                                   *
*    DELX   ---Space descretization step (km)                         *
*    MANN   ---Manning's n for allcross sections                      *
*    WID    ---Width of the channel at each cross section (m)         *
*    ELEV   ---Elevation of the channel at each cross section in      *
*                meters above a datum                                 *
*    WSEBC  ---Water surface elevation for all time at control node   *
*    Q      ---Water discharge (m**3 s**-1)during timestep I          *
*  The internal variables are:                                        *
*    A      ---cross sectional area                                   *
*    P      ---Wetted Perimeter (m)                                   *
*    V      ---Cross sectional mean velocity (m s**-1)                *
*    Sf     ---Friction slope                                         *
*    H      ---Hydraulic head at a section (m)                        *
*    H21    ---Error in Heads (m)                                     *
*    RCH    ---Hydraulic radius (m)                                   *
*    DCH    ---Flow depth (m)                                         *
***********************************************************************
      SUBROUTINE FLDTA
      IMPLICIT NONE
      REAL*8 Q(400),MANN,WSEBC,A,P,V,SF1,H1,H2,H21,RCH,D2,P2
      REAL*8 DELX,SF2,A2,FR,WID(400),ELEV(400),DCH(400),TOX(400)
      REAL*8 STRDATA(6,400),DELT,XPRINT,TPRINT,RHO,G
      REAL*8 SIGMA,VOLCON
      REAL*8 RHOC,RHOM,RHOA,ALPHA,V2,RCH2,POS,DINIT,RTBIS,DY
      INTEGER*4 M,JN,NT,NODE
      CHARACTER*25 RUN
      CHARACTER ANS
C----------------------------------------------------------------------
C  COMMON BLOCKS
C----------------------------------------------------------------------
      COMMON /COM1/M,NT,DELX,STRDATA,WID,ELEV,WSEBC,MANN
      COMMON /COM2/Q,DELT,XPRINT,TPRINT,TOX,RHO,G,RUN
      COMMON /COM4/SIGMA,VOLCON,RHOC,RHOM,RHOA,ALPHA
C----------------------------------------------------------------------
C  DEFINE SOME BOUNDARY CONDITIONS
C----------------------------------------------------------------------
      IF (ELEV(M) .GT. WSEBC) THEN
         WRITE(*,*) 'Pausing...'
         READ(*,*) ANS
      ENDIF
      DINIT = WSEBC - ELEV(M)
C----------------------------------------------------------------------
C  CALCULATE INITIAL CROSS-SECT AREA, THE WETTED
C  PERIMETER(P), THE MEAN VEL(V), THE HYD RAD(RCH),
C  AND THE DEPTH OF THE WATER(DCH).
C----------------------------------------------------------------------
      A = DINIT*WID(M)
      P = WID(M) + 2.*DINIT
      V = Q(M)/A
      RCH = A/P
      DCH(M) = DINIT
      SF1 = ((MANN*V)/RCH**0.666667)**2.
      H1 = ELEV(M) + DINIT + V**2./(2.*9.81)
C----------------------------------------------------------------------
C  SWEEP THE GRID FROM DOWN STREAM TO UPSTREAM
C----------------------------------------------------------------------
      DO 3 JN = M, 2, -1
         NODE = JN - 1
C----------------------------------------------------------------------
C  GUESS A NEW STREAM DEPTH AT THE NEXT UPSTREAM NODE STARTING
C  WITH THE DEPTH AT THE PRESENT NODE. THE ONLY NEW DATA
C  USED IS THE STREAM WIDTH AT THE NEW NODE.
C----------------------------------------------------------------------
         DY = 1000.0
         D2 = 1000.0
         RTBIS = 0.0
    1    CONTINUE
         A2 = D2*WID(NODE)
         P2 = WID(NODE) + 2.*D2
         V2 = Q(NODE)/A2
         RCH2 = A2/P2
         SF2 = ((MANN*V2)/RCH2**0.666667)**2.
         H2 = ELEV(NODE) + D2 + V2**2./(2.*9.81)
         H21 = H2 - H1 - 0.5*DELX*(SF1+SF2)
         IF (H21 .LE. 0.0) RTBIS = D2
         IF (DABS(DY) .GT. 0.01) THEN
            DY = DY*0.5
            D2 = RTBIS + DY
            GO TO 1
         ELSE
    2       CONTINUE
            SF1 = SF2
            H1 = H2
            DCH(NODE) = D2
         ENDIF
         IF (V2**2. .GT. 9.81*D2) DCH(NODE) = ((Q(NODE)/WID(NODE))**2./
     #      9.81)**(1./3.)
    3 CONTINUE
C----------------------------------------------------------------------
C      WRITE TO THE ARRAY STRDATA
C----------------------------------------------------------------------
      DO 4 JN = 1, M
         A = DCH(JN)*WID(JN)
         P = WID(JN) + 2*DCH(JN)
         V = Q(JN)/A
         RCH = A/P
         SF1 = ((MANN*V)/RCH**0.666667)**2.
         POS = FLOAT(JN-1)*DELX
         FR = V/(G*DCH(JN))**.5
         STRDATA(1,JN) = RCH
         STRDATA(2,JN) = FR
         STRDATA(3,JN) = SF1
         STRDATA(4,JN) = POS
         STRDATA(5,JN) = V
         STRDATA(6,JN) = DCH(JN)
         TOX(JN) = RHO*G*RCH*SF1
    4 CONTINUE
      RETURN 
      END
 
 
 
C**************************SUBROUTINE WRFLDTA**************************
C**********************************************************************
C  THIS SUBROUTINE: WRITES INFORMATION ABOUT THE CHANNEL FORM AND
C                   PROFILE TO OUTPUT FILE.
C**********************************************************************
      SUBROUTINE WRFLDTA
      IMPLICIT NONE
      REAL*8    POS,Q(400),DELT,XPRINT,TPRINT
      REAL*8    VV,WID(400),ELEV(400),STRDATA(6,400),TOX(400)
      REAL*8    SF,DELX,RHO,G
      REAL*8    RCH,DCH,WSEBC,MANN
      REAL*8    FR
      REAL*8    POSSTEP
      INTEGER*4 JN
      INTEGER*4 NT,M
      CHARACTER*25 RUN
C----------------------------------------------------------------------
C     COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM1/M,NT,DELX,STRDATA,WID,ELEV,WSEBC,MANN
      COMMON /COM2/Q,DELT,XPRINT,TPRINT,TOX,RHO,G,RUN
C---------------------------------------------------------------------
C   CHECK TO SEE IF ITS IN POSITION TO PRINT
C---------------------------------------------------------------------
      DO 1 JN = 1, M
         POSSTEP = (JN-1)*DELX/XPRINT + 0.00001
         IF (POSSTEP - FLOAT(INT(POSSTEP)) .LE. 0.0001) THEN
            POS = STRDATA(4,JN)
            RCH = STRDATA(1,JN)
            SF = STRDATA(3,JN)
            VV = STRDATA(5,JN)
            FR = STRDATA(2,JN)
            DCH = STRDATA(6,JN)
            WRITE(8,100)POS,ELEV(JN),DCH,WID(JN),RCH,VV,FR,SF,TOX(JN)
         ENDIF
    1 CONTINUE
  100 FORMAT(9E11.5)
      RETURN 
      END
 
C***********************************************************************
C     SUBROUTINE ERODEP90
C     This program calculates the total sediment transport of a river
C     at a node, using Yang's (1973) unit stream power equation.
C     It then solves a conserv. of mass eqn. for the bed, thereby
C     obtaining bed elevations at the next time step.
C     Assume one size and density are present.
C***********************************************************************
C     VARIABLES
C     Y(#space nodes)
C     BLWT(#space nodes)
C     DBDX(#space nodes),XFLUX(#space nodes)
C     WIDM(#space nodes)
C     SIGMA,VOLCON
C     DELT,DELX
C***********************************************************************
      SUBROUTINE ERODEP90(MFLUX,FACTOR,SEDCON,MFEED)
      IMPLICIT NONE
      REAL*8 USTAR(400),TOX(400),RHO,DIMID,OMEGA,A,B,C,D,E,F,VCR,REY
      REAL*8 DLESSU,VISKIN,CONC,STRDATA(6,400),MFLUX(400),MFEED
      REAL*8 WID(400),DELT,SIGMA,VOLCON,DELX,ELEV(400),PHI,DMFLUXDX(400)
      REAL*8 XPRINT,TPRINT,G,Q(400),WSEBC,MANN,FACTOR
      REAL*8 RHOC,RHOM,RHOA,ALPHA,SEDCON,QIN
      INTEGER*4 JN,M,NT
      CHARACTER*25 RUN
C----------------------------------------------------------------------
C     COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM1/M,NT,DELX,STRDATA,WID,ELEV,WSEBC,MANN
      COMMON /COM2/Q,DELT,XPRINT,TPRINT,TOX,RHO,G,RUN
      COMMON /COM3/VISKIN,DIMID,OMEGA,REY
      COMMON /COM4/SIGMA,VOLCON,RHOC,RHOM,RHOA,ALPHA
C***********************************************************************
C     SET SOME CONSTANTS AND CALCULATE SOME VALUES
C***********************************************************************
      VOLCON = 0.7
      SIGMA = 2650.0
      PHI = 0.6
      IF (DIMID .LE. 0.004) THEN
         A = 5.435
         B = 0.286
         C = 0.457
         D = 1.799
         E = 0.409
         F = 0.314
      ELSE
         A = 6.681
         B = 0.633
         C = 4.816
         D = 2.784
         E = 0.305
         F = 0.282
      ENDIF
C***********************************************************************
C     YANG'S EQUATION
C***********************************************************************
      DO 1 JN = 1, M
C***********************************************************************
C     CALCULATE THE SHEAR VELOCITY OF THE FLOW AND CRITICAL VELOCITY
C***********************************************************************
         USTAR(JN) = SQRT(TOX(JN)/RHO)
         IF (USTAR(JN)*DIMID/VISKIN .LE. 70.0) THEN
            VCR = OMEGA*(2.5/(DLOG(USTAR(JN)*DIMID/VISKIN)-0.06D0))
         ELSE
            VCR = 2.05*OMEGA
         ENDIF
         IF (STRDATA(5,JN) - VCR .GT. 0.0D0) THEN
            DLESSU = USTAR(JN)/OMEGA
            CONC = A - B*DLOG10(REY) - C*DLOG10(DLESSU) + (D-E*DLOG10(
     #         REY)-F*DLOG10(DLESSU))*DLOG10((STRDATA(3,JN)/OMEGA)*(
     #         STRDATA(5,JN)-VCR))
            CONC = FACTOR*10**CONC
         ELSE
            CONC = 0.0D0
         ENDIF
         MFLUX(JN) = CONC*STRDATA(5,JN)*STRDATA(6,JN)*WID(JN)*RHO/
     #      1000000.0
    1 CONTINUE
C***********************************************************************
C     SOLVE THE CONS. OF MASS EQN. STARTING AT NODE 2,
C     ADD SEDIMENT BY LATERAL INFLOW AT A RATE PROPORTIONAL TO WATER
C     INFLOW AS GIVEN IN SNOW AND SLINGERLAND (1986)
C***********************************************************************
      DO 2 JN = 2, M - 1
         DMFLUXDX(JN) = ((1.0D0-PHI)*(MFLUX(JN+1)-MFLUX(JN))+PHI*(MFLUX(
     #      JN)-MFLUX(JN-1)))/DELX
C***********************************************************************
C      NOW CALCULATE THE ELEVATION CHANGE AT EACH INTERIOR NODE
C***********************************************************************
         IF (ELEV(JN) .GT. WSEBC) THEN
C***********************************************************************
C     CALCULATE THE AMOUNT OF BED ELEVATION CHANGE DUE TO LATERAL
C     INFLOW OF SEDIMENT. SEDIMENT INPUT IS A FUNCTION OF LATERAL INFLOW 
C     OF WATER IN UNITS OF M**3/S PER METER LENGTH OF STREAM.
C***********************************************************************
            QIN = SEDCON*(Q(JN)-Q(JN-1))*DELT/(WID(JN)*DELX)
         ELSE
            QIN = 0.0
         ENDIF
         ELEV(JN) = ELEV(JN) - DMFLUXDX(JN)*DELT/(WID(JN)*SIGMA*VOLCON)
     #       + QIN
    2 CONTINUE
C***********************************************************************
C    NOW DEAL WITH THE BOUNDARY NODES
C***********************************************************************
      DMFLUXDX(1) = ((1.0D0-PHI)*(MFLUX(2)-MFLUX(1))+PHI*(MFLUX(1)-MFEED
     #   ))/DELX
      IF (ELEV(1) .GT. WSEBC) THEN
         QIN = SEDCON*(Q(2)-Q(1))*DELT/(WID(1)*DELX)
      ELSE
         QIN = 0.0
      ENDIF
      ELEV(1) = ELEV(1) - DMFLUXDX(1)*DELT/(WID(1)*SIGMA*VOLCON) + QIN
      QIN = 0.0
      ELEV(M) = ELEV(M) - DELT*(MFLUX(M)-MFLUX(M-1))/(WID(M)*SIGMA*
     #   VOLCON*DELX) + QIN
      RETURN 
      END
 
C************************SUBROUTINE TECTON******************************
C***********************************************************************
C     This subroutine calculates the tectonic elevation changes
C     at each timestep, through a user-supplied function.
C***********************************************************************
      SUBROUTINE TECTON
      IMPLICIT NONE
      REAL*8 DELX,STRDATA(6,400),WID(400),ELEV(400),
     *       WSEBC,MANN,Q(400),DELT,XPRINT,TPRINT,TOX(400),RHO,G,
     *       SIGMA,VOLCON,RHOC,RHOM,RHOA,ALPHA,
     *       RATE
      INTEGER*4 M,NT,I,TURNON1,TURNOFF1,TURNON2,
     *          TURNOFF2,HINGE
      CHARACTER*25 RUN
C----------------------------------------------------------------------
C     COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM1/M,NT,DELX,STRDATA,WID,ELEV,WSEBC,MANN
      COMMON /COM2/Q,DELT,XPRINT,TPRINT,TOX,RHO,G,RUN
      COMMON /COM4/SIGMA,VOLCON,RHOC,RHOM,RHOA,ALPHA
      COMMON /TECTO/RATE,TURNON1,TURNOFF1,TURNON2,TURNOFF2,HINGE
C***********************************************************************
      IF (NT.GE.TURNON1 .AND. NT.LE.TURNOFF1 .OR. NT.GE.TURNON2 .AND. NT
     #   .LE.TURNOFF2) THEN
         DO I = HINGE, M
            ELEV(I) = ELEV(I) + RATE*(M-I+1)
         END DO
      ENDIF
      RETURN 
      END
 
C*********************SUBROUTINE SETTLE********************************
C**********************************************************************
C      CALCULATES THE CONSTANT TERMINAL SETTLING VELOCITY OF SIZE,
C      DIMID USING DIETRICH'S EQUATION
C**********************************************************************
C*********************************************************************
      SUBROUTINE SETTLE
C----------------------------------------------------------------------
C      IMPLICIT STATEMENT
C----------------------------------------------------------------------
       IMPLICIT NONE
      CHARACTER*25 RUN
      REAL*8 DELX,DELT,XPRINT,TPRINT,WID(400),ELEV(400)
      REAL*8 STRDATA(6,400),RHO,VISKIN,DIMID,G
      REAL*8 WSEBC,MANN,Q(400),TOX(400)
      INTEGER*4 NT,M
      REAL*8 DIMDL,RHS,OMEGA,REY
      COMMON /COM1/M,NT,DELX,STRDATA,WID,ELEV,WSEBC,MANN
      COMMON /COM2/Q,DELT,XPRINT,TPRINT,TOX,RHO,G,RUN
      COMMON /COM3/VISKIN,DIMID,OMEGA,REY
C----------------------------------------------------------------------
C      COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      DIMDL = (2650.0D0-RHO)*G*DIMID**3/(RHO*VISKIN*VISKIN)
      DIMDL = DLOG10(DIMDL)
      RHS = (-3.76715D0) + 1.92944D0*DIMDL - 0.09815D0*DIMDL**2
      RHS = RHS - 0.00575D0*DIMDL**3 + 0.00056D0*DIMDL**4
      OMEGA = (((2650.0D0-RHO)*G*VISKIN*10.0D0**RHS)/RHO)**0.333
      RETURN 
      END
