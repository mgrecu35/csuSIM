      module GPM_TbsimulatorV7_rad_procedures

      use fastem
      
      implicit none

      contains
!-----------------------------------------------------------------------
      
      subroutine absorb_clr(freqy, temp, pres, vapor_pres, kabs_clear)

!     Phil Rosenkranz's absorption code provided to C. kummerow in 2003.
!     Subroutine to provide calling sequence to water vapor, oxygen
!     and air absorption as computed by Phil Rosenkrantz.  The individual
!     absorption codes have not been formally published but rather
!     constitute ongoing work. 
!
!     INPUT
!       freqy      :  frequency [GHz]
!       temp       :  layer avg. temperature [K]
!       pres       :  layer avg. pressure [mb]
!       vapor_pres :  layer average water vapor pressure [mb]
!     OUTPUT
!       kabs       :  absorption coefficient [km^-1]
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real  kabs_H2O, kabs_O2, kabs_N2, kabs_clear
      real  freqy, rho, temp, pres, vapor_pres
      
!     Compute the vapor density, rho [g/m^3]      
      rho = vapor_pres*100*18/(8.314*temp)
      

      call abs_H2O(temp, pres, rho, freqy, kabs_H2O)      
      call abs_O2 (temp, pres, rho, freqy, kabs_O2) 
      call abs_N2 (temp, pres, freqy, kabs_N2)
      
      kabs_clear = kabs_H2O + kabs_O2 + kabs_N2
      
      return
      end subroutine absorb_clr
      
!-----------------------------------------------------------------------
            
      Subroutine abs_H2O(T,P,RHO,F,ABH2O)

!     C. Kummerow, 8/2003.  Changed function to subroutine     
!     Copyright (c) 2002 Massachusetts Institute of Technology

!     NAME- ABH2O    LANGUAGE- FORTRAN 77

!     PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
 
!     CALLING SEQUENCE PARAMETERS-
!     SPECIFICATIONS
      REAL T,P,RHO,F,ABH2O
!      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
!      T       KELVIN    I   TEMPERATURE
!      P       MILLIBAR  I   PRESSURE              .1 TO 1000
!      RHO     G/M**3    I   WATER VAPOR DENSITY
!      F       GHZ       I   FREQUENCY             0 TO 800
!      ABH2O   NEPERS/KM O   ABSORPTION COEFFICIENT

!   REFERENCES-
!   P.W. ROSENKRANZ, RADIO SCIENCE V.33, PP.919-928 (1998); V.34, P.1025 (1999).
!
!   LINE INTENSITIES SELECTION THRESHOLD=
!     HALF OF CONTINUUM ABSORPTION AT 1000 MB.
!   WIDTHS MEASURED AT 22,183,380 GHZ, OTHERS CALCULATED.
!     A.BAUER ET AL.ASA WORKSHOP (SEPT. 1989) (380GHz).
!     M. TRETYAKOV et al., J. MOLEC. SPEC. (2003)

!   REVISION HISTORY-
!    DATE- OCT.6, 1988  P.W.ROSENKRANZ - EQS AS PUBL. IN 1993.
!          OCT.4, 1995  PWR- USE CLOUGH'S DEFINITION OF LOCAL LINE
!                   CONTRIBUTION,  HITRAN INTENSITIES, ADD 7 LINES.
!          OCT. 24, 95  PWR -ADD 1 LINE.
!          JULY 7, 97   PWR -SEPARATE COEFF. FOR SELF-BROADENING, 
!                       REVISED CONTINUUM.
!        Aug. 28, 2002  PWR - CORRECTED LINE INTENSITIES
!        Mar. 2, 2003   PWR - LINE SHIFT

!   LOCAL VARIABLES:
      INTEGER NLINES,I,J
      PARAMETER (NLINES=15)
      REAL DF(2),S1(NLINES),B2(NLINES),W3(NLINES),FL(NLINES),X(NLINES),  &
                 WS(NLINES),XS(NLINES),SR(NLINES)                
      REAL PVAP,PDA,DEN,TI,TI2,SUM,WIDTH,WSQ,S,BASE,RES,CON,SHIFT
!     LINE FREQUENCIES:
      DATA FL/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, &
        443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360,      &
        620.7008, 752.0332, 916.1712/
!     LINE INTENSITIES AT 300K:
      DATA S1/ .1314E-13, .2279E-11, .8058E-13, .2701E-11, .2444E-10,    &
       .2185E-11, .4637E-12, .2568E-10, .8392E-12, .3272E-11, .6676E-12, &
       .1535E-08, .1711E-10, .1014E-08, .4238E-10/
!     T COEFF. OF INTENSITIES:
      DATA B2/ 2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405,    &
        3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/
!     AIR-BROADENED WIDTH PARAMETERS AT 300K:
      DATA W3/.00281, .00287, .0023, .00278, .00287, .0021, .00186,      &
        .00263, .00215, .00236, .0026, .00321, .00244, .00306, .00267/
!     T-EXPONENT OF AIR-BROADENING:
      DATA X/.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69, &
        .71, .68, .70/
!     SELF-BROADENED WIDTH PARAMETERS AT 300K:
      DATA WS/.01349, .01491, .0108, .0135, .01541, .0090, .00788,       &
        .01275, .00983, .01095, .01313, .01320, .01140, .01253, .01275/
!     T-EXPONENT OF SELF-BROADENING:
      DATA XS/ .61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72,    &
        1.0, .68, .84, .78/
!     RATIO OF SHIFT TO WIDTH
      DATA SR/ 0., -.017, 13*0./
!
      IF(RHO.LE.0.) THEN
        ABH2O = 0.
        RETURN
      ENDIF
      PVAP = RHO*T/217.
      PDA = P -PVAP
      DEN = 3.335E16*RHO ! const includes isotopic abundance
      TI = 300./T
      TI2 = TI**2.5

!      CONTINUUM TERMS
      CON = (5.43E-10*PDA*TI**3 + 1.8E-8*PVAP*TI**7.5)*PVAP*F*F 

!      ADD RESONANCES
      SUM = 0.
      DO 30 I=1,NLINES
      WIDTH = W3(I)*PDA*TI**X(I) + WS(I)*PVAP*TI**XS(I)
      SHIFT = SR(I)*WIDTH  ! unknown temperature dependence
      WSQ = WIDTH*WIDTH
      S = S1(I)*TI2*EXP(B2(I)*(1.-TI))
      DF(1) = F - FL(I) - SHIFT
      DF(2) = F + FL(I) + SHIFT
!  USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
      BASE = WIDTH/(562500. + WSQ)
!  DO FOR POSITIVE AND NEGATIVE RESONANCES
      RES = 0.
      DO 20 J=1,2
      IF(ABS(DF(J)).LT.750.) RES = RES + WIDTH/(DF(J)**2+WSQ) - BASE
20    CONTINUE
30    SUM = SUM + S*RES*(F/FL(I))**2
      ABH2O = .3183E-4*DEN*SUM + CON
      RETURN
      END Subroutine abs_H2O

!-----------------------------------------------------------------------
      
      subroutine abs_O2(TEMP,PRES,VAPDEN,FREQY,O2ABS)

!     C. Kummerow, 8/2003.  Changed function to subroutine      
!  Copyright (c) 2003 Massachusetts Institute of Technology

!     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
!              IN NEPERS/KM

!      5/1/95  P. Rosenkranz 
!      11/5/97  P. Rosenkranz - 1- line modification.
!      12/16/98 pwr - updated submm freq's and intensities from HITRAN96
!      8/21/02  pwr - revised width at 425
!      3/20/03  pwr - 1- line mixing and width revised

!     ARGUMENTS:
      REAL TEMP,PRES,VAPDEN,FREQY,O2ABS

!     NAME    UNITS    DESCRIPTION        VALID RANGE

!     TEMP    KELVIN   TEMPERATURE        UNCERTAIN, but believed to be
!                                          valid for atmosphere
!     PRES   MILLIBARS PRESSURE           3 TO 1000
!     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
!                       DUE TO GREATER BROADENING EFFICIENCY OF H2O)
!     FREQ    GHZ      FREQUENCY          0 TO 900

!     REFERENCES FOR EQUATIONS AND COEFFICIENTS:
!     P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
!      BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
!     H.J. Liebe et al, JQSRT V.48, pp.629-643 (1992).
!     M.J. Schwartz, Ph.D. thesis, M.I.T. (1998).
!     A.F. Krupnov et al, J. Mol. Spect. v.215, pp.309-311 (2002).
!     M.Yu. Tretyakov et al, J. Mol. Spect. (2003 preprint).
!     SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.

!     This version differs from Liebe's MPM92 in these significant respects:
!     1. The 1- line has the width and mixing coefficient measured by 
!      Tretyakov et al. 
!     2. It modifies the 1- line width temperature dependence to (1/T)**0.9
!     3. It uses the same temperature dependence (X) for submillimeter 
!      line widths as in the 60 GHz band: (1/T)**0.8 
!     4. The 425 GHz line width is from Krupnov et al.

!     Local variables:
      REAL TH,TH1,B,PRESWV,PRESDA,DEN,DENS,DFNR,SUM,STR,Y,SF1,SF2,FCEN
      real DF
      INTEGER K
      REAL X,WB300,W300(40),F(40),Y300(40),S300(40),V(40),BE(40)
      COMMON /O2COM/ X,WB300,W300,F,Y300,S300,V,BE
!      LINES ARE ARRANGED 1-,1+,3-,3+,ETC. IN SPIN-ROTATION SPECTRUM
      DATA F/118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,   &
       59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,    &
       56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,    &
       55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,    &
       53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,    &
       52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7632,  &
       487.2494, 715.3931, 773.8397, 834.1458/
        DATA S300/.2936E-14,.8079E-15, .2480E-14,.2228E-14,  &
        .3351E-14,.3292E-14, .3721E-14,.3891E-14,   &
        .3640E-14,.4005E-14, .3227E-14,.3715E-14,   &
        .2627E-14,.3156E-14, .1982E-14,.2477E-14,   &
        .1391E-14,.1808E-14, .9124E-15,.1230E-14,   &
        .5603E-15,.7842E-15, .3228E-15,.4689E-15,   &
        .1748E-15,.2632E-15, .8898E-16,.1389E-15,   &
        .4264E-16,.6899E-16, .1924E-16,.3229E-16,   &
        .8191E-17,.1423E-16, .6494E-15, .7083E-14, .3025E-14,   &
        .1835E-14, .1158E-13, .3993E-14/        
      DATA BE/.009,.015, .083,.084, 2*.212, 2*.391, 2*.626,     &
       2*.915, 2*1.260, 1.660,1.665, 2.119,2.115, 2.624,2.625,  &
       2*3.194, 2*3.814, 2*4.484, 2*5.224, 2*6.004, 2*6.844,    &
       2*7.744, .048, .044, .049, .145, .141, .145/
!      WIDTHS IN MHZ/MB
      DATA WB300/.56/, X/.8/
      DATA W300/1.67, 1.646, 1.468, 1.449, 1.382, 1.360,        &
       1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171,  &
       1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 2*1.05,        &
       2*1.02,2*1.00,2*.97,2*.94,2*.92,2*.89, 3*1.64, 3*1.81/ 
      DATA Y300/  -0.036,  0.2408, -0.3486,  0.5227,            &
       -0.5430,  0.5877, -0.3970,  0.3237, -0.1348,  0.0311,    &
        0.0725, -0.1663,  0.2832, -0.3629,  0.3970, -0.4599,    &
        0.4695, -0.5199,  0.5187, -0.5597,  0.5903, -0.6246,    &
        0.6656, -0.6942,  0.7086, -0.7325,  0.7348, -0.7546,    &
        0.7702, -0.7864,  0.8083, -0.8210,  0.8439, -0.8529, 6*0./
      DATA V/  0.0079, -0.0978,  0.0844, -0.1273,               &
        0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,    &
        0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,    &
        0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,    &
        0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,    &
        0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545, 6*0./

      TH = 300./TEMP
      TH1 = TH-1.
      B = TH**X
      PRESWV = VAPDEN*TEMP/217.
      PRESDA = PRES -PRESWV
      DEN = .001*(PRESDA*B + 1.1*PRESWV*TH)
      DENS = .001*(PRESDA*TH**.9 + 1.1*PRESWV*TH)
      DFNR = WB300*DEN
      SUM = 1.6E-17*FREQY*FREQY*DFNR/(TH*(FREQY*FREQY + DFNR*DFNR))
      DO 32 K=1,40
      IF(K.EQ.1) THEN !exception for 1- line
        DF = W300(1)*DENS
      ELSE
        DF = W300(K)*DEN
      ENDIF
      FCEN = F(K)
      Y = .001*PRES*B*(Y300(K)+V(K)*TH1)
      STR = S300(K)*EXP(-BE(K)*TH1)
      SF1 = (DF + (FREQY-FCEN)*Y)/((FREQY-FCEN)**2 + DF*DF)
      SF2 = (DF - (FREQY+FCEN)*Y)/((FREQY+FCEN)**2 + DF*DF)
32    SUM = SUM + STR*(SF1+SF2)*(FREQY/F(K))**2
      O2ABS = .5034E12*SUM*PRESDA*TH**3/3.14159
      O2ABS = AMAX1(O2ABS,0.)
      RETURN
      END subroutine abs_O2

!-----------------------------------------------------------------------

      Subroutine ABS_N2(T,P,F,ABSN2)

!     C. Kummerow, 8/2003.  Changed function to subroutine      
!  Copyright (c) 2002 Massachusetts Institute of Technology
!     ABSN2 = COLLISION-INDUCED ABSORPTION COEFFICIENT (NEPER/KM)
!     IN AIR
!     T = TEMPERATURE (K)
!     P = PRESSURE (MB)
!     F = FREQUENCY (GHZ)(valid 0-1000 GHz)

!     5/22/02 P.Rosenkranz

!     Equations based on:
!      Borysow, A, and L. Frommhold, 
!      Astrophysical Journal, v.311, pp.1043-1057 (1986)
!     with modification of 1.29 to account for O2-O2 and O2-N2
!     collisions, as suggested by
!      J.R. Pardo, E.Serabyn, J.Cernicharo, J. Quant. Spectros.
!      Radiat. Trans. v.68, pp.419-433 (2001).

       real :: T,TH,FDEPEN,BF,P,F,ABSN2

      TH = 300./T
      FDEPEN = .5 + .5/(1.+(F/450.)**2)
      BF = 6.5E-14*FDEPEN*P*P*F*F*TH**3.6
      ABSN2 = 1.29*BF
      RETURN
      END Subroutine ABS_N2

!-----------------------------------------------------------------------
   
      subroutine absorb_clw(freq, temp, water, abliq)
      
!     COMPUTES ABSORPTION IN NEPERS/KM BY SUSPENDED WATER DROPLETS
!     ARGUMENTS (INPUT):
!     WATER IN G/M**3
!     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
!     TEMP IN KELVIN

!     REFERENCES:
!     LIEBE, HUFFORD AND MANABE, INT. J. IR & MM WAVES V.12, pp.659-675
!      (1991);  Liebe et al, AGARD Conf. Proc. 542, May 1993.

!     REVISION HISTORY:
!        PWR 8/3/92   original version
!        PWR 12/14/98 temp. dependence of EPS2 eliminated to agree 
!                     with MPM93 
!        pwr 2/27/02  use exponential dep. on T, eq. 2b instead of eq. 4a 

      COMPLEX :: EPS,RE
      real    :: freq,temp,water,abliq,theta1,eps0,eps1,eps2, FP, FS

      IF(WATER.LE.0.) THEN
       ABLIQ = 0.
       RETURN
      ENDIF
      THETA1 = 1.-300./TEMP
      EPS0 = 77.66 - 103.3*THETA1
      EPS1 = .0671*EPS0
      EPS2 = 3.52                 ! from MPM93
      FP = 20.1*EXP(7.88*THETA1)  ! from eq. 2b
      FS = 39.8*FP
      EPS = (EPS0-EPS1)/CMPLX(1.,FREQ/FP) + (EPS1-EPS2)/CMPLX(1.,FREQ/FS) +EPS2
      RE = (EPS-1.)/(EPS+2.)
      abliq = -.06286*AIMAG(RE)*FREQ*WATER
      RETURN
      END subroutine absorb_clw

!-----------------------------------------------------------------------

      !DSD from LWC and D0
       subroutine calc_DSD_W_D0(D0, W, DSD)

       real :: D0, W, DSD(3), N0, mu, w1

       mu = 3.0
       call calc_lwcD(1.0, D0, mu, w1)
       N0 = W/w1
       DSD(1) = N0
       DSD(2) = D0
       DSD(3) = mu

       return
       end subroutine calc_DSD_W_D0
       
!-----------------------------------------------------------------------       

       !LWC from DSD
       subroutine calc_lwcD(N0, D0, mu, rw)

       real	D0, rw, dD, N0, lam, pi, rho_w, d, num, mu
       integer 	i

       rho_w = 1.0e3	!kg/m^3
       pi = acos(-1.)
       dD = 0.1		!mm

       if(N0 .gt. 0.) then
         lam = (3.67+mu)/D0

         !calculate RW
         rw = 0.
         do i=1,60
           d = dD*float(i)
           num = N0*(d**mu)*exp(-lam*d)
	  rw = rw+1/6.*pi*d**3*1e-3*num*dD
         end do
	
       else
         rw = 0.
       endif

       return
       end subroutine calc_lwcD
       
!-----------------------------------------------------------------------       

      !Rain Rate from DSD
       subroutine calc_RRD(N0, D0, mu, rr, tervels)

       real	D0, rr, dD, N0, lam, pi, rho_w, d, num, mu
       real	dummy1, dummy2, dummy3
       real	tervels(60), tv
       integer	i

       rho_w = 1.0e3	!kg/m^3
       pi = acos(-1.)


       lam = (3.67+mu)/D0

       rr = 0.
       dD = 0.10	!mm
       do i=1,60
         d = dD*float(i)
         num = N0*(d**mu)*exp(-lam*d)
	rr = rr+600*pi*d**3*10.*tervels(i)*1e-9*num*dD
       end do
       !print*, N0, D0, mu, rr

       RETURN
       end subroutine calc_RRD
       
!-----------------------------------------------------------------------       

      SUBROUTINE  EMISSIVITY( F, NPOL, TS, W, ANGLE, EMIS, EBAR )

      integer,parameter :: NANG = 21
      REAL ::  ANG(NANG), MMU(NANG), ESUM(NANG), PI
      real :: F,Ts,W,ANGLE,EMIS,EBAR,S,Tc,Er,Ei,F_angle
      real :: EV,EH,eavg,dmu,avmu,sum, angles
      integer :: npol,i

      PI = 2.*ASIN(1.0)
      S = 35.                   ! SALINITY IN PPM
      Tc = Ts - 273.16                                                   
      ! ANGLE = ACOS(UMU)*180./PI

!     Get DIELECTRIC PROPERTIES of SEA WATER
      CALL eps_sea_english(F,Tc,Er,Ei)
      
!     CALCULATE EMIS AT GIVEN ANGLE
      CALL EMISS (F,ANGLE,TS,W,Er,Ei,EV,EH)
      IF ( NPOL .EQ. 0 ) EMIS = EH
      IF ( NPOL .EQ. 1 ) EMIS = EV

!     CALCULATE EMIS AT VARIOUS ANGLES
      DO 58  I = 1,NANG
       ANG(I) = 4.*( I - 1 )
       MMU(I) = COS(ANG(I)*PI/180.)
       ANGLES = ANG(I)
       CALL EMISS (F,ANGLES,TS,W,Er,Ei,EV,EH)
       ESUM(I) =  EV + EH 
  58  CONTINUE

!     CALCULATE EBAR
      SUM = 0.0
      DO 59  I = 1,NANG-1
       EAVG = 0.5*( ESUM(I) + ESUM(I+1) )
       DMU = MMU(I) - MMU(I+1)
       AVMU = 0.5*( MMU(I) + MMU(I+1) )
       SUM = SUM + EAVG*AVMU*DMU
  59  CONTINUE
      EBAR = SUM
      RETURN
      END SUBROUTINE  EMISSIVITY

!---------------------------------------------------------------------------------       

      SUBROUTINE  EMISSIVITY_FEM6( F, TS, Sal, W, ANGLE, azimuth, EMIS, EBAR )

      integer,parameter :: NANG = 21
      REAL ::  ANG(NANG), MMU(NANG), ESUM(NANG), Sal, PI
      real :: F,Ts,W,ANGLE,EMIS(2),EBAR,S,E(4)
      real :: EV,EH,eavg,dmu,avmu,sum, angles, azimuth
      integer :: npol,i
      real :: Er, Ei, Tc

      PI = 2.*ASIN(1.0)
                                          
!     ANGLE = ACOS(UMU)*180./PI
!     S = 35.                   ! SALINITY IN PPM
!     Tc = Ts - 273.16                                                   
!     ANGLE = ACOS(UMU)*180./PI
!   Get DIELECTRIC PROPERTIES of SEA WATER
!      CALL eps_sea_english(F,Tc,Er,Ei)      
!   CALCULATE EMIS AT GIVEN ANGLE
!      CALL EMISS (F,ANGLE,TS,W,Er,Ei,EV,EH)
!      IF ( NPOL .EQ. 0 ) EMIS = EH
!      IF ( NPOL .EQ. 1 ) EMIS = EV

      CALL compute_fastem(6,F,ANGLE,TS,Sal,W,azimuth,E)
      EMIS(1) = E(1)    !V
      EMIS(2) = E(2)    !H
!      IF ( NPOL .EQ. 0 ) EMIS(1) = E(1)    !V
!      IF ( NPOL .EQ. 1 ) EMIS(2) = E(2)    !H


!     CALCULATE EMIS AT VARIOUS ANGLES
      DO 58  I = 1,NANG
       ANG(I) = 4.*( I - 1 )
       MMU(I) = COS(ANG(I)*PI/180.)
!       write(6,*) i,ang(i),mmu(i),pi
       
       ANGLES = ANG(I)
       CALL compute_fastem(6,F,Angles,Ts,Sal,W,azimuth,E)
       ESUM(I) =  E(1) + E(2)
!       CALL EMISS (F,ANGLES,TS,W,Er,Ei,EV,EH)
!       ESUM(I) =  EV + EH 
  58  CONTINUE

!     CALCULATE EBAR
      SUM = 0.0

      DO 59  I = 1,NANG-1
       EAVG = 0.5*( ESUM(I) + ESUM(I+1) )
       DMU = MMU(I) - MMU(I+1)
       AVMU = 0.5*( MMU(I) + MMU(I+1) )
       SUM = SUM + EAVG*AVMU*DMU
       
!       write(6,'(i2,5F10.4)') i,esum(i),eavg,dmu,avmu,sum
       
  59  CONTINUE
      EBAR = SUM
      RETURN
      END SUBROUTINE  EMISSIVITY_FEM6
      
!-----------------------------------------------------------------------      

      subroutine eps_sea_english(f,T,perm_real,perm_imag)
!     This subroutine returns the complex dielectric constant of sea water,
!     using the model of English & Deblonde (2003) from FASTEM2/RTTOVS-7; 
!     added on 06/01/05 by G. Elsaesser.  Comparison with other permittivity
!     G. Elsaesser's MS Thesis (2006).
!------------------------------------------------------------ ------------------

!     Inputs
      REAL f                    ! [GHz]    Frequency (valid from 0 to 1000 GHz)
      REAL T                    ! [°C]     Temperature

!     Local variables
      REAL T2, T3
      REAL tau1, tau2
      REAL del1, del2
      REAL den1, den2
      REAL einf, estar
      REAL fen, fen2
      REAL sig
      REAL perm_real, perm_real1, perm_real2
      REAL perm_imag, perm_imag1, perm_imag2, perm_imag3

!     Define quadratic and cubic functions for later polynomials
      T2 = T * T
      T3 = T2 * T

!     Define relaxation frequencies
      tau1 = 17.535  - 0.617670*T  + 0.008948*T2
      tau2 = 3.18420 + 0.0191890*T - 0.0108730*T2 + 0.000258180*T3

!     estatic = del1 + del2 + einf
      del1 = 68.3960 - 0.40642*T   + 0.0228320*T2 - 0.000530610*T3
      del2 = 4.76290 + 0.154100*T  - 0.0337170*T2 + 0.000844280*T3
      einf = 5.31250 - 0.0114770*T
      sig  = 2.906   + 0.09437*T
      estar = 8.854             ! Permittivity of free space in these units

!     Calculate eps using double-debye formula
      fen = 2.0 * 3.1415927 * f * 0.001
      fen2 = fen*fen
      den1 = 1.0 + fen2*tau1*tau1
      den2 = 1.0 + fen2*tau2*tau2
      perm_real1 = del1 / den1
      perm_real2 = del2 / den2
      perm_imag1 = del1*fen*tau1 / den1
      perm_imag2 = del2*fen*tau2 / den2
      perm_imag3 = sig / (fen*estar)
      perm_real = perm_real1 + perm_real2 + einf
      perm_imag = perm_imag1 + perm_imag2 + perm_imag3

      RETURN
      end subroutine eps_sea_english
      
!-----------------------------------------------------------------------
      
      SUBROUTINE  EMISS(F,ANGLE,TS,W,Er,Ei,EV,EH)
      
      COMPLEX :: EPS,ETAV,CUNITY,CTERM1V,CTERM1H,CTERM2,CTERMV,CTERMH
      REAL ::  Er,Ei, RV, theta, EMISV, EMISH,F,Angle,Ts,W,EH,Y, DTR,T, EV
      real ::  RSH,RSV, FOAM,GH,GV,A1, RFV, RFH, SQRTF,corrV,corrH,RRV, RRH, RH    
      DATA DTR / 0.01745329252 /
      T = TS-273.16
      THETA = ANGLE*DTR
      EPS = cmplx(Er,Ei)
      CUNITY = cmplx(1.0,0.0)

      !Use Fresnel Equations below (G. Elsaesser MS Thesis (2006) )
 
      CTERM1V = EPS * COS(THETA)
      CTERM1H = CUNITY * COS(THETA)
      CTERM2 = CSQRT(EPS - SIN(THETA)*SIN(THETA))
      CTERMV = (CTERM1V - CTERM2)/(CTERM1V + CTERM2)
      CTERMH = (CTERM1H - CTERM2)/(CTERM1H + CTERM2)
      EMISV = 1.0 - CABS(CTERMV)*CABS(CTERMV)
      EMISH = 1.0 - CABS(CTERMH)*CABS(CTERMH)      
      
!     FROM HERE THE EFFECT OF SURFACE ROUGHNESS AND FOAM ARE INCLUDED
!     BASED ON HOLLINGER MODEL AND FOAM MODEL OF STOGRYN

      RSH = 1.-EMISH
      RSV = 1.-EMISV

      FOAM = 7.751E-06 * W ** 3.231

      GH = 1.-1.748E-3*ANGLE-7.336E-5*ANGLE*ANGLE+1.044E-7*ANGLE*ANGLE*ANGLE
      GV = 1.-9.946E-4*ANGLE+3.218E-5*ANGLE*ANGLE-1.187E-6*ANGLE*ANGLE*ANGLE  &
          +7.0E-20*ANGLE*ANGLE*ANGLE*ANGLE*ANGLE*ANGLE*ANGLE*ANGLE*ANGLE*ANGLE

      A1 = (208.0+1.29*F)/TS

      RFV = 1. - A1 * GV
      RFH = 1. - A1 * GH

      Y = 7.32E-02*ANGLE

!     TS SURFACE TEMP IS IN DEGREE KELVIN

      SQRTF = SQRT(F)

      CORRV = (W*(1.17E-01-2.09E-03*EXP(Y))*SQRTF/TS)
      CORRH = (W*(1.15E-01+3.80E-05*ANGLE*ANGLE)*SQRTF/TS)

      RRV = RSV-CORRV
      RRH = RSH-CORRH

      RV = RRV*(1.-FOAM)+RFV*FOAM
      RH = RRH*(1.-FOAM)+RFH*FOAM
     
      EH = 1.-RH
      EV = 1.-RV

      RETURN
      END SUBROUTINE  EMISS
      
!-----------------------------------------------------------------------      
      
      SUBROUTINE LINPAK(NLYR, A, B, RCOND)

      integer,parameter :: mxlyr = 88
      
      REAL     :: A(2*MXLYR, 2*MXLYR), B(2*MXLYR), Z(2*MXLYR)
      INTEGER  :: PLDA, N, IPVT(2*MXLYR), NLYR, LDA
      REAL     :: RCOND
      
      LDA = 2*MXLYR
      N   = 2*NLYR

      CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
      CALL SGESL(A,LDA,N,IPVT,B,0)
      RETURN
      END SUBROUTINE LINPAK
      
!-----------------------------------------------------------------------      
            
      SUBROUTINE SGECO(A,LDA,N,IPVT,RCOND,Z)

!***BEGIN PROLOGUE  SGECO
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  CONDITION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION AND ESTIMATES
!            THE CONDITION NUMBER OF THE MATRIX.
!***DESCRIPTION

!     SGECO FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
!     AND ESTIMATES THE CONDITION OF THE MATRIX.

!     IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW SGECO BY SGESL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW SGECO BY SGEDI.
!     TO COMPUTE  INVERSE(A) , FOLLOW SGECO BY SGEDI.

!     ON ENTRY

!        A       REAL(LDA, N)
!                THE MATRIX TO BE FACTORED.

!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .

!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .

!     ON RETURN

!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.

!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.

!        RCOND   REAL
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.

!        Z       REAL(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.

!     SUBROUTINES AND FUNCTIONS

!     LINPACK SGEFA
!     BLAS SAXPY,SDOT,SSCAL,SASUM
!     FORTRAN ABS,AMAX1,SIGN
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  SASUM,SAXPY,SDOT,SGEFA,SSCAL
!***END PROLOGUE  SGECO

      INTEGER :: LDA,N,IPVT(1)
      REAL    :: A(LDA,1)
      real    :: Z(1)
      REAL    :: RCOND

      REAL EK,T,WK,WKM
      REAL ANORM,S,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L,m
      
!     COMPUTE 1-NORM OF A

!***FIRST EXECUTABLE STATEMENT  SGECO

      ANORM = 0.0E0
      DO 10 J = 1, N
         ANORM = AMAX1(ANORM,SASUM(N,A(1,J),1))
   10 CONTINUE

!     FACTOR

      CALL SGEFA(A,LDA,N,IPVT,INFO)

!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.

!     SOLVE TRANS(U)*W = E

      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK,-Z(K))
         IF (ABS(EK-Z(K)) .LE. ABS(A(K,K))) GO TO 30
            S = ABS(A(K,K))/ABS(EK-Z(K))
            CALL SSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (A(K,K) .EQ. 0.0E0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0E0
            WKM = 1.0E0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + ABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + ABS(Z(J))
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
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)

!     SOLVE TRANS(L)*Y = W

      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + SDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0E0) GO TO 110
            S = 1.0E0/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)

      YNORM = 1.0E0

!     SOLVE L*V = Y

      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL SAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0E0) GO TO 130
            S = 1.0E0/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM

!     SOLVE  U*Z = V

      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .LE. ABS(A(K,K))) GO TO 150
            S = ABS(A(K,K))/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
         T = -Z(K)
         CALL SAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
!     MAKE ZNORM = 1.0
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM

      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END SUBROUTINE SGECO

!-----------------------------------------------------------------------
      
      SUBROUTINE SGESL(A,LDA,N,IPVT,B,JOB)
!***BEGIN PROLOGUE  SGESL
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  SOLVES THE REAL SYSTEM A*X=B OR TRANS(A)*X=B
!            USING THE FACTORS OF SGECO OR SGEFA
!***DESCRIPTION

!     SGESL SOLVES THE REAL SYSTEM
!     A * X = B  OR  TRANS(A) * X = B
!     USING THE FACTORS COMPUTED BY SGECO OR SGEFA.

!     ON ENTRY

!        A       REAL(LDA, N)
!                THE OUTPUT FROM SGECO OR SGEFA.

!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .

!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .

!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM SGECO OR SGEFA.

!        B       REAL(N)
!                THE RIGHT HAND SIDE VECTOR.

!        JOB     INTEGER
!                = 0         TO SOLVE  A*X = B ,
!                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
!                            TRANS(A)  IS THE TRANSPOSE.

!     ON RETURN

!        B       THE SOLUTION VECTOR  X .

!     ERROR CONDITION

!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
!        OR SGEFA HAS SET INFO .EQ. 0 .

!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND IS TOO SMALL) GO TO ...
!           DO 10 J = 1, P
!              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE

!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.

!     SUBROUTINES AND FUNCTIONS
!
!     BLAS SAXPY,SDOT
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  SAXPY,SDOT
!***END PROLOGUE  SGESL

      INTEGER LDA,N,IPVT(1),JOB
      REAL A(LDA,1),B(1)
!
      REAL T
      INTEGER K,KB,L,NM1
      
!***FIRST EXECUTABLE STATEMENT  SGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50

!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE  L*Y = B
!
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL SAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE

!        NOW SOLVE  U*X = Y

         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL SAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE

!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B

         DO 60 K = 1, N
            T = SDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE

!        NOW SOLVE TRANS(L)*X = Y

         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + SDOT(N-K,A(K+1,K),1,B(K+1),1)
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
      END SUBROUTINE SGESL
      
!-----------------------------------------------------------------------
      
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)

!***PURPOSE  S.P. COMPUTATION Y = A*X + Y

!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       SA  SINGLE PRECISION SCALAR MULTIPLIER
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
!       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY

!     --OUTPUT--
!       SY  SINGLE PRECISION RESULT (UNCHANGED IF N .LE. 0)

!     OVERWRITE SINGLE PRECISION SY WITH SINGLE PRECISION SA*SX +SY.
!     FOR I = 0 TO N-1, REPLACE  SY(LY+I*INCY) WITH SA*SX(LX+I*INCX) +
!       SY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N
!       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.

      REAL    :: SX(1),SY(1),SA
      integer :: n,incx,incy,ix,iy,m,mp1,ns,i
      
!-------------------------------------------------------------      

      IF(N.LE.0.OR.SA.EQ.0.E0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE

!        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.

      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1

!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.

   20 M = IMOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN

!        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.

   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SA*SX(I) + SY(I)
   70     CONTINUE
      RETURN
      END SUBROUTINE SAXPY
      
!-----------------------------------------------------------------------      
      
      SUBROUTINE SSCAL(N,SA,SX,INCX)

!***PURPOSE  S.P. VECTOR SCALE X = A*X

!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       SA  SINGLE PRECISION SCALE FACTOR
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX

!     --OUTPUT--
!       SX  SINGLE PRECISION RESULT (UNCHANGED IF N .LE. 0)

!     REPLACE SINGLE PRECISION SX BY SINGLE PRECISION SA*SX.
!     FOR I = 0 TO N-1, REPLACE SX(1+I*INCX) WITH  SA * SX(1+I*INCX)
!
      REAL    :: SA,SX(1)
      Integer :: N, INCX,NS,M,MP1,i

      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20

!        CODE FOR INCREMENTS NOT EQUAL TO 1.

      NS = N*INCX
          DO 10 I = 1,NS,INCX
          SX(I) = SA*SX(I)
   10     CONTINUE
      RETURN

!        CODE FOR INCREMENTS EQUAL TO 1.

!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.

   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I + 1) = SA*SX(I + 1)
        SX(I + 2) = SA*SX(I + 2)
        SX(I + 3) = SA*SX(I + 3)
        SX(I + 4) = SA*SX(I + 4)
   50 CONTINUE
      RETURN
      END SUBROUTINE SSCAL
      
!-----------------------------------------------------------------------
      
      SUBROUTINE SGEFA(A,LDA,N,IPVT,INFO)
!***BEGIN PROLOGUE  SGEFA
!***DATE WRITTEN   780814   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  D2A1
!***KEYWORDS  FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX
!***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
!***PURPOSE  FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
!***DESCRIPTION

!     SGEFA FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.

!     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .

!     ON ENTRY

!        A       REAL(LDA, N)
!                THE MATRIX TO BE FACTORED.

!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .

!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .

!     ON RETURN

!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.

!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.

!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
!                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.

!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.

!     SUBROUTINES AND FUNCTIONS

!     BLAS SAXPY,SSCAL,ISAMAX
!***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
!                 *LINPACK USERS  GUIDE*, SIAM, 1979.
!***ROUTINES CALLED  ISAMAX,SAXPY,SSCAL
!***END PROLOGUE  SGEFA

      INTEGER :: LDA,N,IPVT(1),INFO
      REAL    :: A(LDA,N)

      REAL    :: T
      INTEGER :: J,K,KP1,L,NM1

!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING

!***FIRST EXECUTABLE STATEMENT  SGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1

!        FIND L = PIVOT INDEX

         L = ISAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L

!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED

         IF (A(L,K) .EQ. 0.0E0) GO TO 40

!           INTERCHANGE IF NECESSARY

            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE

!           COMPUTE MULTIPLIERS

            T = -1.0E0/A(K,K)
            CALL SSCAL(N-K,T,A(K+1,K),1)

!           ROW ELIMINATION WITH COLUMN INDEXING

            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL SAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0E0) INFO = N
      RETURN
      END SUBROUTINE SGEFA
      
!----------------------------------------------------------------------- 
      
      FUNCTION SASUM(N,SX,INCX)

!***PURPOSE  SUM OF MAGNITUDES OF S.P VECTOR COMPONENTS

!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX

!     --OUTPUT--
!    SASUM  SINGLE PRECISION RESULT (ZERO IF N .LE. 0)

!     RETURNS SUM OF MAGNITUDES OF SINGLE PRECISION SX.
!     SASUM = SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))

      integer,intent(in) :: n
      real,intent(in)    :: sx(1)
      integer,intent(in) :: incx
      real :: sasum

      integer :: NS,M,MP1,I

      SASUM = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20

!        CODE FOR INCREMENTS NOT EQUAL TO 1.

      NS = N*INCX
          DO 10 I=1,NS,INCX
          SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
      RETURN

!        CODE FOR INCREMENTS EQUAL TO 1.

!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.

   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SASUM = SASUM + ABS(SX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        SASUM = SASUM + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2))   &
       + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))
   50 CONTINUE
      RETURN
      END FUNCTION SASUM
            
!-----------------------------------------------------------------------      
      
      FUNCTION SDOT(N,SX,INCX,SY,INCY)

!***PURPOSE  S.P. INNER PRODUCT OF S.P. VECTORS

!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
!       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY

!     --OUTPUT--
!     SDOT  SINGLE PRECISION DOT PRODUCT (ZERO IF N .LE. 0)

!     RETURNS THE DOT PRODUCT OF SINGLE PRECISION SX AND SY.
!     SDOT = SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
!     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
!     DEFINED IN A SIMILAR WAY USING INCY.

      integer :: ix,iy,m,mp1,ns,I

      integer, intent(in) :: N
      real,    intent(in) :: SX(1)
      integer, intent(in) :: INCX
      real,    intent(in) :: SY(1)
      integer, intent(in) :: INCY
      real :: SDOT

      SDOT = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1)5,20,60
    5 CONTINUE

!        CODE FOR UNEQUAL INCREMENTS OR NONPOSITIVE INCREMENTS.

      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SDOT = SDOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1

!        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.

   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SDOT = SDOT + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SDOT = SDOT + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +           &
        SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE
      RETURN

!        CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.

   60 CONTINUE
      NS=N*INCX
      DO 70 I=1,NS,INCX
        SDOT = SDOT + SX(I)*SY(I)
   70   CONTINUE
      RETURN
      END FUNCTION SDOT
      
!-----------------------------------------------------------------------      
      
      FUNCTION ISAMAX(N,SX,INCX)

!***PURPOSE  FIND LARGEST COMPONENT OF S.P. VECTOR

!     --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
!       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
!     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX

!     --OUTPUT--
!   ISAMAX  SMALLEST INDEX (ZERO IF N .LE. 0)

!     FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF SINGLE PRECISION SX.
!     ISAMAX =  FIRST I, I = 1 TO N, TO MINIMIZE  ABS(SX(1-INCX+I*INCX)

      real :: SMAX,XMAG
      integer :: NS,II,I

      integer, intent(in) :: N
      real,intent(in)     :: SX(1)
      integer,intent(in)  :: INCX
      integer :: isamax

      ISAMAX = 0
      IF(N.LE.0) RETURN
      ISAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20

!        CODE FOR INCREMENTS NOT EQUAL TO 1.

      SMAX = ABS(SX(1))
      NS = N*INCX
      II = 1
          DO 10 I=1,NS,INCX
          XMAG = ABS(SX(I))
          IF(XMAG.LE.SMAX) GO TO 5
          ISAMAX = II
          SMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN

!        CODE FOR INCREMENTS EQUAL TO 1.

   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N
         XMAG = ABS(SX(I))
         IF(XMAG.LE.SMAX) GO TO 30
         ISAMAX = I
         SMAX = XMAG
   30 CONTINUE
      RETURN
      END  FUNCTION ISAMAX
            
!----------------------------------------------------------------------- 
     

      end module GPM_TbsimulatorV7_rad_procedures
