module GPM_TbsimulatorV7_mie_code
 
 use GPM_TbsimulatorV7_variables 
! use GPM_TbsimulatorV7_procedures     
 
 implicit none
       
 contains

!-------------------------------------------------------------------------------------------

 subroutine mie_ci_lut(freqy, temp, lwc, ksca, asca, gsca)
      

! Compute the extinction, absorption, asymmetry parameter and
! backscatter for a given water content of cloud ice in [g/m^3].

! Input:
! freqy	     frequency of radiation [GHz]
! temp 	     temperature of particles [K]
! lwc  	     water content of cloud ice distribution [g/m**3]
! reff_ciw	     effective radius of particles [mm]

! Output:
! ksca 	     extinction coefficient [1/km]
! asca 	     single-scatter albedo []
! gsca 	     asymmetry factor []
! pbck 	     backscatter phase function/(4*pi) []
!

  real :: freqy, temp, lwc, reff_ciw
  real :: ksca, asca, gsca, pbck

  real :: wavel, x
  real :: density_ice, density_cloudice
  real :: rad, dropmass, num
  real :: qext, qsca, asym ,qbsca
  real :: bext, bsca, bsym, bq11
 
  real :: epsreal, epsimag
  real :: fincl      
  complex :: eice, eair, emg, cref

!  data pi/3.14159265/
  data density_ice/0.917e+3/
  data density_cloudice/0.917e+3/
  data reff_ciw  / 0.20 /	       ! [mm]
      
!---Assign some useful constants
  wavel = 300./freqy
   

!---  Begin by checking if hydrometeors of this species are present.
!---  If not, set scattering parameters to zero and return.

  if(lwc .lt. 0.001) then
    ksca=0.
    asca=0.
    gsca=0.
    pbck=0.
    return
  endif

!--- Initialize the scattering parameters

  bext=0.
  bsca=0.
  bsym=0.
  bq11=0.

!--- Compute scattering properties

  rad = reff_ciw
  dropmass = 4/3.*pi*density_ice*rad*rad*rad*1.E-09
  num = lwc/dropmass

!--- Get complex refractive index of cloud ice

  call iceoptic(freqy,temp,epsreal,epsimag)
  eice = cmplx(epsreal,epsimag)
  cref = csqrt(eice)
  x = 2.*pi*rad/wavel

! call mie_sphere(x,cref,qsca,qext,asym,qbsca)

!--- call mie_ice with lookup table

  call mie_LJBtable_all(x,cref,qsca,qext,asym,qbsca)

  bext = num*qext*pi*rad*rad*1.e-6
  bsca = num*qsca*pi*rad*rad*1.e-6
  bsym = num*qsca*asym*pi*rad*rad*1.e-6
  bq11 = num*qbsca*pi*rad*rad*1.e-6

!--- check for distribution with very small extinction;
!--- set parameters to zero to avoid numerical problems
     
  if( bext .gt. 1.e-6) then
    ksca=bext
    asca=bsca/bext
    gsca=bsym/bsca
    pbck=bq11/bsca
  else
    ksca=0.
    asca=0.
    gsca=0.
    pbck=0.
  endif

  return
  end subroutine mie_ci_lut
      
!-------------------------------------------------------------------------

  subroutine mie_graup_lut(freqy, temp, lwc, ksca, asca, gsca)
      

!--- Compute the extinction, absorption, asymmetry parameter and
!--- backscatter for a given water content of graupel in [g/m^3], and
!--- a particle size distribution n(D) with intercept n0g.

!     Input:
!     freqy		frequency of radiation [GHz]
!     temp		temperature of particles [K]
!     lwc		water content of graupel distribution [g/m**3]
!     density_snow      density of graupel [gm/m^3]
!     N0                Intercept of expeonential DSD [1/m^4]

!     Output:
!     ksca		extinction coefficient [1/km]
!     asca		single-scatter albedo []
!     gsca		asymmetry factor []
!     pbck		backscatter phase function/(4*pi) []

  integer ::  i      
  real :: freqy, temp, lwc
  real :: ksca, asca, gsca, pbck
  real :: wavel, density_ice, density_graupel, rad, fincl
  real :: N0, lam, num, x
  real :: qext, qsca, asym, qbsca
  real :: bext, bsca, bsym, bq11      
  real :: epsreal, epsimag	
  complex :: eice, eair, emg, cref

!  data pi/3.14159265/
  data density_ice/0.917e+3/
  data density_graupel/400./
  data N0 / 4.0E+6 /  
      
!--- Assign some useful constants

  wavel = 300./freqy

!--- Begin by checking if hydrometeors of this species are present.
!--- If not, set scattering parameters to zero and return.

  if(lwc .lt. 0.0025) then
    ksca=0.
    asca=0.
    gsca=0.
    pbck=0.
    return
  endif

!--- If hydrometeors are present, initialize the scattering parameters

  bext=0.
  bsca=0.
  bsym=0.
  bq11=0.

!--- Loop over particle sizes:

!--- increments of particle radius are 0.05 mm; the particle
!--- size distribution is expressed as a particle number density,
!--- num, per radius increment; the original psd is
!--- n(D) = n0g * exp(-lam * D) where n is the number density
!--- per diameter increment, n0g is the distribution intercept,
!--- lam is the slope of the distribution of ln(n(D)), and D is
!--- the particle diameter. It is assumed here that n0g is
!--- and lwc are prescribed, and lam is therefore constrained 
!--- to yield the prescribed water content:
!--- lwc = integral {n(D) * pi * (D**3) * density(D) * dD/6}
!--- therefore:
!--- lam=(n0g*pi*density/lwc)**(0.25)
      
  do i=0,200
    rad=0.025+0.05*float(i)
    x =  2.*pi*rad/wavel 
    lam = (N0*pi*density_graupel/(lwc*(1.e-3)))**(0.25)

!--- num is in terms of radius - factors of 2 account for conversion
    
    num = 2.*N0*exp(-lam*2.*rad*(1.e-3))

!--- complex refractive index of graupel

    call iceoptic(freqy,temp,epsreal,epsimag)
    eice = cmplx(epsreal,epsimag)
    eair = cmplx(1.0006,0.0)

!--- calculate dielectric constant of grauple as
!--- ice matrix with air inclusions, using Maxwell-Garnett 
!--- mixing for 2-component media w. elliptical inclusions 

    fincl = 1.-(density_graupel/density_ice)
    call mg_ellips(fincl, eice, eair, emg)
    cref = csqrt(emg)

! call mie_sphere(x,cref,qsca,qext,asym,qbsca)

    call mie_LJBtable_all(x,cref,qsca,qext,asym,qbsca)

    qext = qext 
    qsca = qsca 
    asym = asym 
    qbsca = qbsca         

!--- integrate over particle size distribution;

    bext=bext+num*qext*pi*rad*rad*(0.05)*1.e-6
    bsca=bsca+num*qsca*pi*rad*rad*(0.05)*1.e-6
    bsym=bsym+num*qsca*asym*pi*rad*rad*(0.05)*1.e-6
    bq11=bq11+num*qbsca*pi*rad*rad*(0.05)*1.e-6 

  end do

!--- check for distribution with very small extinction;
!--- set parameters to zero to avoid numerical problems
    
  if( bext .gt. 1.e-6) then
    ksca=bext
    asca=bsca/bext
    gsca=bsym/bsca
    pbck=bq11/bsca
  else
    ksca=0.
    asca=0.
    gsca=0.
    pbck=0.
  end if

 return
 end subroutine mie_graup_lut
      
!----------------------------------------------------------------------------      

 subroutine mie_melt_lut(freqy,temp,lwc,Dm1,Nw1,Mu1,f_melt,ksca,asca,gsca)
                         

!--- Compute the extinction, absorption, asymmetry parameter and
!--- backscatter for a given partially melted snow water content.
!--- Refractive index and density of aprticles are computed based 
!--- upon the melted fraction of snow.   The original DSD of snow
!--- is assumed for melting particles as well.

!--- Input:
!--- freqy	    frequency of radiation [GHz]
!--- temp	    temperature of particles [K]
!--- lwc	    water content of snow distribution [g/m**3]
!--- melt_frac      melted fraction of snow water content [0-1]

!--- Output:
!---  ksca	     extinction coefficient [1/km]
!---  asca	     single-scatter albedo []
!---  gsca	     asymmetry factor []
!---  pbck	     backscatter phase function/(4*pi) []

  integer ::  i
  real  ::   freqy, temp, lwc, f_melt
  real  ::   ksca, asca, gsca, pbck
  real  ::   density_snow, density_ice, density_liq, density_melt 
  real  ::   wavel, fincl, N0_snow, N0, N0_rain, Dm1,Nw1,Mu1
  real  ::   Nw2, rad,  lam, num, x
  real  ::   num_rain, num_snow
  real  ::   qext, qsca, asym, qbsca
  real  ::   bext, bsca, bsym, bq11      
  real  ::   epsreal, epsimag	 
  real  ::   gamma,f1,f2,f3,f4
  complex ::  eice, ewat, eair, ei, ewi, cref, cref_snow, cref_rain

!  data pi/3.14159265/
  data  density_liq /1.0e+3 /
  data density_ice /0.917e+3/
  data density_snow/100./ 
  data    N0_rain / 8.E+06 /		   ! [1/m^4]
  data    N0_snow / 1.0e+8 /		   ! [1/m^4]
      
!---  Assign some useful constants
     
  wavel = 300./freqy

!--- Begin by checking if hydrometeors of this species are present.
!--- If not, set scattering parameters to zero and return.

  if(lwc .lt. 0.0025) then
    ksca=0.
    asca=0.
    gsca=0.
    pbck=0.
    return
  endif

!--- If hydrometeors are present, initialize the scattering parameters
!--- for each sublayer

  bext=0.
  bsca=0.
  bsym=0.
  bq11=0.

!--- increments of particle radius are 0.05 mm; the particle
!--- size distribution is expressed as a particle number density,
!--- num, per radius increment; the original psd is
!--- n(D) = Nws * exp(-lam * D) where n is the number density
!--- per diameter increment, Nws is the distribution intercept,
!--- lam is the slope of the distribution of ln(n(D)), and D is
!--- the particle diameter. It is assumed here that Nws is
!--- and lwc are prescribed, and lam is therefore constrained 
!--- to yield the prescribed water content:
!--- lwc = integral {n(D) * pi * (D**3) * density(D) * dD/6}
!--- therefore:
!--- lam=(Nws*pi*density_snow/lwc)**(0.25)

!--- prescribe N0 [1/m**4] for melt water
     
  N0 = (1.-f_melt)*N0_snow  + f_melt*N0_rain
  density_melt = (1.-f_melt)*density_snow + f_melt*density_liq
      
! write(*,*) 'melt',n_sublyr, f_melt, Nw, density_melt 



  do i=0,200
  
    rad = 0.025 + 0.05*float(i)        
    x =  2.*pi*rad/wavel 

!- Snow DSD

    lam = (N0_snow*pi*density_snow/(lwc*(1.e-3)))**(0.25)

!- num is in terms of radius - factors of 2 account for conversion
       
    num_snow = 2.*N0_snow*exp(-lam*2.*rad*(1.e-3))

!- Rain DSD

!- Convert Nw1 (given log10Nw)

    Nw2 = (10.**Nw1)*(1.0e-3)       
 
    gamma = 6.*5.*4.*3.*2*1.
    IF(Mu1.eq.2)gamma = 5.*4.*3.*2*1.
    f1 = (6./4.**4.)
    f2 = ((Mu1+4.0)**(Mu1+4))/(gamma)
    f3 = ((2.*rad)/Dm1)**Mu1 
    f4 = exp(-(4.+Mu1)*((2.*rad)/Dm1))
    num_rain = Nw2*f1*f2*f3*f4
	
!- num_rain = 2.*Nw*exp(-lam*2.*rad*(1.e-3))

!- complex refractive index melting snow constituents

    call iceoptic(freqy,temp,epsreal,epsimag)
    eice = cmplx(epsreal,epsimag)
    eair = cmplx(1.0006,0.0)

!- calculate dielectric constant of snow as an ice matrix  
!- with air inclusions, using Maxwell-Garnett mixing for  
!- 2-component media w. elliptical inclusions 

    fincl = 1. - (density_snow/density_ice)

    call mg_ellips(fincl, eice, eair, ei)
    cref_snow = csqrt(ei)
        
!- compute dielectric constant of water
       
    call watoptic(freqy,temp,0.0,epsreal,epsimag)
    ewat = cmplx(epsreal,epsimag)
    cref_rain = csqrt(ewat)

!- define number of dielectric constant of melting particles

    num = (1-f_melt)*num_snow + f_melt*num_rain
    cref = 0.5d0*cref_snow + 0.5d0*cref_rain

!- call Mie program
 
!call mie_sphere(x,cref,qsca,qext,asym,qbsca)	
!write(6,*)x,cref_snow, cref_rain, cref	

    call mie_LJBtable_all(x,cref,qsca,qext,asym,qbsca)

!-integrate over particle size distribution;
    
    bext=bext+num*qext*pi*rad*rad*(0.05)*1.e-6
    bsca=bsca+num*qsca*pi*rad*rad*(0.05)*1.e-6
    bsym=bsym+num*qsca*asym*pi*rad*rad*(0.05)*1.e-6
    bq11=bq11+num*qbsca*pi*rad*rad*(0.05)*1.e-6
  
  end do

!check for distribution with very small extinction;
!set parameters to zero to avoid numerical problems
 
  if( bext .gt. 1.e-6) then
      ksca=bext
      asca=bsca/bext
      gsca=bsym/bsca
      pbck=bq11/bsca
  else
      ksca=0.
      asca=0.
      gsca=0.
      pbck=0.
  endif

  return      
 end subroutine mie_melt_lut
      
!----------------------------------------------------------------------------      

 subroutine mie_rain_lut(freqy,temp,lwc,Dm1,Nw1,Mu1,ksca,asca,gsca)
      
!--- Changed 2015 by Sarah Ringerud to use normalized gamma DSD with 
!--- input Nw and Dm (mass weighted mean diameter) following 
!--- Bringi et al. 2003

!--- Compute the extinction, absorption, asymmetry parameter and
!--- backscatter for a given water content of rain in [g/m^3], and
!--- a particle size distribution n(D) with intercept N0 (N0 rain).

!--- Input:
!--- freqy  	   frequency of radiation [GHz]
!--- temp		   temperature of particles [K]
!--- lwc		   water content of rain distribution [g/m^3]
!--- Dm		   Intercept of exponential DSD [1/m^4]

!--- Output:
!--- ksca		   extinction coefficient [1/km]
!--- asca		   single-scatter albedo []
!--- gsca		   asymmetry factor []
!--- pbck		   backscatter phase function/(4*pi) []
!

  integer :: i      
  real    :: freqy, temp, lwc
  real    :: f1, f2,f3,f4
  real    :: ksca, asca, gsca, pbck
  real    :: wavel, density_liq, rad
  real    :: Dm1, Mu1, Nw1, lam, num, xm,  Nw2
  real    :: mom6, gamma, dr, dD, D,x
  real    :: qext, qsca, asym, qbsca
  real    :: bext, bsca, bsym, bq11	
  real    :: epsreal, epsimag	  
  complex :: ewat, cref

!  data pi/ 3.14159265 /
  data density_liq / 1.0 /   ! in g/cm^3
!      data Mu1   / 3. /           ! gamma DSD exponent
      
!-Assign some useful constants
      
  wavel=300./freqy
      

!-Begin by checking if hydrometeors of this species are present.
!-If not, set scattering parameters to zero and return.

  if(lwc .lt. 0.001) then
    ksca=0.
    asca=0.
    gsca=0.
    pbck=0.
    return
  endif

!-If hydrometeors are present, initialize the scattering parameters

  bext=0.
  bsca=0.
  bsym=0.
  bq11=0.
  mom6=0.

!- Loop over particle sizes:

!-increments of particle radius are 0.05 mm; the particle
!-size distribution is expressed as a particle number density,
!-num, per radius increment; the original psd is
!-n(D) = Nw * exp(-lam * D) where n is the number density
!-per diameter increment, Nw is the distribution intercept,
!-lam is the slope of the distribution of ln(n(D)), and D is
!-the particle diameter. It is assumed here that Nw and lwc  
!-are prescribed, and lam is therefore constrained to yield the  
!-prescribed water content:
!-lwc = integral {n(D) * pi * (D**3) * density(D) * dD/6}
!-therefore:
!-lam=(Nw*pi*density/lwc)**(0.25)

!-gamma (4+Mu1) = gamma(7) = 6!

!- Convert Nw (given log10Nw)
    
  Nw2 = (10.**Nw1)*(1.0e-3) 

  gamma = 6.*5.*4.*3.*2*1.
  IF(Mu1.eq.2)gamma = 5.*4.*3.*2*1.
  
  dr = 0.05
  dD = dr*2.
  do i = 0, 200
    rad = 0.5*dr + 0.05*float(i)
    D = 2.*rad
    x = 2.*pi*rad/wavel
    ! num is in terms of radius - factors of 2 account for conversion
    f1 = (6./4.**4.)
    f2 = ((Mu1+4.0)**(Mu1+4))/(gamma)
    f3 = (D/Dm1)**Mu1 
    f4 = exp(-(4.+Mu1)*(D/Dm1))
    num = Nw2*f1*f2*f3*f4

!- complex refractive index of liquid water

    call watoptic(freqy,temp,0.0,epsreal,epsimag)
    ewat = cmplx(epsreal,epsimag)
    cref = csqrt(ewat)

!-call Mie program 

! call mie_sphere(x,cref,qsca,qext,asym,qbsca)

!-call mie_water with lookup table

    call mie_LJBtable_all(x,cref,qsca,qext,asym,qbsca)

!- integrate over particle size distribution;

    bext=bext+num*qext*pi*rad*rad*dD*1e-3
    bsca=bsca+num*qsca*pi*rad*rad*dD*1e-3
    bsym=bsym+num*qsca*asym*pi*rad*rad*dD*1e-3
    bq11=bq11+num*qbsca*pi*rad*rad*dD*1e-3
    mom6=mom6*D**6*dD*1e-3

  end do


!- check for distribution with very small extinction;
!- set parameters to zero to avoid numerical problems
  if( bext .gt. 1.e-6) then
      ksca=bext
      asca=bsca/bext
      gsca=bsym/bsca
      pbck=bq11/bsca 
  else
      ksca=0.
      asca=0.
      gsca=0.
      pbck=0.
  end if

  return
 end subroutine mie_rain_lut
           
!-----------------------------------------------------------------------------------------------

 subroutine mie_snow_lut(freqy, temp, lwc, ksca, asca, gsca)
      
!- Compute the extinction, absorption, asymmetry parameter and
!- backscatter for a given water content of snow in [g/m^3], and
!- a particle size distribution n(D) with intercept n0s.

!- Input:
!- freqy	    frequency of radiation [GHz]
!- temp 	    temperature of particles [K]
!- lwc  	    water content of snow distribution [g/m**3]
!- density_snow      density of snow [kg/m^3]
!- N0		     Intercept of expeonential DSD [1/m^4]

!- Output:
!- ksca 	    extinction coefficient [1/km]
!- asca 	    single-scatter albedo []
!- gsca 	    asymmetry factor []
!- pbck 	    backscatter phase function/(4*pi) []

  integer :: i	  
  real :: freqy, temp, lwc
  real :: ksca, asca, gsca, pbck
  real :: wavel, density_ice, density_snow, rad, fincl
  real :: N0, lam, num, x
  real :: qext, qsca, asym, qbsca
  real :: bext, bsca, bsym, bq11      
  real :: epsreal, epsimag	 
  complex :: eice, eair, emg, cref

!  data pi/3.14159265/
  data density_ice/0.917e+3/
  data density_snow/100./  
  data N0 / 1.0e+8 /		   ! [1/m^4]
      
!- Assign some useful constants
      
  wavel = 300./freqy
      
!- Begin by checking if hydrometeors of this species are present.
!- If not, set scattering parameters to zero and return.
  if(lwc .lt. 0.0025) then
      ksca=0.
      asca=0.
      gsca=0.
      pbck=0.
      return
  endif    

!- If hydrometeors are present, initialize the scattering parameters

  bext=0.
  bsca=0.
  bsym=0.
  bq11=0.

!- Loop over particle sizes:

!- increments of particle radius are 0.05 mm; the particle
!- size distribution is expressed as a particle number density,
!- num, per radius increment; the original psd is
!- n(D) = n0s * exp(-lam * D) where n is the number density
!- per diameter increment, n0s is the distribution intercept,
!- lam is the slope of the distribution of ln(n(D)), and D is
!- the particle diameter. It is assumed here that n0s is
!- and lwc are prescribed, and lam is therefore constrained 
!- to yield the prescribed water content:
!- lwc = integral {n(D) * pi * (D**3) * density(D) * dD/6}
!- therefore:
!- lam=(n0s*pi*density_snow/lwc)**(0.25)

  do i=0,200
    rad=0.025+0.05*float(i)
    x =  2.*pi*rad/wavel 
    lam = (N0*pi*density_snow/(lwc*(1.e-3)))**(0.25)
 
!-- num is in terms of radius - factors of 2 account for conversion
       
    num = 2.*N0*exp(-lam*2.*rad*(1.e-3))

!-- complex refractive index of snow

    call iceoptic(freqy,temp,epsreal,epsimag)
    eice = cmplx(epsreal,epsimag)
    eair = cmplx(1.0006,0.0)

!-- calculate dielectric constant of snow as an ice matrix  
!-- with air inclusions, using Maxwell-Garnett mixing for  
!-- 2-component media w. elliptical inclusions 
    
    fincl = 1. - (density_snow/density_ice)
    call mg_ellips(fincl, eice, eair, emg)
    cref = csqrt(emg)

!-- call Mie program
 
! call mie_sphere(x,cref,qsca,qext,asym,qbsca)

    call mie_LJBtable_all(x,cref,qsca,qext,asym,qbsca)

    qext = qext 
    qsca = qsca 
    asym = asym 
    qbsca = qbsca         

!-- integrate over particle size distribution;

    bext=bext+num*qext*pi*rad*rad*(0.05)*1.e-6
    bsca=bsca+num*qsca*pi*rad*rad*(0.05)*1.e-6
    bsym=bsym+num*qsca*asym*pi*rad*rad*(0.05)*1.e-6
    bq11=bq11+num*qbsca*pi*rad*rad*(0.05)*1.e-6

  end do

!- check for distribution with very small extinction;
!- set parameters to zero to avoid numerical problems
     
  if( bext .gt. 1.e-6) then
      ksca=bext
      asca=bsca/bext
      gsca=bsym/bsca
      pbck=bq11/bsca
  else
      ksca=0.
      asca=0.
      gsca=0.
      pbck=0.
  end if

  return
 end subroutine mie_snow_lut
      
!----------------------------------------------------------------------------      
      
 SUBROUTINE MG_ELLIPS (FINCL, EMATRIX, EINCL, EMG)

  SAVE

!** Maxwell-Garnett formula for effective permittivity of 2-component media 
!** (elliptical inclusions) P. Bauer 1996
!**
!** FINCL     volume fraction of inclusions
!** EMATRIX   permittivity of matrix
!** EINCL     permittivity of inclusions
!**
!** EMG       effective permittivity

  REAL ::    FINCL
  COMPLEX  :: EMATRIX, EINCL, EMG, GAMMA, Q

  Q = (EINCL / (EINCL - EMATRIX)) * CLOG (EINCL / EMATRIX) - 1.0
  GAMMA = 2.0 * EMATRIX * Q / (EINCL - EMATRIX)

  EMG = ((1.0 - FINCL) * EMATRIX + FINCL * GAMMA * EINCL) / (1.0 - FINCL + FINCL * GAMMA) 

  RETURN
 END subroutine MG_ELLIPS

!---------------------------------------------------------------------------      
      
 SUBROUTINE MIE_SPHERE (X, MIN, QSCAT, QEXTI, ASYM, QBSCAT)

  SAVE

!-- Mie Routine P. Bauer 

  integer,parameter   ::   limitx = 1500

  REAL ::	 X
  REAL ::	MR, MI, N1, N2
  REAL ::	QSCAT, QEXTI, QABSO, ASYM, QBSCAT

  REAL ::	RFAC1, RFAC2
  REAL ::	RHELP1(2), RHELP2(2)

  COMPLEX ::	M, MX, MIN
  COMPLEX ::	CHELP1, CHELP2, CFAC1, CFAC2, CBSCAT

  COMPLEX ::	DN(0:LIMITX), WN(-1:LIMITX)
  COMPLEX ::	AN(LIMITX), BN(LIMITX)

  INTEGER ::	NEND
  INTEGER ::	I100, I101

  EQUIVALENCE (CHELP1, RHELP1(1))
  EQUIVALENCE (CHELP2, RHELP2(1))

!************************************************************************

  M	 = CONJG (MIN)
  CHELP1 = M
  MR	 =	  RHELP1(1)
  MI	 = -1.0 * RHELP1(2)
  
  MX   = M  * X
  N1   = MR * X
  N2   = MI * X

  IF (X .LE. 20000.0) NEND = X + 4.00 * X ** (1.0 / 3.0) + 2.0
  IF (X .LE.  4200.0) NEND = X + 4.05 * X ** (1.0 / 3.0) + 2.0
  IF (X .LE.	 8.0) NEND = X + 4.00 * X ** (1.0 / 3.0) + 1.0
  IF (NEND .LE.    5) NEND = 5
  IF (NEND .GT. LIMITX) NEND = LIMITX

  RFAC1     = SIN(N1) * SIN(N1) + SINH(N2) * SINH(N2)
  RHELP1(1) = SIN(N1) * COS(N1) / RFAC1
  RHELP1(2) = SINH(N2) * COSH(N2) / RFAC1

  DN (0) = CHELP1

  RHELP1(1) = COS(X)
  RHELP1(2) = -1.0E+00 * SIN(X)
  RHELP2(1) = SIN(X)
  RHELP2(2) = COS(X)

  WN (-1) = CHELP1
  WN ( 0) = CHELP2

  QEXTI  = 0.0
  QSCAT  = 0.0
  QBSCAT = 0.0
  QABSO  = 0.0
  ASYM   = 0.0 
  CBSCAT = CMPLX (0.0,0.0)

  DO I100 = 1, NEND
    DN (I100) = -1.0 * I100 / MX +  1.0 / (I100 / MX - DN (I100 - 1))
    WN (I100) = WN (I100 - 1) * (2.0 * I100 - 1.0) / X - WN (I100 - 2)

    CFAC1 = DN (I100) / M + I100 / X
    CFAC2 = M * DN (I100) + I100 / X

    CHELP1 = WN (I100)
    CHELP2 = WN (I100 - 1)

    AN (I100) = (CFAC1 * RHELP1 (1) - RHELP2 (1)) / (CFAC1 * CHELP1 - CHELP2)
    BN (I100) = (CFAC2 * RHELP1 (1) - RHELP2 (1)) / (CFAC2 * CHELP1 - CHELP2)

    CHELP1 = AN (I100)
    CHELP2 = BN (I100)

    RFAC1 = RHELP1 (1) + RHELP2 (1)
    RFAC2 = CABS(AN(I100)) * CABS(AN (I100)) + CABS(BN(I100)) * CABS(BN(I100))

    QEXTI  = QEXTI  + (2.0 * I100 + 1.0) * RFAC1
    QSCAT  = QSCAT  + (2.0 * I100 + 1.0) * RFAC2
    CBSCAT = CBSCAT + (2.0 * I100 + 1.0) * (-1.0) ** I100 * (AN(I100) - BN(I100))

    IF (I100 .EQ. 1) cycle

    CHELP1 = AN(I100 - 1) * CONJG(AN(I100)) + BN(I100 - 1) * CONJG(BN(I100))
    CHELP2 = AN(I100 - 1) * CONJG(BN (I100 - 1))

    I101 = I100 - 1
    RFAC1  = I101 * (I101 + 2) / (I101 + 1.0)
    RFAC2  = (2.0 * I101 + 1.0) / (I101 * (I101 + 1.0))
    
    ASYM = ASYM + RFAC1 * RHELP1 (1) + RFAC2 * RHELP2 (1)

  ENDDO

  QEXTI  = QEXTI * 2.0 / (X * X)
  QSCAT  = QSCAT * 2.0 / (X * X)
  ASYM   = ASYM  * 4.0 / (X * X * QSCAT)
  QBSCAT = CABS (CBSCAT) * CABS (CBSCAT) / (X * X)
  IF (QSCAT .GT. QEXTI) QSCAT = QEXTI

  RETURN
 END SUBROUTINE MIE_SPHERE

!-----------------------------------------------------------------------------

 subroutine watoptic(freqy,temp,salinity,epsreal,epsimag)

  implicit none

  real :: freqy,temp,salinity
  real :: epsreal,epsimag

!- internal variables      
  real :: freqhz, ctemp
  real :: omega, epsstat, trelax, fac1
  real, parameter :: epshigh = 4.90
  real, parameter :: pi = 3.141592654

  freqhz = freqy*1.0E+09
  ctemp  = temp - 273.16
  omega  = 2.*pi*freqhz

  epsstat = (87.134E+00 - 1.949E-01 * ctemp              &               
             - 1.276E-02 * ctemp * ctemp                 &            
             + 2.491E-04 * ctemp * ctemp * ctemp)        &        
             * (1.0E+00 + 1.613E-05 * salinity * ctemp   &  
             - 3.656E-03 * salinity                      &                     
             + 3.210E-05 * salinity * salinity           &
             - 4.232E-07 * salinity * salinity * salinity)

  trelax = (1.768E-11 - 6.086E-13 * ctemp		 &
            + 1.104E-14 * ctemp * ctemp		         &
            - 8.111E-17 * ctemp * ctemp * ctemp) 	 &   
            * (1.0E+00 + 2.282E-05 * salinity * ctemp    &
            - 7.638E-04  * salinity   			 &  
            - 7.760E-06  * salinity * salinity		 & 
            + 1.105E-08  * salinity * salinity * salinity)

  fac1    = 1.0 + omega*omega*trelax*trelax
  epsreal = epshigh + (epsstat - epshigh) / fac1
  epsimag = ((epsstat - epshigh) * omega * trelax) / fac1
  
  return
 end subroutine watoptic

!-----------------------------------------------------------------------------

 SUBROUTINE ICEOPTIC (freqy, temp, epsreal, epsimag)

!- Hufford (1991), see Brussard and Watson (1995), p.297
 
!- Input & output variables; eps the dielectric constant, epsilon
  real ::  freqy, temp
  real :: epsreal, epsimag
      
!- internal variables 
  real :: t_ice, theta, A , B    

  epsreal = 3.15

  if (temp .gt. 273.16) then
    t_ice = 273.16
  else
    t_ice = temp
  endif

  theta  = 300.0 / t_ice
  A  = 1.0E-04 * (50.4 + 62.0 * (theta - 1.0)) * EXP (-22.1 * (theta - 1.0))
  B  = 1.0E-04 * (0.633 / theta - 0.131) + (7.36E-04 * theta / (theta - 0.9927))   &
                * (7.36E-04 * theta / (theta - 0.9927))
  epsimag = A / freqy + B * freqy

  return
 end subroutine iceoptic

!-----------------------------------------------------------------------------

 subroutine mie_LJBtable_ALL(xin,crefin,qscaout,qextout,asymout,qbscaout)      
      
  real  	      :: xin
  complex	      :: crefin   
  real  	      :: qscaout,qextout,asymout,qbscaout      
  real, parameter     :: missing_mie_value=-99.99

!-- variables will be needed to apply 1D linear interpolation
!-- dimensions are assigned as follows
   
  real  	      :: ra, rb, rc, rd, re	 
  real  	      :: crefrcnt, creficnt, xcnt
  logical	      :: crefrfound, crefifound, xfound      
  real  	      :: r1,r2,i1,i2,x1,x2
  	
  type  	      :: interpolation_structure
      real	      :: qsca
      real	      :: qext
      real	      :: asym
      real	      :: qbsca
  end type interpolation_structure	
  type(interpolation_structure) :: bestcref_x1,bestcref_x2
      
!-- location closest cref real locations 
   
  crefrcnt=0
  crefrfound=.false.
 
  ra=1
  rb=ndielec_rlist_ALL
  rc=int((ra+rb)/2)
 
  do while ((crefrcnt.lt.10).and.(.not.crefrfound))
    if (real(crefin).ge.dielec_rlist_ALL(rc)) then
!--on the right-hand side
        ra=rc
        rb=rb
        if (rc.eq.int((ra+rb)/2)) crefrfound=.true.
        rc=int((ra+rb)/2)
    else
!--on the left-hand side
        ra=ra
        rb=rc-1
        if (rc.eq.int((ra+rb)/2)) crefrfound=.true.
        rc=int((ra+rb)/2)
    endif
    crefrcnt=crefrcnt+1.0
  enddo

! write(*,*) crefrcnt,rc
! write(*,*) dielec_rlist_ALL((rc-1):(rc+1))
      
  ra=rc-3
  rb=rc+3
  if (ra.lt.1) ra=1
  if (rb.gt.ndielec_rlist_ALL) rb=ndielec_rlist_ALL
  rd=10000.0
  re=0.0
  do rc=ra, rb, 1
    if (abs(real(crefin)-dielec_rlist_ALL(rc)).lt.rd) then
	rd=abs(real(crefin)-dielec_rlist_ALL(rc))
	re=rc
    endif
  enddo
    
  if (real(crefin).ge.dielec_rlist_ALL(re)) then
      ra=re
      rb=re+1
      if (rb.gt.ndielec_rlist_ALL) then
	  ra=ndielec_rlist_ALL-1
	  rb=ndielec_rlist_ALL
	  crefin=cmplx(dielec_rlist_ALL(rb)-1.0e-6,aimag(crefin))
      endif
  else
      ra=re-1
      rb=re
      if (ra.lt.1) then
	  ra=1
	  rb=2
	  crefin=cmplx(dielec_rlist_ALL(ra)+1.0e-6,aimag(crefin))
      endif
  endif
  r1=ra
  r2=rb
      
!      write(*,*) 'Neighbourng CrefR Values: '   &, 
!                 dielec_rlist_ALL(r1),dielec_rlist_ALL(r2)
 
 
!-- location closest cref imaginary locations      
  
  creficnt=0
  crefifound=.false.
      
  ra=1
  rb=ndielec_ilist_ALL
  rc=int((ra+rb)/2)
      
  do while ((creficnt.lt.15).and.(.not.crefifound))
    if (aimag(crefin).ge.dielec_ilist_ALL(rc)) then
!-- on the right-hand side
	  ra=rc
	  rb=rb
	  if (rc.eq.int((ra+rb)/2)) crefifound=.true.
	  rc=int((ra+rb)/2)
    else
!-- on the left-hand side
	  ra=ra
	  rb=rc-1
	  if (rc.eq.int((ra+rb)/2)) crefifound=.true.
	  rc=int((ra+rb)/2)
    endif
    creficnt=creficnt+1.0
  enddo

!  write(*,*) creficnt,rc
!  write(*,*) dielec_ilist_ALL((rc-1):(rc+1))

  ra=rc-3
  rb=rc+3
  if (ra.lt.1) ra=1
  if (rb.gt.ndielec_ilist_ALL) rb=ndielec_ilist_ALL
  rd=10000.0
  re=0.0
  do rc=ra, rb, 1
    if (abs(aimag(crefin)-dielec_ilist_ALL(rc)).lt.rd) then
       rd=abs(aimag(crefin)-dielec_ilist_ALL(rc))
       re=rc
    endif
  enddo
  if (aimag(crefin).ge.dielec_ilist_ALL(re)) then
      ra=re
      rb=re+1
      if (rb.gt.ndielec_ilist_ALL) then
         ra=ndielec_ilist_ALL-1
         rb=ndielec_ilist_ALL
         crefin=cmplx(real(crefin),dielec_ilist_ALL(rb)-1.0e-6)
      endif
  else
      ra=re-1
      rb=re
      if (ra.lt.1) then
          ra=1
          rb=2
          crefin=cmplx(real(crefin),dielec_ilist_ALL(ra)+1.0e-6)
      endif
  endif
  i1=ra
  i2=rb
  
!-- write(*,*) 'Neighbourng CrefI Values: ',    &
!            dielec_ilist_ALL(i1),dielec_ilist_ALL(i2)
          
!-- location closest xsize locations      
  xcnt=0
  xfound=.false.
      
  ra=1
  rb=nxlist_ALL
  rc=int((ra+rb)/2)
  
  do while ((xcnt.lt.10).and.(.not.xfound))
    if (xin.ge.xlist_ALL(rc)) then
      !on the right-hand side
      ra=rc
      rb=rb
      if (rc.eq.int((ra+rb)/2)) xfound=.true.
      rc=int((ra+rb)/2)
    else
      !on the left-hand side
      ra=ra
      rb=rc-1
      if (rc.eq.int((ra+rb)/2)) xfound=.true.
      rc=int((ra+rb)/2)
    endif
    xcnt=xcnt+1.0
  enddo!

!-- write(*,*) xcnt,rc
!-- write(*,*) xlist_ALL((rc-1):(rc+1))
 
  ra=rc-3
  rb=rc+3
  if (ra.lt.1) ra=1
  if (rb.gt.nxlist_ALL) rb=nxlist_ALL
  rd=10000.0
  re=0.0
  do rc=ra, rb, 1
    if (abs(xin-xlist_ALL(rc)).lt.rd) then
      rd=abs(xin-xlist_ALL(rc))
      re=rc
    endif
  enddo
  if (xin.ge.xlist_ALL(re)) then
      ra=re
      rb=re+1
      if (rb.gt.nxlist_ALL) then
         ra=nxlist_ALL-1
         rb=nxlist_ALL
         xin=xlist_ALL(rb)-1.0e-6
      endif
  else
      ra=re-1
      rb=re
      if (ra.lt.1) then
          ra=1
          rb=2
          xin=xlist_ALL(ra)+1.0e-6
      endif
  endif
  x1=ra
  x2=rb
      
!    write(*,*) 'Neighbourng Xsize Values: ', 
!     >            xlist_ALL(x1),xlist_ALL(x2)

!-APPLY 2D (for crefr+crefi) + 1D(for x) Linear Interpolation!
!-For all variables at x=x1
   
  bestcref_x1%qsca=1.0/( (crefindex_ALL(r2,i1)%real_cref       &
    			  -crefindex_ALL(r1,i1)%real_cref)     &
    			*(crefindex_ALL(r1,i2)%imag_cref       &
    			  -crefindex_ALL(r1,i1)%imag_cref))    &
   * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x1)%qsca &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))      &
    	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))    &
      + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x1)%qsca &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)      &
  	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))    &
      + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x1)%qsca &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))      &
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref))    &
      + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x1)%qsca &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)      &
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))


  bestcref_x1%qext=1.0/( (crefindex_ALL(r2,i1)%real_cref       &
    			  -crefindex_ALL(r1,i1)%real_cref)     &
    			*(crefindex_ALL(r1,i2)%imag_cref       &
    			  -crefindex_ALL(r1,i1)%imag_cref))    &
   * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x1)%qext &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))      &
    	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))    &
      + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x1)%qext &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)      &
  	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))    &
      + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x1)%qext &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))      &
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref))    &
      + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x1)%qext &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)      &
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))
  
  bestcref_x1%asym=1.0/( (crefindex_ALL(r2,i1)%real_cref       &
    			  -crefindex_ALL(r1,i1)%real_cref)     &
    			*(crefindex_ALL(r1,i2)%imag_cref       &
    			  -crefindex_ALL(r1,i1)%imag_cref))    &
   * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x1)%asym &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))      &
    	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))    &
      + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x1)%asym &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)      &
  	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))    &
      + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x1)%asym &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))      &
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref))    &
      + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x1)%asym &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)      &
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))   
  
  bestcref_x1%qbsca=1.0/( (crefindex_ALL(r2,i1)%real_cref	&
    			  -crefindex_ALL(r1,i1)%real_cref)	&
    			*(crefindex_ALL(r1,i2)%imag_cref	&
    			  -crefindex_ALL(r1,i1)%imag_cref))	&
   * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x1)%qbsca &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))	&
    	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))	&
      + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x1)%qbsca &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)	&
  	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))	&
      + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x1)%qbsca &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))	&
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref))	&
      + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x1)%qbsca &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)	&
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))	      
  

!-- For all variables at x=x2

  bestcref_x2%qsca=1.0/( (crefindex_ALL(r2,i1)%real_cref	&
    			  -crefindex_ALL(r1,i1)%real_cref)	&
    			*(crefindex_ALL(r1,i2)%imag_cref	&
    			  -crefindex_ALL(r1,i1)%imag_cref))	&
   * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x2)%qsca  &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))	&
    	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))	&
      + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x2)%qsca  &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)	&
  	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))	&
      + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x2)%qsca  &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))	&
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref))	&
      + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x2)%qsca  &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)	&
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))


  bestcref_x2%qext=1.0/( (crefindex_ALL(r2,i1)%real_cref	&
    			  -crefindex_ALL(r1,i1)%real_cref)	&
    			*(crefindex_ALL(r1,i2)%imag_cref	&
    			  -crefindex_ALL(r1,i1)%imag_cref))	&
   * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x2)%qext  &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))	&
    	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))	&
      + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x2)%qext  &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)	&
  	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))	&
      + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x2)%qext  &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))	&
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref))	&
      + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x2)%qext  &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)	&
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))
  
  bestcref_x2%asym=1.0/( (crefindex_ALL(r2,i1)%real_cref	&
    			  -crefindex_ALL(r1,i1)%real_cref)	&
    			*(crefindex_ALL(r1,i2)%imag_cref	&
    			  -crefindex_ALL(r1,i1)%imag_cref))	&
   * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x2)%asym  &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))	&
    	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))	&
      + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x2)%asym  &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)	&
  	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))	&
      + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x2)%asym  &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))	&
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref))	&
      + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x2)%asym  &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)	&
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))
  
  bestcref_x2%qbsca=1.0/( (crefindex_ALL(r2,i1)%real_cref	&
    			  -crefindex_ALL(r1,i1)%real_cref)	&
    			*(crefindex_ALL(r1,i2)%imag_cref	&
    			  -crefindex_ALL(r1,i1)%imag_cref))	&
   * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x2)%qbsca &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))	&
    	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))	&
      + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x2)%qbsca &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)	&
  	   *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin)))	&
      + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x2)%qbsca &
    	 *( (crefindex_ALL(r2,i1)%real_cref-real(crefin))	&
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref))	&
      + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x2)%qbsca &
    	 *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref)	&
    	   *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))	      
  
!--Interpolate between x1 and x2

  qscaout=bestcref_x1%qsca+(bestcref_x2%qsca-bestcref_x1%qsca)     &
    	 *(xin-xlist_ALL(x1))/(xlist_ALL(x2)-xlist_ALL(x1))

  qextout=bestcref_x1%qext+(bestcref_x2%qext-bestcref_x1%qext)     &
    	 *(xin-xlist_ALL(x1))/(xlist_ALL(x2)-xlist_ALL(x1))

  asymout=bestcref_x1%asym+(bestcref_x2%asym-bestcref_x1%asym)     &
    	 *(xin-xlist_ALL(x1))/(xlist_ALL(x2)-xlist_ALL(x1))
  
  qbscaout=bestcref_x1%qbsca+(bestcref_x2%qbsca-bestcref_x1%qbsca) &
     	     *(xin-xlist_ALL(x1))/(xlist_ALL(x2)-xlist_ALL(x1))          


  if ((qscaout.lt.0).or.(qextout.lt.0).or.(qbscaout.lt.0)) then
      if (min(qscaout,qextout,qbscaout).lt. (missing_mie_value/2.0)) then

!---the value need to be interpolated is closer to the NAN, instead of any valid numbers
	
	  qscaout=0.0/0.0
	  qextout=0.0/0.0
	  asymout=0.0/0.0
	  qbscaout=0.0/0.0		
      
      else
	  
!--- the value need to be interpolated is closer to a nearby valid point, instead of NAN	  
	  
          if ((bestcref_x1%qsca.gt.0) .and.   &
     	      (bestcref_x1%qext.gt.0) .and.   &
     	      (bestcref_x1%qbsca.gt.0)) then      !--if interpolation at x1 is good
	       qscaout=bestcref_x1%qsca
	       qextout=bestcref_x1%qext
	       asymout=bestcref_x1%asym
	       qbscaout=bestcref_x1%qbsca
	  else
	    if ((bestcref_x2%qsca.gt.0) .and.  &
     	        (bestcref_x2%qext.gt.0) .and.  &
     		(bestcref_x2%qbsca.gt.0))then     !--- if interpolation at x2 is good		
	           qscaout=bestcref_x2%qsca
	           qextout=bestcref_x2%qext
	           asymout=bestcref_x2%asym
	           qbscaout=bestcref_x2%qbsca	        
	    else                                  !---if interpolation is bad at both x1 and x2
	           qscaout=0.0/0.0
	           qextout=0.0/0.0
	           asymout=0.0/0.0
	           qbscaout=0.0/0.0	      
	    endif  !for if interpolation at x2	  
	  endif!for if interpolation at x1
	
      endif!for if relative location of interpolated result	

  endif!for if results are unreasonable
      
  return
 end subroutine mie_LJBtable_ALL   
      
!-----------------------------------------------------------------------  

      end module GPM_TbsimulatorV7_mie_code
