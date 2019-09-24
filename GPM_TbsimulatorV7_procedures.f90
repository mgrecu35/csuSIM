module GPM_TbsimulatorV7_procedures

 Use GPM_TbsimulatorV7_variables
 Use GPM_TbsimulatorV7_rad_procedures 
 Use GPM_TbsimulatorV7_mie_code
 Use GPM_TbsimulatorV7_IO
 Use MONOrtm

 implicit none

 contains 

!---------------------------------------------------------------------- 

 subroutine GPM_compute_Tbs(ichan,iCMB)
  
 !--- this routine will take the combined/mirs profiles and run the Radiative
 !--- transfer forward model over these profiles

 !--- general program variables
   include 'fordef.for'
  character(len=256) :: mietable_all_filenamein
  character(len=100) :: lut_file
  character(len=1)   :: cpol(0:1)=(/'V','H'/)
  
  integer :: iscan, ipix, nf, ilev, jchan, verbose, i, z, k 
  real    :: ci
 
  real    :: tavg,pavg,umu,kextcw,kextsn,salbsn,asymsn
  real    :: kextci,salbci,asymci   
  real    :: atm_ext,kextrr,albrr,salbrr,asymrr
  real    :: emisdif, oceanfrac
 
  real,allocatable :: diff(:), freqs(:), pol(:), ipol(:), buf(:,:)  
  real,parameter   :: fastemfreq(8) = (/10.65,18.7,23.8,36.64, &
                                        89.0,166.0,180.31,190.31/)
  integer,parameter :: fastempol(8)  = (/2,2,0,2,2,2,1,1/)           !0=V,1=H,2=both				
  real,parameter    :: fastemEIA(8)  =(/52.80,52.80,52.80,52.80,  &
                                       52.80,49.19,49.19,49.19/) 
  
  integer, parameter :: Kuofreqindex(13) = (/1,1,2,2,3,4,4,5,5,6,6,7,7/)
  real  :: ebartemp, emistemp(2)
 
  real    :: mp_kexttot(mplyrs),mp_salbtot(mplyrs),mp_asymtot(mplyrs)
  real    :: sum_kext, sum_salb, sum_asym
  integer :: n_sublyr, ichan
  real :: kext1(8),salb1(8),asym1(8)
  integer:: nfreq1, iCMB

  nfreq1=8
  
 !--- allocate the variables.  (npixs,nscans,nlevs) come from read_strmfile
  if(ichan==0) then
     allocate(kexttot(npixs,nscans,nlevs))
     allocate(salbtot(npixs,nscans,nlevs))
     allocate(asymtot(npixs,nscans,nlevs))
     
     allocate(emis(npixs,nscans))
     allocate(ebar(npixs,nscans))
     allocate(tb(npixs,nscans))
     
     allocate(dTbssim(nchans,npixs,nscans))
     allocate(deltadTbssim(nchans,npixs,nscans))
     allocate(emissall(nchans,npixs,nscans))
     allocate(ebarall(nchans,npixs,nscans))
     allocate(lambert(npixs,nscans))
     
     allocate(diff(nchans), freqs(nchans), pol(nchans), ipol(nchans), buf(nlevs,5))
     
     !--- start the radiative transfer computation
     
     write(*,*)' starting to compute Tbs'
     
     !--- read in the mie scattering tables for ice and water
     
     write(*,*)'  read in the mie scattering tables for ice and water'  
     mietable_all_filenamein = 'external_files/binary_MIE_original_ALL_mis-99_2pct.tbl'
     
     call read_mietable_all(mietable_all_filenamein)
     
     !--- Read in the Kuo Scattering Table
     
     write(*,*)'  reading Kuo scattering table'
     Kuoscattable = 'external_files/Kuo_Scattering_Table.txt'
     call read_Kuo_scattering_table
     
     !---  read the MONO RTM look-up table for absorption
     
     write(*,*)'  reading MonoRTM-GMI file'
     verbose = 0
     lut_file = 'external_files/MonoRTM-GMI.tbl'
     call read_lut(lut_file, verbose) 
     
     !--- start loop over all frequencies
     print*, 'reading CMB scattering tables'
     call readTablesLiang2(5,8)  !5 mus (2,2,2,2,2) 8 frequencies (10,19,22,36,89,166,183.3+3,183.3+7)
     print*, 'done reading tables and allocating memory'
     return
  end if
  write(*,*)'  starting  freq calculations'
  do nf =   ichan, ichan    !total number of frequencies to run the RT over (freq(1:nchan)
     write(*,*)
     write(*,'(a,i2,a,i2,a,f6.2,a1)')'    starting frequency: ',nf,'/', &
                                  nchans,' ',gmifreqs(nf),cpol(gmipol(nf))

    write(*,*)'    starting Mie Assignment from lookup table, and other RT coeff'
    kexttot = dflt
    salbtot = dflt
    asymtot = dflt
    atm_ext = 0.0

    do iscan = 1801,2100!1,nscans
      do ipix = 1,npixs
      
!---  loop over all layers for absorption and mie scattering calculations

        do z = nlevs,1,-1
	  
          tavg = 0.5*(tempprof(z,ipix,iscan) + tempprof(z-1,ipix,iscan))

          pavg = (pressprof(z,ipix,iscan)-pressprof(z-1,ipix,iscan))/ &
	      log(pressprof(z,ipix,iscan)/pressprof(z-1,ipix,iscan))
          if(pressprof(z,ipix,iscan) .eq. pressprof(z-1,ipix,iscan)) then
	        pavg=pressprof(z,ipix,iscan)
          endif
	     
          call monortm_lut(nf,pavg,tavg,mixr(z,ipix,iscan), atm_ext) !use MONORTM lookup 

          if(cloudwater(z,ipix,iscan) .lt. 0.0) cloudwater(z,ipix,iscan) = 0.0
          call absorb_clw(gmifreqs(nf),tavg,cloudwater(z,ipix,iscan),kextcw)   !cloud H20 absorption 
                                                                               !(scattering negligible)
          if(fp_class(snowwater(z,ipix,iscan)).eq.FOR_K_FP_POS_INF) then
             snowwater(z,ipix,iscan)=0
          endif
          if(fp_class(rainwater(z,ipix,iscan)).eq.FOR_K_FP_POS_INF) then
             rainwater(z,ipix,iscan)=0
          endif
          call mie_rain_lut(gmifreqs(nf),tavg,rainwater(z,ipix,iscan),   &        !rain water
                          dm(z,ipix,iscan), nw(z,ipix,iscan),   &
                          mu(z,ipix,iscan), kextrr, salbrr, asymrr) 
          if(rainwater(z,ipix,iscan)>1e-3.and.iCMB==1) then
             call getrainwp(nw(z,ipix,iscan)-6.90309,rainwater(z,ipix,iscan),kext1,salb1,asym1,nfreq1)
             if(rainwater(z,ipix,iscan)>1.5) then
                print*, nw(z,ipix,iscan)-6.90309,rainwater(z,ipix,iscan)
                print*, kextrr,salbrr,asymrr
                print*, kext1(Kuofreqindex(nf)), salb1(Kuofreqindex(nf)), asym1(Kuofreqindex(nf))
             endif
             kextrr=kext1(Kuofreqindex(nf))
             salbrr=salb1(Kuofreqindex(nf))
             asymrr=asym1(Kuofreqindex(nf))
             
          endif

	  call Kuo_intplte_snow(kuofreqindex(nf),snowwater(z,ipix,iscan),kextsn,salbsn,asymsn)
          if(snowwater(z,ipix,iscan)>1e-3.and.iCMB==1) then
             call getsnowp2(nw(z,ipix,iscan)-6.90309,snowwater(z,ipix,iscan),kext1,salb1,asym1,nfreq1)
             !print*, kextsn,salbsn,asymsn
             !print*, kext1(Kuofreqindex(nf)), salb1(Kuofreqindex(nf)), asym1(Kuofreqindex(nf))
             if(isnan(kextsn)) then
                print*, snowwater(z,ipix,iscan)
                stop
             endif
              
             kextsn=kext1(Kuofreqindex(nf))
             salbsn=salb1(Kuofreqindex(nf))
             asymsn=asym1(Kuofreqindex(nf))
          endif
!---     Total up the Kext
          
          kexttot(ipix,iscan,z) = atm_ext + kextcw + kextrr + kextsn !+ kextci   !sum all scattering params
	      
          if(kexttot(ipix,iscan,z) .gt. 1.e-06 )then
             salbtot(ipix,iscan,z) = (salbrr*kextrr + salbsn*kextsn) /  &
	                              kexttot(ipix,iscan,z)				      	     
	  else
             salbtot(ipix,iscan,z) = 0.
          endif
              
          if(salbtot(ipix,iscan,z) .gt. 1.e-06 )then
             asymtot(ipix,iscan,z)  =     &
	      (asymrr*salbrr*kextrr + asymsn*salbsn*kextsn)  /   &
	      (salbtot(ipix,iscan,z) * kexttot(ipix,iscan,z))     
          else
             asymtot(ipix,iscan,z) = 0.
          endif		    	     
 
         if(ipix .eq. tx .and. iscan .eq. -ty) then
	   if(z .eq. nlevs) then
	        write(6,'(a,a,a)')'level   mixr    atm_ext | cloudwater kextCLW | rainprof kextrain salbrain  asymrain |',   &
	                          '  snowprof  kextsn   salbsnow  asymsn  |  kexttot   salbtot   asymtot'		      
                write(6,'(a,a,a)')'-------------------------------------------------------------------------------------',   &
		                  '---------------------------------------------------------------------'               	   
	   endif
	   write(6,'(i3, 15F10.4)') z,mixr(z,ipix,iscan),atm_ext,     &
	                     cloudwater(z,ipix,iscan),kextcw, &
		             rainwater(z,ipix,iscan), kextrr,salbrr,asymrr,  &
	                     snowwater(z,ipix,iscan), kextsn,salbsn,asymsn,  &
	                     kexttot(ipix,iscan,z), salbtot(ipix,iscan,z), asymtot(ipix,iscan,z)
         endif     

        enddo    !nylr (z)    
	 
      enddo   !ipix
    enddo    !iscan 
    
!--- set the emiss and ebar for this frequency that's used inside the eddington routine
   
    emis(:,:) = emiss_strm(nf,:,:)
    ebar(:,:) = emiss_strm(nf,:,:)
    
    umu = cos(gmieia(nf) * 3.14159/180.)
    goto 10
    write(*,*)'    freq              = ',nf   
    write(*,*)'    umu               = ',umu
    write(*,*)'    tempprof(0,tx,ty) = ',tempprof(0,tx,ty)
    write(*,*)'    pressprof(0,tx,ty)= ',pressprof(0,tx,ty)
    write(*,*)'    fov_dpr           = ',fov_dpr
    write(*,*)'    retr_flag         = ',retr_flag(tx,ty)
    write(*,*)'    emis,ebar(tx,ty)  = ',emis(tx,ty),ebar(tx,ty)
    write(*,*)'    sfctype(tx,ty)    = ',sfctype(tx,ty)
    write(*,*)' CMB precip(tx,ty)    = ',csfcprcp(tx,ty)
    write(*,*)' mirs precip(tx,ty)   = ',msfcprcp(tx,ty)
    write(6,*)' wind speed           = ',wind(tx,ty)
    
    write(*,*)
    write(*,'(a)')'lev   plevel   airtemp   cldwater  snow/graupel rainwater  mixr   kexttot   salbtot   asymptot'     
    10 continue 
    do i = -88,1,-1   
      write(6,'(i3,9F10.4)')i,pressprof(i,tx,ty),tempprof(i,tx,ty),cloudwater(i,tx,ty),  &
                snowwater(i,tx,ty), rainwater(i,tx,ty), mixr(i,tx,ty),kexttot(tx,ty,i),  &  
	        salbtot(tx,ty,i),asymtot(tx,ty,i)
    enddo
    i=0
    write(6,'(i3,2F10.4)')i,pressprof(i,tx,ty),tempprof(i,tx,ty)
    

!--- compute brightness temperatures using eddington routine
!--- with slant path (at dpr resolution)

    write(*,*)
!    write(*,*)'    calling eddington '
    write(*,*)'    eddington Tb routine commented out for now'
          
!    lambert = .false.  
    call eddington_sp(umu)       
    dTbssim(nf,:,:) = Tb(:,:)    
!    write(6,*)' Tbs for tx,ty = ', Tb(tx,ty)
    
    
!---- need to assign deltaTb here ******************************************************    
    
    deltadTbssim(nf,:,:) = Tb(:,:)
        
 
  enddo  ! nf  - channel loop

  !write(*,*)
  !write(6,*)' tx,ty pixel values'
  !write(*,*)' freq/pol    DRsimTbs    GMITbobs   CombsimTbs     emiss       ebar'
  !do nf = 1,nchans
  !  write(6,'(F8.2,1x,a1,5F12.3)') gmifreqs(nf),cpol(gmipol(nf)), dTbssim(nf,tx,ty),   &
  !           gTbsobs(nf,tx,ty),cTbssim(nf,tx,ty),emissall(nf,tx,ty),ebarall(nf,tx,ty)
  !enddo

 end subroutine GPM_compute_Tbs

!----------------------------------------------------------------------

 subroutine GPM_compute_fpa
  
 !--- this routine will take the combined/mirs Tb retrieval and calculate
 !--- the sensor footprint averages
 
 
 !--- general subroutine variables
 
  integer,parameter :: outchans = 13
  real,parameter    :: minTb = 45, maxTb = 340
  integer,parameter :: fpfov = 3
  integer           :: icnt(5) = 0
  real,parameter    :: fov_dpr = 5.0  
  real              :: xa, ya, p, wgt      
 
  real :: fov_crosstrack(13), fov_downtrack(13)
 
  integer :: badchan(outchans), pix,scan,i,j,k, n_xt,n_dt, nf, nc, nl
  integer :: nx,ny
  logical :: firstpixres = .true.
  real    :: sum_wgt, sum_wgt_mirs, sum_wgt_gTbsobs, sum_wgt_cTbssim
  real    :: sum_wgt_mTbssim, sum_wgt_dTbssim

!--- results printing
  
  character(len=128) :: out_label='        10V     10H     19V     19H   23.8V      '//     &
            '     36.5V   36.5H     89V     89H    166V    166H' //   &
            '         183+/-3 183+/-7'
 
!--- sums and cnts for footprint averages  

  real(8) :: sum_datasource,sum_emiss(outchans)                     !surface variables
  real(8) :: sum_csfcprcp, sum_msfcprcp, sum_ccnvprcp		  
  real(8) :: sum_cTbssim, sum_mTbssim           !Tbs variables
  real(8) :: sum_dTbssim, sum_deltadTbssim, sum_gTbsobs
  real(8) :: sum_rainwater(nlevs)		          	    !layered variables
  real(8) :: sum_snowwater(nlevs)
  real(8) :: sum_cloudwater(nlevs)
  real(8) :: sum_SLH(nlevs)
  
!  real(8) :: cnt_datasource,cnt_emiss(outchans)                     !cnt surface variables
!  real(8) :: cnt_csfcprcp, cnt_msfcprcp, cnt_ccnvprcp		    
!  real(8) :: cnt_cTbssim(outchans), cnt_mTbssim(outchans)           !cnt Tbs variables
!  real(8) :: cnt_dTbssim(outchans), cnt_deltadTbssim(outchans)
!  real(8) :: cnt_rainwater(nlevs)		          	    !cnt layered variables
!  real(8) :: cnt_snowwater(nlevs)
!  real(8) :: cnt_cloudwater(nlevs)
!  real(8) :: cnt_SLH(nlevs)

!-- diagnostic variables  - sums, cnts, and averages of all pixels in single orbit

  real(8) :: sum_datasource_orb, sum_emiss_orb(13)                  !surface variables
  real(8) :: sum_csfcprcp_orb, sum_msfcprcp_orb, sum_ccnvprcp_orb
  real(8) :: sum_cTbssim_orb(13), sum_mTbssim_orb(13)               !Tbs variables
  real(8) :: sum_dTbssim_orb(13), sum_deltadTbssim_orb(13)
  real(8) :: sum_gTbsobs_orb(13)
  real(8) :: sum_rainwater_orb(nlevs)		          	    !layered variables
  real(8) :: sum_snowwater_orb(nlevs)
  real(8) :: sum_cloudwater_orb(nlevs)
  real(8) :: sum_SLH_orb(nlevs)

  real(8) :: cnt_datasource_orb, cnt_emiss_orb(13)                  !surface variables
  real(8) :: cnt_csfcprcp_orb, cnt_msfcprcp_orb, cnt_ccnvprcp_orb
  real(8) :: cnt_cTbssim_orb(13), cnt_mTbssim_orb(13)               !Tbs variables
  real(8) :: cnt_dTbssim_orb(13), cnt_deltadTbssim_orb(13)
  real(8) :: cnt_gTbsobs_orb(13)
  real(8) :: cnt_rainwater_orb(nlevs)		          	    !layered variables
  real(8) :: cnt_snowwater_orb(nlevs)
  real(8) :: cnt_cloudwater_orb(nlevs)
  real(8) :: cnt_SLH_orb(nlevs)

  real :: ave_datasource_orb, ave_emiss_orb(13) 		    !surface variables
  real :: ave_csfcprcp_orb, ave_msfcprcp_orb, ave_ccnvprcp_orb
  real :: ave_cTbssim_orb(13), ave_mTbssim_orb(13) 		    !Tbs variables
  real :: ave_dTbssim_orb(13), ave_deltadTbssim_orb(13)
  real :: ave_gTbsobs_orb(13)
  real :: ave_rainwater_orb(nlevs)				    !layered variables
  real :: ave_snowwater_orb(nlevs)
  real :: ave_cloudwater_orb(nlevs)
  real :: ave_SLH_orb(nlevs)
  
 !--- Allocate memory for footprint average arrays
       
  allocate (datasource_fpa(npixs,nscans), emiss_fpa(13,npixs,nscans))   !surface fpa variables     
  allocate (csfcprcp_fpa(npixs,nscans),msfcprcp_fpa(npixs,nscans))   
  allocate (ccnvprcp_fpa(npixs,nscans))

  allocate (cTbssim_fpa(outchans,npixs,nscans))                          !footprint Tbs
  allocate (mTbssim_fpa(outchans,npixs,nscans))
  allocate (dTbssim_fpa(outchans,npixs,nscans))
  allocate (deltadTbssim_fpa(outchans,npixs,nscans))
!  allocate (gTbsobs_fpa(outchans,npixs,nscans))

  allocate (rainwater_fpa(nlevs,npixs,nscans))	                        !footprint layered variable       
  allocate (cloudwater_fpa(nlevs,npixs,nscans))
  allocate (snowwater_fpa(nlevs,npixs,nscans))	  
  allocate (SLH_fpa(nlevs,npixs,nscans))

!--- set GMI footprint pixel sizes

 fov_crosstrack(1:outchans) = (/19.4, 19.4, 10.9, 10.9, 9.7,  &     !1/2 power 'crosst' (km)
                                 9.4,  9.4,  4.4,  4.4, 4.1,  &
                                 4.1,  3.8,  3.8 /)
 fov_downtrack(1:outchans)  = (/32.1, 32.1, 18.1, 18.1, 16.0, &     !1/2 power 'downt' (km)
                                15.6, 15.6,  7.2,  7.2, 6.3,  & 
                                6.3,  5.8,  5.8/)
 
!--- Error check the dr simulated Tbs
    
  do nf = 1, outchans
  
    if(gmifreqs(nf) .eq. dflt) cycle
         
       do scan = 1,nscans
         do pix = 1,npixs	   
	   if(dTbssim(nf,pix,scan) .ge. mintb .and. dTbssim(nf,pix,scan) .le. maxtb) then     
!---                Tb is good  (doing IF statement this way identifies NaN as bad)     
           else 	     
!	    write(6,*)'badtbs', nf,scan,pix,dprsim_tb(nf,pix,scan,2)	     
                    badchan(nf) = badchan(nf) + 1
		    retr_flag(pix,scan) = 0
           endif
	 enddo

       enddo
  enddo
  write(*,'(a,<outchans>I6)')'    number unexpected bad Tbs:',badchan(1:outchans)

!---  calculate footprint averages for all hydrometeors at the selected frequency as
!---  defined in variable 'fpfov'. Also creates footprint ave Tbs at the sensor frequencies.

!---  loop over all frequencies for this sensor

  do nf = 1, outchans

    n_xt = nint(fov_crosstrack(nf) / fov_dpr / 2.)        !define num of pixs in ftprint
    n_dt = nint(fov_downtrack(nf) / fov_dpr / 2.)         !   in both directions

!--- 1st good freq has largest footprint. Therefore limit pixels in scan to having all
!--- channels at the full footprint  (generally eliminates DPR center pixel 1,2,24,25 
    
    if(nf .eq. 1) then
       pixstart = 1  + nint(fov_crosstrack(nf) / fov_dpr / 2.)   !only calc footprints
       pixend   = 25 - nint(fov_crosstrack(nf) / fov_dpr / 2.)   !for these pixels

       scanstart = 1  + nint(fov_downtrack(nf) / fov_dpr / 2.)	!only calc footprints
       scanend   = nscans - nint(fov_downtrack(nf) / fov_dpr / 2.)  !for these pixels       
       write(*,'(a,7i5)')'  selected ftprnt in dpr pixels: nf,start,end:',    &
                          nf, n_xt, n_dt, pixstart, pixend, scanstart,scanend
    endif

!--- loop over all DPR pixels, sum up and keep counter within the footprint

    do pix = pixstart,pixend	     !  middle of the scan which has PR data filling all fpa
      do scan = scanstart,scanend    !  eliminates a few  scans/pixs at the beginning and
    				     !  end that are not full		
        sum_wgt = 0.    	
        sum_wgt_mirs = 0.
	sum_wgt_gTbsobs = 0.
	sum_wgt_cTbssim = 0.
	sum_wgt_mTbssim = 0.
	sum_wgt_dTbssim = 0.
	sum_datasource=0.; sum_emiss=0.;sum_csfcprcp=0.     !surface variables
        sum_msfcprcp=0.; sum_ccnvprcp=0.		  
        sum_cTbssim=0.;  sum_mTbssim=0.                      !Tbs variables
        sum_dTbssim=0. ; sum_deltadTbssim=0.
	sum_gTbsobs=0.
        sum_rainwater=0.; sum_snowwater=0.;                 !layered variables
        sum_cloudwater=0.; sum_SLH=0.
	
    	do i = pix-n_xt, pix+n_xt
    	  if(i .lt. 1 .or. i .gt. npixs) then
!    	     write(*,*)' skipping i: ', i
	     cycle
	  endif
	   
	  do j = scan-n_dt, scan+n_dt
	    if(j .lt. 1 .or. j .gt. nscans) then
!	       write(*,*)' skipping j: ', j
	       cycle
	    endif
	       
	    if(retr_flag(i,j) .le. 0) cycle          !only do when good Combined pixels
	       
!---       calculate the gaussian pixel distance weight
	       
	    xa = fov_dpr * (i-pix)                          
            ya = fov_dpr * (j-scan)
            p = (xa**2 / fov_crosstrack(nf)**2) + (ya**2) / (fov_downtrack(nf)**2)
            wgt = exp(-2.7726 * p)            
	    sum_wgt = sum_wgt + wgt
	       
!---       calculate the total weight = sum(parameters * distance weight)
            
            sum_datasource = sum_datasource + (wgt * datasource(i,j))
	    sum_csfcprcp   = sum_csfcprcp   + (wgt * csfcprcp(i,j))
            sum_ccnvprcp   = sum_ccnvprcp    + (wgt * ccnvprcp(i,j))

	    if(msfcprcp(i,j) .ge. 0.0) then
	        sum_msfcprcp   = sum_msfcprcp   + (wgt * msfcprcp(i,j))
                sum_wgt_mirs = sum_wgt_mirs + wgt
            endif
	    
	    sum_emiss(:) = sum_emiss(:) + (wgt * emiss_strm(:,i,j))          !emissivity
		 
            sum_rainwater(:) = sum_rainwater(:) + (wgt * rainwater(:,i,j))       !layered variables                             	       
            sum_cloudwater(:)= sum_cloudwater(:)+ (wgt * cloudwater(:,i,j))                 
            sum_snowwater(:) = sum_snowwater(:) + (wgt * snowwater(:,i,j))	 
            sum_SLH(:)       = sum_SLH(:)       + (wgt * SLH(:,i,j))             		     
	    
	    if(cTbssim(nf,i,j) .ge. minTb .and. cTbssim(nf,i,j).le. maxTb) then
	        sum_cTbssim     = sum_cTbssim     + (wgt * cTbssim(nf,i,j))      !CMB simulated Tbs
		sum_wgt_cTbssim = sum_wgt_cTbssim + wgt
            endif 
	    if(mTbssim(nf,i,j) .ge. minTb .and. mTbssim(nf,i,j).le. maxTb) then
	        sum_mTbssim     = sum_mTbssim     + (wgt * mTbssim(nf,i,j))      !MIRS simulated Tbs
		sum_wgt_mTbssim = sum_wgt_mTbssim + wgt
	    endif
	    if(dTbssim(nf,i,j) .ge. minTb .and. dTbssim(nf,i,j).le. maxTb) then
	        sum_dTbssim     = sum_dTbssim     + (wgt * dTbssim(nf,i,j))      !DR simulated Tbs
	        sum_wgt_dTbssim = sum_wgt_dTbssim + wgt
	    endif
!	    if(gTbsobs(nf,i,j) .ge. minTb .and. gTbsobs(nf,i,j).le. maxTb) then
!	        sum_gTbsobs     = sum_gTbsobs     + (wgt * gTbsobs(nf,i,j))      !GMI observed Tbs
!                sum_wgt_gTbsobs = sum_wgt_gTbsobs + wgt
!            endif
	    
	    sum_deltadTbssim = sum_deltadTbssim + (wgt * deltadTbssim(nf,i,j)) !DeltaTbs 	
     
	  enddo    !scan-n_dt
	enddo     !pix-n_nt


!---   calculate footprint averages for all hydrometeors if freq = selected freq (fpfov)
  
        if(nf .eq. fpfov) then	
            if(sum_wgt_mirs .gt. 0.0) then
	        msfcprcp_fpa(pix,scan)   = sum_msfcprcp   / sum_wgt_mirs
            else
	        msfcprcp_fpa(pix,scan) = dflt
            endif

            if(sum_wgt .gt. 0.0) then          
                datasource_fpa(pix,scan) = sum_datasource / sum_wgt            !surface variables
      	        csfcprcp_fpa(pix,scan)   = sum_csfcprcp   / sum_wgt
	        ccnvprcp_fpa(pix,scan)   = sum_ccnvprcp   / sum_wgt
	        emiss_fpa(:,pix,scan)    = sum_emiss(:)   / sum_wgt
	      
	        rainwater_fpa(:,pix,scan)  = sum_rainwater(:)  / sum_wgt       !layer variables
	        cloudwater_fpa(:,pix,scan) = sum_cloudwater(:)  / sum_wgt
	        snowwater_fpa(:,pix,scan)  = sum_snowwater(:)  / sum_wgt
	        SLH_fpa(:,pix,scan)        = sum_SLH(:)        / sum_wgt				           
	    else
                datasource_fpa(pix,scan) = dflt                                !surface variables
	        csfcprcp_fpa(pix,scan)   = dflt 
	        ccnvprcp_fpa(pix,scan)    = dflt
	        emiss_fpa(:,pix,scan)    = dflt
	      
	        rainwater_fpa(:,pix,scan)   = dflt	                       !layer variables
	        cloudwater_fpa(:,pix,scan)  = dflt 
	        snowwater_fpa(:,pix,scan)   = dflt 
	        SLH_fpa(:,pix,scan)         = dflt 
            endif
        endif
	
!---   calculate Tb, and deltaTb frequencies footprint averages for each native footprint size

        if(sum_wgt .gt. 0.0 .and. sum_cTbssim .gt. 0.0) then
           cTbssim_fpa(nf,pix,scan) = sum_cTbssim / sum_wgt_cTbssim     !Combined Simulated Tbs
	else
	   cTbssim_fpa(nf,pix,scan) = dflt
	endif	    
        if(sum_wgt .gt. 0.0 .and. sum_mTbssim .gt. 0.0) then            !mirs Simulated Tbs
           mTbssim_fpa(nf,pix,scan) = sum_mTbssim / sum_wgt_mTbssim
	else
	   mTbssim_fpa(nf,pix,scan) = dflt
	endif  
        if(sum_wgt .gt. 0.0 .and. sum_dTbssim .gt. 0.0) then            !GPROF Simulated Tbs
           dTbssim_fpa(nf,pix,scan) = sum_dTbssim / sum_wgt_dTbssim
	else
	   dTbssim_fpa(nf,pix,scan) = dflt
	endif  
!        if(sum_wgt_gTbsobs .gt. 0.0 .and. sum_gTbsobs .gt. 0.0) then    !GMI obs Tbs
!           gTbsobs_fpa(nf,pix,scan) = sum_gTbsobs / sum_wgt_gTbsobs
!	else
!	   gTbsobs_fpa(nf,pix,scan) = dflt
!	endif
	
        if(sum_wgt .gt. 0.0) then
	   deltadTbssim_fpa(nf,pix,scan) = sum_deltadTbssim / sum_wgt   !GPROF Delta Tbs
	else
           deltadTbssim_fpa(nf,pix,scan) = dflt       
	endif
                  
	   	   
      enddo   !scan
    enddo   !pix
	 
  enddo   !nf (frequency)

  write(6,*)' starting averaging'


!--- calculate average of all the Tbs other parameteres for all the pixels in the orbit
      
  sum_datasource_orb=0.0; sum_emiss_orb=0.0 		 !surface variables
  sum_csfcprcp_orb=0.0;   sum_msfcprcp_orb=0.0; sum_ccnvprcp_orb=0.0
  sum_cTbssim_orb=0.0;    sum_mTbssim_orb=0.0		 !Tbs variables
  sum_dTbssim_orb=0.0;    sum_deltadTbssim_orb=0.0
  sum_gTbsobs_orb = 0
  sum_rainwater_orb = 0; sum_snowwater_orb = 0.0      !layered variables
  sum_cloudwater_orb= 0; sum_SLH_orb = 0.0

  cnt_datasource_orb=0.0; cnt_emiss_orb=0.0 		 !surface variables
  cnt_csfcprcp_orb=0.0;   cnt_msfcprcp_orb=0.0; cnt_ccnvprcp_orb=0.0
  cnt_cTbssim_orb=0.0;    cnt_mTbssim_orb=0.0		 !Tbs variables
  cnt_dTbssim_orb=0.0;    cnt_deltadTbssim_orb=0.0
  cnt_gTbsobs_orb = 0
  cnt_rainwater_orb = 0; cnt_snowwater_orb = 0.0      !layered variables
  cnt_cloudwater_orb= 0; cnt_SLH_orb = 0.0

  do nf = 1, outchans
    do scan = scanstart, scanend
      do pix = pixstart, pixend
   	      
   	if(retr_flag(pix,scan) .le. 0)cycle
  
   	icnt(1) = icnt(1) + 1

   	if(gTbsobs(nf,pix,scan) .gt. 0.0) then   !only use the GMI pixel location
  	    icnt(2) = icnt(2) + 1 

            if(gTbsobs(nf,pix,scan).ge.minTb .and. gTbsobs(nf,pix,scan).le.maxTb) then	    
	        sum_gTbsobs_orb(nf) = sum_gTbsobs_orb(nf) + gTbsobs(nf,pix,scan)
	        cnt_gTbsobs_orb(nf) = cnt_gTbsobs_orb(nf) + 1
	    endif
	    if(dTbssim_fpa(nf,pix,scan).ge.minTb .and. dTbssim_fpa(nf,pix,scan).le.maxTb) then	    
	        sum_dTbssim_orb(nf) = sum_dTbssim_orb(nf) + dTbssim_fpa(nf,pix,scan)
	        cnt_dTbssim_orb(nf) = cnt_dTbssim_orb(nf) + 1    
	    endif
   	    if(cTbssim_fpa(nf,pix,scan).ge.minTb .and. cTbssim_fpa(nf,pix,scan).le.maxTb) then
	        sum_cTbssim_orb(nf) = sum_cTbssim_orb(nf) + cTbssim_fpa(nf,pix,scan)
                cnt_cTbssim_orb(nf) = cnt_cTbssim_orb(nf) + 1
   	    endif
	    if(mTbssim_fpa(nf,pix,scan).ge.minTb .and. mTbssim_fpa(nf,pix,scan).le.maxTb) then
	        sum_mTbssim_orb(nf) = sum_mTbssim_orb(nf) + mTbssim_fpa(nf,pix,scan)
                cnt_mTbssim_orb(nf) = cnt_mTbssim_orb(nf) + 1
	    endif
	    if(deltadTbssim_fpa(nf,pix,scan) .ge. minTb .and.   &
	       deltadTbssim_fpa(nf,pix,scan) .le. maxTb) then	    
	        sum_deltadTbssim_orb(nf) = sum_deltadTbssim_orb(nf) +  &
	                               deltadTbssim_fpa(nf,pix,scan)
                cnt_deltadTbssim_orb(nf) = cnt_deltadTbssim_orb(nf) + 1
	    endif		       
                

	    if(emiss_fpa(nf,pix,scan) .gt. 0.0) then
	        sum_emiss_orb(nf)  = sum_emiss_orb(nf)  + emiss_fpa(nf,pix,scan)
	    	cnt_emiss_orb(nf)  = cnt_emiss_orb(nf)  + 1
            endif

		  
            if(nf .eq. fpfov) then                !do these only one time

	        if(csfcprcp_fpa(pix,scan) .ge. 0.0) then
	            sum_csfcprcp_orb = sum_csfcprcp_orb + csfcprcp_fpa(pix,scan)
	            cnt_csfcprcp_orb = cnt_csfcprcp_orb + 1
	        endif
	        if(msfcprcp_fpa(pix,scan) .ge. 0.0) then
	            sum_msfcprcp_orb = sum_msfcprcp_orb + msfcprcp_fpa(pix,scan)
	            cnt_msfcprcp_orb = cnt_msfcprcp_orb + 1
	        endif
	        if(ccnvprcp_fpa(pix,scan) .ge. 0.0) then
	            sum_ccnvprcp_orb = sum_ccnvprcp_orb + ccnvprcp_fpa(pix,scan)
	            cnt_ccnvprcp_orb = cnt_ccnvprcp_orb + 1
	        endif		
		if(datasource_fpa(pix,scan) .ge. 1) then
	            sum_datasource_orb = sum_datasource_orb + datasource_fpa(pix,scan)
	            cnt_datasource_orb = cnt_datasource_orb + 1
	        endif
		
		do nl = 1, nlevs
		  if(rainwater_fpa(nl,pix,scan) .ge. 0.0) then
		     sum_rainwater_orb(nl)= sum_rainwater_orb(nl) + rainwater_fpa(nl,pix,scan)
		     cnt_rainwater_orb(nl)= cnt_rainwater_orb(nl) + 1
		  endif		
		  if(cloudwater_fpa(nl,pix,scan) .ge. 0.0) then
		     sum_cloudwater_orb(nl)= sum_cloudwater_orb(nl) + cloudwater_fpa(nl,pix,scan)
		     cnt_cloudwater_orb(nl)= cnt_cloudwater_orb(nl) + 1
	  	  endif		
		  if(snowwater_fpa(nl,pix,scan) .ge. 0.0) then
		     sum_snowwater_orb(nl)= sum_snowwater_orb(nl) + snowwater_fpa(nl,pix,scan)
		     cnt_snowwater_orb(nl)= cnt_snowwater_orb(nl) + 1
		  endif		
		  if(rainwater_fpa(nl,pix,scan) .ge. 0.0) then
		     sum_SLH_orb(nl)= sum_SLH_orb(nl) + SLH_fpa(nl,pix,scan)
		     cnt_SLH_orb(nl)= cnt_SLH_orb(nl) + 1
		  endif
		enddo

            endif   !nf = fpfov
		  	  
	endif   !gTbsobs > 0.0
     
      enddo !pix
    enddo  !scan
  enddo   !nf
  
  
  write(6,*)' number pix = ', icnt(1)
  write(6,*)' number pix good gTBS = ',icnt(2)
	     	
!---- calc overall orbit averages

  write(6,*)' calc orbit averages'

!--- initial average values to missing

  ave_datasource_orb = dflt; ave_csfcprcp_orb=dflt; ave_msfcprcp_orb=dflt
  ave_ccnvprcp_orb = dflt; ave_emiss_orb=dflt; ave_cTbssim_orb=dflt
  ave_dTbssim_orb=dflt; ave_mTbssim_orb=dflt;ave_deltadTbssim_orb=dflt
  ave_rainwater_orb=dflt;ave_cloudwater_orb=dflt;ave_snowwater_orb=dflt
  ave_SLH_orb=dflt

  if(cnt_datasource_orb .gt. 0.0) then                                     !surface variables
      ave_datasource_orb = sum_datasource_orb / cnt_datasource_orb
  endif  
  if(cnt_csfcprcp_orb .gt. 0.0) then
       ave_csfcprcp_orb   = sum_csfcprcp_orb  / cnt_csfcprcp_orb
  endif
  if(cnt_msfcprcp_orb .gt. 0.0) then
       ave_msfcprcp_orb   = sum_msfcprcp_orb  / cnt_msfcprcp_orb
  endif
  if(cnt_ccnvprcp_orb .gt. 0.0) then    
       ave_ccnvprcp_orb    = sum_ccnvprcp_orb   / cnt_ccnvprcp_orb  
  endif
  
  do nf = 1, outchans
    if(cnt_emiss_orb(nf) .gt. 0.0) then
        ave_emiss_orb(nf) = sum_emiss_orb(nf) / cnt_emiss_orb(nf)         !Emiss
    endif
    if(cnt_cTbssim_orb(nf) .gt. 0.0) then
        ave_cTbssim_orb(nf) = sum_cTbssim_orb(nf) / cnt_cTbssim_orb(nf)   !cmb Tbsim
    endif
    if(cnt_mTbssim_orb(nf) .gt. 0.0) then
        ave_mTbssim_orb(nf) = sum_mTbssim_orb(nf) / cnt_mTbssim_orb(nf)   !mirs Tbsim
    endif
    if(cnt_dTbssim_orb(nf) .gt. 0.0) then
        ave_dTbssim_orb(nf) = sum_dTbssim_orb(nf) / cnt_dTbssim_orb(nf)   !GPROF Tbsim
    endif
    if(cnt_deltadTbssim_orb(nf) .gt. 0.0) then                             !GPROF Delta Tbsim
        ave_deltadTbssim_orb(nf) = sum_deltadTbssim_orb(nf) / cnt_deltadTbssim_orb(nf)
    endif
    if(cnt_gTbsobs_orb(nf) .gt. 0.0) then                                  !GMI Obs Tbs
        ave_gTbsobs_orb(nf) = sum_gTbsobs_orb(nf) / cnt_gTbsobs_orb(nf)
    endif
  enddo
  
  do nl = 1, nlevs
    if(cnt_rainwater_orb(nl) .gt. 0.0) then
       ave_rainwater_orb(nl)  = sum_rainwater_orb(nl) / cnt_rainwater_orb(nl)
    endif  
    if(cnt_cloudwater_orb(nl) .gt. 0.0) then
       ave_cloudwater_orb(nl) = sum_cloudwater_orb(nl) / cnt_cloudwater_orb(nl)
    endif
    if(cnt_snowwater_orb(nl) .gt. 0.0) then
       ave_snowwater_orb(nl) = sum_snowwater_orb(nl) / cnt_snowwater_orb(nl)
    endif
    if(cnt_SLH_orb(nl) .gt. 0.0) then
       ave_SLH_orb(nl) = sum_SLH_orb(nl) / cnt_SLH_orb(nl)
    endif
  enddo
  
!--- print out the results

  write(6,*)
  write(*,'(9x,a5,a,a5,a)')satcode,'  ',sensor, out_label 
  write(*,'(a,a,a)')'------------------------------------------------',   &
                    '------------------------------------------------',   &
                    '-------------------------------------------------'
!--- list tb averages

  write(*,*)' MIRS,CMB, GPROF simulated Tbs on GMI footprints'       
  write(*,'(5x,a,5F8.2,8x,6F8.2,8x,2F8.2)')'GMI   obsTbs   :   ',ave_gTbsobs_orb(:)
  write(*,'(5x,a,5F8.2,8x,6F8.2,8x,2F8.2)')'GPROF simTbs   :   ',ave_dTbssim_orb(:)
  write(*,'(5x,a,5F8.2,8x,6F8.2,8x,2F8.2)')'CMB   simTbs   :   ',ave_cTbssim_orb(:)
  write(*,'(5x,a,5F8.2,8x,6F8.2,8x,2F8.2)')'MIRS  simTbs   :   ',ave_mTbssim_orb(:)
  write(*,*)
  write(*,'(5x,a,5F8.2,8x,6F8.2,8x,2F8.2)')'GPROF DeltaTbs :   ',ave_deltadTbssim_orb(:)
  write(*,*)
  write(*,'(a,13F8.2)')  '     Frequencies out  : ',gmifreqs(:)
  write(*,'(a,13F8.2)')  '        Nominal EIAs  : ',gmieia(:)
  write(*,'(5x,a,13F8.2)')'   Blended emiss : ',ave_emiss_orb(:)   
  write(*,*)
  write(*,*)' AVE PRCP VARIABLES   (mm/hr),  (mm/day): '
  write(*,'(a,2F12.2)')'     CMB precip : ',ave_csfcprcp_orb,ave_csfcprcp_orb*24
  write(*,'(a,2F12.2)')'    MIRS precip : ',ave_msfcprcp_orb,ave_msfcprcp_orb*24
  write(*,'(a,2F12.2)')' CMB CNV precip : ',ave_ccnvprcp_orb,ave_ccnvprcp_orb*24
  write(*,'(a,F12.2)') ' Average source : ',ave_datasource_orb
  write(*,*)  
  write(*,*)'LAYERED PRODUCT AVES  RAIN     CLOUD     SNOW      SLH'
  write(*,*)'-------------------------------------------------------------'
  do nl = nlevs,1,-1
    write(*,'(15x,i3,4F10.3)') nl, ave_rainwater_orb(nl),ave_cloudwater_orb(nl),   &
                      ave_snowwater_orb(nl),ave_SLH_orb(nl)
  enddo

 end subroutine GPM_compute_fpa

!------------------------------------------------------------------------

 subroutine GPM_compute_aves
  
 !--- this routine will take the combined/mirs Tb retrieval and calculate
 !--- the averages of all parameters at the DPR/CMB pixel resolution
 
 
 !--- general subroutine variables
 
  integer,parameter :: outchans = 13
  real,parameter    :: minTb = 45, maxTb = 340
  integer,parameter :: fpfov = 3
  integer           :: icnt(5) = 0
  real,parameter    :: fov_dpr = 5.0  
  real              :: xa, ya, p, wgt      
 
  integer :: badchan(outchans), pix,scan,i,j,k, n_xt,n_dt, nf, nc, nl
  integer :: nx,ny
  logical :: firstpixres = .true.
  real    :: sum_wgt, sum_wgt_mirs, sum_wgt_gTbsobs, sum_wgt_cTbssim
  real    :: sum_wgt_mTbssim, sum_wgt_dTbssim

!--- results printing
  
  character(len=128) :: out_label='        10V     10H     19V     19H   23.8V      '//     &
            '     36.5V   36.5H     89V     89H    166V    166H' //   &
            '         183+/-3 183+/-7'

!-- diagnostic variables  - sums, cnts, and averages of all pixels in single orbit

  real(8) :: sum_datasource_orb, sum_emiss_orb(13)                  !surface variables
  real(8) :: sum_csfcprcp_orb, sum_msfcprcp_orb, sum_ccnvprcp_orb
  real(8) :: sum_cTbssim_orb(13), sum_mTbssim_orb(13)               !Tbs variables
  real(8) :: sum_dTbssim_orb(13), sum_deltadTbssim_orb(13)
  real(8) :: sum_gTbsobs_orb(13)
  real(8) :: sum_rainwater_orb(nlevs)		          	    !layered variables
  real(8) :: sum_snowwater_orb(nlevs)
  real(8) :: sum_cloudwater_orb(nlevs)
  real(8) :: sum_SLH_orb(nlevs)

  real(8) :: cnt_datasource_orb, cnt_emiss_orb(13)                  !surface variables
  real(8) :: cnt_csfcprcp_orb, cnt_msfcprcp_orb, cnt_ccnvprcp_orb
  real(8) :: cnt_cTbssim_orb(13), cnt_mTbssim_orb(13)               !Tbs variables
  real(8) :: cnt_dTbssim_orb(13), cnt_deltadTbssim_orb(13)
  real(8) :: cnt_gTbsobs_orb(13)
  real(8) :: cnt_rainwater_orb(nlevs)		          	    !layered variables
  real(8) :: cnt_snowwater_orb(nlevs)
  real(8) :: cnt_cloudwater_orb(nlevs)
  real(8) :: cnt_SLH_orb(nlevs)

  real :: ave_datasource_orb, ave_emiss_orb(13) 		    !surface variables
  real :: ave_csfcprcp_orb, ave_msfcprcp_orb, ave_ccnvprcp_orb
  real :: ave_cTbssim_orb(13), ave_mTbssim_orb(13) 		    !Tbs variables
  real :: ave_dTbssim_orb(13), ave_deltadTbssim_orb(13)
  real :: ave_gTbsobs_orb(13)
  real :: ave_rainwater_orb(nlevs)				    !layered variables
  real :: ave_snowwater_orb(nlevs)
  real :: ave_cloudwater_orb(nlevs)
  real :: ave_SLH_orb(nlevs)
  
!--- Error check the dr simulated Tbs
    
  do nf = 1, outchans
  
    if(gmifreqs(nf) .eq. dflt) cycle
         
       do scan = 1,nscans
         do pix = 1,npixs	   
	   if(dTbssim(nf,pix,scan) .ge. mintb .and. dTbssim(nf,pix,scan) .le. maxtb) then     
!---                Tb is good  (doing IF statement this way identifies NaN as bad)     
           else 	     
!	    write(6,*)'badtbs', nf,scan,pix,dprsim_tb(nf,pix,scan,2)	     
                    badchan(nf) = badchan(nf) + 1
		    retr_flag(pix,scan) = 0
           endif
	 enddo

       enddo
  enddo
  write(*,'(a,<outchans>I6)')'    number unexpected bad Tbs:',badchan(1:outchans)

!--- define pixel and scan start and end
   
  pixstart  = 1
  pixend    = npixs
  scanstart = 1
  scanend   = nscans

!---  calculate averages for all hydrometeors at the selected frequency as
!---  defined in variable 'fpfov'. Also creates Tbs at the sensor frequencies.

!--- calculate average of all the Tbs other parameteres for all the pixels in the orbit
      
  sum_datasource_orb=0.0; sum_emiss_orb=0.0 		 !surface variables
  sum_csfcprcp_orb=0.0;   sum_msfcprcp_orb=0.0; sum_ccnvprcp_orb=0.0
  sum_cTbssim_orb=0.0;    sum_mTbssim_orb=0.0		 !Tbs variables
  sum_dTbssim_orb=0.0;    sum_deltadTbssim_orb=0.0
  sum_gTbsobs_orb = 0
  sum_rainwater_orb = 0; sum_snowwater_orb = 0.0      !layered variables
  sum_cloudwater_orb= 0; sum_SLH_orb = 0.0

  cnt_datasource_orb=0.0; cnt_emiss_orb=0.0 		 !surface variables
  cnt_csfcprcp_orb=0.0;   cnt_msfcprcp_orb=0.0; cnt_ccnvprcp_orb=0.0
  cnt_cTbssim_orb=0.0;    cnt_mTbssim_orb=0.0		 !Tbs variables
  cnt_dTbssim_orb=0.0;    cnt_deltadTbssim_orb=0.0
  cnt_gTbsobs_orb = 0
  cnt_rainwater_orb = 0; cnt_snowwater_orb = 0.0      !layered variables
  cnt_cloudwater_orb= 0; cnt_SLH_orb = 0.0

  do nf = 1, outchans
    do scan = scanstart, scanend
      do pix = pixstart, pixend
   	      
   	if(retr_flag(pix,scan) .le. 0)cycle
  
   	icnt(1) = icnt(1) + 1

        if(gTbsobs(nf,pix,scan).ge.minTb .and. gTbsobs(nf,pix,scan).le.maxTb) then	
	    sum_gTbsobs_orb(nf) = sum_gTbsobs_orb(nf) + gTbsobs(nf,pix,scan)
	    cnt_gTbsobs_orb(nf) = cnt_gTbsobs_orb(nf) + 1
	endif
	if(dTbssim(nf,pix,scan).ge.minTb .and. dTbssim(nf,pix,scan).le.maxTb) then	
	    sum_dTbssim_orb(nf) = sum_dTbssim_orb(nf) + dTbssim(nf,pix,scan)
	    cnt_dTbssim_orb(nf) = cnt_dTbssim_orb(nf) + 1    
	endif
   	if(cTbssim(nf,pix,scan).ge.minTb .and. cTbssim(nf,pix,scan).le.maxTb) then
	    sum_cTbssim_orb(nf) = sum_cTbssim_orb(nf) + cTbssim(nf,pix,scan)
            cnt_cTbssim_orb(nf) = cnt_cTbssim_orb(nf) + 1
   	endif
	if(mTbssim(nf,pix,scan).ge.minTb .and. mTbssim(nf,pix,scan).le.maxTb) then
	    sum_mTbssim_orb(nf) = sum_mTbssim_orb(nf) + mTbssim(nf,pix,scan)
            cnt_mTbssim_orb(nf) = cnt_mTbssim_orb(nf) + 1
	endif
	if(deltadTbssim(nf,pix,scan) .ge. minTb .and.   &
	   deltadTbssim(nf,pix,scan) .le. maxTb) then	
	    sum_deltadTbssim_orb(nf) = sum_deltadTbssim_orb(nf) +  &
				   deltadTbssim(nf,pix,scan)
            cnt_deltadTbssim_orb(nf) = cnt_deltadTbssim_orb(nf) + 1
	endif			   
            

	if(emiss_strm(nf,pix,scan) .gt. 0.0) then
	    sum_emiss_orb(nf)  = sum_emiss_orb(nf)  + emiss_strm(nf,pix,scan)
	    cnt_emiss_orb(nf)  = cnt_emiss_orb(nf)  + 1
        endif

	      
        if(nf .eq. fpfov) then  	      !do these only one time

	    if(csfcprcp(pix,scan) .ge. 0.0) then
		sum_csfcprcp_orb = sum_csfcprcp_orb + csfcprcp(pix,scan)
		cnt_csfcprcp_orb = cnt_csfcprcp_orb + 1
	    endif
	    if(msfcprcp(pix,scan) .ge. 0.0) then
		sum_msfcprcp_orb = sum_msfcprcp_orb + msfcprcp(pix,scan)
		cnt_msfcprcp_orb = cnt_msfcprcp_orb + 1
	    endif
	    if(ccnvprcp(pix,scan) .ge. 0.0) then
		sum_ccnvprcp_orb = sum_ccnvprcp_orb + ccnvprcp(pix,scan)
		cnt_ccnvprcp_orb = cnt_ccnvprcp_orb + 1
	    endif	    
	    if(datasource(pix,scan) .ge. 1) then
		sum_datasource_orb = sum_datasource_orb + datasource(pix,scan)
		cnt_datasource_orb = cnt_datasource_orb + 1
	    endif
	
	    do nl = 1, nlevs
	      if(rainwater(nl,pix,scan) .ge. 0.0) then
	    	 sum_rainwater_orb(nl)= sum_rainwater_orb(nl) + rainwater(nl,pix,scan)
	    	 cnt_rainwater_orb(nl)= cnt_rainwater_orb(nl) + 1
	      endif	    
	      if(cloudwater(nl,pix,scan) .ge. 0.0) then
	    	 sum_cloudwater_orb(nl)= sum_cloudwater_orb(nl) + cloudwater(nl,pix,scan)
	    	 cnt_cloudwater_orb(nl)= cnt_cloudwater_orb(nl) + 1
	      endif	    
	      if(snowwater(nl,pix,scan) .ge. 0.0) then
	    	 sum_snowwater_orb(nl)= sum_snowwater_orb(nl) + snowwater(nl,pix,scan)
	    	 cnt_snowwater_orb(nl)= cnt_snowwater_orb(nl) + 1
	      endif	    
	      if(rainwater(nl,pix,scan) .ge. 0.0) then
	    	 sum_SLH_orb(nl)= sum_SLH_orb(nl) + SLH(nl,pix,scan)
	    	 cnt_SLH_orb(nl)= cnt_SLH_orb(nl) + 1
	      endif
	    enddo

        endif	!nf = fpfov
	
      enddo !pix
    enddo  !scan
  enddo   !nf
  
  
  write(6,*)' number pix = ', icnt(1)
  write(6,*)' number pix good gTBS = ',icnt(2)
	     	
!---- calc overall orbit averages

  write(6,*)' calc orbit averages'

!--- initial average values to missing

  ave_datasource_orb = dflt; ave_csfcprcp_orb=dflt; ave_msfcprcp_orb=dflt
  ave_ccnvprcp_orb = dflt; ave_emiss_orb=dflt; ave_cTbssim_orb=dflt
  ave_dTbssim_orb=dflt; ave_mTbssim_orb=dflt;ave_deltadTbssim_orb=dflt
  ave_rainwater_orb=dflt;ave_cloudwater_orb=dflt;ave_snowwater_orb=dflt
  ave_SLH_orb=dflt

  if(cnt_datasource_orb .gt. 0.0) then                                     !surface variables
      ave_datasource_orb = sum_datasource_orb / cnt_datasource_orb
  endif  
  if(cnt_csfcprcp_orb .gt. 0.0) then
       ave_csfcprcp_orb   = sum_csfcprcp_orb  / cnt_csfcprcp_orb
  endif
  if(cnt_msfcprcp_orb .gt. 0.0) then
       ave_msfcprcp_orb   = sum_msfcprcp_orb  / cnt_msfcprcp_orb
  endif
  if(cnt_ccnvprcp_orb .gt. 0.0) then    
       ave_ccnvprcp_orb    = sum_ccnvprcp_orb   / cnt_ccnvprcp_orb  
  endif
  
  do nf = 1, outchans
    if(cnt_emiss_orb(nf) .gt. 0.0) then
        ave_emiss_orb(nf) = sum_emiss_orb(nf) / cnt_emiss_orb(nf)         !Emiss
    endif
    if(cnt_cTbssim_orb(nf) .gt. 0.0) then
        ave_cTbssim_orb(nf) = sum_cTbssim_orb(nf) / cnt_cTbssim_orb(nf)   !cmb Tbsim
    endif
    if(cnt_mTbssim_orb(nf) .gt. 0.0) then
        ave_mTbssim_orb(nf) = sum_mTbssim_orb(nf) / cnt_mTbssim_orb(nf)   !mirs Tbsim
    endif
    if(cnt_dTbssim_orb(nf) .gt. 0.0) then
        ave_dTbssim_orb(nf) = sum_dTbssim_orb(nf) / cnt_dTbssim_orb(nf)   !GPROF Tbsim
    endif
    if(cnt_deltadTbssim_orb(nf) .gt. 0.0) then                             !GPROF Delta Tbsim
        ave_deltadTbssim_orb(nf) = sum_deltadTbssim_orb(nf) / cnt_deltadTbssim_orb(nf)
    endif
    if(cnt_gTbsobs_orb(nf) .gt. 0.0) then                                  !GMI Obs Tbs
        ave_gTbsobs_orb(nf) = sum_gTbsobs_orb(nf) / cnt_gTbsobs_orb(nf)
    endif
  enddo
  
  do nl = 1, nlevs
    if(cnt_rainwater_orb(nl) .gt. 0.0) then
       ave_rainwater_orb(nl)  = sum_rainwater_orb(nl) / cnt_rainwater_orb(nl)
    endif  
    if(cnt_cloudwater_orb(nl) .gt. 0.0) then
       ave_cloudwater_orb(nl) = sum_cloudwater_orb(nl) / cnt_cloudwater_orb(nl)
    endif
    if(cnt_snowwater_orb(nl) .gt. 0.0) then
       ave_snowwater_orb(nl) = sum_snowwater_orb(nl) / cnt_snowwater_orb(nl)
    endif
    if(cnt_SLH_orb(nl) .gt. 0.0) then
       ave_SLH_orb(nl) = sum_SLH_orb(nl) / cnt_SLH_orb(nl)
    endif
  enddo
  
!--- print out the results

  write(6,*)
  write(*,'(9x,a5,a,a5,a)')satcode,'  ',sensor, out_label 
  write(*,'(a,a,a)')'------------------------------------------------',   &
                    '------------------------------------------------',   &
                    '-------------------------------------------------'
!--- list tb averages

  write(*,*)' MIRS,CMB, GPROF simulated Tbs - original DPR pixels'       
  write(*,'(5x,a,5F8.2,8x,6F8.2,8x,2F8.2)')'GMI   obsTbs   :   ',ave_gTbsobs_orb(:)
  write(*,'(5x,a,5F8.2,8x,6F8.2,8x,2F8.2)')'GPROF simTbs   :   ',ave_dTbssim_orb(:)
  write(*,'(5x,a,5F8.2,8x,6F8.2,8x,2F8.2)')'CMB   simTbs   :   ',ave_cTbssim_orb(:)
  write(*,'(5x,a,5F8.2,8x,6F8.2,8x,2F8.2)')'MIRS  simTbs   :   ',ave_mTbssim_orb(:)
  write(*,*)
  write(*,'(5x,a,5F8.2,8x,6F8.2,8x,2F8.2)')'GPROF DeltaTbs :   ',ave_deltadTbssim_orb(:)
  write(*,*)
  write(*,'(a,13F8.2)')  '     Frequencies out  : ',gmifreqs(:)
  write(*,'(a,13F8.2)')  '        Nominal EIAs  : ',gmieia(:)
  write(*,'(5x,a,13F8.2)')'   Blended emiss : ',ave_emiss_orb(:)   
  write(*,*)
  write(*,*)' AVE PRCP VARIABLES   (mm/hr),  (mm/day): '
  write(*,'(a,2F12.2)')'     CMB precip : ',ave_csfcprcp_orb,ave_csfcprcp_orb*24
  write(*,'(a,2F12.2)')'    MIRS precip : ',ave_msfcprcp_orb,ave_msfcprcp_orb*24
  write(*,'(a,2F12.2)')' CMB CNV precip : ',ave_ccnvprcp_orb,ave_ccnvprcp_orb*24
  write(*,'(a,F12.2)') ' Average source : ',ave_datasource_orb
  write(*,*)  
  write(*,*)'LAYERED PRODUCT AVES  RAIN     CLOUD     SNOW      SLH'
  write(*,*)'-------------------------------------------------------------'
  do nl = nlevs,1,-1
    write(*,'(15x,i3,4F10.3)') nl, ave_rainwater_orb(nl),ave_cloudwater_orb(nl),   &
                      ave_snowwater_orb(nl),ave_SLH_orb(nl)
  enddo

 end subroutine GPM_compute_aves

!------------------------------------------------------------------------
 
 subroutine read_Kuo_scattering_table
  
  integer :: i,j, klun
  real    :: Kuoprate, Kuodm, Zku
!       real    :: pwc(253), prate(253), Kuodm(253), Zku(253)
  character(len=298) :: headerbuff
       
!--- the 7 frequencies are the GMI primary freqs:
!---  10.65,18.70,23.80,35.5,89.0,165.5,183.3
       
!--- open Kuo Table
       
  call GPM_lun(klun)
  open(unit=klun,file=Kuoscattable,status='old',access='sequential',readonly)
       
!--- read the header (just the column identifier)       
  read(klun,'(a)') headerbuff
       
!--- loop over all the Kext, salbedo, gfactor values based on the PWC (0.001-35.49)
       
  do i = 1,253
    read(klun,'(f11.8,f13.8,f12.8,f8.2,15F12.8,2(F13.8,F12.8,F12.8))') &
       Kuo_pwc(i),Kuoprate,Kuodm,Zku,	&
      (Kuo_kext(j,i),Kuo_salbedo(j,i),Kuo_gfactor(j,i),j=1,7)	    
  enddo        
  close(klun)

!       do i = 1,253
!         write(6,'(F11.8,15F12.8,2(F13.8,F12.8,F12.8))') Kuo_pwc(i),              &
!                    (Kuo_kext(j,i),Kuo_salbedo(j,i),Kuo_gfactor(j,i),j=1,7)
!       enddo
       
  return
 end subroutine read_Kuo_scattering_table

!-------------------------------------------------------------------- 

 SUBROUTINE Kuo_intplte_snow(nf, Snow, kext, salb, asym)
 IMPLICIT NONE
  
 REAL ::    Snow, dX, dY,i     
 INTEGER :: nf, idxS
 INTEGER,PARAMETER :: Kuo_nSWC=253
 REAL ::    kext, salb, asym
 LOGICAL :: check_bdry

!      write(6,*)'k',nf,snow, kuo_pwc(1)

           
 IF( Snow .lt. Kuo_pwc(1))THEN
   kext = 0.
   salb = 0.
   asym = 0.
   RETURN
 END IF
      
!--Determine the index (entry just above the variable of
!--interest in the Kuo LUT).  LUT(index-1) < Variable < LUT(index)    
      check_bdry = .true.
 i = 1
 DO WHILE (check_bdry)
   i = i + 1
   if(i .gt. Kuo_nSWC) then
      i = Kuo_nSWC
      exit
   endif
   IF (Snow .le. Kuo_pwc(i)) check_bdry = .false.
 END DO
 idxS =  i
!
!--Simple linear interpolation
!

 dX = (Snow - Kuo_pwc(idxS-1))/(Kuo_pwc(idxS) - Kuo_pwc(idxS-1))    
 kext	= (1-dX) * Kuo_kext(nf,idxS-1) + dX * Kuo_kext(nf,idxS)
 salb	= (1-dX) * Kuo_salbedo(nf,idxS-1) + dX * Kuo_salbedo(nf,idxS)
 asym	= (1-dX) * Kuo_gfactor(nf,idxS-1) + dX * Kuo_gfactor(nf,idxS)

!      write(6,*)'k',kext,salb,asym
     
 RETURN 
 END SUBROUTINE Kuo_intplte_snow
      
!-----------------------------------------------------------------------
      
 subroutine read_mietable_all(mietable_all_filenamein)  	  
  implicit none
 
  character(len=256)  :: mietable_all_filenamein
  integer	      :: alun,ios

!---  read in lookup table for ice

  call GPM_lun(alun)
  open(unit=alun, file=trim(mietable_all_filenamein),	&
      access='stream', status='old',iostat = ios)
  if(ios .ne. 0) then
      write(*,*)' error opening mie all file:',trim(mietable_all_filenamein)
      stop    
  endif
 
  read(alun) dielec_rlist_ALL
  read(alun) dielec_ilist_ALL
  read(alun) crefindex_ALL
  read(alun) xlist_ALL
  read(alun) mietable_ALL_db
  close(alun)
  
  return
 end subroutine read_mietable_all

!----------------------------------------------------------------------
 subroutine eddington_sp(umu)

!---   chris kummerow
!---   includes asymptotic expressions for term3, term4, and term5 in
!---   the limit of small effective optical depths; bill olson feb, 1995.
     
  integer::  nx,ny,nz,j,i, k
!  real:: lyrtemp(npixs,nscans,0:nglevs) 
!  real:: kext(nglevs),salb(nglevs), asym(nglevs), l(nglevs), h(nglevs)
!  real:: b0(nglevs), b1(nglevs), w(2*nglevs,2*nglevs)
!  real:: bb(2*nglevs), dp(nglevs), dm1(nglevs), z(0:nglevs)
!  real ::iout(0:nglevs), i_in(nglevs+1,100)
  real:: umu, numu, rcond

  real:: lyrtemp(npixs,nscans,0:nlevs) 
  real:: kext(nlevs),salb(nlevs), asym(nlevs), l(nlevs), h(nlevs)
  real:: b0(nlevs), b1(nlevs), w(2*nlevs,2*nlevs)
  real:: bb(2*nlevs), dp(nlevs), dm1(nlevs), z(0:nlevs)
  real ::iout(0:nlevs), i_in(nlevs+1,100)



  real:: xa,xb,xc,xd,ya,yb,yc,dz
  real:: term1,term2,term3,term4,term5

  integer   :: nymid
  real      :: zmid
  integer   :: jin(0:nlevs),jout(0:nlevs)
  real      :: fisot, btemp, xnu, xiup
  integer   :: nang, nn

  fisot = 2.7	

  do nx = 1, npixs
    do ny = 1, nscans
      do k = 0, nlevs
 	lyrtemp(nx,ny,k)=tempprof(k,nx,ny)
      enddo  
    enddo
  enddo
    
!  do i= 88,0,-1
!    write(6,*) i,lyrtemp(tx,ty,i),glev(i)
!  enddo
       
              
!---  slant-path eddington computes tb for one oblique radiance
!---  path at a time.  the coordinate nx,ny is the gridpoint where
!---  the radiance path intersects the earth's surface.

!---  calculate the mid-layer coordinates of the slant path for
!---  both incoming and outgoing paths.

  do nz=0,nlevs
    if(nz .eq. 0) then  			!find altitude of midlayer (or surface)
      zmid = clev(nz)
    else
      zmid = (clev(nz-1)+clev(nz))/2.
    endif

    if (sc_orient .eq. 0) then  		! find nearest pixel to slant path intersection of midlayer.
       jin(nz)  = nint( zmid*sqrt(1.-umu*umu)/(umu*fov_dpr))
       jout(nz) = nint(-zmid*sqrt(1.-umu*umu)/(umu*fov_dpr))
    else
       jin(nz)  = nint(-zmid*sqrt(1.-umu*umu)/(umu*fov_dpr))
       jout(nz) = nint( zmid*sqrt(1.-umu*umu)/(umu*fov_dpr))
    endif
  enddo

!---  loop over pixels  

!  write(6,*)'      nymid      stlev       endlev         Kexttot        Salbtot'

  do nx = 1, npixs
    do ny = 1800, 2100
	       
!---    construct 1-d "column" along incoming ray path

       z(0) = clev(0)
       btemp = tempprof(0,nx,ny)

!---    skip a column with the default value

      if (retr_flag(nx,ny) .eq. 0.) cycle

      do j = 1,nlevs           
	nymid = mod(2*nscans + ny + jin(j) - 1, nscans) + 1     !index of pather intersection of mid-layer
                                                                   ! grid is assumed periodic in y
        z(j)  = clev(j)
        b0(j) = lyrtemp(nx,nymid,j-1)

!      if(ny .eq. tx .and. nx .eq. ty) then
!        write(6,*) nymid,z(j-1),z(j),kexttot(nx,nymid,j),salbtot(nx,nymid,j),jin(j)
!     endif	
	
        b1(j) = (lyrtemp(nx,nymid,j) -lyrtemp(nx,nymid,j-1)) / (z(j) - z(j-1))
        kext(j) = kexttot(nx,nymid,j)
        salb(j) = salbtot(nx,nymid,j)
        asym(j) = asymtot(nx,nymid,j)
        
        if ( nint(kext(j)) .eq. nint(dflt) )cycle    
        l(j) = sqrt(3.*kext(j)*kext(j)*(1. - salb(j))*(1. - salb(j)*asym(j)))
        h(j) = 1.5*kext(j)*(1. - salb(j)*asym(j))
	
      enddo 		 	 

!--- fill in the matrix elements which form the boundary conditions
!--- at the top, bottom and layer interfaces of the cloud.  there are
!--- two quantities, d+ "dp", and d- "dm1" at each boundary, so the 
!--- matrix has dimensions  2*nlevs by 2*nglev.
!--- order of d's:  d+(1),d-(1),d+(2),d-(2), .... , d+(nglev),d-(nglev)

!--  set all matriz elements to zero	

      w = 0.0

!--- fill in the non-zero matrix elements
     
      w(1,1)   = ((ebar(nx,ny) - 2.)*l(1)/h(1)) + ebar(nx,ny)
      w(1,2)   = ((2. - ebar(nx,ny))*l(1)/h(1)) + ebar(nx,ny)
      do  i = 2,2*(nlevs-1),2
        w(i,i-1)   =  (1. - l(i/2)/h(i/2))*exp(+l(i/2)*(z(i/2) - z(i/2-1)))
        w(i,i  )   =  (1. + l(i/2)/h(i/2))*exp(-l(i/2)*(z(i/2) - z(i/2-1)))
        w(i,i+1)   = -(1. - l(i/2+1)/h(i/2+1))
        w(i,i+2)   = -(1. + l(i/2+1)/h(i/2+1))
    
        w(i+1,i-1) =  (1. + l(i/2)/h(i/2))*exp(+l(i/2)*(z(i/2) - z(i/2-1)))
        w(i+1,i)   =  (1. - l(i/2)/h(i/2))*exp(-l(i/2)*(z(i/2) - z(i/2-1)))
        w(i+1,i+1) = -(1. + l(i/2+1)/h(i/2+1))
        w(i+1,i+2) = -(1. - l(i/2+1)/h(i/2+1))
      enddo
      w(2*nlevs,2*nlevs-1) =  (1. + l(nlevs)/h(nlevs))*exp(+l(nlevs) * (z(nlevs)-z(nlevs-1)))
      w(2*nlevs,2*nlevs)   =  (1. - l(nlevs)/h(nlevs))*exp(-l(nlevs) * (z(nlevs)-z(nlevs-1)))

!--- fill in the row of constants in the linear equations
     
      bb(1)    = ebar(nx,ny)*btemp - ebar(nx,ny)*b0(1) - (ebar(nx,ny) - 2.)*b1(1)/h(1)
     
!      write(6,*) btemp,b0(1),ebar(nx,ny)*btemp - ebar(nx,ny)*b0(1)
!      write(6,*) b1(1), h(1), b1(1)/h(1)
!      write(6,*) (ebar(nx,ny) - 2.)*b1(1)/h(1)
!      if(pau) pause 4
     
      do  i = 2,2*(nlevs-1),2
        bb(i)   =  + b1(i/2)/h(i/2) - b1(i/2+1)/h(i/2+1)
        bb(i+1) =  - b1(i/2)/h(i/2) + b1(i/2+1)/h(i/2+1)
      enddo
      bb(2*nlevs)  =  fisot - b0(nlevs) - b1(nlevs)*(z(nlevs) - z(nlevs-1) + 1/h(nlevs))

!      do i = 1,nlevs*2    
!        write(6,*) i,bb(i)
!      enddo
!      if(pau) pause 5

!--- matrix inversion in done in subroutine linpak
     
!      write(6,*)' ready to call linpak'
      
      call linpak(nlevs, w, bb, rcond)
      do i = 1,nlevs
        dp(i) = bb(2*i-1)
        dm1(i) = bb(2*i)
      enddo

!--- after d's are known, calculate surface radiance

      numu = -umu     !set the negative of umu
      
!--- for the following calculations, refer to appendix b of thesis
      
      if (lambert(nx,ny) ) then 
	  	  
!---   calculate the downwelling flux at 81 angles
         nang = 81
         do  nn = 1,nang
           xnu = -(2.*nn - 1.)/(nang*2.)
           i_in(nlevs+1,nn) = fisot

!---     loop through the remaining layers
           do  j = nlevs,1,-1

!---      calculate radiance from top of layer "j"
             xa = b0(j) - 1.5*salb(j)*asym(j)*xnu*b1(j)/h(j)
             xb = b1(j)
             xc = salb(j)*dp(j)*(1. - 1.5*asym(j)*xnu*l(j)/h(j)) 	      
             xd = salb(j)*dm1(j)*(1. + 1.5*asym(j)*xnu*l(j)/h(j))
             ya = kext(j)/xnu
             yb = ya + l(j)
             yc = ya - l(j)
             dz = z(j) - z(j-1) 

             term1 = i_in(j+1,nn)*exp(ya*dz)
             term2 = xa*(1. - exp(ya*dz))
            if(abs(ya*dz) .lt. 1.e-5) then
                term3=-xb*ya*dz*dz
             else
                term3 = xb/ya*(exp(ya*dz)*(1. - ya*dz) - 1.)
             endif
             if(abs(yb*dz) .lt. 1.e-5) then
                term4=-xc*ya*dz
             else
                term4 = xc*ya/yb*(1. - exp(yb*dz))
             endif
             if(abs(yc*dz) .lt. 1.e-5) then
                term5=-xd*ya*dz
             else
                term5 = xd*ya/yc*(1. - exp(yc*dz))
             endif
             i_in(j,nn) = term1 + term2 + term3 + term4 + term5
           enddo
         enddo

!--- calculate the total downwelling flux reaching the surface
         xiup = 0.
         do  nn = 1,nang
           xiup = xiup + i_in(1,nn)*(1./nang)*(2.*nn-1.)/(2.*nang)
         end do
         xiup = 2.*xiup
           
      else     !if lambert
	  
!--- calculate the downwelling flux at angle umu only
            
	 nn = 22            ! this is a dummy index for i_in
         i_in(nlevs+1,nn) = fisot

!--- loop through the remaining layers
         do  j = nlevs,1,-1

!---    calculate radiance from top of layer "j"
           xa = b0(j) - 1.5*salb(j)*asym(j)*numu*b1(j)/h(j)
           xb = b1(j)
           xc = salb(j)*dp(j)*(1. - 1.5*asym(j)*numu*l(j)/h(j))		  
           xd = salb(j)*dm1(j)*(1. + 1.5*asym(j)*numu*l(j)/h(j))
           ya = kext(j)/numu
           yb = ya + l(j)
           yc = ya - l(j)
           dz = z(j) - z(j-1)

           term1 = i_in(j+1,nn)*exp(ya*dz)
           term2 = xa*(1. - exp(ya*dz))
           if(abs(ya*dz) .lt. 1.e-5) then
              term3=-xb*ya*dz*dz
           else
              term3 = xb/ya*(exp(ya*dz)*(1. - ya*dz) - 1.)
           endif
           if(abs(yb*dz) .lt. 1.e-5) then
              term4=-xc*ya*dz
           else
              term4 = xc*ya/yb*(1. - exp(yb*dz))
           endif
           if(abs(yc*dz) .lt. 1.e-5) then
              term5=-xd*ya*dz
           else
              term5 = xd*ya/yc*(1. - exp(yc*dz))
           endif
           i_in(j,nn) = term1 + term2 + term3 + term4 + term5
        
	 enddo
         xiup = i_in(1,22)

      endif

      iout(0) = emis(nx,ny)*btemp + (1. - emis(nx,ny))*xiup

!---     recalculate flux conditions along outgoing path.
!---      construct 1-d "column" along outgoing ray path

      z(0) = clev(0)     !0
      btemp = tempprof(0,nx,ny)

      do j = 1,nlevs

!--- index of path intersection of midlayer; grid is assumed
!--- periodic in y.
    
        nymid = mod(2*nscans + ny + jout(j) - 1, nscans) + 1
            
        z(j)  = clev(j)
        b0(j) = lyrtemp(nx,nymid,j-1)
        b1(j) = (lyrtemp(nx,nymid,j) - lyrtemp(nx,nymid,j-1)) / (z(j) - z(j-1))
        kext(j) = kexttot(nx,nymid,j)
        salb(j) = salbtot(nx,nymid,j)
        asym(j) = asymtot(nx,nymid,j)

	if ( nint(kext(j)) .eq. nint(dflt) )cycle

        l(j) = sqrt(3.*kext(j)*kext(j)*(1. - salb(j))* (1. - salb(j)*asym(j)))
        h(j) = 1.5*kext(j)*(1. - salb(j)*asym(j))

      enddo

!--- fill in the matrix elements which form the boundary conditions
!--- at the top, bottom and layer interfaces of the cloud.  there are
!--- two quantities, d+ "dp", and d- "dm1" at each boundary, so the 
!--- matrix has dimensions  2*nlevs by 2*nlevs.
!--- order of d's:  d+(1),d-(1),d+(2),d-(2), .... , d+(nlevs),d-(nlevs)

!--- set all matriz elements to zero	
      
      w = 0.0
      
!--- fill in the non-zero matrix elements
         
      w(1,1)   = ((ebar(nx,ny) - 2.)*l(1)/h(1)) + ebar(nx,ny)
      w(1,2)   = ((2. - ebar(nx,ny))*l(1)/h(1)) + ebar(nx,ny)
      do  i = 2,2*(nlevs-1),2
        w(i,i-1)   =  (1. - l(i/2)/h(i/2))*exp(+l(i/2)*(z(i/2) - z(i/2-1)))
        w(i,i  )   =  (1. + l(i/2)/h(i/2))*exp(-l(i/2)*(z(i/2) - z(i/2-1)))
        w(i,i+1)   = -(1. - l(i/2+1)/h(i/2+1))
        w(i,i+2)   = -(1. + l(i/2+1)/h(i/2+1))

        w(i+1,i-1) =  (1. + l(i/2)/h(i/2))*exp(+l(i/2)*(z(i/2) - z(i/2-1)))
        w(i+1,i)   =  (1. - l(i/2)/h(i/2))*exp(-l(i/2)*(z(i/2) - z(i/2-1)))
        w(i+1,i+1) = -(1. + l(i/2+1)/h(i/2+1))
        w(i+1,i+2) = -(1. - l(i/2+1)/h(i/2+1))
      enddo
      w(2*nlevs,2*nlevs-1) =  (1. + l(nlevs)/h(nlevs))*exp(+l(nlevs) * (z(nlevs)-z(nlevs-1)))
      w(2*nlevs,2*nlevs)   =  (1. - l(nlevs)/h(nlevs))*exp(-l(nlevs) * (z(nlevs)-z(nlevs-1)))

!--- fill in the row of constants in the linear equations

      bb(1) = ebar(nx,ny)*btemp - ebar(nx,ny)*b0(1) - (ebar(nx,ny) - 2.)*b1(1)/h(1)
      do  i = 2,2*(nlevs-1),2
        bb(i)	=  + b1(i/2)/h(i/2) - b1(i/2+1)/h(i/2+1)
        bb(i+1) =  - b1(i/2)/h(i/2) + b1(i/2+1)/h(i/2+1)
      enddo
      bb(2*nlevs)  =  fisot - b0(nlevs) - b1(nlevs)*(z(nlevs) - z(nlevs-1) + 1/h(nlevs))

!--- matrix inversion in done in subroutine linpak
         
      call linpak(nlevs, w, bb, rcond)

      do  i = 1,nlevs
        dp(i) = bb(2*i-1)
        dm1(i) = bb(2*i)
      enddo

      do  j = 1,nlevs
!--- calculate the upwelling radiances at the top of each layer j
        xa = b0(j) - 1.5*salb(j)*asym(j)*umu*b1(j)/h(j)
        xb = b1(j)
        xc = salb(j)*dp(j)*(1. - 1.5*asym(j)*umu*l(j)/h(j))		 
        xd = salb(j)*dm1(j)*(1. + 1.5*asym(j)*umu*l(j)/h(j))
        ya = kext(j)/umu
        yb = ya + l(j)
        yc = ya - l(j)
        dz = z(j) - z(j-1)
        term1 = iout(j-1)*exp(-ya*dz)
        term2 = xa*(1. - exp(-ya*dz))
        if(abs(ya*dz) .lt. 1.e-5) then
            term3=0.
        else
            term3 = xb/ya*(ya*dz - 1. + exp(-ya*dz))
        endif
        if(abs(yb*dz) .lt. 1.e-5) then
            term4=xc*ya*dz*exp(-ya*dz)
        else
            term4 = xc*ya/yb*(exp( (yb-ya)*dz ) - exp(-ya*dz) )
        endif
        if(abs(yc*dz) .lt. 1.e-5) then
            term5=xd*ya*dz*exp(-ya*dz)
        else
            term5 = xd*ya/yc*exp(-ya*dz)*(exp(yc*dz) - 1.)
        endif
        iout(j) = term1 + term2 + term3 + term4 + term5
      end do

      tb(nx,ny) = iout(nlevs)
	  
    enddo    !nscan
  enddo   !npixs

  return
 end subroutine eddington_sp

!------------------------------------------------------------------------

end module GPM_TbsimulatorV7_procedures
