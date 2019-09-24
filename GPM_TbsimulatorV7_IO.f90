module GPM_TbsimulatorV7_IO
  
use GPM_TbsimulatorV7_variables      

implicit none
       
contains

!-----------------------------------------------------------------------

subroutine GPM_read_strmfile(input_file)
 
 integer :: rlun,ios,i,j,k, ilen, intran, glun1, pix,scan
 character(len=256) :: input_file, infile, command
 logical            :: comp
 real               :: rannum
 character(len=5)   :: cran  
 character(len=16)  :: cfile
 
!--- define input data structure for reading file

 type :: strmstructure
  integer :: pix1, scan1, sc_orient
  real    :: lat,lon,elev
  integer :: scntime(6)
  integer :: retr_flag
  integer :: sfctype
  integer :: datasource
    
  real    :: emissivity(13)
  real    :: wind
  real    :: msfcprcp
  real    :: csfcprcp
  integer :: prcpflag
  integer :: prcptype
  real    :: tempprof(0:88)
  real    :: pressprof(0:88)
  real    :: mixr(88)
  real    :: cloudwater(88)
  real    :: rainwater(88)
  real    :: snowwater(88)
  real    :: SLH(80)

  integer :: closestgmi_dpr(2)
  real    :: dm(88),mu(88),nw(88)

  real    :: cTbssim(13)
  real    :: gTbsobs(13)
  real    :: mTbssim(13)

 end type strmstructure
 type(strmstructure) :: in 
 
!--- check if file is compressed, if so then uncompress it
 
 comp = .false.
 ilen  = len_trim(input_file)
 do i = 1,ilen-2
   if(input_file(i:i+2) .eq. '.gz') comp = .true.
 enddo
 ilen=14
 if(comp) then
       call random_seed
       call random_number(rannum)
       intran = nint(rannum*100000)
       write(cran,'(I5.5)') intran
       
       cfile(1:16) = 'gunziptemp/' // cran
       write(6,'(a,a)')'  cp & decmprss:',trim(input_file)
       command = 'gunzip -c ' // trim(input_file) // '>' // cfile
       call system(trim(command))
       infile = cfile
 else
       infile = input_file
       write(6,'(a,a)')'  reading file : ',trim(input_file(1:ilen))
 endif

!--- open input file

 write(*,*)'  opening input file: ', trim(infile(1:ilen))      
 open(unit=glun1,file=trim(infile(1:14)),access='stream',iostat=ios,status='old')
 if(ios .ne. 0) then
       write(*,*)'  error opening input file : ',trim(input_file)
       stop 80
 endif

!--- read header info on nscans (dpr scans), gscans (GMI scans)

 read(glun1,iostat=ios) nscans,gscans
 if(ios .ne. 0) then
       write(6,*)'  error reading nscans and gscans header'
       stop 81
 endif   
 write(6,*)'   nscans, gscans = ', nscans,gscans

!--- Allocate all program arrays to hold input_file data

 allocate(lat(npixs,nscans), lon(npixs,nscans), elev(npixs,nscans))
 allocate(scntime(nscans,6),retr_flag(npixs,nscans),sfctype(npixs,nscans))
 allocate(datasource(npixs,nscans))
 allocate(emiss_strm(nchans,npixs,nscans),wind(npixs,nscans))
 allocate(msfcprcp(npixs,nscans), csfcprcp(npixs,nscans), ccnvprcp(npixs,nscans))
 allocate(prcpflag(npixs,nscans), prcptype(npixs,nscans))
 allocate(tempprof(0:nlevs,npixs,nscans),pressprof(0:nlevs,npixs,nscans))
 allocate(mixr(nlevs,npixs,nscans), cloudwater(nlevs,npixs,nscans))
 allocate(rainwater(nlevs,npixs,nscans),snowwater(nlevs,npixs,nscans),SLH(nlevs,npixs,nscans))
 allocate(closestGMI_DPR(2,npixs,nscans))
 allocate(dm(nlevs,npixs,nscans),mu(nlevs,npixs,nscans),nw(nlevs,npixs,nscans))
 allocate(cTbssim(nchans,npixs,nscans),gTbsobs(nchans,npixs,nscans))
 allocate(mTbssim(nchans,npixs,nscans))
  	       
!---  read input CMB profile data file

 do scan = 1,nscans
   do pix = 1,npixs
      read(glun1,iostat=ios) in
      if(ios .ne. 0) then
         if(in%scan1 .eq. 1 .and. in%pix1 .eq. 1) then
            write(*,*)' Error reading datafile at first pixel'
	    stop 81
          endif
      endif
 
!--- assign input values to program arrays
     
      sc_orient                 = in%sc_orient
      lat(pix,scan)             = in%lat
      lon(pix,scan)             = in%lon      
      elev(pix,scan)            = in%elev  
      scntime(scan,:)           = in%scntime(:)
      retr_flag(pix,scan)       = in%retr_flag
      sfctype(pix,scan)         = in%sfctype
      datasource(pix,scan)      = in%datasource
      emiss_strm(:,pix,scan)    = in%emissivity(:)
      wind(pix,scan)            = in%wind
      msfcprcp(pix,scan)        = in%msfcprcp
      csfcprcp(pix,scan)        = in%csfcprcp
      prcpflag(pix,scan)        = in%prcpflag
      prcptype(pix,scan)        = in%prcptype
      tempprof(:,pix,scan)      = in%tempprof(:)
      pressprof(:,pix,scan)     = in%pressprof(:)
      mixr(:,pix,scan)          = in%mixr(:)
      cloudwater(:,pix,scan)    = in%cloudwater(:)
      rainwater(:,pix,scan)     = in%rainwater
      snowwater(:,pix,scan)     = in%snowwater
      SLH(:,pix,scan)           = in%SLH
      closestgmi_dpr(:,pix,scan)= in%closestgmi_dpr(:)
      dm(:,pix,scan)            = in%dm(:)
      mu(:,pix,scan)            = in%mu(:)
      nw(:,pix,scan)            = in%nw(:)
      cTbssim(:,pix,scan)       = in%cTbssim(:)
      gTbsobs(:,pix,scan)       = in%gTbsobs(:)
      mTbssim(:,pix,scan)       = in%mTbssim(:)
      
      if(prcptype(pix,scan) .eq. 2) ccnvprcp(pix,scan) = csfcprcp(pix,scan)   !set cnvprcp precip



      if(pix .eq. tx .and. scan .eq. ty) then    !these are set in *_variables.f90 routine
         write(6,*)' Profiles directly from .strm file - these match CMB profiles for pixel tx,ty'
	 write(6,*)' level  Mixing Ratio   Cloudwater  Snowwater   Rainwater '
	 do i = 88,1,-1
	   write(6,'(i5,5F12.3)')i, mixr(i,pix,scan),cloudwater(i,pix,scan), snowwater(i,pix,scan),rainwater(i,pix,scan)
	   
	 enddo
	 
	 write(6,*)
	 write(6,*)' Observed Tbs   CMB Tbs   CMB Emissitity'
	 do i = 1,13
	   write(6,'(3x,3F10.3)') gTbsobs(i,pix,scan),cTbssim(i,pix,scan),emiss_strm(i,pix,scan)
	 enddo
      endif
      
      
   enddo !pix
 enddo  !scan
 
 return
 end subroutine GPM_read_strmfile

!-----------------------------------------------------------------------

 subroutine GPM_write_simfile(outfile, doFPA)
  integer olun, ios
  character(len=256) :: outfile
!  character(len=5) satcode, sensor
  integer :: outchans=13, scan, pix, icnt(4)=0
  logical :: doFPA
  
!--- open output file
   
  call GPM_lun(olun)
   
  open(unit=olun,file=outfile,access='stream',status='unknown',iostat=ios)
  if(ios .ne. 0) then
      write(*,*)' error opening output_file : ',trim(outfile)
      stop 90
  endif
  write(*,*)' writing output to : ', trim(outfile)

!--- write the header information for this file

  write(olun,iostat=ios) satcode,sensor
  if(ios .ne. 0) then
      write(*,*)' error writing header : ',satcode,' ',sensor
      stop 91
  endif    
	
  write(olun,iostat=ios) gmifreqs(1:outchans)
  if(ios .ne. 0) then
      write(*,*)' error writing header : freq_sat'
      stop 92
  endif
  write(olun,iostat=ios) gmieia(1:outchans)
  if(ios .ne. 0) then
      write(*,*)' error writing header : nominal_EIA'
      stop 93
  endif

  write(olun,iostat=ios) pixstart, pixend, scanstart, scanend
  if(ios .ne. 0) then
      write(*,*)' error writing npixs,nscans '
      stop 94
  endif


!--- output the data

  write(*,*)' writing output scans, doFPA = ',doFPA
  do scan = scanstart, scanend              !these are the ranges of the largest
    do pix = pixstart, pixend               !footprint frequency
	   
!--- write out the CMB (DPR) pixel locations and All simulated Tbs
     
      if(.not. doFPA) then
        if(retr_flag(pix,scan) .eq. 1  .and. gTbsobs(1,pix,scan) .ne. -1) then         !only writes pixels at GMI pixels locations
           write(olun,iostat=ios) pix,scan, datasource(pix,scan), lat(pix,scan),    & 
            lon(pix,scan), elev(pix,scan), scntime(scan,1:6), sfctype(pix,scan),  &
	    csfcprcp(pix,scan), msfcprcp(pix,scan), ccnvprcp(pix,scan),           &
            rainwater(:,pix,scan), cloudwater(:,pix,scan),                        &
	    snowwater(:,pix,scan), SLH(:,pix,scan), emiss_strm(:,pix,scan),     &                             
	    cTbssim(:,pix,scan), mTbssim(:,pix,scan), gTbsobs(:,pix,scan),        &
	    dTbssim(:,pix,scan)
!	    ,deltadTbssim_fpa(:,pix,scan), Tbs_bias(:,pix,scan)

           if(ios .ne. 0) then	      
              write(*,*)' error writing pixel data ',pix,scan
	      stop 95
	   endif
	      	      
           if(sfctype(pix,scan) .eq. 1) then 
              icnt(1) = icnt(1) + 1            !ocean pixels
	   else 
              icnt(2) = icnt(2) + 1          !all other pixels
           endif	   
        endif
      
      else      

!--- write out the footprint averaged parameters and Tbs

        if(retr_flag(pix,scan) .eq. 1  .and. gTbsobs(1,pix,scan) .ne. -1) then         !only writes pixels at GMI pixels locations
           write(olun,iostat=ios) pix,scan, datasource_fpa(pix,scan), lat(pix,scan),     & 
            lon(pix,scan), elev(pix,scan), scntime(scan,1:6), sfctype(pix,scan),   &
	    csfcprcp_fpa(pix,scan), msfcprcp_fpa(pix,scan), ccnvprcp_fpa(pix,scan), &
            rainwater_fpa(:,pix,scan), cloudwater_fpa(:,pix,scan),                 &
	    snowwater_fpa(:,pix,scan), SLH_fpa(:,pix,scan), emiss_fpa(:,pix,scan), &                             
	    cTbssim_fpa(:,pix,scan), mTbssim_fpa(:,pix,scan), gTbsobs(:,pix,scan), &
	    dTbssim_fpa(:,pix,scan)
!	    ,deltadTbssim_fpa(:,pix,scan), Tbs_bias(:,pix,scan)

           if(ios .ne. 0) then	      
              write(*,*)' error writing pixel data ',pix,scan
	      stop 95
	   endif
	      	      
           if(sfctype(pix,scan) .eq. 1) then 
              icnt(1) = icnt(1) + 1            !ocean pixels
	   else 
              icnt(2) = icnt(2) + 1          !all other pixels
           endif	   
        endif      

      endif  ! doFPA


  
    enddo   !pix
  enddo   !scan

  close(olun)


  write(*,*)' number of ocean sensor pixels written : ', icnt(1)
  write(*,*)' number of  land sensor pixels written : ', icnt(2)
  write(*,*)' total number sensor pixels written    : ', icnt(1) + icnt(2)
  write(*,*)


end subroutine GPM_write_simfile

!-----------------------------------------------------------------------
	
 subroutine GPM_lun(ilun)
!    
!--- This routine gets an open logical unit number starting with 100
!
  integer  :: ilun
  logical  :: llun
  character(len=128) :: blank=' '
!
  do ilun = 100,201
    inquire(unit=ilun, opened=llun)	    
    if(.not. llun) exit
  enddo
  if(ilun .eq. 201) then
      write(*,*)' error assigning LUN '
      stop 'lun error'
  endif
  return
 end subroutine GPM_lun  

!-----------------------------------------------------------------------


end module GPM_TbsimulatorV7_IO 
