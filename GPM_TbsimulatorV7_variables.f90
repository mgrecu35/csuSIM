MODULE GPM_TbsimulatorV7_variables

implicit none
save

!--------------------------------------------------------------------------------

 integer	      :: LUT_nfreq
 integer,allocatable  :: snowf(:)

 real,parameter    :: salinity = 35.0

 integer,parameter :: tx = 6, ty = 1968
 logical,parameter :: pau = .false.
 real,parameter    :: pi = 3.14159265

!-- variables 

 integer,parameter :: npixs=25, nchans=13, nlevs=88, mplyrs=46, nglevs=28
 integer           :: nscans, gscans
 real              :: sc_orient
 character(len=5)  :: sensor ,satcode

 integer :: pixstart, pixend, scanstart, scanend

!-- define memory space for input file profile variables
  
real,allocatable    :: lat(:,:), lon(:,:), elev(:,:)
integer,allocatable :: scntime(:,:),retr_flag(:,:),sfctype(:,:)
real,allocatable    :: emiss_strm(:,:,:),wind(:,:),datasource(:,:)
real,allocatable    :: msfcprcp(:,:), csfcprcp(:,:), ccnvprcp(:,:)
integer,allocatable :: prcpflag(:,:), prcptype(:,:)
real,allocatable    :: tempprof(:,:,:),pressprof(:,:,:)
real,allocatable    :: mixr(:,:,:), cloudwater(:,:,:)
real,allocatable    :: rainwater(:,:,:),snowwater(:,:,:),SLH(:,:,:)
real,allocatable    :: closestGMI_DPR(:,:,:)
real,allocatable    :: dm(:,:,:),mu(:,:,:),nw(:,:,:)
real,allocatable    :: cTbssim(:,:,:),gTbsobs(:,:,:)
real,allocatable    :: mTbssim(:,:,:) 
  
!-- define program  variables
 
 integer,parameter :: gmipol(13) = (/0,1,0,1,0,0,1,0,1,0,1,1,1/) 
 real,parameter    :: gmifreqs(13)=(/10.65,10.65,18.7,18.7,23.8,36.64,   &   !GMI frequencies
                           36.64,89.0,89.0,166.0,166.0,180.31,190.31/)
 
 real,parameter    :: gmieia(13)=(/52.80, 52.80, 52.80, 52.80, 52.80,   &
                                   52.80, 52.80, 52.80, 52.80, 49.19,   &
                                   49.19, 49.19, 49.19/)
 
 	 
 real,allocatable     :: kexttot(:,:,:)
 real,allocatable     :: salbtot(:,:,:)
 real,allocatable     :: asymtot(:,:,:)
 real,allocatable     :: emis(:,:)
 real,allocatable     :: ebar(:,:)
 real,allocatable     :: tb(:,:) 
 real,allocatable     :: emissall(:,:,:)
 real,allocatable     :: ebarall(:,:,:)
! real,allocatable     :: emiss_stream(:,:,:)

 real,allocatable     :: dTbssim(:,:,:)
 real,allocatable     :: deltadTbssim(:,:,:)
 
 logical,allocatable :: lambert(:,:)
 integer,parameter   :: fov_dpr=5.0
 real,parameter      :: dflt = -999.
 real,parameter      :: iflt = -999.

!--- variables for fpa

 real,allocatable :: csfcprcp_fpa(:,:)           ! surface variables (npix,nscans)
 real,allocatable :: msfcprcp_fpa(:,:)          
 real,allocatable :: ccnvprcp_fpa(:,:)
 real,allocatable :: datasource_fpa(:,:)
 real,allocatable :: emiss_fpa(:,:,:)            !emissivity(13,npixs,nscans)

 real,allocatable :: rainwater_fpa(:,:,:)        !hydro profiles (nlev,npix,nscans)
 real,allocatable :: cloudwater_fpa(:,:,:)
 real,allocatable :: snowwater_fpa(:,:,:)
 real,allocatable :: SLH_fpa(:,:,:)

 real,allocatable :: dTbssim_fpa(:,:,:)            ! Tbs (nchans,npix,nscans)
 real,allocatable :: mTbssim_fpa(:,:,:)
 real,allocatable :: cTbssim_fpa(:,:,:)
 real,allocatable :: deltadTbssim_fpa(:,:,:)
 real,allocatable :: gTbsobs_fpa(:,:,:)

!--- used in the re-layer routine 
 
 real    :: glev(0:28)= (/0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,  & 
                          5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,  &
			 10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0/)			 
			 
 real    :: clev(0:88)= (/ 0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75,  &
                           2.00, 2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75,  &
                           4.00, 4.25, 4.50, 4.75, 5.00, 5.25, 5.50, 5.75,  &
                           6.00, 6.25, 6.50, 6.75, 7.00, 7.25, 7.50, 7.75,  &
                           8.00, 8.25, 8.50, 8.75, 9.00, 9.25, 9.50, 9.75,  &
                          10.00,10.25,10.50,10.75,11.00,11.25,11.50,11.75,  &
                          12.00,12.25,12.50,12.75,13.00,13.25,13.50,13.75,  &
                          14.00,14.25,14.50,14.75,15.00,15.25,15.50,15.75,  &
                          16.00,16.25,16.50,16.75,17.00,17.25,17.50,17.75,  &
                          18.00,18.25,18.50,18.75,19.00,19.25,19.50,19.75,  &
                          20.00,20.25,20.50,20.75,21.00,21.25,21.50,21.75,22.00/) 

 
!-- variables for the Kuo Scattering Table

 character(len=128)   :: Kuoscattable		
 real :: Kuo_pwc(253), Kuo_kext(7,253)
 real :: Kuo_salbedo(7,253), Kuo_gfactor(7,253)
 
!--- arrays for re-layering near freezing level

 real, allocatable :: mp_rainprof(:,:,:),mp_cloudprof(:,:,:)
 real, allocatable :: mp_snowprof(:,:,:)
 real, allocatable :: mp_cloud_ice(:,:,:),mp_graupprof(:,:,:)
 real, allocatable :: mp_mixr(:,:,:), mp_reflyr(:,:,:), mp_dm(:,:,:)
 real, allocatable :: mp_nw(:,:,:),mp_mu(:,:,:)
 real, allocatable :: mp_melt_water(:,:,:), mp_meltfrac(:,:,:)
 real, allocatable :: mp_tprof(:,:,:), mpht(:,:,:)
 real, allocatable :: mp_pressprof(:,:,:) 
  
!--- variables for MonoRTM lookup table

 integer	:: nfreq_lut
 integer	:: nchan_lut
 integer	:: npres_lut
 integer	:: ntemp_lut
 integer	:: nrmix_lut

 integer, allocatable :: ifreq_lut(:)
 integer, allocatable :: ipol_lut(:)
 real, allocatable    :: freq_lut(:)
 real, allocatable    :: pres_lut(:)
 real, allocatable    :: temp_lut(:)
 real, allocatable    :: rmix_lut(:)
 real, allocatable    :: kabs_lut(:,:,:,:)

!--- variables for the mie scattering lookup table access -----------------
      
 real,parameter      :: ndensityALL_list=10	 
 real,parameter      :: tmin_ALL=223.15, tmax_ALL=303.15
 real,parameter      :: tinc_ALL=0.1
 
 real,parameter      :: xmin_ALL=1.0e-3, xmax_ALL=2.0e2
 real		     :: xinc_ALL
 real,parameter      :: nxlist_ALL=184
 real		     :: xlist_ALL(nxlist_ALL)
 
 real,parameter      :: dielec_rmin_ALL=1.05   !1.085046
 real,parameter      :: dielec_rmax_ALL=8.90   !8.722480
 real		     :: dielec_rinc_ALL
 real,parameter      :: ndielec_rlist_ALL= 108 !5%=44  1%=206, 2%=108 increments
 real		     :: dielec_rlist_ALL(ndielec_rlist_ALL)
        
 real,parameter      :: dielec_imin_ALL=6.0e-6  !6.8959e-6
 real,parameter      :: dielec_imax_ALL=3.15e+0  !3.044035
 real		     :: dielec_iinc_ALL 
 real,parameter      :: ndielec_ilist_ALL= 666    !5%=270 1%=1324, 2%=666 increments
 real		     :: dielec_ilist_ALL(ndielec_ilist_ALL)    
       
 type		     :: crefindex_structure_ALL
    real	     :: real_cref
    real	     :: imag_cref
    real	     :: loc_cref      ! This is the sorting INDEX
 end type crefindex_structure_ALL
 type(crefindex_structure_ALL) :: crefindex_ALL(ndielec_rlist_ALL,ndielec_ilist_ALL)  
  
 real,parameter      :: ncrefindex_ALL = 16496.  !5%=3387.0  2%=16496  increments
 type		     :: mietable_ALL_structure
    real	     :: real_cref
    real	     :: imag_cref
    real	     :: xsize
    real	     :: qsca
    real	     :: qext
    real	     :: asym
    real	     :: qbsca
 end type mietable_ALL_structure
 type(mietable_ALL_structure) :: mietable_ALL_db(ncrefindex_ALL,nxlist_ALL)

!----------------------------------------------------------------------------------------------

end module GPM_TbsimulatorV7_variables
