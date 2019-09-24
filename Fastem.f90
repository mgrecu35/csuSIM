  MODULE FASTEM

  USE Type_Kinds, ONLY: fp
  USE CRTM_FastemX
  USE CRTM_MWwaterCoeff
  USE Message_Handler

  IMPLICIT NONE

  INTEGER :: FASTEM_MODEL=0

  CONTAINS

!--------------------------------------------------------------------------------
! NAME:
!       Compute_Fastem
!
! PURPOSE:
!       Subroutine to compute the Fastem4, Fastem5, or Fastem6 microwave sea surface
!       emissivity and reflectivity.
!
! CALLING SEQUENCE:
!       CALL Compute_Fastem(
!              Frequency    , &  ! Input
!              Zenith_Angle , &  ! Input
!              Temperature  , &  ! Input
!              Salinity     , &  ! Input
!              Wind_Speed   , &  ! Input
!              Emissivity   , &  ! Output
!              Reflectivity , &  ! Output
!              Azimuth_Angle)
!
!
! INPUTS:
!       Frequency:      Microwave frequency.
!                       UNITS:      GHz
!                       TYPE:       REAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Zenith_Angle:   Sensor zenith angle at the sea surface
!                       UNITS:      Degrees
!                       TYPE:       REAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Temperature:    Sea surface temperature
!                       UNITS:      Kelvin, K
!                       TYPE:       REAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Salinity:       Water salinity
!                       UNITS:      ppt (parts per thousand)
!                       TYPE:       REAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Wind_Speed:     Sea surface wind speed
!                       UNITS:      m/s
!                       TYPE:       REAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Azimuth_Angle:  Relative azimuth angle (wind direction - sensor azimuth)
!                       UNITS:      Degrees
!                       TYPE:       REAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
! OUTPUTS:
!       Emissivity:     The surface emissivity
!                       UNITS:      N/A
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Rank-1, 4-elements (n_Stokes)
!                       ATTRIBUTES: INTENT(OUT)
!
!       Reflectivity:   The surface reflectivity.
!                       UNITS:      N/A
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Rank-1, 4-elements (n_Stokes)
!                       ATTRIBUTES: INTENT(OUT)
!
!--------------------------------------------------------------------------------

  SUBROUTINE Compute_Fastem(model,freq,eia,sst,salinity,wsp,azimuth,emis,verbose)

  integer            :: model
  real               :: freq
  real               :: eia
  real               :: sst
  real               :: salinity
  real               :: wsp
  real               :: azimuth
  real               :: emis(4)
  integer, optional  :: verbose

  character(len=200) :: File_Path
  character(len=200) :: MWWaterCoeff_File
  character(len=200) :: msg,pid_msg
  integer            :: err_stat
  logical            :: Quiet
  integer            :: Process_ID
  integer            :: Output_Process_ID
  real(fp)           :: frequency
  real(fp)           :: view_angle
  real(fp)           :: tsfc
  real(fp)           :: salin
  real(fp)           :: wndspd
  real(fp)           :: azimuth_angle
  real(fp)           :: emissivity(4)
  real(fp)           :: reflectivity(4)
  TYPE(iVar_type)    :: iVar

  if (PRESENT(verbose)) then
    Quiet=0
  else
    Quiet=1
  endif
  File_Path = '../fastem/Little_Endian/'
  if (model .eq. 4) then
    MWwaterCoeff_File = 'FASTEM4.MWwater.EmisCoeff.bin'
  else
    if (model .eq. 5) then
      MWwaterCoeff_File = 'FASTEM5.MWwater.EmisCoeff.bin'
    else
      MWwaterCoeff_File = 'FASTEM6.MWwater.EmisCoeff.bin'
    endif
  endif
  MWwaterCoeff_File  = TRIM(ADJUSTL(File_Path)) // TRIM(MWwaterCoeff_File)

  ! Load coefficient file for MW water

  if (FASTEM_MODEL .eq. 0 .or. model .ne. FASTEM_MODEL) then
    err_stat = CRTM_MWwaterCoeff_Load( MWwaterCoeff_File, &
                 Quiet             = Quiet            , &
                 Process_ID        = Process_ID       , &
                 Output_Process_ID = Output_Process_ID  )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error loading MWwaterCoeff data from '//TRIM(MWwaterCoeff_File)
      CALL Display_Message("FASTEM",TRIM(msg)//TRIM(pid_msg),err_stat )
      STOP
    END IF
    FASTEM_MODEL = model
  endif

  ! FastemX model

  frequency     = freq
  view_angle    = eia
  tsfc          = sst
  salin         = salinity
  wndspd        = wsp
  azimuth_angle = azimuth

  CALL Compute_FastemX(           &
         MWwaterC,                &  ! Input model coefficients
         frequency,               &  ! Input
         view_angle,              &  ! Input
         tsfc,                    &  ! Input
         salin,                   &  ! Input
         wndspd,                  &  ! Input
         iVar,                    &  ! Internal variable output
         emissivity,              &  ! Output
         reflectivity,            &  ! Output
         Azimuth_Angle = azimuth_angle)  ! Optional input (wind direction - sensor azimuth)

  emis(:) = emissivity(:)

  END SUBROUTINE Compute_Fastem
  END MODULE Fastem
