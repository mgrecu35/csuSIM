  program cmpr_emis

!*---------------------------------------------------------------------------*/
!* cmpr_emis: Sample code to run and compare output emissivity values using  */
!*            FASTEM4/5/6 and RSS ocean emissivity models. Currently the     */
!*            values for output frequency, eia, and salinity are hardcoded   */
!*            in the variable declaration section. Values for SST and wind   */
!*            speed are specified in do loops to create an output array of   */
!*            emissivity values.                                             */
!*                                                                           */
!* subroutines:                                                              */
!*    find_surface_tb: The RSS emissivity code (Meisner and Wentz, 2012).    */
!*                     Output values include specular sea surface emis (e0)  */
!*                     and wind induced isotropic emis. (ewind) for both     */
!*                     v-pol and h-pol as well as and wind direction emis    */
!*                     (edirstokes) for all four stokes vectors.             */
!*                                                                           */
!*    compute_fastem:  The Fastem code. Note that the same code is used for  */
!*                     FASTEM4, FASTEM5, and FASTEM6. The first input        */
!*                     parameter specifies which model to use, which selects */
!*                     the appropriate input emissivity file. The output     */
!*                     emissivity values are provided for each of the four   */
!*                     Stokes parameters (i.e. vert, horiz, U, V).           */
!*---------------------------------------------------------------------------*/

  USE Fastem
  USE RSS_RTM
  IMPLICIT NONE

  integer :: i
  integer :: model=6
  real    :: freq(6)=(/10.65,18.7,23.8,36.64,89.0,166.0/)
  real    :: eia=53.0
  real    :: salinity=35.0
  real    :: sst
  real    :: wsp
  real    :: azimuth=0.0
  real    :: emis_rss(4)
  real    :: emis_fm4(4)
  real    :: emis_fm5(4)
  real    :: emis_fm6(4)
  real    :: e0(2),ewind(2),edir(4)
  integer :: ifreq,isst,iwsp,ios

  type :: output_struct
    real :: freq
    real :: sst
    real :: wsp
    real :: emis_rss(2)
    real :: emis_fm4(2)
    real :: emis_fm5(2)
    real :: emis_fm6(2)
  end type output_struct
  type(output_struct) :: drec

  open(unit=20,file='cmpr_emis.dat',form='binary',status='unknown',iostat=ios)

  do ifreq=1,6
    do isst=0,15
      sst = 275.0 + isst*2.0
      do iwsp=0,15
        wsp = iwsp*2.0
	
        call find_surface_tb(freq=freq(ifreq), surtep=sst,ssws=wsp, phir=azimuth, tht=eia, &
                             sal=salinity, e0=e0, ewind=ewind, edirstokes=edir)
        emis_rss(1) = e0(1) + ewind(1) + edir(1)
        emis_rss(2) = e0(2) + ewind(2) + edir(2)
  
        call compute_fastem(4,freq(ifreq),eia,sst,salinity,wsp,azimuth,emis_fm4)
        call compute_fastem(5,freq(ifreq),eia,sst,salinity,wsp,azimuth,emis_fm5)
        call compute_fastem(6,freq(ifreq),eia,sst,salinity,wsp,azimuth,emis_fm6)


        drec.freq = freq(ifreq)
        drec.sst  = sst
        drec.wsp  = wsp
        drec.emis_rss = emis_rss(1:2)
        drec.emis_fm4 = emis_fm4(1:2)
        drec.emis_fm5 = emis_fm5(1:2)
        drec.emis_fm6 = emis_fm6(1:2)
        write(20,iostat=ios) drec
      enddo
    enddo
  enddo
  close(20)

  stop
  end
