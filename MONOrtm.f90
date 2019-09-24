MODULE MonoRTM

 use GPM_TbsimulatorV7_variables
       
 implicit  none

 contains

!-----------------------------------------------------------------------
     
      subroutine read_lut(lut_file, verbose)

!--- Read routine for MonoRTM look-up table

      character(len=100) :: lut_file
      integer, intent(in), optional :: verbose

      open(unit=10, file=trim(lut_file), access='stream')

      read(10) nfreq_lut
      allocate(freq_lut(nfreq_lut))
      read(10) freq_lut

      read(10) nchan_lut
      allocate(ifreq_lut(nchan_lut))
      allocate(ipol_lut(nchan_lut))
      read(10) ifreq_lut
      read(10) ipol_lut

      read(10) npres_lut
      allocate(pres_lut(npres_lut))
      read(10) pres_lut

      read(10) ntemp_lut
      allocate(temp_lut(ntemp_lut))
      read(10) temp_lut

      read(10) nrmix_lut
      allocate(rmix_lut(nrmix_lut))
      read(10) rmix_lut

!      write(6,*) nfreq_lut, nchan_lut, npres_lut, ntemp_lut, nrmix_lut
!      pause

      allocate(kabs_lut(npres_lut,ntemp_lut,nrmix_lut,nfreq_lut))
      read(10) kabs_lut
      close(10)
      
      if (present(verbose)) then
        if (verbose) then
          write(6,'("nchan = ",i3)') nchan_lut
          write(6,'("ifreq = ",20(i3))') ifreq_lut
          write(6,'("ipol  = ",20(i3))') ipol_lut
          write(6,'("npres = ",i5)') npres_lut
          write(6,'("ntemp = ",i5)') ntemp_lut
          write(6,'("nrmix = ",i5)') nrmix_lut
        endif
      endif

      return
      end subroutine read_lut

!-----------------------------------------------------------------------
     
      subroutine monortm_lut(ichan, pavg, tavg, ravg, kext)

!     Look-up table version of MonoRTM

      integer :: ichan
      real    :: pavg
      real    :: tavg
      real    :: ravg
      real    :: kext

      integer :: freq_index
      integer :: i,j,k,n
      integer :: p1,p2
      integer :: t1,t2
      integer :: r1,r2
      real    :: pw1,pw2
      real    :: tw1,tw2
      real    :: rw1,rw2

      p1 = 1
      do n=1,npres_lut
        if (pres_lut(n) .gt. pavg) then
          p2 = n
          exit
        else
          p1 = n
          p2 = n
        endif
      enddo
      if (p1 .eq. p2) then
        pw1 = 0.5
        pw2 = 0.5
      else
        pw1 = (pres_lut(p2) - pavg) / (pres_lut(p2) - pres_lut(p1))
        pw2 = (pavg - pres_lut(p1)) / (pres_lut(p2) - pres_lut(p1))
      endif

      t1 = 1
      do n=1,ntemp_lut
        if (temp_lut(n) .gt. tavg) then
          t2 = n
          exit
        else
          t1 = n
          t2 = n
        endif
      enddo
      if (t1 .eq. t2) then
        tw1 = 0.5
        tw2 = 0.5
      else
        tw1 = (temp_lut(t2) - tavg) / (temp_lut(t2) - temp_lut(t1))
        tw2 = (tavg - temp_lut(t1)) / (temp_lut(t2) - temp_lut(t1))
      endif

      r1 = 1
      do n=1,nrmix_lut
        if (rmix_lut(n) .gt. ravg) then
          r2 = n
          exit
        else
          r1 = n
          r2 = n
        endif
      enddo
      if (r1 .eq. r2) then
        rw1 = 0.5
        rw2 = 0.5
      else
        rw1 = (rmix_lut(r2) - ravg) / (rmix_lut(r2) - rmix_lut(r1))
        rw2 = (ravg - rmix_lut(r1)) / (rmix_lut(r2) - rmix_lut(r1))
      endif

      freq_index = ifreq_lut(ichan)
      kext = (rw1 * tw1 * pw1 * kabs_lut(p1,t1,r1,freq_index)) +   &
             (rw1 * tw1 * pw2 * kabs_lut(p2,t1,r1,freq_index)) +   & 
             (rw1 * tw2 * pw1 * kabs_lut(p1,t2,r1,freq_index)) +   &
             (rw1 * tw2 * pw2 * kabs_lut(p2,t2,r1,freq_index)) +   &
             (rw2 * tw1 * pw1 * kabs_lut(p1,t1,r2,freq_index)) +   &
             (rw2 * tw1 * pw2 * kabs_lut(p2,t1,r2,freq_index)) +   &
             (rw2 * tw2 * pw1 * kabs_lut(p1,t2,r2,freq_index)) +   &
             (rw2 * tw2 * pw2 * kabs_lut(p2,t2,r2,freq_index))

      return
      end subroutine monortm_lut

!-----------------------------------------------------------------------
     
END MODULE MonoRTM
