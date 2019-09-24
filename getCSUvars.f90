subroutine getcsuvars(temp,pressj,kexttotj,salbtotj,asymtotj,tbj,ctbj,rainj,snowj,nwj,&
     nlevsj,nscansj,npixsj,nchansj)
 Use GPM_TbsimulatorV7_variables
 integer :: nlevsj, nscansj, npixsj, nchansj
 real :: temp(0:nlevsj,npixsj,nscansj),pressj(0:nlevsj,npixsj,nscansj)
 real :: kexttotj(npixsj,nscansj,nlevsj),salbtotj(npixsj,nscansj,nlevsj),&
      asymtotj(npixsj,nscansj,nlevsj)
 real :: tbj(nchansj,npixsj,nscansj),ctbj(nchansj,npixsj,nscansj),&
      rainj(nlevsj,npixsj,nscansj),snowj(nlevsj,npixsj,nscansj),&
      nwj(nlevsj,npixsj,nscansj)
 temp=tempprof
 pressj=pressprof
 kexttotj=kexttot
 salbtotj=salbtot
 asymtotj=asymtot
 nwj=nw
 rainj=rainwater
 snowj=snowwater 
 tbj=dtbssim
 ctbj=ctbssim
end subroutine getcsuvars
