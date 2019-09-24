program GPM_TbsimulatorV7

use GPM_TbsimulatorV7_variables
use GPM_TbsimulatorV7_IO
use GPM_TbsimulatorV7_procedures


implicit none

!--- this routine will read a GPM V6 profile file (from MIRS and CMB) output from the 
!--- mkstrm routine. It'll run the Kummerow Radiative transfer on all the pixels, both 
!--- raining and non-raining.
!---
!--- written by Dave Randel CSU/ATMOS
!--- 2/2019 - 9/2019
!
!--- misc variables

 integer           :: i,j,k,z

! character(len=5)  :: sat, sensor

!--- misc variables

 integer            :: n_args, iargc
 character(len=256) :: input_file, output_file
 logical            :: doFPA
  
!--- get input parameters

 n_args = iargc()
 if (n_args .eq. 4) then
    call getarg(1, satcode)
    call getarg(2, sensor)
    call getarg(3, input_file)
    call getarg(4, output_file)
 else
    write(*,*)'ERROR reading in arguments - 4 args needed'
    stop
 endif
 write(*,*)'input arguments'
 write(*,*)' satellite   = ',trim(satcode)
 write(*,*)' sensor      = ',trim(sensor)
 write(*,*)' input file  = ',trim(input_file)
 write(*,*)' output file = ',trim(output_file)


!--- Do FPA or just print out pixel aves

 doFPA = .false.

!--- read in the file which holds the combined profiles

 write(*,*)
 write(*,*)'calling read_strmfile'
 call GPM_read_strmfile(input_file)
 write(*,*)' finished reading strmfile'

!--- call the Tb simulator routine 

 write(*,*)'  calling GPM_compute_Tbs'
 call GPM_compute_Tbs

!--- Footprint average the Tbs, calculated at each DPR pixel into Sensor footprint

 if(doFPA) then
    write(*,*)
    write(*,*)' calling GPM_compute_fpa (footprint ave routine) '
    call GPM_compute_fpa
 endif

!--- compute simple pixel aves (if NO FPA) over the whole orbit

 if( .not. doFPA) then
    write(*,*)'  calling GPM_compute_aves'
    call GPM_compute_aves
 endif
 
!--- write out the results

 write(*,*)' writing simulated Tbs to : ', trim(output_file)
 write(*,*)' write out FPA? : ',doFPA
 
 call GPM_write_simfile(output_file, doFPA)


write(*,*)'successful completion'
stop
end program GPM_TbsimulatorV7

