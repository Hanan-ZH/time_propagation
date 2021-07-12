!>>>> The main program.
program main
  use mod_types
  use mod_math
  use mod_imRK4_parameters
  use mod_imRK4
  implicit none
  character(len=6)               :: output_file = "output"
  integer                        :: i, counter, k
  real(kind=dp)                  :: t, h_new, ERR, h_old, ERR_old, diff
  real(kind=dp)                  :: exp_x,  exp_y,  exp_z
  complex(kind=dp), dimension(2) :: Yn, Yn_new

  call initialize()
  Yn = Yn_initial
  open(unit=11,file=trim(output_file), status= 'replace')
  write(unit=11,fmt=*) '#      Time      ', '        M_x       ', '        M_y       ', '        M_z       ', '      M      ', '      norm_psi     '
  
  t = 0.d0
  i = 0
  do while (i<2000) !(t <= integration_time)
    t = t + step
    counter = 0
    do
      call iterate_Zki(Yn, t, Yn_new, h_new, ERR, h_old, ERR_old, counter, diff, k)
      ! the condition (h_new < delta* step) is equivalent to the condition (ERR > 1)
      ! if ( h_new < 0.9 * step) then
         ! do while ( h_new < step) 
      ! save ERR from iterate_Zki to ERR_old
      ERR_old = ERR
      if ( ERR > 1.d0) then
        ! repeat the calculation using h_new
        t = t - step + h_new
        h_old = step
        step = h_new
        write(*,*) t, step, ERR
      else
        Yn = Yn_new
        step = h_new
        exit
      end if
    ! this condition seems to mantain a small step size
      !else if (10*step <= h_new <= 1.2*step) then
        !step = step
    end do

    write(*,*)  "Accepted", t, step, ERR, diff, k

    call spin_exp(Yn, exp_x, exp_y, exp_z)

    write(unit=11,fmt="(6(es16.9,2x))") t*6.582d-7 , exp_x, exp_y, exp_z, (exp_x**2 + exp_y**2 +exp_z**2),vec_norm(Yn, hdim)

    !write(unit=11,fmt="(6(es16.9,2x))") h_new, ERR
    !write(*,*) 't = ', t*6.582d-7, 'm = ', 'norm = ', (exp_x**2 + exp_y**2 +exp_z**2), vec_norm(Yn, hdim)
    counter = counter + 1
    i = i + 1
  end do 

  close(unit=11)

  
end program main

