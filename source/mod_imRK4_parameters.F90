!> module for system parameters.
module mod_imRK4_parameters
  use mod_types
  use mod_constants
  implicit none
  character(len=*), parameter    :: input_file = "imRK4_input_parameters"
  real(kind=dp)                  :: w0, w1, w, lambda_p, lambda_o, omega, integration_time, step, sc_tol, delta, ATOL, RTOL, alpha, beta
  integer                        :: N, time , dim, hdim
  complex(kind=dp), dimension(2) :: Yn_initial ! the initial eigenstate of the system.
  logical                        :: linear

  contains
  
  subroutine initialize()
    use mod_matrices, only: A, build_identity
    implicit none
    integer              :: dim_I
    integer              :: ios

    !! Error variable (for I/O)
    ! read parameters w0, w1, w, integration_time, Yn_initial.
    namelist /parameters/ w0 , w1, w, lambda_p, lambda_o, N, integration_time, Yn_initial, linear,sc_tol, delta, ATOL, RTOL, alpha, beta
    open(unit=10, file=input_file, status='old')
    read(unit=10, nml=parameters,  iostat=ios)
    if(ios /= 0) then
      write(*,"('Error reading ""parameter"" list from ',a)") input_file
      stop
    end if

    close(unit=10)

  ! Writing the parameters to output file
    open(  unit=11,  file="output_file" )
    write( unit=11,  fmt="('Input parameters:')" )
    write( unit=11,  nml=parameters )
    close( unit=11 )

    omega= sqrt(w1**2+ (w0-0.5*w)**2) 
    step= (2*pi)/(N*max(w, omega))
    time= int(integration_time/step)
write(*,*) "w0 = ", w0
write(*,*) "w1 = ", w1
write(*,*) "w  = ", w
write(*,*) "omega = ", omega
write(*,*) "step = ", step
    ! get the dimensions
    hdim = size(Yn_initial)
    dim  = 2*hdim

    ! Building identity
    dim_I = size(A,1)*size(Yn_initial)
    call build_identity(dim_I)

    

  end subroutine initialize

end module mod_imRK4_parameters
