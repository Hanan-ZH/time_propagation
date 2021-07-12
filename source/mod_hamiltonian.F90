!>  module contains variables and subroutine to build the hamiltonian.
module mod_hamiltonian
  use mod_types
  implicit none
  complex(kind=dp), dimension(2,2)  :: hamilt, Jacobian
   
contains


!-------------*** this is to build the linear Hamiltonian and Jacobian ***----------------------!


  ! subroutine to build the Jacobian
subroutine build_linear_Jacobian(Yn, t)
  use mod_types
  use mod_constants,  only: ci
  use mod_imRK4_parameters, only: w0, w1, w
  implicit none    
  complex(kind=dp), dimension(2), intent(in)  :: Yn 
  real(kind=dp),                  intent(in)  :: t

  call build_linear_hamiltonian(Yn, t)

  Jacobian= -ci*hamilt
  
end subroutine build_linear_Jacobian


  !  This subroutine builds the spin hamiltonian(t, Z_k)
  !  Hamiltonian(t, Z_k) = h_bar * (   w0 * sigma_z + w1 * ( cos(wt) * sigma_x + sin(wt) * sigma_y )   )
  subroutine build_linear_hamiltonian(Yn, t)
    use mod_constants
    use mod_matrices,                          only: sigma_x, sigma_y, sigma_z
    use mod_imRK4_parameters,                        only: hdim, w0, w1, w
    implicit none
    real(kind=dp),    intent(in)                  :: t
    complex(kind=dp), intent(in), dimension(hdim) :: Yn
    ! Building the Hamiltonian
    hamilt = w0*sigma_z + w1*(cos(w*t)*sigma_x + sin(w*t)*sigma_y)
    !hamilt = w0

    ! Checking if Hamiltonian is hermitian
    if( sum(abs(conjg(transpose(hamilt))-hamilt)) > 1.d-12 ) then
      write(*,"('Hamiltonian is not hermitian!')")
      stop
    end if
  end subroutine build_linear_hamiltonian


!-------------*** this is to build the non-linear Hamiltonian and non-linear Jacobian ***-------------------!



  ! Jacobian(t,Z_k) =  2 * non_linear_Hamiltonian - linear_Hamiltonian  
  !                 =  2 * hamilt_nl - h_bar * ( w0 * sigma_z + w1 * ( cos(wt) * sigma_x + sin(wt) * sigma_y ) )

  subroutine build_non_linear_Jacobian(Yn, t)
    use mod_types
    use mod_constants,                       only: ci
    use mod_matrices,                        only: sigma_x, sigma_y, sigma_z
    use mod_imRK4_parameters,                      only: hdim, w0, w1, w, lambda_p, lambda_o
    implicit none    
    complex(kind=dp), dimension(2), intent(in)  :: Yn 
    real(kind=dp),                  intent(in)  :: t

    call build_non_linear_hamiltonian(Yn, t)
  
    Jacobian= -ci * ( 2.d0 * hamilt - (   w0 * sigma_z + w1 * (cos(w*t) * sigma_x + sin(w*t) * sigma_y)  )  )


  end subroutine build_non_linear_Jacobian
  
  
    !> Hamiltonian(t, Z_k) = H1(t) + H2(t, Z_k)+ H3(t, Z_k)
    !> H1(t)      = h_bar ( w0 * sigma_z + w1 * ( cos(wt) * sigma_x + sin(wt) * sigma_y ) )
     
    !> H2(t, Z_k) = lambda_p * ( <sigma_x> * sigma_x + <sigma_y> * sigma_y + <sigma_z> * sigma_z )

    !> H3(t, Z_k) = lambda_o * {    ( w0 * <sigma_y>         - w1 * sinwt * <sigma_z> ) * sigma_x      
    !                            -  ( w0 * <sigma_x>         - w1 * coswt * <sigma_z> ) * sigma_y 
    !                            +  ( w1 * sinwt * <sigma_x> - w1 * coswt * <sigma_y> ) * sigma_z   }
                                                 
    subroutine build_non_linear_hamiltonian(Yn, t)
      use mod_constants
      use mod_matrices,                          only: sigma_x, sigma_y, sigma_z
      use mod_imRK4_parameters,                        only: hdim, w0, w1, w, lambda_p, lambda_o
      use mod_math,                              only: spin_exp
      implicit none
      real(kind=dp),    intent(in)                  :: t
      complex(kind=dp), intent(in), dimension(hdim) :: Yn
      complex(kind=dp), dimension(hdim, hdim)       :: hamilt_1, hamilt_2, hamilt_3
      real(kind=dp)                                 :: exp_x, exp_y, exp_z

      ! Building the Hamiltonian
      hamilt_1  = w0 * sigma_z + w1 * (cos(w*t) * sigma_x + sin(w*t) * sigma_y)

      call spin_exp(Yn, exp_x, exp_y, exp_z)

      hamilt_2  = lambda_p * ( exp_x * sigma_x + exp_y * sigma_y + exp_z * sigma_z )

      hamilt_3  = lambda_o * (  ( w0 * exp_y - w1 * sin(w*t) *  exp_z ) * sigma_x -  ( w0 * exp_x - w1 * cos(w*t) *  exp_z ) * sigma_y + ( w1 * sin(w*t) * exp_x - w1 * cos(w*t) *  exp_y ) * sigma_z )

      hamilt = hamilt_1 + hamilt_2 + hamilt_3

      ! Checking if Hamiltonian is hermitian
      if( sum(abs(conjg(transpose(hamilt))-hamilt)) > 1.d-12 ) then
        write(*,"('Hamiltonian is not hermitian!')")
        stop
      end if
    end subroutine build_non_linear_hamiltonian


end module mod_hamiltonian
