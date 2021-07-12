!> module for imRK4 routines
module mod_imRK4
  use mod_types
  implicit none
  complex(kind=dp), dimension(4)    :: Z_k 

contains
  ! subroutine to find the vectors Z_ki
  subroutine iterate_Zki(Yn, t, Yn_new, h_new, ERR, h_old, ERR_old, counter, diff, k)
    use mod_types
    use mod_constants
    use mod_imRK4_parameters, only: step, dim, hdim, sc_tol, ATOL, RTOL, delta, alpha, beta
    use mod_math
    use mod_matrices
    implicit none
    ! define the variables: 
    real(kind=dp)                     , intent(in)    :: t
    complex(kind=dp), dimension(hdim) , intent(inout) :: Yn
    complex(kind=dp), dimension(hdim) , intent(out)   :: Yn_new
    real(kind=dp)   ,                   intent(out)   :: h_new,ERR, diff
    real(kind=dp)   ,                   intent(in)    :: h_old,ERR_old
    integer,                            intent(out)   :: k

    integer                                           :: i, j, s, counter
    real(kind=dp)                                     :: h_new_1, h_new_2, error, norm_k_old, norm_k_new, theta_k, eta_k, tol, ERR_i
    complex(kind=dp), dimension(hdim)                 :: F, Yn_hat
    complex(kind=dp), dimension(dim)                  :: deltaZ_k_old, deltaZ_k_new, deltaZ_k  
    complex(kind=dp), dimension(dim,dim)              :: M2n
    complex(kind=dp), dimension(dim)                  :: Fun 
    real(kind=dp),    dimension(hdim)                 :: Eta

    ! z = (g - Y)
    ! z = 0.0 is the solution
    ! Initial value of z_k (Taylor expansion point) 
    Z_k   = 0.d0
    k     = 0
    error = 1.d0
    !
    ! Solve the linear system of equations of size 2 * n ,
    ! where n is the dimension of the hamiltonian:
    !   (I- hA\otimes J)  \Delta z^k = -z^k+ h(A\otimes I)F(z^k)
    !   ----------------               -------------------------
    !        M2n                               Fun
    !
    sc_loop: do while (error >= sc_tol) 
      ! Build M2n matrix
      call M_2n(Yn, t, step, M2n)
      ! Build Fun vector
      call Fsystem(Yn, t, step, Z_k, Fun)
      ! Solve the Ax=b linear equation
      call LS_solver(dim,M2n,Fun)

      deltaZ_k  = Fun
      
      ! Save old deltaZ(k) to deltaZ(k-1)
      deltaZ_k_old= Z_k
      ! Update Z_k
      Z_k  = Z_k + deltaZ_k

      ! Obtaining new deltaZ(k)
      ! deltaZ_k_new = Z_k   
      if (k >= 1) then
        ! call vec_norm(deltaZ_k_old, dim,norm)
        norm_k_old = vec_norm(deltaZ_k_old, dim)
        ! call vec_norm(deltaZ_k_new, dim,norm)
        norm_k_new = vec_norm(Z_k, dim)
        theta_k    = norm_k_new/norm_k_old
        if (theta_k == 1.d0) then 
          error = 0.d0
        else
          eta_k = theta_k/(1.d0-theta_k)
          tol   = abs(norm_k_new - norm_k_old)/min(norm_k_new, norm_k_old)
          error = eta_k*norm_k_new/tol
        end if 
      else
        error = 1.d0
        k     = k+1
      end if 
    end do sc_loop
    ! write(*,*) k


    ! Find Yn_hat
    Yn_hat = Yn + d1_hat*Z_k(1:hdim)+ d2_hat*Z_k(hdim+1:dim)
    
    ! Find Yn
    Yn_new = Yn + d1*Z_k(1:hdim) + d2*Z_k(hdim+1:dim)
     
    ! Find the difference Eta in the two solutions Yn and Yn_hat (error)
    Eta = abs(Yn_new - Yn_hat)
    diff= abs( vec_norm(Yn_new - Yn_hat,hdim) )
    ! sum over the error components to get the scaled norm (ERR)
    ! ERR= sqrt[ 1\n Sum_i{ abs(Yn_new - Yn_hat)/( ATOL + RTOL * Yn(i) ) }^2 ]
    ERR = 0.d0
    ERR_loop: do i=1, hdim
      ! Add the change to TITAN
      ERR_i = ( Eta(i)/(ATOL+ RTOL* max(abs(Yn(i)), abs(Yn_new(i)) ) ) )**2 
      ERR = ERR_i + ERR
    end do ERR_loop
    ERR = sqrt(ERR/hdim)
    
    ! Find the new step size h_new, step= h_used
    ! h_new = delta * h_used / (ERR)^(1/p+1) where, 
    ! p = 2*s, delta is some saftey factor (e.g. 0.9d0)
    ! delta = 0.9 * (2*K_max + 1)/ (2*K_max + NEWT), where NEWT is the number of newton iterations
    ! In our case NEWT and K_max are equal to 1! 
    s = 2.0 ! get s from the method 
    ! h_new = delta * step * (ERR)**(-1.0/(2.0*s)) !! simpler formula
    ! set a value of ERR_old at the first step 

    if (counter /= 0) then
      !! Gustafsson (1994) formula: 
      h_new_1 = delta * step * (ERR)**(-1.0/(2.0*s)) * (step/h_old) * (ERR_old/ERR)**(1.0/(2.0*s))
      h_new_2 = delta * step * (ERR)**(-1.0/(2.0*s))
      h_new = min(h_new_1, h_new_2)
    else 
      h_new = delta * step * (ERR)**(-1/(2.0*s))
      
    end if 

  end subroutine iterate_Zki


 ! Subroutine to build the matrix M_2n
  subroutine M_2n(Yn, t, step, M2n)
    use mod_types
    use mod_constants
    use mod_imRK4_parameters, only: dim, hdim, linear 
    use mod_hamiltonian
    use mod_matrices
    use mod_math
    implicit none
    integer                                            :: i,j
    complex(kind=dp), dimension(hdim),    intent(in)   :: Yn
    real(kind=dp),                        intent(in)   :: t, step
    complex(kind=dp), dimension(dim,dim), intent(out)  :: M2n

    integer                                            :: nax, nay, nbx, nby
    complex(kind=dp), dimension(dim,dim)               :: Kprod(dim,dim)

    ! choose linear or non-linear systems
    if (linear) then
      call build_linear_Jacobian(Yn, t)
    else
      call build_non_linear_Jacobian(Yn, t)
    end if

    call KronProd(hdim,hdim,hdim,hdim,A,Jacobian,Kprod)
    M2n = id - step * Kprod
  end subroutine M_2n

  ! subroutine f(Z_k) that builds the right side of the linear system.
  subroutine Fsystem(Yn, t, step, Z_k, Fun) 
    use mod_types
    use mod_matrices
    use mod_imRK4_parameters, only: dim, hdim
    implicit none
    real(kind=dp),                     intent(in)  :: t, step    
    complex(kind=dp), dimension(hdim), intent(in)  :: Yn
    complex(kind=dp), dimension(dim),   intent(in)  :: Z_k
    complex(kind=dp), dimension(dim),  intent(out) :: Fun 
    complex(kind=dp), dimension(hdim)              :: F1, F2   

    call Build_Function(Yn, t, step, Z_k, F1, F2)  
    Fun = [ F1, F2 ]

    Fun = -Z_k + step * matmul(M1,Fun)
  end subroutine Fsystem

  ! subroutine to build the right side of the differential equation H*psi.
  subroutine Build_Function(Yn, t, step, Z_k, F1, F2) 
    use mod_types
    use mod_constants,   only: ci
    use mod_hamiltonian
    use mod_imRK4_parameters, only: hdim, dim, linear 
    use mod_matrices,   only: c1, c2
    implicit none
    real(kind=dp),                     intent(in)  :: t, step
    complex(kind=dp), dimension(hdim), intent(in)  :: Yn
    complex(kind=dp), dimension(dim), intent(in)   :: Z_k
    complex(kind=dp), dimension(hdim), intent(out) :: F1, F2   

    real(kind=dp)                                  :: t1, t2

    t1 =  t + step * c1
    ! choose linear or non-linear systems
    if (linear) then
      call build_linear_hamiltonian(Yn, t1)
    else
      call build_non_linear_hamiltonian(Yn, t1)
    end if
    F1 = -ci * matmul(hamilt,Z_k(     1:hdim)+Yn)


    t2=  t + step * c2
    ! choose linear or non-linear systems
    if (linear) then
      call build_linear_hamiltonian(Yn, t2)
    else
      call build_non_linear_hamiltonian(Yn, t2)
    end if
    F2= -ci * matmul(hamilt,Z_k(hdim+1: dim)+Yn)


  end subroutine Build_Function


end module mod_imRK4

