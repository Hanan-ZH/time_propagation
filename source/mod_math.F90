!> module for mathematical operations(kronecker product, norm of a vector,
!> expec_vals of pauli matricies,  eigenvalue/vector solver).
module mod_math
  use mod_types                                           

contains
  ! subroutine for kronecker product.
  subroutine KronProd(nax,nay,nbx,nby,B1,B2,Kprod)
    use mod_types
    implicit none
    integer,          intent(in)  :: nax, nay, nbx, nby
    integer                       :: i, j, p, q, l, m                                              
    complex(kind=dp), intent(in)  :: B1(nax,nay), B2(nbx,nby) 
    complex(kind=dp), intent(out) :: Kprod(nax*nbx,nay*nby) 
    do i = 1,nax
      do j = 1,nay
        l=(i-1)*nbx + 1
        m=l+nbx-1
        p=(j-1)*nby + 1
        q=p+nby-1
        Kprod(l:m,p:q) = B1(i,j)*B2
      end do
    end do
  end subroutine KronProd

  ! subroutine for vector norm.
  function vec_norm(v, dim_v)
    use mod_types
    implicit none
    integer,                            intent(in)  :: dim_v ! vector dimension
    complex(kind=dp), dimension(dim_v), intent(in)  :: v ! vector v
    real(kind=dp)                                   :: vec_norm
    real(kind=dp)                                   :: sum
    integer                                         :: i
    sum= 0.d0
    do i = 1, dim_v
      sum= sum+ v(i)*conjg(v(i))
    end do
    vec_norm= sqrt(sum)
  end function vec_norm

  ! subroutine for expectation value of a vector in spin matricies
  subroutine spin_exp(Yn, exp_x, exp_y, exp_z)
    use mod_types                         
    use mod_matrices                       
    implicit none
    complex(kind=dp), intent(in), dimension(2)  :: Yn
    real(kind=dp)   , intent(out)               :: exp_x,  exp_y,  exp_z
    real(kind=dp)                               :: exp_x_0, exp_y_0, exp_z_0
    integer                                     :: j, k, n
    n= size(Yn)
    exp_x= 0.d0
    exp_y= 0.d0
    exp_z= 0.d0
    do j=1, n
      do k=1, n
        exp_x_0= real(CONJG(Yn(j))* sigma_x(j,k)*Yn(k))
        exp_x= exp_x + exp_x_0

        exp_y_0= real(CONJG(Yn(j))* sigma_y(j,k)*Yn(k))
        exp_y= exp_y + exp_y_0

        exp_z_0= real(CONJG(Yn(j))* sigma_z(j,k)*Yn(k))
        exp_z= exp_z + exp_z_0
      end do
    end do
  end subroutine spin_exp

     
! interface to the LAPACK eigensolver
  subroutine eigensolver(n,A,w)

    implicit none
    
    integer,          intent(in)    :: n
    double complex,   intent(inout) :: A(n,n)
    double precision, intent(out)   :: w(n)
    ! Workspace variables
    integer                       :: lwork, info
    double complex, allocatable   :: work(:)
    double precision, allocatable :: rwork(:)
    
    lwork = 3*n
    allocate(work(lwork),rwork(lwork))
    call zheev('V','U',n,A,n,w,work,lwork,rwork,info)
    if (info /= 0) write(*,'("eigensolver: failure in zheev")')
    deallocate(work,rwork)
    
    end subroutine eigensolver



    ! interface to the LAPACK linear system solver A*X=B
  subroutine LS_solver(n,A,b)
  
    implicit none
    
    integer,          intent(in)    :: n
    double complex,   intent(in)    :: A(n,n)
    double complex, intent(out)     :: b(n)
    ! Workspace variables
    integer                         :: nrhs, lda, ldb, ldx, info, iter
    integer, allocatable            :: ipiv(:)
    complex, allocatable            :: swork(:)
    double complex, allocatable     :: work(:), X(:)
    double precision, allocatable   :: rwork(:)

    lda= n
    ldb=n
    ldx= n
    nrhs=1

    allocate(rwork(n),swork(n*(n+1)),work(n), ipiv(n), X(n) )
    ! call zcgesv( n, nrhs, a, lda, ipiv, b, ldb, X, ldx, work, swork, rwork, iter, info )
    ! b = X
    call zgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
    
    if (info /= 0) write(*,'("eigensolver: failure in zheev")')
    deallocate(work,swork,rwork,ipiv,X)
    
    end subroutine LS_solver

end module mod_math