!Code for normal coordinates
MODULE nco
  USE linal

CONTAINS
!------------------------------------------------------------
! nco_cHO
!   - compute normal coordinates of a system of linearly 
!     coupled harmonic oscillators
!------------------------------------------------------------
! ndim          : int, number of dimensions
! k             : 1D real*8, force constants
! l             : 2D real*8, couplings
! w             : 1D real*8, eigenvalues
! q             : 2D real*8, array of normal coordinates
! qi            : 2D real*8, inverse array of normal coordinates
! error         : int, error code

SUBROUTINE nco_cHO(ndim,k,l,w,q,qi,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1), INTENT(INOUT) :: q,qi
  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1), INTENT(IN) :: l
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(INOUT) :: w
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(IN) :: k
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim

  REAL(KIND=8) :: symtol
  INTEGER :: i,j
  LOGICAL :: sym

  error = 0
  symtol = 1.0D-15
  sym = .TRUE.
  WRITE(*,*) "nco_cHO : called"
  WRITE(*,*) "Determining renormalized coordinates"

  !test if matrix is symmetric
  DO j=0,ndim-1
    DO i=0,j-1
      IF (ABS(l(i,j) - l(j,i)) .GT. symtol ) THEN
        sym = .FALSE. 
      END IF
    END DO
  END DO

  IF (sym) THEN
    WRITE(*,*) "Matrix is symmetric"
  ELSE
    WRITE(*,*) "Matrix is not symmetric"
  END IF
  WRITE(*,*) 

  IF (sym) THEN
    q = 0
    qi = 0
    DO j=0,ndim-1
      q(j,j) = k(j) + l(j,j) 
      DO i=j+1,ndim-1
        q(i,j) = l(i,j)
      END DO
    END DO
    WRITE(*,*) "The undiagonalized matrix is"
    DO i=0,ndim-1
      WRITE(*,*) q(i,0:ndim-1)
    END DO
    CALL linal_sym_diag(ndim,q(0:ndim-1,0:ndim-1),w(0:ndim-1),error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "nco_cHO : Error out of linal_sym_diag",error
      RETURN
    END IF
    !check this
    WRITE(*,*) 
    WRITE(*,*) "This should return the original matrix"
    DO i=0,ndim-1 
      qi(i,i) = w(i)
    END DO
    qi = MATMUL(q,MATMUL(qi,TRANSPOSE(q))) 
    DO i=0,ndim-1
      WRITE(*,*) qi(i,0:ndim-1) 
    END DO
    !CALL linal_invert(ndim,q(0:ndim-1,0:ndim-1),qi(0:ndim-1,0:ndim-1),error) 
    qi = TRANSPOSE(q(0:ndim-1,0:ndim-1))
    IF (error .NE. 0) THEN
      WRITE(*,*) "nco_cHO : Error out of linal_invert",error
      RETURN
    END IF
    w = SQRT(w)
  ELSE
    WRITE(*,*) "Sorry, not implemented yet"
    error = 1
    RETURN
  END IF


  WRITE(*,*)
  WRITE(*,*) "Renormalized Frequencies"
  WRITE(*,*) w(0:ndim-1)
  WRITE(*,*) 
  WRITE(*,*) "Renormalized coordinate matrix"
  DO i=0,ndim-1
    WRITE(*,*) q(i,0:ndim-1)
  END DO
  WRITE(*,*) 
  OPEN(file='normco.dat',unit=103,status='replace')
  DO i=0,ndim-1
    WRITE(103,*) "Renormalized Frequency", w(i)
    DO j=0,ndim-1
      WRITE(103,*) q(j,i)
    END DO
    WRITE(103,*)
  END DO
  CLOSE(unit=103)
  

END SUBROUTINE nco_cHO


!------------------------------------------------------------
! nco_general
!   - computes renormalized coordinates for the general case
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nHO           : int, number of harmonic oscillators
! qHO           : 1D int, ids of harmonic oscillators 
! HO            : 1D real*8, force constants of harmonic oscillators
! nMO           : int, number of Morse oscillators
! qMO           : 1D int, ids of morse oscillators
! MO            : 1D real*8, Morse oscillator [{k,De}:i]
! nl1           : int, number of order 1 terms
! ql1           : 1D int, ids of order 1 terms 
! l1            : 1D real*8, order 1 terms
! nl2           : int, nubmer of order 2 terms
! ql2           : 1D int, ids of order 2 terms
! l2            : 1D real*8, order 2 terms
! nl3           : int, number of order 3 terms
! ql3           : 1D int, ids of order 3 terms
! l3            : 1D real*8, order 3 terms 
! nl4           : int, number of order 4 terms
! ql4           : 1D int, ids of order 4 terms
! l4            : 1D real*8, order 4 terms 
! w             : 1D real*8, eigenvalues
! q             : 2D real*8, array of normal coordinates
! qi            : 2D real*8, inverse array of normal coordinates
! error         : int, exit code

SUBROUTINE nco_general(ndim,nHO,qHO,HO,nMO,qMO,MO,nl1,ql1,l1,nl2,ql2,l2,&
                       nl3,ql3,l3,nl4,ql4,l4,w,q,qi,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1), INTENT(INOUT) :: q,qi
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(INOUT) :: w
  REAL(KIND=8), DIMENSION(0:2*nMO-1), INTENT(IN) :: MO
  REAL(KIND=8), DIMENSION(0:nHO-1), INTENT(IN) :: HO
  REAL(KIND=8), DIMENSION(0:nl1-1), INTENT(IN) :: l1
  REAL(KIND=8), DIMENSION(0:nl2-1), INTENT(IN) :: l2
  REAL(KIND=8), DIMENSION(0:nl3-1), INTENT(IN) :: l3
  REAL(KIND=8), DIMENSION(0:nl4-1), INTENT(IN) :: l4
  INTEGER, DIMENSION(0:4*nl4-1), INTENT(IN) :: ql4
  INTEGER, DIMENSION(0:3*nl3-1), INTENT(IN) :: ql3
  INTEGER, DIMENSION(0:2*nl2-1), INTENT(IN) :: ql2
  INTEGER, DIMENSION(0:1*nl1-1), INTENT(IN) :: ql1
  INTEGER, DIMENSION(0:nMO-1), INTENT(IN) :: qMO
  INTEGER, DIMENSION(0:nHO-1), INTENT(IN) :: qHO
  INTEGER, INTENT(IN) :: ndim,nHO,nMO,nl1,nl2,nl3,nl4
  INTEGER, INTENT(INOUT) :: error

  INTEGER :: n,i,j,k,l

  error = 0

  WRITE(*,*)
  WRITE(*,*) "nco_general : called"
  WRITE(*,*) "Generating Renormalized Coordinates" 
  WRITE(*,*) "Terms higher than 2nd order will be evaluated at q = 0"
  WRITE(*,*) 
  WRITE(*,*) "Constructing Hessian Matrix at q = 0"

  q = 0.0D0
  qi = 0.0D0
  w = 0.0D0
  
  !Harmonic terms
  DO n=0,nHO-1
    i = qHO(n)
    q(i,i) = q(i,i) + HO(n)  
  END DO 

  !Morse terms, at q  = 0
  DO n=0,nMO-1
    i = qMO(n)
    q(i,i) = q(i,i) + 2.0D0*MO(2*n)**2.0D0*MO(2*n+1) !1st is alpha, 2nd is De
  END DO

  !Order 1 terms do not contribute
  
  !Order 2 terms contribute to i>j
  DO n=0,nl2-1
    i = MAX(ql2(2*n),ql2(2*n+1))
    j = MIN(ql2(2*n),ql2(2*n+1))
    q(i,j) = q(i,j) + l2(n)
  END DO
  
  !Order 3 terms do not contribute at q = 0

  !Order 4 terms do not contriute at q = 0

  WRITE(*,*) 
  WRITE(*,*) "The Hessian is"
  DO i=0,ndim-1
    WRITE(*,*) q(i,0:i) 
  END DO
  WRITE(*,*) 
  !diagonalize the hessian
  WRITE(*,*) "Diagonalizing Hessian" 
  CALL linal_sym_diag(ndim,q(0:ndim-1,0:ndim-1),w(0:ndim-1),error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "nco_general : Error out of linal_sym_diag",error
    RETURN
  END IF
  WRITE(*,*) 
  !print out the check
  WRITE(*,*) "This should return the original matrix"
  DO i=0,ndim-1 
    qi(i,i) = w(i)
  END DO
  qi = MATMUL(q,MATMUL(qi,TRANSPOSE(q))) 
  DO i=0,ndim-1
    WRITE(*,*) qi(i,0:ndim-1) 
  END DO
  !get the projection matrix
  qi = TRANSPOSE(q(0:ndim-1,0:ndim-1))
  IF (error .NE. 0) THEN
    WRITE(*,*) "nco_cHO : Error out of linal_invert",error
    RETURN
  END IF
  w = SQRT(w)

  WRITE(*,*)
  WRITE(*,*) "Renormalized Frequencies"
  WRITE(*,*) w(0:ndim-1)
  WRITE(*,*) 
  WRITE(*,*) "Renormalized coordinate matrix"
  DO i=0,ndim-1
    WRITE(*,*) q(i,0:ndim-1)
  END DO
  WRITE(*,*) 
  OPEN(file='normco.dat',unit=103,status='replace')
  DO i=0,ndim-1
    WRITE(103,*) "Renormalized Frequency", w(i)
    DO j=0,ndim-1
      WRITE(103,*) q(j,i)
    END DO
    WRITE(103,*)
  END DO

END SUBROUTINE nco_general

!------------------------------------------------------------

END MODULE nco
