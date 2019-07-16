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

  !test if matrix is symmetric
  DO j=0,ndim-1
    DO i=0,j-1
      IF (ABS(l(i,j) - l(j,i)) .GT. symtol ) THEN
        sym = .FALSE. 
      END IF
    END DO
  END DO

  WRITE(*,*) "Matrix is sym: ", sym

  IF (sym) THEN
    q = 0
    qi = 0
    DO j=0,ndim-1
      q(j,j) = k(j) + l(j,j) 
      DO i=j+1,ndim-1
        q(i,j) = l(i,j)
      END DO
    END DO
    CALL linal_sym_diag(ndim,q(0:ndim-1,0:ndim-1),w(0:ndim-1),error)
    IF (error .NE. 0) THEN
      WRITE(*,*) "nco_cHO : Error out of linal_sym_diag",error
      RETURN
    END IF
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

END MODULE nco
