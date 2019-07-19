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
! R             : 1D real*8, positions to evalute at
! w             : 1D real*8, eigenvalues
! q             : 2D real*8, array of normal coordinates
! qi            : 2D real*8, inverse array of normal coordinates
! error         : int, exit code

SUBROUTINE nco_general(ndim,nHO,qHO,HO,nMO,qMO,MO,nl1,ql1,l1,nl2,ql2,l2,&
                       nl3,ql3,l3,nl4,ql4,l4,R,w,q,qi,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1), INTENT(INOUT) :: q,qi
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(INOUT) :: w
  REAL(KIND=8), DIMENSION(0:2*nMO-1), INTENT(IN) :: MO
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(IN) :: R
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

  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1) :: hHO,hMO,hl1,hl2,hl3,hl4
  INTEGER :: n,i,j,k,l

  error = 0
  q = 0.0D0
  qi = 0.0D0
  w = 0.0D0
  hHO = 0.0D0
  hMO = 0.0D0
  
  !Harmonic terms
  DO n=0,nHO-1
    i = qHO(n)
    hHO(i,i) = hHO(i,i) + HO(n)  
  END DO 

  !Morse terms
  ! MO(2*i) is ith alpha
  ! MO(2*i+1) is ith De
  DO n=0,nMO-1
    i = qMO(n)
    hMO(i,i) = hMO(i,i) + 2.0D0*MO(2*n)**2.0D0*MO(2*n+1)&
             *(2.0D0*EXP(-2.0D0*MO(2*n)*R(i)) - EXP(-1.0D0*MO(2*n)*R(i))) 
  END DO

  !Order 1 terms do not contribute

  !Order 2 terms
  CALL nco_hl2(ndim,nl2,ql2,l2,R,hl2)

  !Order 3 terms
  CALL nco_hl3(ndim,nl3,ql3,l3,R,hl3)
  
  !Order 4 terms
  CALL nco_hl4(ndim,nl4,ql4,l4,R,hl4)

  q = hHO + hMO + hl2 + hl3 + hl4

  !WRITE(*,*) "Before diag"
  !DO i=0,ndim-1
  !  WRITE(*,*) q(i,0:ndim-1)
  !END DO
  CALL linal_sym_diag(ndim,q(0:ndim-1,0:ndim-1),w(0:ndim-1),error)
  !CALL linal_diag_inplace(ndim,q,w,error)
  !WRITE(*,*) "After diag"
  !DO i=0,ndim-1
  !  WRITE(*,*) q(i,0:ndim-1)
  !END DO
  !WRITE(*,*) "Eigenvalues"
  !WRITE(*,*) w
  !STOP

  
  IF (error .NE. 0) THEN
    WRITE(*,*) "nco_general : Error out of linal_sym_diag",error
    RETURN
  END IF
  qi = TRANSPOSE(q(0:ndim-1,0:ndim-1))
  w = SQRT(w)

END SUBROUTINE nco_general

!------------------------------------------------------------
! nco_hl2
!    - calculates order 2 contributions to hessian in a 
!      complicated manner that reduces excess movement 
!      of memory
!    - the code would be more easy to read/implement if 
!      we did not care about dealing with the symmetry of
!      (i,j) and (j,i)
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nl2           : int, number of l2 terms
! ql2           : 1D int, ids of l2 terms
! l2            : 1D real*8, constants of l2 terms
! R             : 1D real*8, positions of dimensions
! hl2           : 2D real*8, 2D matrix of hessian to return

SUBROUTINE nco_hl2(ndim,nl2,ql2,l2,R,hl2)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1), INTENT(INOUT) :: hl2
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(IN) :: R
  REAL(KIND=8), DIMENSION(0:nl2-1), INTENT(IN) :: l2
  INTEGER, DIMENSION(0:2*nl2-1), INTENT(IN) :: ql2
  INTEGER, INTENT(IN) :: ndim,nl2

  INTEGER :: n,i,j,q1,q2

  hl2 = 0.0D0

  !Order 2 terms contribute to i>j
  DO n=0,nl2-1
    ! λ * q1 * q2 case
    q1 = ql2(2*n)
    q2 = ql2(2*n+1)
    IF (ql2(2*n).NE.ql2(2*n+1)) THEN
      i = MAX(q1,q2)
      j = MIN(q1,q2)
      hl2(i,j) = hl2(i,j) + l2(n)
    ! λ * q1 * q1 case
    ELSE
      i = q1 
      hl2(i,i) = hl2(i,i) + 2.0D0*l2(n)
    END IF
  END DO

END SUBROUTINE nco_hl2

!------------------------------------------------------------
! nco_hl3
!    - calculates order 3 contributions to hessian in a 
!      complicated manner that reduces excess movement 
!      of memory
!    - the code would be more easy to read/implement if 
!      we did not care about dealing with the symmetry of
!      (i,j) and (j,i)
!    - this code could actually be made even more complicated,
!      and even faster, by branching if statements...
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nl3           : int, number of l3 terms
! ql3           : 1D int, ids of l3 terms
! l3            : 1D real*8, constants of l3 terms
! R             : 1D real*8, positions of dimensions
! hl3           : 2D real*8, 2D matrix of hessian to return

SUBROUTINE nco_hl3(ndim,nl3,ql3,l3,R,hl3)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1), INTENT(INOUT) :: hl3
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(IN) :: R
  REAL(KIND=8), DIMENSION(0:nl3-1), INTENT(IN) :: l3
  INTEGER, DIMENSION(0:3*nl3-1), INTENT(IN) :: ql3
  INTEGER, INTENT(IN) :: ndim,nl3

  INTEGER :: n,i,j,k,q1,q2,q3

  hl3 = 0.0D0

  !Order 3 terms λ * q1 * q2 * q3  
  DO n=0,nl3-1
    ! λ * q1 * q2 * q3 case
    IF (ql3(3*n) .NE. ql3(3*n+1) .AND. &
        ql3(3*n+1) .NE. ql3(3*n+2) .AND. &
        ql3(3*n) .NE. ql3(3*n+2) ) THEN
      !q1,q2
      i = MAX(ql3(3*n),ql3(3*n+1))
      j = MIN(ql3(3*n),ql3(3*n+1))
      k = ql3(3*n+2)
      hl3(i,j) = hl3(i,j) + l3(n)*R(k) !dq1,dq2
      !q1,q3
      i = MAX(ql3(3*n),ql3(3*n+2))
      j = MIN(ql3(3*n),ql3(3*n+2))
      k = ql3(3*n+1)
      hl3(i,j) = hl3(i,j) + l3(n)*R(k) !dq1,dq3
      !q2,q3
      i = MAX(ql3(3*n+1),ql3(3*n+2))
      j = MIN(ql3(3*n+1),ql3(3*n+2))
      k = ql3(3*n)
      hl3(i,j) = hl3(i,j) + l3(n)*R(k) !dq2,dq3
    ! λ * q1 * q1 * q2
    ELSE IF (ql3(3*n) .EQ. ql3(3*n+1) .AND. &
             ql3(3*n+1) .NE. ql3(3*n+2)) THEN
      i = MAX(ql3(3*n),ql3(3*n+2))
      j = MIN(ql3(3*n),ql3(3*n+2))
      hl3(i,i) = hl3(i,i) * 2.0D0*l3(n)*R(j) !dq1,dq1
      hl3(i,j) = hl3(i,j) * 2.0D0*l3(n)*R(i) !dq1,dq2
    ! λ * q1 * q2 * q2
    ELSE IF (ql3(3*n) .NE. ql3(3*n+1) .AND. &
             ql3(3*n+1) .EQ. ql3(3*n+2)) THEN
      i = MAX(ql3(3*n),ql3(3*n+2))
      j = MIN(ql3(3*n),ql3(3*n+2))
      hl3(i,i) = hl3(i,i) * 2.0D0*l3(n)*R(j) !dq2,dq2
      hl3(i,j) = hl3(i,j) * 2.0D0*l3(n)*R(i) !dq1,dq2
    ! λ * q1 * q2 * q1
    ELSE IF (ql3(3*n+1) .NE. ql3(3*n+2) .AND. &
             ql3(3*n) .EQ. ql3(3*n+2)) THEN
      i = MAX(ql3(3*n),ql3(3*n+1))
      j = MIN(ql3(3*n),ql3(3*n+1))
      hl3(i,i) = hl3(i,i) * 2.0D0*l3(n)*R(j) !dq3,dq3
      hl3(i,j) = hl3(i,j) * 2.0D0*l3(n)*R(i) !dq1,dq2
    ! λ * q1 * q1 * q case
    ELSE 
      i = ql3(3*n) 
      hl3(i,i) = hl3(i,i) + 6.0D0*l3(n)  
    END IF


    !Alternate method
    !q1 = ql3(3*n)
    !q2 = ql3(3*n+1)
    !q3 = ql3(3*n+2)
    !
    !IF (q1 .GE. q2) hl3(q1,q2) = hl3(q1,q2) + l3(n)*R(q3)
    !IF (q1 .LE. q2) hl3(q2,q1) = hl3(q2,q1) + l3(n)*R(q3)
    !
    !IF (q1 .GE. q3) hl3(q1,q3) = hl3(q1,q3) + l3(n)*R(q2)
    !IF (q1 .LE. q3) hl3(q3,q1) = hl3(q3,q1) + l3(n)*R(q2)
    ! 
    !IF (q2 .GE. q3) hl3(q2,q3) = hl3(q2,q3) + l3(n)*R(q1)
    !IF (q2 .LE. q3) hl3(q3,q2) = hl3(q3,q2) + l3(n)*R(q1)
  END DO

END SUBROUTINE nco_hl3

!------------------------------------------------------------
! nco_hl4
!    - calculates order 4 contributions to hessian
!    - this code could actually be made even faster with a 
!       branching if statement structure
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nl4           : int, number of l4 terms
! ql4           : 1D int, ids of l4 terms
! l4            : 1D real*8, constants of l4 terms
! R             : 1D real*8, positions of dimensions
! hl4           : 2D real*8, 2D matrix of hessian to return

SUBROUTINE nco_hl4(ndim,nl4,ql4,l4,R,hl4)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1), INTENT(INOUT) :: hl4
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(IN) :: R
  REAL(KIND=8), DIMENSION(0:nl4-1), INTENT(IN) :: l4
  INTEGER, DIMENSION(0:4*nl4-1), INTENT(IN) :: ql4
  INTEGER, INTENT(IN) :: ndim,nl4

  INTEGER :: n,i,j,q1,q2,q3,q4

  hl4 = 0.0D0

  ! λ * q1 * q2 * q3 * q4
  DO n=0,nl4-1
    q1 = ql4(4*n)
    q2 = ql4(4*n+1)
    q3 = ql4(4*n+2)
    q4 = ql4(4*n+3)

    IF (q1 .GE. q2) hl4(q1,q2) = hl4(q1,q2) + l4(n)*R(q3)*R(q4)  
    IF (q1 .LE. q2) hl4(q2,q1) = hl4(q2,q1) + l4(n)*R(q3)*R(q4)  

    IF (q1 .GE. q3) hl4(q1,q3) = hl4(q1,q3) + l4(n)*R(q2)*R(q4)
    IF (q1 .LE. q3) hl4(q3,q1) = hl4(q3,q1) + l4(n)*R(q2)*R(q4)

    IF (q1 .GE. q4) hl4(q1,q4) = hl4(q1,q4) + l4(n)*R(q2)*R(q3)
    IF (q1 .LE. q4) hl4(q4,q1) = hl4(q4,q1) + l4(n)*R(q2)*R(q3)

    IF (q2 .GE. q3) hl4(q2,q3) = hl4(q2,q3) + l4(n)*R(q1)*R(q4)
    IF (q2 .LE. q3) hl4(q3,q2) = hl4(q3,q2) + l4(n)*R(q1)*R(q4)

    IF (q2 .GE. q4) hl4(q2,q4) = hl4(q2,q4) + l4(n)*R(q1)*R(q3)
    IF (q2 .LE. q4) hl4(q4,q2) = hl4(q4,q2) + l4(n)*R(q1)*R(q3)

    IF (q3 .GE. q4) hl4(q3,q4) = hl4(q3,q4) + l4(n)*R(q1)*R(q2)
    IF (q3 .LE. q4) hl4(q4,q3) = hl4(q4,q3) + l4(n)*R(q1)*R(q2)


    !IF (q1 .NE. q2 .AND. q1 .NE. q3 .AND. q1 .NE. q4 .AND. &
    !    q2 .NE. q3 .AND. q2 .NE. q4 .AND. q3 .NE. q4) THEN
    !  i = MAX(q1,q2)
    !  j = MIN(q1,q2)
    !  q(i,j) = q(i,j) + l4(n)*R(q3)*R(q4) !dq1,dq2
    !  i = MAX(q1,q3)
    !  j = MIN(q1,q3)
    !  q(i,j) = q(i,j) + l4(n)*R(q2)*R(q4) !dq1,dq3
    !  i = MAX(q1,q4)
    !  j = MIN(q1,q4)
    !  q(i,j) = q(i,j) + l4(n)*R(q2)*R(q3) !dq1,dq4
    !  i = MAX(q2,q3)
    !  j = MIN(q2,q3)
    !  q(i,j) = q(i,j) + l4(n)*R(q1)*R(q4) !dq2,dq3
    !  i = MAX(q2,q4)
    !  j = MIN(q2,q4)
    !  q(i,j) = q(i,j) + l4(n)*R(q1)*R(q3) !dq2,dq4
    !  i = MAX(q3,q4)
    !  j = MIN(q3,q4)
    !  q(i,j) = q(i,j) + l4(n)*R(q1)*R(q2) !dq3,dq4
    ! λ * q1 * q1 * q3 * q4
    !ELSE IF (q1 .EQ. q2 .AND. q1 .NE. q3 .AND. q1 .NE. q4 .AND. &
    !    q3 .NE. q4) THEN
    !  q(q1,q1) = q(q1,q1) + 4.0D0*l4(n)*R(q3)*R(q4) !dq1,qd1
    !  i = MAX(q1,q3)
    !  j = MIN(q1,q3)
    !  q(i,j) = q(i,i) + 2.0D0*l4(n)*R(q1)*R(q4)     !dq1,dq3
    !  i = MAX(q1,q4)
    !  j = MIN(q1,q4)
    !  q(i,j) = q(i,i) + 2.0D0*l4(n)*R(q1)*R(q3)     !dq1,dq4
    !  i = MAX(q3,q4)
    !  j = MIN(q3,q4)
    !  q(i,j) = q(i,j) + l4(n)*R(q1)**2.0D0          !dq3,dq4
    !λ * q1 * q2 * q1 * q4
    !ELSE IF (q1 .EQ. q3 .AND. q1 .NE. q2 .AND. q1 .NE. q4 .AND. &
    !    q3 .NE. q4) THEN
    !  q(q1,q1) = q(q1,q1) + 4.0D0*l4(n)*R(q3)*R(q4) !dq1,qd1
    !  i = MAX(q1,q3)
    !  j = MIN(q1,q3)
    !  q(i,j) = q(i,i) + 2.0D0*l4(n)*R(q1)*R(q4)     !dq1,dq3
    !  i = MAX(q1,q4)
    !  j = MIN(q1,q4)
    !  q(i,j) = q(i,i) + 2.0D0*l4(n)*R(q1)*R(q3)     !dq1,dq4
    !  i = MAX(q3,q4)
    ! NOT FINISHED, MANY MORE TERMS
    !END IF
  END DO

END SUBROUTINE nco_hl4

!------------------------------------------------------------
! nco_analysis_incore
!   - performs analysis of Normal Coordinates
!   - assumes all of R data is held in core
!------------------------------------------------------------
! job           : int, job type
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
! R0            : 1D real*8, position at R0 
! R             : 2D real*8, positions over time 
! w             : 1D real*8, eigenvalues
! q             : 2D real*8, array of normal coordinates
! qi            : 2D real*8, inverse array of normal coordinates
! dt            : real*8, timestep
! tsteps        : int, number of timesteps
! error         : int, exit code

SUBROUTINE nco_analysis_incore(job,ndim,nHO,qHO,HO,nMO,qMO,MO,nl1,&
                    ql1,l1,nl2,ql2,l2,nl3,ql3,l3,nl4,ql4,l4,&
                    dt,tsteps,R0,R,w,q,qi,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:ndim-1,0:tsteps-1), INTENT(INOUT) :: R
  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1), INTENT(INOUT) :: q,qi
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(INOUT) :: w
  REAL(KIND=8), DIMENSION(0:2*nMO-1), INTENT(IN) :: MO
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(IN) :: R0
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
  REAL(KIND=8), INTENT(IN) :: dt
  INTEGER, INTENT(IN) :: job,ndim,nHO,nMO,nl1,nl2,nl3,nl4,tsteps
  INTEGER, INTENT(INOUT) :: error

  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1) :: q0,qo
  REAL(KIND=8), DIMENSION(0:ndim-1) :: zero,Rn
  REAL(KIND=8) :: ti,tf
  INTEGER :: n,i,j

  CALL CPU_TIME(ti)
  error = 0
  q0 = 0.0D0
  zero = 0.0D0

  WRITE(*,*) "-------------------------------------------------------------"
  WRITE(*,*) "Beginning Normal Coordinate Analysis"
  WRITE(*,*)

  !Analysis of q = R0
  CALL nco_general(ndim,nHO,qHO,HO,nMO,qMO,MO,nl1,ql1,l1,nl2,ql2,l2,&
                   nl3,ql3,l3,nl4,ql4,l4,R0,w,q,qi,error)
  WRITE(*,*) "@ q = R0"
  WRITE(*,*) "Renormalized Frequencies @ q = R0"
  WRITE(*,*) w(0:ndim-1)
  WRITE(*,*)
  WRITE(*,*) "Renormalized Coordinates @ q = R0"
  DO n=0,ndim-1
    WRITE(*,*) q(n,0:ndim-1)
  END DO
  WRITE(*,*)
  qo = q

  !Analysis of q = 0 
  CALL nco_general(ndim,nHO,qHO,HO,nMO,qMO,MO,nl1,ql1,l1,nl2,ql2,l2,&
                   nl3,ql3,l3,nl4,ql4,l4,zero,w,q0,qi,error)
  WRITE(*,*) "@ q = 0"
  WRITE(*,*) "Renormalized Frequencies @ q = 0"
  WRITE(*,*) w(0:ndim-1)
  WRITE(*,*)
  WRITE(*,*) "Renormalized Coordinates @ q = 0"
  DO n=0,ndim-1
    WRITE(*,*) q0(n,0:ndim-1)
  END DO
  WRITE(*,*)

  OPEN(file='nco.dat',unit=103,status='replace')
  ! just compare everything to q=0
  IF (job .EQ. 0) THEN
    WRITE(*,*) "Transforming to q=0 normal coordinates"
    R = MATMUL(qi(0:ndim-1,0:ndim-1),R(0:ndim-1,0:tsteps-1))
    DO n=0,tsteps-1
      WRITE(103,*) n*dt,R(0:ndim-1,n)
    END DO

  ! Adapt the normal coordinates for each new position
  ELSE IF (job .EQ. 1) THEN
    WRITE(*,*) "Transforming to adaptive normal coordinates"
    OPEN(unit=104,file='w.dat',status='replace')
    DO n=0,tsteps-1
      CALL nco_general(ndim,nHO,qHO,HO,nMO,qMO,MO,nl1,ql1,l1,nl2,ql2,l2,&
                       nl3,ql3,l3,nl4,ql4,l4,R(0:ndim-1,n),w,q,qi,error)
      !reorder w,q, and wi to overlap of last q
      !CALL nco_order(ndim,q0,q,w,error)
      !CALL nco_angle(ndim,q0,q,error)
      !CALL nco_order(ndim,qo,q,w,error)
      !qi = TRANSPOSE(q)
      !qo = q
      CALL linal_matvec(ndim,qi,R(0:ndim-1,n),Rn)
      WRITE(103,*) n*dt,Rn(0:ndim-1)
      WRITE(104,*) n*dt,w(0:ndim-1)
    END DO
    CLOSE(unit=104)
  ELSE
    WRITE(*,*) "Normal Coordinates will not be written"
  END IF
  CLOSE(unit=103)

  CALL CPU_TIME(tf)
  WRITE(*,*) 
  WRITE(*,*) "Normal coordinate analysis finished in ",tf-ti,"(s)"
END SUBROUTINE nco_analysis_incore 

!------------------------------------------------------------
! nco_order
!    - order normal coordinates to be closest to the ones in 
!      the previous evaluation
!------------------------------------------------------------
! ndim          : int, number of dimensions
! qo            : 2D real*8, last normal coordinates
! q             : 2D real*8, new normal coordinates
! w             : 2D real*8, new frequencies 
! error         : int, exit code

SUBROUTINE nco_order(ndim,qo,q,w,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1), INTENT(INOUT) :: q
  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1), INTENT(IN) :: qo
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(INOUT) :: w
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim

  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1) :: qn
  REAL(KIND=8), DIMENSION(0:ndim-1) :: wn,theta
  INTEGER :: i,j
  error = 0
  wn = 0.D0
  qn = 0.0D0

  ! go through vectors of q
  DO i=0,ndim-1
    !find angles with qo
    DO j=0,ndim-1
      theta(j) = ACOS(ABS(SUM(q(0:ndim-1,i)*qo(0:ndim-1,j))))
    END DO
    j = MINLOC(theta(0:ndim-1),1)-1 !this is the vector we want 
    qn(0:ndim-1,i) = qo(0:ndim-1,j)  
    wn(i) = w(j)
  END DO 

  q = qn
  w = wn

END SUBROUTINE

!------------------------------------------------------------

END MODULE nco
