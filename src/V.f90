!Module containing subroutines to calculate potential energy and
!gradients

MODULE V 

CONTAINS
!------------------------------------------------------------
! V_cHO
!   - evaluate potential  
!   - currently, this is hard coded, but could be generalized
!------------------------------------------------------------
! ndim          : int, number of dimensions
! R             : 1D real*8, coordinates [0:ndim-1]
! k             : 1D real*8, force constants
! l             : 2D real*8, coupling matrix
! V             : real*8, potential 
! dV            : 1D real*8, gradient    [0:ndim-1]
! error         : int, exit code

SUBROUTINE V_cHO(ndim,R,k,l,V,dV,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1),INTENT(IN) :: l
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(INOUT) :: dV
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(IN) :: R,k
  REAL(KIND=8), INTENT(INOUT) :: V
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim

  REAL(KIND=8), DIMENSION(0:ndim-1) :: dV0,dV1,dV2,dL0,dL1,dL2 
  REAL(KIND=8) :: V0,V1,V2,L0,L1,L2 
  INTEGER :: i,j

  error = 0

  !Scalar terms
  V0 = 0.0D0
  dV0 = 0.0D0
  L0 = 0.0D0
  dL0 = 0.0D0

  !Linear terms
  V1 = 0.0D0
  L1 = 0.0D0
  dV1 = 0.0D0
  dL1 = 0.0D0
  DO j=0,ndim-1
    dV1(j) = dV1(j) + k(j)*R(j)
    DO i=0,j-1 
      L1 = L1 + l(i,j)*R(i)*R(j)
      dL1(j) = dL1(j) + l(i,j)*R(i)
    END DO
    DO i=j,ndim-1
      dL1(j) = dL1(j) + l(i,j)*R(i)
    END DO
  END DO

  !Quadratic terms
  V2 = 0.0D0
  L2 = 0.0D0
  dV2 = 0.0D0
  dL2 = 0.0D0
  DO i=0,ndim-1
    V2 = V2 + k(i)*R(i)**2.0D0 
  END DO
  V2 = V2 * 0.5D0 

  V = V0 + V1 + V2 + L0 + L1 + L2
  dV = dV0 + dV1 + dV2 + dL0 + dL1 + dL2

!  WRITE(*,*) "V is:", V
!  WRITE(*,*) "dV is:", dV
  

END SUBROUTINE V_cHO


!------------------------------------------------------------
! V_general
!    - calculates potential and gradient for the general case
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
! R             : 1D real*8, coordinates [0:ndim-1]
! V             : real*8, potential 
! dV            : 1D real*8, gradient    [0:ndim-1]
! error         : int, exit code

SUBROUTINE V_general(ndim,nHO,qHO,HO,nMO,qMO,MO,nl1,ql1,l1,nl2,ql2,l2,&
                     nl3,ql3,l3,nl4,ql4,l4,R,V,dV,error)

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(INOUT) :: dV
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
  REAL(KIND=8), INTENT(INOUT) :: V
  INTEGER, INTENT(IN) :: ndim,nHO,nMO,nl1,nl2,nl3,nl4
  INTEGER, INTENT(INOUT) :: error

  REAL(KIND=8), DIMENSION(0:ndim-1) :: dvHO,dvMO,dvl1,dvl2,dvl3,dvl4
  REAL(KIND=8) :: vHO,vMO,vl1,vl2,vl3,vl4
  INTEGER :: n,i,j,k,l

  error = 0 
  vHO = 0.0D0
  vMO = 0.0D0
  vl1 = 0.0D0
  vl2 = 0.0D0
  vl3 = 0.0D0
  vl4 = 0.0D0
  dvHO = 0.0D0
  dvMO = 0.0D0
  dvl1 = 0.0D0
  dvl2 = 0.0D0
  dvl3 = 0.0D0
  dvl4 = 0.0D0

  !WRITE(*,*) "R ", R

  !Harmonic Terms
  DO n=0,nHO-1
    i = qHO(n)
    vHO = vHO + HO(n)*R(i)**2.0D0 
    dvHO(i) = dvHO(i) + HO(n)*R(i)
  END DO
  vHO = vHO * 0.5D0

  !Morse Terms
  ! MO(2*n) is alpha
  ! MO(2*n+1) is De
  DO n=0,nMO-1
    i = qMO(n)
    vMO = vMO + MO(2*n+1) + MO(2*n+1)*&
                (EXP(-2.0D0*MO(2*n)*R(i)) - 2.0D0*EXP(-1.0D0*MO(2*n)*R(i))) 
    dvMO(i) = dvMO(i) + MO(2*n+1)*(-2.0D0*MO(2*n)*EXP(-2.0D0*MO(2*n)*R(i)) &
                   + 2.0D0*MO(2*n)*EXP(-1.0D0*MO(2*n)*R(i)))
  END DO

  !Order 1 terms
  ! λ * qi 
  DO n=0,nl1-1
    i = ql1(n)
    vl1 = vl1 + l1(n)*R(i)
    dvl1(i) = dvl1(i) + l1(n) 
  END DO

  !Order 2 terms
  ! λ * qi * qj 
  DO n=0,nl2-1
    i = ql2(2*n)
    j = ql2(2*n+1)
    vl2 = vl2 + l2(n)*R(i)*R(j)
    dvl2(i) = dvl2(i) + l2(n)*R(j)
    dvl2(j) = dvl2(j) + l2(n)*R(i)
  END DO

  !Order 3 terms
  ! λ * qi * qj  * qk
  DO n=0,nl3-1
    i = ql3(3*n)
    j = ql3(3*n+1)
    k = ql3(3*n+2)
    vl3 = vl3 + l3(n)*R(i)*R(j)*R(k)
    dvl3(i) = dvl3(i) + l3(n)*R(j)*R(k)
    dvl3(j) = dvl3(j) + l3(n)*R(i)*R(k)
    dvl3(k) = dvl3(k) + l3(n)*R(i)*R(j)
  END DO

  !Order 4 terms
  ! λ * qi * qj  * qk * ql
  DO n=0,nl4-1
    i = ql4(4*n)
    j = ql4(4*n+2)
    k = ql4(4*n+3)
    l = ql4(4*n+4)
    vl4 = vl4 + l4(n)*R(i)*R(j)*R(k)*R(l)
    dvl4(i) = dvl4(i) + l4(n)*R(j)*R(k)*R(l)
    dvl4(j) = dvl4(j) + l4(n)*R(i)*R(k)*R(l)
    dvl4(k) = dvl4(k) + l4(n)*R(i)*R(j)*R(l)
    dvl4(l) = dvl4(l) + l4(n)*R(i)*R(j)*R(k)
  END DO

  V = vHO + vMO + vl1 + vl2 + vl3 + vl4
  dV = dvHO + dvMO + dvl1 + dvl2 + dvl3 + dvl4

  !WRITE(*,*) 
  !WRITE(*,*) "vHO", vHO
  !WRITE(*,*) "vMO", vMO
  !WRITE(*,*) "vl1", vl1
  !WRITE(*,*) "vl2", vl2
  !WRITE(*,*) "vl3", vl3
  !WRITE(*,*) "vl4", vl4
  !WRITE(*,*) "V is", V
  !WRITE(*,*)  
  !WRITE(*,*) "dvHO", dvHO
  !WRITE(*,*) "dvMO", dvMO
  !WRITE(*,*) "dvl1", dvl1
  !WRITE(*,*) "dvl2", dvl2
  !WRITE(*,*) "dvl3", dvl3
  !WRITE(*,*) "dvl4", dvl4
  !WRITE(*,*) "dV is", dV

END SUBROUTINE V_general

!------------------------------------------------------------

END MODULE V 
