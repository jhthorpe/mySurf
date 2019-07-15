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
! error         : int, error code
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

END MODULE V 
