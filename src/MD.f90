!Module for MD simulations using velocity verlet
MODULE MD
  USE input
  USE V
  USE nco
  
CONTAINS
!------------------------------------------------------------
! MD_start
!------------------------------------------------------------
! job           : int, jobtype
! ndim          : int, number of dimensions
! tmax          : real*8, max time to run simulation
! dt            : real*8, time increment
! tsteps        : real*8, number of timesteps
! mem           : int*8, memory in MB
! error         : int, error code
! R             : 2D real*8, coordinates as a function of time 
! R0            : 1D real*8, initial coordinates

SUBROUTINE MD_start(job,ndim,tmax,dt,tsteps,mem,error)
  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(INOUT) :: mem
  REAL(KIND=8), INTENT(IN) :: tmax,dt
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,tsteps,job
  REAL(KIND=8), DIMENSION(0:ndim-1) :: R0  
  INTEGER :: i
  LOGICAL :: ex
  
  error = 0
  WRITE(*,*) "MD_start : called"

  !Read in initial position
  INQUIRE(file='R0.in',EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "You need the R0.in input file"
    error = 1
    RETURN
  END IF
  OPEN(file='R0.in',unit=100,status='old')
  DO i=0,ndim-1
    READ(100,*) R0(i)
  END DO 
  CLOSE(unit=100)
  WRITE(*,*) "R0 is", R0

  CALL MD_run(job,ndim,tmax,dt,tsteps,R0,mem,error) 


END SUBROUTINE MD_start

!------------------------------------------------------------
! MD_run
!  -runs MD
!------------------------------------------------------------
! job           : int, jobtype
! ndim          : int, number of dimensions
! tmax          : real*8, max time to run simulation
! dt            : real*8, time increment
! tsteps        : real*8, number of timesteps
! mem           : int*8, memory in MB
! error         : int, error code
! R0            : 1D real*8, initial coordinates

SUBROUTINE MD_run(job,ndim,tmax,dt,tsteps,R0,mem,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(IN) :: R0
  INTEGER(KIND=8), INTENT(IN) :: mem
  REAL(KIND=8), INTENT(IN) :: tmax,dt
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,tsteps,job

  INTEGER(KIND=8) :: qword

  error = 0
  WRITE(*,*) "MD_run : called"
  WRITE(*,*) "MD_run : starting memory analysis"

  !Memory Analysis
  qword = 1000000*mem/8
  qword = qword - 10000
  IF (qword .LT. 0) THEN
    WRITE(*,*) "MD_run : We ran out of memory"
    error = 1 
    RETURN
  END IF

  IF (qword - tsteps*ndim - tsteps - 10*ndim .GT. 100) THEN
    WRITE(*,*) "MD_run : Holding everything in memory" 
    CALL MD_incore(job,ndim,tmax,dt,tsteps,R0,error)
  ELSE
    WRITE(*,*) "MD_run : We ran out of memory"
    error = 1
    RETURN
  END IF
  
END SUBROUTINE MD_run

!------------------------------------------------------------
! MD_incore
!  - runs MD with everything held in core
!------------------------------------------------------------
! job           : int, jobtype
! ndim          : int, number of dimensions
! tmax          : real*8, max time to run simulation
! dt            : real*8, time increment
! tsteps        : real*8, number of timesteps
! mem           : int*8, memory in MB
! error         : int, error code
! R             : 2D real*8, coordinates as a function of time 
! R0            : 1D real*8, initial coordinates
! dR            : 1D real*8, velocities
! ddR           : 1D real*8, acceleration
! w             : 1D real*8, eigenvalues of Hessian
! q             : 2D real*8, normal coordinate vectors
! qi            : 2D real*8, inverse of q

SUBROUTINE MD_incore(job,ndim,tmax,dt,tsteps,R0,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(IN) :: R0
  REAL(KIND=8), INTENT(IN) :: tmax,dt
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,tsteps,job

  REAL(KIND=8), DIMENSION(0:ndim-1,0:tsteps-1) :: R
  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1) :: q,qi
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: l
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: k,HO,MO,l1,l2,l3,l4
  INTEGER, DIMENSION(:), ALLOCATABLE :: qHO,qMO,ql1,ql2,ql3,ql4
  REAL(KIND=8), DIMENSION(0:tsteps-1) :: energy
  REAL(KIND=8), DIMENSION(0:ndim-1) :: dR,ddR,ddRo,dV,w
  REAL(KIND=8) :: V,dtdt,ti,tf
  INTEGER :: n,nHO,nMO,nl1,nl2,nl3,nl4

  CALL CPU_TIME(ti)
  WRITE(*,*) "MD_incore : Starting simulation"
  WRITE(*,*) "WARNING -- assuming mass is 1" 

  error = 0
  dtdt = dt**2.0D0

  IF (job .EQ. 0) THEN !general case
    CALL input_general(ndim,nHO,qHO,HO,nMO,qMO,MO,nl1,ql1,l1,nl2,ql2,l2,&
                       nl3,ql3,l3,nl4,ql4,l4,error)
    IF (error .NE. 0) RETURN
    CALL nco_general(ndim,nHO,qHO,HO,nMO,qMO,MO,nl1,ql1,l1,nl2,ql2,l2,&
                     nl3,ql3,l3,nl4,ql4,l4,w,q,qi,error)
    IF (error .NE. 0) RETURN
    CALL V_general(ndim,nHO,qHO,HO,nMO,qMO,MO,nl1,ql1,l1,nl2,ql2,l2,&
                   nl3,ql3,l3,nl4,ql4,l4,R0,V,dV,error)
    IF (error .NE. 0) RETURN
  ELSE IF (job .EQ. 1) THEN
    ALLOCATE(k(0:ndim-1))
    ALLOCATE(l(0:ndim-1,0:ndim-1))
    CALL input_cHO(ndim,k(0:ndim-1),l(0:ndim-1,0:ndim-1),error)
    IF (error .NE. 0) RETURN
    CALL nco_cHO(ndim,k(0:ndim-1),l(0:ndim-1,0:ndim-1),w(0:ndim-1),&
                 q(0:ndim-1,0:ndim-1),qi(0:ndim-1,0:ndim-1),error)
    IF (error .NE. 0) RETURN
    CALL V_cHO(ndim,R0(0:ndim-1),k(0:ndim-1),l(0:ndim-1,0:ndim-1),V,dV(0:ndim-1),error) 
    IF (error .NE. 0 ) RETURN
  ELSE
    WRITE(*,*) "Sorry, that jobtype is not supported"
    error = 1
    RETURN
  END IF

  R(0:ndim-1,0) = R0(0:ndim-1)
  dR = 0
  ddR = -dV  !ignoring mass for now, it should be built into potential...
  energy(0) = V
  WRITE(*,*) "Initial energy is", V

  DO n=1,tsteps-1
    R(0:ndim-1,n) = R(0:ndim-1,n-1) + dR*dt &
                    + 0.5D0*ddR*dtdt
    IF (job .EQ. 0) THEN
    CALL V_general(ndim,nHO,qHO,HO,nMO,qMO,MO,nl1,ql1,l1,nl2,ql2,l2,&
                   nl3,ql3,l3,nl4,ql4,l4,R(0:ndim-1,n),V,dV,error)
    ELSE IF (job .EQ. 1) THEN
      call V_cHO(ndim,R(0:ndim-1,n),k(0:ndim-1),l(0:ndim-1,0:ndim-1),V,dV(0:ndim-1),error)
    ELSE
      WRITE(*,*) "Sorry, that job type is not supported"
      error = 1
      RETURN
    END IF
    IF (error .NE. 0) RETURN
    ddRo = ddR
    ddR = -dV
    dR = dR + 0.5D0*(ddRo + ddR)*dt 
    energy(n) = V + 0.5D0*SUM(dR**2.0D0)

  END DO

  CALL CPU_TIME(tf)
  WRITE(*,*) "MD_incore : simulation finished in ",tf-ti,"(s)"
  WRITE(*,*) "Writing output to R.dat, E.dat"
  OPEN(file='R.dat',unit=101,status='replace')
  DO n=0,tsteps-1
    WRITE(101,*) n*dt,R(0:ndim-1,n)
  END DO
  CLOSE(unit=101)
  OPEN(file='E.dat',unit=102,status='replace')
  DO n=0,tsteps-1
    WRITE(102,*) n*dt,energy(n)
  END DO
  CLOSE(unit=102)
  !transform into normal coordinates
  R = MATMUL(qi(0:ndim-1,0:ndim-1),R(0:ndim-1,0:tsteps-1))
  OPEN(file='nco.dat',unit=103,status='replace')
  DO n=0,tsteps-1
    WRITE(103,*) n*dt,R(0:ndim-1,n)
  END DO
  CLOSE(unit=103)

  IF (ALLOCATED(k)) DEALLOCATE(k)
  IF (ALLOCATED(l)) DEALLOCATE(l)


END SUBROUTINE MD_incore

!------------------------------------------------------------
END MODULE MD
