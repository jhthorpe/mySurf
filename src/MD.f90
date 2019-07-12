!Module for MD simulations using velocity verlet
MODULE MD
  USE hamiltonian
  
CONTAINS
!------------------------------------------------------------
! MD_start
!------------------------------------------------------------
! ndim          : int, number of dimensions
! tmax          : real*8, max time to run simulation
! dt            : real*8, time increment
! tsteps        : real*8, number of timesteps
! mem           : int*8, memory in MB
! error         : int, error code
! R             : 2D real*8, coordinates as a function of time 
! R0            : 1D real*8, initial coordinates

SUBROUTINE MD_start(ndim,tmax,dt,tsteps,mem,error)
  IMPLICIT NONE

  INTEGER(KIND=8), INTENT(INOUT) :: mem
  REAL(KIND=8), INTENT(IN) :: tmax,dt
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim,tsteps
  REAL(KIND=8), DIMENSION(0:ndim-1) :: R0  
  INTEGER :: i
  LOGICAL :: ex
  
  error = 0
  WRITE(*,*) "MD_start : called"

  INQUIRE(file='R0.dat',EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "You need the R0.dat input file"
    error = 1
    RETURN
  END IF

  OPEN(file='R0.dat',unit=100,status='old')
  DO i=0,ndim-1
    READ(100,*) R0(i)
  END DO 
  CLOSE(unit=100)

  CALL MD_run(ndim,tmax,dt,tsteps,R0,mem,error) 


END SUBROUTINE MD_start

!------------------------------------------------------------
! MD_run
!  -runs MD
!------------------------------------------------------------
! ndim          : int, number of dimensions
! tmax          : real*8, max time to run simulation
! dt            : real*8, time increment
! tsteps        : real*8, number of timesteps
! mem           : int*8, memory in MB
! error         : int, error code
! R0            : 1D real*8, initial coordinates

SUBROUTINE MD_run(ndim,tmax,dt,tsteps,R0,mem,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(IN) :: R0
  INTEGER(KIND=8), INTENT(IN) :: mem
  REAL(KIND=8), INTENT(IN) :: tmax,dt
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim,tsteps

  INTEGER(KIND=8) :: qword

  error = 0
  WRITE(*,*) "MD_run : called"
  WRITE(*,*) "MD_run : starting memory analysis"

  !Memory Analysis
  qword = 1000000*FLOOR(mem/8)
  qword = qword - 10000
  IF (qword .LT. 0) THEN
    WRITE(*,*) "MD_run : We ran out of memory"
    error = 1 
    RETURN
  END IF

  IF (qword - tsteps*ndim - tsteps -4*ndim .GT. 0)  
    WRITE(*,*) "MD_run : Holding everything in memory" 
    CALL MD_incore(ndim,tmax,dt,tsteps,R0,error)
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
! ndim          : int, number of dimensions
! tmax          : real*8, max time to run simulation
! dt            : real*8, time increment
! tsteps        : real*8, number of timesteps
! mem           : int*8, memory in MB
! error         : int, error code
! R             : 2D real*8, coordinates as a function of time 
! R0            : 1D real*8, initial coordinates

SUBROUTINE MD_incore(ndim,tmax,dt,tsteps,R0,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(IN) :: R0
  INTEGER(KIND=8), INTENT(IN) :: mem
  REAL(KIND=8), INTENT(IN) :: tmax,dt
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: job,ndim,tsteps

  REAL(KIND=8), DIMENSION(0:ndim-1,0:tmax-1) :: R
  REAL(KIND=8), DIMENSION(0:ndim-1) :: Vo,Vn,Ao,An
  INTEGER :: i,n

  error = 0
  WRITE(*,*) "MD_incore : Starting simulation"
  R(0:ndim-1,0) = R0(0:ndim-1)
  Vo = 0
  

  DO n=0,tsteps-1
    
  

  END DO


END SUBROUTINE MD_incore

!------------------------------------------------------------
END MODULE MD
