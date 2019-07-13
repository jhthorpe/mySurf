!Module for MD simulations using velocity verlet
MODULE MD
  USE V
  
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
  INTEGER, INTENT(IN) :: ndim,tsteps
  REAL(KIND=8), DIMENSION(0:ndim-1) :: R0  
  INTEGER :: i
  LOGICAL :: ex
  
  error = 0
  WRITE(*,*) "MD_start : called"

  !Read in initial position
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
  WRITE(*,*) "R0 is", R0

  !Read in parameters -- implement later

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
  INTEGER, INTENT(IN) :: ndim,tsteps

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
! dR            : 1D real*8, velocities
! ddR           : 1D real*8, acceleration

SUBROUTINE MD_incore(ndim,tmax,dt,tsteps,R0,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(IN) :: R0
  REAL(KIND=8), INTENT(IN) :: tmax,dt
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim,tsteps

  REAL(KIND=8), DIMENSION(0:ndim-1,0:tsteps-1) :: R
  REAL(KIND=8), DIMENSION(0:tsteps-1) :: energy
  REAL(KIND=8), DIMENSION(0:ndim-1) :: dR,ddR,ddRo,dV
  REAL(KIND=8) :: V,dtdt,ti,tf
  INTEGER :: n

  CALL CPU_TIME(ti)
  WRITE(*,*) "MD_incore : Starting simulation"
  WRITE(*,*) "WARNING -- assuming mass is 1" 

  error = 0
  dtdt = dt**2.0D0

  CALL V_eval(ndim,R0,V,dV,error) 
  IF (error .NE. 0 ) THEN
    WRITE(*,*) "MD_incore : error out of V_eval"
    RETURN
  END IF

  R(0:ndim-1,0) = R0(0:ndim-1)
  dR = 0
  ddR = -dV  !ignoring mass for now, it should be built into potential...
  energy(0) = V

  DO n=1,tsteps-1
    R(0:ndim-1,n) = R(0:ndim-1,n-1) + dR*dt &
                    + 0.5D0*ddR*dtdt
    call V_eval(ndim,R(0:ndim-1,n),V,dV(0:ndim-1),error)
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


END SUBROUTINE MD_incore

!------------------------------------------------------------
END MODULE MD
