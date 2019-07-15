!Module for input parsing 

MODULE input

CONTAINS
!------------------------------------------------------------
! input_read
!   - reads input
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of dimensions
! tmax          : real*8, max time to run simulation
! dt            : real*8, time increment
! tsteps        : real*8, number of timesteps
! mem           : int*8, memory in megabytes
! error         : int, error code

SUBROUTINE input_get(job,ndim,tmax,dt,tsteps,mem,error)
  IMPLICIT NONE
  INTEGER(KIND=8), INTENT(INOUT) :: mem
  REAL(KIND=8), INTENT(INOUT) :: tmax,dt
  INTEGER, INTENT(INOUT) :: job,ndim,error,tsteps
  LOGICAL :: ex

  error = 0
  CALL input_read(job,ndim,tmax,dt,mem,error)
  IF (error .NE. 0) RETURN
  CALL input_check(job,ndim,tmax,dt,mem,error)
  IF (error .NE. 0) RETURN
  tsteps = CEILING(tmax/dt)
  CALL input_write(job,ndim,tmax,dt,tsteps,mem,error) 
  IF (error .NE. 0) RETURN

END SUBROUTINE input_get

!------------------------------------------------------------
! input_read
!   -gets the input from surf.in
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of dimensions
! tmax          : real*8, max time to run simulation
! dt            : real*8, time increment
! mem           : int, mem in MB
! error         : int, error code

SUBROUTINE input_read(job,ndim,tmax,dt,mem,error)
  IMPLICIT NONE
  INTEGER(KIND=8), INTENT(INOUT) :: mem
  REAL(KIND=8), INTENT(INOUT) :: tmax,dt
  INTEGER, INTENT(INOUT) :: job,ndim,error
  LOGICAL :: ex
  error = 0
  INQUIRE(file='surf.in',EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "You must create the input file, surf.in"
    error = 1
    RETURN
  END IF
  OPEN(file='surf.in',unit=100,status='old')
  READ(100,*) job
  READ(100,*) ndim
  READ(100,*) tmax
  READ(100,*) dt
  READ(100,*) mem
  CLOSE(unit=100)
END SUBROUTINE input_read

!------------------------------------------------------------
! input_check
!  -checks the input from surf.in
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of dimensions
! tmax          : real*8, max time to run simulation
! dt            : real*8, time increment
! error         : int, error code

SUBROUTINE input_check(job,ndim,tmax,dt,mem,error)
  IMPLICIT NONE
  INTEGER(KIND=8), INTENT(IN) :: mem
  REAL(KIND=8), INTENT(IN) :: tmax,dt
  INTEGER, INTENT(INOUT) :: error 
  INTEGER, INTENT(IN) :: job,ndim
  error = 0
  IF (job .LT. 0 .OR. job .GT. 0) THEN
    WRITE(*,*) "That jobtype is not supported. Options are..."
    WRITE(*,*) " 0 : basic integration from one starting position"
    error = 1
  END IF
  IF (ndim .LT. 1) THEN
    WRITE(*,*) "Must have at least one dimension"
    error = 2
  END IF
  IF (tmax .LT. 0.D0) THEN
    WRITE(*,*) "Time must be positive"
    error = 3
  END IF
  IF (dt .LT. 0.0D0) THEN
    WRITE(*,*) "dt must be positive"
    error = 4
  END IF
  IF (mem .LT. 1) THEN
    WRITE(*,*) "mem must be more integer more than 1"
    error = 5
  END IF
END SUBROUTINE input_check

!------------------------------------------------------------
! input_write
!  -writes out input back to user
!------------------------------------------------------------
! job           : int, job type
! ndim          : int, number of dimensions
! tmax          : real*8, max time to run simulation
! dt            : real*8, time increment
! tsteps        : int, number of timesteps
! error         : int, error code

SUBROUTINE input_write(job,ndim,tmax,dt,tsteps,mem,error)
  IMPLICIT NONE
  INTEGER(KIND=8), INTENT(IN) :: mem
  REAL(KIND=8), INTENT(IN) :: tmax,dt
  INTEGER, INTENT(INOUT) :: error 
  INTEGER, INTENT(IN) :: job,ndim,tsteps
  error = 0
  WRITE(*,*) "======================================" 
  WRITE(*,*) " mySurf input was...                  "
  IF (job .EQ. 0) THEN
    WRITE(*,*) " job    :  0 - basic MD from coord.dat" 
  END IF
  WRITE(*,*) " ndim   :", ndim
  WRITE(*,*) " tmax   :", tmax, " s"
  WRITE(*,*) " dt     :", dt, " s"
  WRITE(*,*) " tsteps :", tsteps
  WRITE(*,*) " mem    :", mem, " MB"
END SUBROUTINE input_write

!------------------------------------------------------------
! input_cHO
!  - reads in parameters for coupled Harmonic Oscillator 
!    jobs
!------------------------------------------------------------
! ndim          : int, number of dimensions
! k             : 1D real*8, force constants
! l             : 2D real*8, coupling constants
! error         : int, error codes

SUBROUTINE input_cHO(ndim,k,l,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:ndim-1,0:ndim-1), INTENT(INOUT) :: l
  REAL(KIND=8), DIMENSION(0:ndim-1), INTENT(INOUT) :: k
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: ndim
  REAL(KIND=8) :: val
  INTEGER :: n,i,j
  LOGICAL :: ex
  k = 0.0D0
  l = 0.0D0
  error = 0

  WRITE(*,*) "input_cHO : reading from k.in and l.in"
  
  INQUIRE(file='k.in',EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "You must supply force constants in 'k.in'"
    error = 1
    RETURN
  ELSE
    OPEN(file='k.in',unit=100,status='old')
    DO n=0,ndim-1
      READ(100,*) i,val
      k(i-1) = val
    END DO
    CLOSE(unit=100)
  END IF

  INQUIRE(file='l.in',EXIST=ex)
  IF (.NOT. ex) THEN
    WRITE(*,*) "You must supply force constants in 'l.in'"
    error = 1
    RETURN
  ELSE
    OPEN(file='l.in',unit=101,status='old')
    DO n=0,(ndim)**2-1
      READ(101,*) i,j,val
      l(i-1,j-1) = val
    END DO
    CLOSE(unit=101)
  END IF

  WRITE(*,*) "k vector was" 
  DO i=0,ndim-1
    WRITE(*,*) k(i)
  END DO

  WRITE(*,*) "l matrix is"
  DO i=0,ndim-1
    WRITE(*,*) l(i,0:ndim-1)
  END DO

END SUBROUTINE input_cHO
!------------------------------------------------------------

END MODULE input
