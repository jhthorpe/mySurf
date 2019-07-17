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
    WRITE(*,*) " 0 : general potential - HO,MO, order 1-4 terms" 
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
    WRITE(*,*) " job    :  0 - Harmonic, Morse, terms up to order 4"
  END IF
  WRITE(*,*) " ndim   :", ndim
  WRITE(*,*) " tmax   :", tmax, " s"
  WRITE(*,*) " dt     :", dt, " s"
  WRITE(*,*) " tsteps :", tsteps
  WRITE(*,*) " mem    :", mem, " MB"
  WRITE(*,*) " ======================================"
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

  WRITE(*,*)

END SUBROUTINE input_cHO

!------------------------------------------------------------
! input_general
!       - read in input for general simulation
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
! error         : int, exit code

SUBROUTINE input_general(ndim,nHO,qHO,HO,nMO,qMO,MO,nl1,ql1,l1,nl2,ql2,l2,&
                         nl3,ql3,l3,nl4,ql4,l4,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: HO,MO,l1,l2,l3,l4 
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qHO,qMO,ql1,ql2,ql3,ql4
  INTEGER, INTENT(INOUT) :: nHO,nMO,nl1,nl2,nl3,nl4,error
  INTEGER, INTENT(IN) :: ndim
  
  INTEGER :: i,j
  LOGICAL :: ex,flag

  nHO = 0
  nMO = 0
  nl1 = 0
  nl2 = 0
  nl3 = 0
  nl4 = 0

  CALL input_HO(ndim,nHO,qHO,HO,error)
  IF (error .NE. 0) RETURN
  CALL input_MO(ndim,nMO,qMO,MO,error)
  IF (error .NE. 0) RETURN
  CALL input_l1(ndim,nl1,ql1,l1,error)
  IF (error .NE. 0) RETURN
  CALL input_l2(ndim,nl2,ql2,l2,error)
  IF (error .NE. 0) RETURN
  CALL input_l3(ndim,nl3,ql3,l3,error)
  IF (error .NE. 0) RETURN
  CALL input_l4(ndim,nl4,ql4,l4,error)
  IF (error .NE. 0) RETURN
 
  IF (nHO + nMO + nl1 + nl2 + nl3 + nl4 .LE. 0) THEN
    error = 1
    WRITE(*,*) "input_general : no input files found"
    RETURN
  END IF

END SUBROUTINE input_general

!------------------------------------------------------------
! input_HO
!       - reads input from HO file
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nHO           : int, number of harmonic oscillators
! qHO           : 1D int, ids of harmonic oscillators
! HO            : 1D real*8, force constants of HOs
! error         : int, exit code

SUBROUTINE input_HO(ndim,nHO,qHO,HO,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: HO
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qHO
  INTEGER, INTENT(INOUT) :: nHO,error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,q1,q2,q3,q4,uid
  LOGICAL :: ex
  error = 0
  nHO = 0
  uid = 200
  INQUIRE(file='HO.in',EXIST=ex)
  IF (.NOT. ex) RETURN
  OPEN(file='HO.in',unit=uid,status='old')
  READ(uid,*) nHO
  ALLOCATE(qHO(0:nHO-1))
  ALLOCATE(HO(0:nHO-1))
  qHO = 0
  HO = 0.0D0
  DO i=0,nHO-1
    READ(uid,*) q1,HO(i) 
    IF (q1 .LT. 1 .OR. q1 .GT. ndim) THEN
      error = 1
      WRITE(*,*) "In HO.in, input ", i," is outside range [1:ndim]"
    END IF 
    qHO(i) = q1 - 1
  END DO
  CLOSE(unit=uid)
END SUBROUTINE input_HO

!------------------------------------------------------------
! input_MO
!       - reads input from MO file
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nMO           : int, number of morse oscillators
! qMO           : 1D int, ids of morse oscillators
! MO            : 1D real*8, constants of MOs [{k@0,De}]
! error         : int, exit code

SUBROUTINE input_MO(ndim,nMO,qMO,MO,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: MO
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: qMO
  INTEGER, INTENT(INOUT) :: nMO,error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,q1,q2,q3,q4,uid
  LOGICAL :: ex
  error = 0
  nMO = 0
  uid = 201
  INQUIRE(file='MO.in',EXIST=ex)
  IF (.NOT. ex) RETURN
  OPEN(file='MO.in',unit=uid,status='old')
  READ(uid,*) nMO
  ALLOCATE(qMO(0:nMO-1))
  ALLOCATE(MO(0:2*nMO-1))
  DO i=0,nMO-1
    READ(uid,*) q1,MO(2*i),MO(2*i+1) !id, k, De 
    IF (q1 .LT. 1 .OR. q1 .GT. ndim) THEN
      error = 1
      WRITE(*,*) "In MO.in, input ", i," is outside range [1:ndim]"
    END IF 
    qMO(i) = q1 - 1
  END DO
  CLOSE(unit=uid)
  !transform k's to alphas
  DO i=0,nMO-1
    MO(2*i) = SQRT(MO(2*i)/(2.0D0*MO(2*i+1))) 
  END DO
  WRITE(*,*) "qMO is",qMO
  WRITE(*,*) "MO is", MO
END SUBROUTINE input_MO

!------------------------------------------------------------
! input_l1
!       - reads input from l1 file
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nl1           : int, number of order 1 terms 
! ql1           : 1D int, ids of order 1 terms 
! l1            : 1D real*8, constants of order 1 terms 
! error         : int, exit code

SUBROUTINE input_l1(ndim,nl1,ql1,l1,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: l1
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ql1
  INTEGER, INTENT(INOUT) :: nl1,error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,q1,q2,q3,q4,uid
  LOGICAL :: ex
  error = 0
  nl1 = 0
  uid = 202
  INQUIRE(file='l1.in',EXIST=ex)
  IF (.NOT. ex) RETURN
  OPEN(file='l1.in',unit=uid,status='old')
  READ(uid,*) nl1
  ALLOCATE(ql1(0:nl1-1))
  ALLOCATE(l1(0:2*nl1-1))
  ql1 = 0
  l1 = 0.0D0
  DO i=0,nl1-1
    READ(uid,*) q1,l1(i) 
    IF (q1 .LT. 1 .OR. q1 .GT. ndim) THEN
      error = 1
      WRITE(*,*) "In l1.in, input ", i," is outside range [1:ndim]"
    END IF 
    ql1(i) = q1 - 1
  END DO
  CLOSE(unit=uid)
END SUBROUTINE input_l1

!------------------------------------------------------------
! input_l2
!       - reads input from l2 file
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nl2           : int, number of order 2 terms 
! ql2           : 1D int, ids of order 2 terms 
! l2            : 1D real*8, constants of order 2 terms 
! error         : int, exit code

SUBROUTINE input_l2(ndim,nl2,ql2,l2,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: l2
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ql2
  INTEGER, INTENT(INOUT) :: nl2,error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,q1,q2,q3,q4,uid
  LOGICAL :: ex
  error = 0
  nl2 = 0
  uid = 203
  INQUIRE(file='l2.in',EXIST=ex)
  IF (.NOT. ex) RETURN
  OPEN(file='l2.in',unit=uid,status='old')
  READ(uid,*) nl2
  ALLOCATE(ql2(0:2*nl2-1))
  ALLOCATE(l2(0:nl2-1))
  ql2 = 0
  l2 = 0.0D0
  DO i=0,nl2-1
    READ(uid,*) q1,q2,l2(i) 
    IF (q1 .LT. 1 .OR. q2 .LT. 1 .OR. q1 .GT. ndim .OR. q2 .GT. ndim) THEN
      error = 1
      WRITE(*,*) "In l2.in, input ", i," is outside range [1:ndim]"
    END IF
    ql2(2*i) = q1-1
    ql2(2*i+1) = q2 - 1
  END DO
  CLOSE(unit=uid)
END SUBROUTINE input_l2
!------------------------------------------------------------
! input_l3
!       - reads input from l3 file
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nl3           : int, number of order 3 terms 
! ql3           : 1D int, ids of order 3 terms 
! l3            : 1D real*8, constants of order 3 terms 
! error         : int, exit code

SUBROUTINE input_l3(ndim,nl3,ql3,l3,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: l3
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ql3
  INTEGER, INTENT(INOUT) :: nl3,error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,q1,q2,q3,q4,uid
  LOGICAL :: ex
  error = 0
  nl3 = 0
  uid = 204
  INQUIRE(file='l3.in',EXIST=ex)
  IF (.NOT. ex) RETURN
  OPEN(file='l3.in',unit=uid,status='old')
  READ(uid,*) nl3
  ALLOCATE(ql3(0:3*nl3-1))
  ALLOCATE(l3(0:nl3-1))
  DO i=0,nl3-1
    READ(uid,*) q1,q2,q3,l3(i) 
    IF (q1 .LT. 1 .OR. q2 .LT. 1 .OR. q1 .GT. ndim .OR. q2 .GT. ndim &
        .OR. q3 .LT. 1 .OR. q3 .GT. ndim) THEN
      error = 1
      WRITE(*,*) "In l3.in, input ", i," is outside range [1:ndim]"
    END IF
    ql3(3*i) = q1-1
    ql3(3*i+1) = q2-1
    ql3(3*i+2) = q3-1
  END DO
  CLOSE(unit=uid)
END SUBROUTINE input_l3
!------------------------------------------------------------
! input_l4
!       - reads input from l4 file
!------------------------------------------------------------
! ndim          : int, number of dimensions
! nl4           : int, number of order 4 terms 
! ql4           : 1D int, ids of order 4 terms 
! l4            : 1D real*8, constants of order 4 terms 
! error         : int, exit code

SUBROUTINE input_l4(ndim,nl4,ql4,l4,error)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: l4
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: ql4
  INTEGER, INTENT(INOUT) :: nl4,error
  INTEGER, INTENT(IN) :: ndim
  INTEGER :: i,q1,q2,q3,q4,uid
  LOGICAL :: ex
  error = 0
  nl4 = 0
  uid = 205
  INQUIRE(file='l4.in',EXIST=ex)
  IF (.NOT. ex) RETURN
  OPEN(file='l4.in',unit=uid,status='old')
  READ(uid,*) nl4
  ALLOCATE(ql4(0:4*nl4-1))
  ALLOCATE(l4(0:nl4-1))
  DO i=0,nl4-1
    READ(uid,*) q1,q2,q3,q4,l4(i) 
    IF (q1 .LT. 1 .OR. q2 .LT. 1 .OR. q1 .GT. ndim .OR. q2 .GT. ndim &
        .OR. q3 .LT. 1 .OR. q3 .GT. ndim .OR. q4 .LT. 1 .OR. q4 .GT. ndim) THEN
      error = 1
      WRITE(*,*) "In l4.in, input ", i," is outside range [1:ndim]"
    END IF
    ql4(4*i) = q1 - 1
    ql4(4*i+1) = q2 - 1
    ql4(4*i+2) = q3 - 1
    ql4(4*i+3) = q4 - 1
  END DO
  CLOSE(unit=uid)
END SUBROUTINE input_l4
!------------------------------------------------------------


END MODULE input
