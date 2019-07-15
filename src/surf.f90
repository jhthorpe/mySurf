!Program for modeling dynamics on a surface

PROGRAM surf
  USE input
  USE MD
  IMPLICIT NONE
  INTEGER(KIND=8) :: mem
  REAL(KIND=8) :: tmax,dt
  INTEGER :: job,ndim,error,tsteps
  REAL(KIND=8) :: ti

  CALL CPU_TIME(ti)
  WRITE(*,*) "Starting mySurf"

  CALL input_get(job,ndim,tmax,dt,tsteps,mem,error)
  IF (error .NE. 0) CALL finish(error,ti)

  IF (job .EQ. 0) THEN
    CALL MD_start(job,ndim,tmax,dt,tsteps,mem,error)
  END IF
  IF (error .NE. 0) CALL finish(error,ti)

  CALL finish(error,ti) 

CONTAINS
SUBROUTINE finish(error,ti)
  IMPLICIT NONE
  REAL(KIND=8), INTENT(IN) :: ti
  INTEGER, INTENT(IN) :: error
  REAL(KIND=8) :: tf
  CALL CPU_TIME(tf)
  WRITE(*,*)
  WRITE(*,*) "mySurf finished with status", error," in ",tf-ti," (s)"
  IF (error .NE. 0) THEN
    STOP 1 
  END IF 
   
END SUBROUTINE finish

END PROGRAM surf


