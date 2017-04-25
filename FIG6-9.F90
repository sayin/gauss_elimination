PROGRAM test_simul

IMPLICIT NONE

! Declare parameters:
INTEGER, PARAMETER :: max_size = 10    ! Max number of eqns 

! Declare variables:
INTEGER :: i, j, n, istat, er
REAL    :: ran
REAL, DIMENSION(max_size,max_size) :: a ,b1
REAL, DIMENSION(max_size) :: b


READ (*,*) ran
! Get the name of the disk file containing the equations.
!WRITE (*,1000) 
!!READ (*,'(A20)') file_name
 
! Open input data file.  Status is OLD because the input data must 
! already exist.
OPEN ( UNIT=1, FILE='OUTPUT.DAT', STATUS='REPLACE', ACTION='READWRITE', &
       IOSTAT=istat )
READ(*,*) n
WRITE(1,1000) n
1000 FORMAT (1X,I4)

!RANDOM GENERATOR
DO i=1,n
 DO j=1,n+1
    CALL RANDOM_NUMBER(ran)
    b1(i,j)=ran*10  
END DO
END DO

 DO i=1,n
WRITE(1,100) (b1(i,j), j=1,n+1)
     100 FORMAT (1X,7F11.4)
END DO
CLOSE (UNIT=1)


! Open input data file.  Status is OLD because the input data must 
! already exist.
OPEN ( UNIT=1, FILE='OUTPUT.DAT', STATUS='OLD',ACTION='READ',IOSTAT=istat )

! Was the OPEN successful? 
fileopen: IF ( istat == 0 ) THEN
   ! The file was opened successfully, so read the number of 
   ! equations in the system.
   READ (1,*) n
 
   ! If the number of equations is <= max_size, read them in
   ! and process them.
   size_ok: IF ( n <= max_size ) THEN
      DO i = 1, n
         READ (1,*) (a(i,j), j=1,n) ,b(i)
      END DO
 
      ! Display coefficients.
      WRITE (*,1020)
      1020 FORMAT (/,1X,'Coefficients before call:')
      DO i = 1, n
         WRITE (*,1030) (a(i,j), j=1,n), b(i)
         1030 FORMAT (1X,7F11.4)
      END DO
 
      ! Solve equations.
      CALL simul (a, b, max_size, n, er )
 
      ! Check for error.
      error_check: IF ( er /= 0 ) THEN

         WRITE (*,1040) 
         1040 FORMAT (/1X,'Zero pivot encountered!', &
                     //1X,'There is no unique solution to this system.')

      ELSE error_check
 
         ! No errors. Display coefficients.
         WRITE (*,1050)
         1050 FORMAT (/,1X,'Coefficients after call:')
         DO  i = 1, n
             WRITE (*,1030) (a(i,j), j=1,n), b(i)
         END DO
 
         ! Write final answer.
         WRITE (*,1060)
         1060 FORMAT (/,1X,'The solutions are:')
         DO i = 1, n
            WRITE (*,1070) i, b(i)
            1070 FORMAT (3X,'X(',i2,') = ',F16.6)
         END DO

      END IF error_check
   END IF size_ok
ELSE fileopen

   ! Else file open failed.  Tell user.
   WRITE (*,1080) istat
   1080 FORMAT (1X,'File open failed--status = ', I6)

END IF fileopen
END PROGRAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE simul ( a, b, ndim, n, er )

IMPLICIT NONE

! Declare calling arguments:
INTEGER, INTENT(IN) :: ndim          ! Dimension of arrays a and b
REAL, INTENT(INOUT), DIMENSION(ndim,ndim) :: a 
                                     ! Array of coefficients (n x n).
                                     ! This array is of size ndim x 
                                     ! ndim, but only n x n of the 
                                     ! coefficients are being used.
                                     ! The declared dimension ndim 
                                     ! must be passed to the sub, or
                                     ! it won't be able to interpret
                                     ! subscripts correctly.  (This
                                     ! array is destroyed during
                                     ! processing.)
REAL, INTENT(INOUT), DIMENSION(ndim) :: b 
                                     ! Input: Right-hand side of eqns.
                                     ! Output: Solution vector.
INTEGER, INTENT(IN) :: n             ! Number of equations to solve.
INTEGER, INTENT(OUT) :: er           ! Error flag:
                                     !   0 -- No error
                                     !   1 -- Singular equations


! Declare local parameters
REAL, PARAMETER :: epsilon = 1.0E-6  ! A "small" number for comparison
                                     ! when determining singular eqns 

! Declare local variables:
REAL :: factor                       ! Factor to multiply eqn irow by
                                     ! before adding to eqn jrow
INTEGER :: irow                      ! Number of the equation currently
                                     ! being processed
INTEGER :: ipeak                     ! Pointer to equation containing
                                     ! maximum pivot value
INTEGER :: jrow                      ! Number of the equation compared
                                     ! to the current equation
INTEGER :: kcol                      ! Index over all columns of eqn
REAL :: temp                         ! Scratch value

! Process n times to get all equations...
mainloop: DO irow = 1, n
 
   ! Find peak pivot for column irow in rows irow to n
   ipeak = irow
max_pivot: DO jrow = irow+1, n
      IF (ABS(a(jrow,irow)) > ABS(a(ipeak,irow))) THEN
         ipeak = jrow
      END IF
   END DO max_pivot
 
   ! Check for singular equations.  
   singular: IF ( ABS(a(ipeak,irow)) < epsilon ) THEN
      er = 1
      RETURN
   END IF singular
 
   ! Otherwise, if ipeak /= irow, swap equations irow & ipeak
   swap_eqn: IF ( ipeak /= irow ) THEN
      DO kcol = 1, n
         temp          = a(ipeak,kcol)
         a(ipeak,kcol) = a(irow,kcol)
         a(irow,kcol)  = temp 
      END DO
      temp     = b(ipeak)
      b(ipeak) = b(irow)
      b(irow)  = temp 
   END IF swap_eqn
 
   ! Multiply equation irow by -a(jrow,irow)/a(irow,irow),  
   ! and add it to Eqn jrow (for all eqns except irow itself).
eliminate: DO jrow = 1, n
      IF ( jrow /= irow ) THEN
         factor = -a(jrow,irow)/a(irow,irow)
         DO kcol = 1, n
            a(jrow,kcol) = a(irow,kcol)*factor + a(jrow,kcol)
         END DO
         b(jrow) = b(irow)*factor + b(jrow)
      END IF
   END DO eliminate
END DO mainloop
   
! End of main loop over all equations.  All off-diagonal
! terms are now zero.  To get the final answer, we must
! divide each equation by the coefficient of its on-diagonal
! term.
divide: DO irow = 1, n
   b(irow)      = b(irow) / a(irow,irow)
   a(irow,irow) = 1.
END DO divide
 
! Set error flag to 0 and return.
er = 0
END SUBROUTINE simul






