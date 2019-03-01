      SUBROUTINE DIFF(
     &                 NLx, NLy, u, u_ex, MaxNorm, MPI_COMM
     &               )

C ==========================================================================C
C  --------------------
C   DESCRIPTION:
C  --------------------
C     
C     Compute the MaxNorm of the error.  Since this is rigged to
C     to give the exact solution, the MaxNorm should be on the same
C     order as the tolerance that we chose for the iterative solve.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C ---------------------------
C    Includes:
C
      INCLUDE 'mpif.h'
      INCLUDE 'BMG_constants.h'

C ---------------------------
C    Argument Declarations:
C     
      INTEGER  NLx, NLy, MPI_COMM

      REAL*8   MaxNorm, u_ex(0:NLx+1,0:NLy+1), u(0:NLx+1,0:NLy+1) 

C ---------------------------
C    Local Declarataions:
C     
      REAL*8   AbsErr, TMP_Buffer
      INTEGER  i, j,  MPI_IERR
      
C ==========================================================================

      MaxNorm = rZERO
      DO  j = 1, NLy
         DO i = 1, NLx
            AbsErr = ABS(u(i,j)-u_ex(i,j))
            IF ( MaxNorm .LT. AbsErr ) THEN
               MaxNorm = AbsErr
            ENDIF
         ENDDO
      ENDDO

      CALL MPI_Allreduce(
     &         MaxNorm, TMP_Buffer, 1, MPI_DOUBLE_PRECISION, 
     &         MPI_MAX, MPI_COMM, MPI_IERR
     &         )

      MaxNorm = TMP_Buffer

C ==========================================================================

      RETURN
      END
