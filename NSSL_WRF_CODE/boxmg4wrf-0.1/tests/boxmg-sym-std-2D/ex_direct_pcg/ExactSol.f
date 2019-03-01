      SUBROUTINE ExactSol( u_ex, NLx, NLy, iGs, jGs, x, y )

C ==========================================================================C
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Create the array containing the exact solution.  When things
C     are going badly, setting this to just a constant is a good way
C     to track down where things are getting derailed.
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
      INCLUDE  'BMG_stencils.h'
      INCLUDE  'BMG_constants.h'

C ---------------------------
C    Argument Declarations:
C     
      INTEGER   NLx, NLy, iGs, jGs
      REAL*8    u_ex(0:NLx+1,0:NLy+1), x(0:NLx+1), y(0:NLy+1)

C ---------------------------
C    Local Declarations:
C     
      INTEGER   i, j

C ==========================================================================

      DO j=0, NLy+1
         DO i=0, NLx+1
            u_ex(i,j) =  x(i)**2 + y(j)**2
C           u_ex(i,j) =  (x(i) + y(j))/2
C           u_ex(i,j) = 1
         ENDDO
      ENDDO

C ==========================================================================

      RETURN
      END
