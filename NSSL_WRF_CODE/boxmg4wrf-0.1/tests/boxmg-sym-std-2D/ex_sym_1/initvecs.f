      subroutine initvecs ( u, v, L, M )

C =========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_constants.h'

C -----------------------------
C     Argument Declarations
C 
      integer  L, M
      real*8   u(0:L+1,0:M+1), v(0:L+1,0:M+1)

C -----------------------------
C     Local Declarations
C
      INTEGER  i,j

C -----------------------------
C     External Functions
C
      REAL*8   rand
      EXTERNAL rand

C =========================================================================

      DO i=0,L+1
         DO j=0,M+1
            IF (i.eq.0.or.i.eq.L+1.or.j.eq.0.or.j.eq.M+1) THEN
               u(i,j) = 0.0
               v(i,j) = 0.0
            ELSE
               u(i,j) = rand(rZERO)
               v(i,j) = rand(rZERO)
            ENDIF
         ENDDO
      ENDDO

C =========================================================================

      return
      end
