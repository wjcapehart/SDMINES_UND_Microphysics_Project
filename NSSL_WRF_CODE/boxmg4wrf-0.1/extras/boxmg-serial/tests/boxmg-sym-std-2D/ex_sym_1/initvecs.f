      subroutine initvecs ( u, v, L, M )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Initalize vectors to use in the inner product check of symmetry.
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
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
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

C -----------------------------
C     Includes
C
      INCLUDE 'BMG_SER_constants.h'

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
      REAL*8   BMG_rand
      EXTERNAL BMG_rand

C ==========================================================================

      DO i=0,L+1
         DO j=0,M+1
            IF (i.eq.0.or.i.eq.L+1.or.j.eq.0.or.j.eq.M+1) THEN
               u(i,j) = 0.0
               v(i,j) = 0.0
            ELSE
               u(i,j) = BMG_rand(rZERO)
               v(i,j) = BMG_rand(rZERO)
            ENDIF
         ENDDO
      ENDDO

C ==========================================================================

      return
      end