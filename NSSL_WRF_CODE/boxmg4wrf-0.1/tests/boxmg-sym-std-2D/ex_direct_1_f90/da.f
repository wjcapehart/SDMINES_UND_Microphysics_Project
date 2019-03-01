      REAL*8 FUNCTION da( i, j, hx, hy )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Sets the (2,2) entry of (diagonal) diffusion tensor.
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

C ----------------
C     Arguments
C
      INTEGER  i, j
      REAL*8   hx, hy

C ----------------
C     Local
C
      INTEGER  n
      REAL*8   dxl, dxr, eps, xl, xr, y

C -----------------
C     Includes
C
      INCLUDE 'common2.h'

C ==========================================================================

      DATA eps /1.e-08/

      da = diy(1)

      IF ( ireg.EQ.1 ) RETURN

      xl = (i-2.5)*hx
      xr = xl+hx

      y = (j-2.5)*hy

      dxl =diy(1)

      DO 10 n=2,ireg
         IF ( x1(n)-eps.GT.xl .OR. xl.GT.x2(n)+eps  .OR.
     &        y1(n)-eps.GT.y  .OR. y.GT.y2(n)+eps ) GO TO 10
         dxl = diy(n)
         GO TO 20
 10   CONTINUE

 20   dxr=diy(1)

      DO 30 n=2,ireg
         IF ( x1(n)-eps.GT.xr .OR. xr.GT.x2(n)+eps .or. 
     &        y1(n)-eps.GT.y  .OR. y.GT.y2(n)+eps ) GO TO 30
         dxr = diy(n)
         GO TO 40
 30   CONTINUE

 40   da = 0.5*( dxl+dxr )

C ==========================================================================

      RETURN
      END
