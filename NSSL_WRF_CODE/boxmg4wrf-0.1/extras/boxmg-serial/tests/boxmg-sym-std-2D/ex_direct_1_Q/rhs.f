      SUBROUTINE RHS( i, j, hx, hy, fs, s )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Sets the source and removal terms.
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

C ==========================================================================

C ----------------
C     Arguments
C
      INTEGER  i, j
      REAL*8   hx, hy, fs, s

C ----------------
C     Local
C
      INTEGER  n
      REAL*8   eps,  
     &         flb, flt, frb, frt, 
     &         slb, slt, srb, srt,
     &         xl, xr, yb, yt

C -----------------
C     Includes
C
      INCLUDE 'common2.h'

C ==========================================================================

      DATA eps /1.e-08/

      fs =fi(1)
      s = si(1)

      IF ( ireg.EQ.1 ) RETURN

      xl = (i-2.5)*hx
      xr = xl+hx
      yb = (j-2.5)*hy
      yt = yb+hy

      flb=fi(1)
      slb=si(1)

      DO 10 n=2,ireg
         IF ( x1(n)-eps.GT.xl .OR. xl.GT.x2(n)+eps .OR.
     &        y1(n)-eps.GT.yb .OR. yb.GT.y2(n)+eps) GO TO 10
         flb = fi(n)
         slb = si(n)
         GO TO 20
 10   CONTINUE

   20 flt = fi(1)
      slt = si(1)

      DO 30 n=2,ireg
         IF ( x1(n)-eps.GT.xl .OR. xl.GT.x2(n)+eps  .OR.
     &        y1(n)-eps.GT.yt .OR. yt.GT.y2(n)+eps ) go to 30
         flt = fi(n)
         slt = si(n)
         GO TO 40
 30   CONTINUE

   40 frb = fi(1)
      srb = si(1)

      DO 50 n=2,ireg
         IF ( x1(n)-eps.GT.xr .OR. xr.GT.x2(n)+eps  .OR.
     &        y1(n)-eps.GT.yb .OR. yb.GT.y2(n)+eps) GO TO 50
         frb = fi(n)
         srb = si(n)
         GO TO 60
 50   CONTINUE

   60 frt = fi(1)
      srt = si(1)

      DO 70 n=2,ireg
         IF ( x1(n)-eps.GT.xr .OR. xr.GT.x2(n)+eps  .or.
     &        y1(n)-eps.GT.yt .OR. yt.gt.y2(n)+eps) GO TO 70
         frt = fi(n)
         srt = si(n)
         GO TO 80
 70   CONTINUE

 80   fs = 0.25*( flb + flt + frb + frt )
      s  = 0.25*( slb + slt + srb + srt )

C ==========================================================================

      RETURN
      END
