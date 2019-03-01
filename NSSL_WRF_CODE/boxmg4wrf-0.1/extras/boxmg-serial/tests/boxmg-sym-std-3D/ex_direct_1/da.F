      REAL*RKIND FUNCTION da( i, j, k, hx, hy, hz)

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

C ==========================================================================

C ----------------
C     Arguments
C
      INTEGER  i, j, k
      REAL*RKIND   hx, hy, hz

C ----------------
C     Local
C
      INTEGER  n
      REAL*RKIND dxlzb, dxlzt, dxrzb, dxrzt, eps, xl, xr, y, zb, zt 

C -----------------
C     Includes
#include     'common2.h'

C ==========================================================================

      DATA eps /1.e-08/

      da=diy(1)
      if (ireg.eq.1) return
      xl=(i-2.5)*hx
      xr=xl+hx
      y=(j-2.5)*hy
      zb=(k-2.5)*hz
      zt=zb+hz
      dxlzb=diy(1)
      do 10 n=2,ireg
      if (x1(n)-eps.gt.xl.or.xl.gt.x2(n)+eps.or.y1(n)-eps.gt.y.or.y.gt
     1.y2(n)+eps.or.z1(n)-eps.gt.zb.or.zb.gt.z2(n)+eps)go to 10
      dxlzb=diy(n)
      go to 20
   10 continue
   20 dxrzb=diy(1)
      do 30 n=2,ireg
      if (x1(n)-eps.gt.xr.or.xr.gt.x2(n)+eps.or.y1(n)-eps.gt.y.or.y.gt
     1.y2(n)+eps.or.z1(n)-eps.gt.zb.or.zb.gt.z2(n)+eps)go to 30
      dxrzb=diy(n)
      go to 40
   30 continue
   40 dxlzt=diy(1)
      do 50 n=2,ireg
      if(x1(n)-eps.gt.xl.or.xl.gt.x2(n)+eps.or.y1(n)-eps.gt.y.or.
     1y.gt.y2(n)+eps.or.z1(n)-eps.gt.zt.or.zt.gt.z2(n)+eps)go to 50
      dxlzt=diy(n)
      go to 60
   50 continue
   60 dxrzt=diy(1)
      do 70 n=2,ireg
      if(x1(n)-eps.gt.xr.or.xr.gt.x2(n)+eps.or.y1(n)-eps.gt.y.or.
     1y.gt.y2(n)+eps.or.z1(n)-eps.gt.zt.or.zt.gt.z2(n)+eps)go to 70
      dxrzt=diy(n)
      go to 80
   70 continue
   80 da=.25*(dxlzb+dxrzb+dxlzt+dxrzt)

C ==========================================================================

      return
      end
