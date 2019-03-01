      REAL*RKIND FUNCTION dz( i, j, k, hx, hy, hz )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Sets the (3,3) entry of (diagonal) diffusion tensor.
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
      INTEGER  i, j, k
      REAL*RKIND   hx, hy, hz

C ----------------
C     Local
C
      INTEGER  n
      REAL*RKIND   dxlyb, dxlyt, dxryb, dxryt, eps, xl, xr, yb, yt, z

C -----------------
C     Includes
C
#include     'common2.h'

C ==========================================================================

      DATA eps /1.e-08/
      dz=diz(1)
      if(ireg.eq.1)return
      xl=(i-2.5)*hx
      xr=xl+hx
      yb=(j-2.5)*hy
      yt=yb+hy
      z=(k-2.5)*hz
      dxlyb=diz(1)
      do 10 n=2,ireg
      if(x1(n)-eps.gt.xl.or.xl.gt.x2(n)+eps.or.y1(n)-eps.gt.yb.or.
     1yb.gt.y2(n)+eps.or.z1(n)-eps.gt.z.or.z.gt.z2(n)+eps)go to 10
      dxlyb=diz(n)
      go to 20
   10 continue
   20 dxlyt=diz(1)
      do 30 n=2,ireg
      if(x1(n)-eps.gt.xl.or.xl.gt.x2(n)+eps.or.y1(n)-eps.gt.yt.or.
     1yt.gt.y2(n)+eps.or.z1(n)-eps.gt.z.or.z.gt.z2(n)+eps)go to 30
      dxlyt=diz(n)
      go to 40
   30 continue
   40 dxryb=diz(1)
      do 50 n=2,ireg
      if(x1(n)-eps.gt.xr.or.xr.gt.x2(n)+eps.or.y1(n)-eps.gt.yb.or.
     1yb.gt.y2(n)+eps.or.z1(n)-eps.gt.z.or.z.gt.z2(n)+eps)go to 50
      dxryb=diz(n)
      go to 60
   50 continue
   60 dxryt=diz(1)
      do 70 n=2,ireg
      if(x1(n)-eps.gt.xr.or.xr.gt.x2(n)+eps.or.y1(n)-eps.gt.yt.or.
     2yt.gt.y2(n)+eps.or.z1(n)-eps.gt.z.or.z.gt.z2(n)+eps)go to 70
      dxryt=diz(n)
      go to 80
   70 continue
   80 dz=.25*(dxlyb+dxlyt+dxryb+dxryt)

C ==========================================================================

      RETURN
      END
