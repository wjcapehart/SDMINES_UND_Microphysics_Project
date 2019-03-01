      REAL*8 FUNCTION dr( i, j, k, hx, hy, hz )

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
      REAL*8   hx, hy, hz

C ----------------
C     Local
C
      INTEGER  n
      REAL*8   dybzb, dybzt, dytzb, dytzt, eps, xo, yb, yt, zb, zt 

C -----------------
C     Includes
C
      INCLUDE 'common2.h'

C ==========================================================================

      DATA eps /1.e-08/

      dr=dix(1)
      if (ireg.eq.1) return
      xo=(i-2.5)*hx
      yb=(j-2.5)*hy
      yt=yb+hy
      zb=(k-2.5)*hz
      zt=zb+hz
      dybzb=dix(1)
      do 10 n=2,ireg
      if (x1(n)-eps.gt.xo.or.xo.gt.x2(n)+eps.or.y1(n)-eps.gt.yb.or.yb.gt
     1.y2(n)+eps.or.z1(n)-eps.gt.zb.or.zb.gt.z2(n)+eps)go to 10
      dybzb=dix(n)
      go to 20
   10 continue
   20 dytzb=dix(1)
      do 30 n=2,ireg
      if (x1(n)-eps.gt.xo.or.xo.gt.x2(n)+eps.or.y1(n)-eps.gt.yt.or.yt.gt
     1.y2(n)+eps.or.z1(n)-eps.gt.zb.or.zb.gt.z2(n)+eps)go to 30
      dytzb=dix(n)
      go to 40
   30 continue
   40 dybzt=dix(1)
      do 50 n=2,ireg
      if(x1(n)-eps.gt.xo.or.xo.gt.x2(n)+eps.or.y1(n)-eps.gt.yb.or.
     1yb.gt.y2(n)+eps.or.z1(n)-eps.gt.zt.or.zt.gt.z2(n)+eps)go to 50
      dybzt=dix(n)
      go to 60
   50 continue
   60 dytzt=dix(1)
      do 70 n=2,ireg
      if(x1(n)-eps.gt.xo.or.xo.gt.x2(n)+eps.or.y1(n)-eps.gt.yt.or.
     1yt.gt.y2(n)+eps.or.z1(n)-eps.gt.zt.or.zt.gt.z2(n)+eps)go to 70
      dytzt=dix(n)
   70 continue
   80 dr=.25*(dybzb+dytzb+dybzt+dytzt)

C ==========================================================================

      return
      end
