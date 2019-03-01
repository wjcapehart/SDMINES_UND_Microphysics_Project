      SUBROUTINE RHS ( i, j, k, hx, hy, hz, fs, s )

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
      INTEGER  i, j, k
      REAL*8   hx, hy, hz, fs, s

C ----------------
C     Local
C
      INTEGER  n
      REAL*8   eps,  
     &         flbb, flbt, fltb, fltt, frbb, frtb, frbt, frtt,
     &         slbb, slbt, sltb, sltt, srbb, srbt, srtb, srtt, 
     &         xl, xr, yb, yt, zb, zt

C -----------------
C     Includes
C
      INCLUDE 'common2.h'

C ==========================================================================

      DATA eps /1.e-08/

      fs=fi(1)
      s=si(1)
      if (ireg.eq.1) return
      xl=(i-2.5)*hx
      xr=xl+hx
      yb=(j-2.5)*hy
      yt=yb+hy
      flbb=fi(1)
      slbb=si(1)
      zb=(k-2.5)*hz
      zt=zb+hz
      do 10 n=2,ireg
      if (x1(n)-eps.gt.xl.or.xl.gt.x2(n)+eps.or.y1(n)-eps.gt.yb.or.yb.gt
     1.y2(n)+eps.or.z1(n)-eps.gt.zb.or.zb.gt.z2(n)+eps)go to 10
      flbb=fi(n)
      slbb=si(n)
      go to 20
   10 continue
   20 fltb=fi(1)
      sltb=si(1)
      do 30 n=2,ireg
      if (x1(n)-eps.gt.xl.or.xl.gt.x2(n)+eps.or.y1(n)-eps.gt.yt.or.yt.gt
     1.y2(n)+eps.or.z1(n)-eps.gt.zb.or.zb.gt.z2(n)+eps)go to 30
      fltb=fi(n)
      sltb=si(n)
      go to 40
   30 continue
   40 frbb=fi(1)
      srbb=si(1)
      do 50 n=2,ireg
      if (x1(n)-eps.gt.xr.or.xr.gt.x2(n)+eps.or.y1(n)-eps.gt.yb.or.yb.gt
     1.y2(n)+eps.or.z1(n)-eps.gt.zb.or.zb.gt.z2(n)+eps)go to 50
      frbb=fi(n)
      srbb=si(n)
      go to 60
   50 continue
   60 frtb=fi(1)
      srtb=si(1)
      do 70 n=2,ireg
      if (x1(n)-eps.gt.xr.or.xr.gt.x2(n)+eps.or.y1(n)-eps.gt.yt.or.yt.gt
     1.y2(n)+eps.or.z1(n)-eps.gt.zb.or.zb.gt.z2(n)+eps)go to 70
      frtb=fi(n)
      srtb=si(n)
      go to 80
   70 continue
   80 flbt=fi(1)
      slbt=si(1)
      do 90 n=2,ireg
      if(x1(n)-eps.gt.xl.or.xl.gt.x2(n)+eps.or.y1(n)-eps.gt.yb.or.
     1yb.gt.y2(n)+eps.or.z1(n)-eps.gt.zt.or.zt.gt.z2(n)+eps)go to 100
      flbt=fi(n)
      slbt=si(n)
      go to 100
   90 continue
  100 fltt=fi(1)
      sltt=si(1)
      do 110 n=2,ireg
      if(x1(n)-eps.gt.xl.or.xl.gt.x2(n)+eps.or.y1(n)-eps.gt.yt.or.yt.
     1gt.y2(n)+eps.or.z1(n)-eps.gt.zt.or.zt.gt.z2(n)+eps)go to 110
      fltt=fi(n)
      sltt=si(n)
      go to 120
  110 continue
  120 frbt=fi(1)
      srbt=si(1)
      do 130 n=2,ireg
      if(x1(n)-eps.gt.xr.or.xr.gt.x2(n)+eps.or.y1(n)-eps.gt.yb.or.yb.
     1gt.y2(n)+eps.or.z1(n)-eps.gt.zt.or.zt.gt.z2(n)+eps)go to 130
      frbt=fi(n)
      srbt=si(n)
      go to 140
  130 continue
  140 frtt=fi(1)
      srtt=si(1)
      do 150 n=2,ireg
      if(x1(n)-eps.gt.xr.or.xr.gt.x2(n)+eps.or.y1(n)-eps.gt.yt.or.yt.
     1gt.y2(n)+eps.or.z1(n)-eps.gt.zt.or.zt.gt.z2(n)+eps)go to 150
      frtt=fi(n)
      srtt=si(n)
      go to 160
  150 continue
  160 fs=.125*(flbb+fltb+frbb+frtb+flbt+fltt+frbt+frtt)
      s=.125*(slbb+sltb+srbb+srtb+slbt+sltt+srbt+srtt)
      return
      end
