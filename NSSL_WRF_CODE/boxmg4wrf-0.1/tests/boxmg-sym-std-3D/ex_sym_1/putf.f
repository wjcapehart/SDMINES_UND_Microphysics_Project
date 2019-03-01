      SUBROUTINE PUTF( SO, QF, Q, ii, jj, kk, hx, hy, hz, kg )

      IMPLICIT NONE

C ==========================================================================

C ----------------
C     Parameters
C
      INTEGER   kb, kp, kpw, kps
      PARAMETER ( kp=1, kpw=2, kps=3, kb=4 )

C ----------------
C     Arguments
C
      INTEGER  ii, jj, kk, kg
      REAL*8   hx, hy, hz, SO(ii,jj,kk,4), QF(ii,jj,kk), Q(ii,jj,kk)

C ----------------
C     Local
C
      INTEGER  i, i1, i2, j, j1, j2, k, k1
      REAL*8   cm, fs, h2, s, xh, yh, zh

C -----------------
C     External
C
      REAL*8   da, dr, dz
      EXTERNAL da, dr, dz

C -----------------
C     Includes
C
      INCLUDE 'common2.h'

C ==========================================================================

C
C   Zero everything
C
      DO k=1,kk
         DO j=1,jj
            DO i=1,ii
               SO(i,j,k,kp) = 0
               SO(i,j,k,kpw) = 0
               SO(i,j,k,kps) = 0
               SO(i,j,k,kb) = 0
               Q(i,j,k)=0
               QF(i,j,k)=0
            END DO
         END DO
      END DO
c
c   initialize right hand side and difference operator on grid k.
c
      h2=hx*hy*hz
      xh=hy*hz/hx
      yh=hx*hz/hy
      zh=hx*hy/hz
      i1=ii-1
      j1=jj-1
      i2=ii-2
      j2=jj-2
      k1=kk-1
      do 10 k=2,k1
      do 10 j=3,j1
      do 10 i=2,i1
      so(i,j,k,kps)=da(i,j,k,hx,hy,hz)*yh
   10 continue
      do 20 k=2,k1
      do 20 j=2,j1
      do 20 i=3,i1
      so(i,j,k,kpw)=dr(i,j,k,hx,hy,hz)*xh
   20 continue
      do 25 k=3,k1
      do 25 i=2,i1
      do 25 j=2,j1
      so(i,j,k,kb)=dz(i,j,k,hx,hy,hz)*zh
   25 continue
      do 30 k=2,k1
      do 30 j=2,j1
      so(2,j,k,kps)=.5*so(2,j,k,kps)
      so(i1,j,k,kps)=.5*so(i1,j,k,kps)
      so(2,j,k,kb)=.5*so(2,j,k,kb)
      so(i1,j,k,kb)=.5*so(i1,j,k,kb)
   30 continue
      do 40 i=2,i1
      do 40 k=2,k1
      so(i,2,k,kpw)=.5*so(i,2,k,kpw)
      so(i,j1,k,kpw)=.5*so(i,j1,k,kpw)
      so(i,2,k,kb)=.5*so(i,2,k,kb)
      so(i,j1,k,kb)=.5*so(i,j1,k,kb)
   40 continue
      do 45 j=2,j1
      do 45 i=2,i1
      so(i,j,2,kpw)=.5*so(i,j,2,kpw)
      so(i,j,k1,kpw)=.5*so(i,j,k1,kpw)
   45 continue
      do 46 j=2,j1
      do 46 i=2,i1
      so(i,j,2,kps)=.5*so(i,j,2,kps)
      so(i,j,k1,kps)=.5*so(i,j,k1,kps)
   46 continue
      do 50 k=2,k1
      do 50 j=2,j1
      do 50 i=2,i1
      call rhs (i,j,k,hx,hy,hz,fs,s)
      qf(i,j,k)=fs*h2
      so(i,j,k,kp)=s*h2
c     q(i,j,k)=1.
   50 continue
      if (ibl.eq.0) go to 70
      do 60 j=2,j1
      do 60 k=2,k1
      so(2,j,k,kp)=so(2,j,k,kp)+hy*hz
   60 continue
   70 if (ibr.eq.0) go to 90
      do 80 k=2,k1
      do 80 j=2,j1
      so(i1,j,k,kp)=so(i1,j,k,kp)+hy*hz
   80 continue
   90 if (ibyb.eq.0) go to 110
      do 100 k=2,k1
      do 100 i=2,i1
      so(i,2,k,kp)=so(i,2,k,kp)+hx*hz
  100 continue
  110 if (ibyt.eq.0) go to 130
      do 120 k=2,k1
      do 120 i=2,i1
      so(i,j1,k,kp)=so(i,j1,k,kp)+hx*hz
  120 continue
  130 if(ibzb.eq.0)go to 150
      do 140 j=2,j1
      do 140 i=2,i1
      so(i,j,2,kp)=so(i,j,2,kp)+hx*hy
  140 continue
  150 if(ibzt.eq.0)go to 170
      do 160 j=2,j1
      do 160 i=2,i1
      so(i,j,k1,kp)=so(i,j,k1,kp)+hx*hy
  160 continue
  170 do 180 k=2,k1
      do 180 j=2,j1
      do 180 i=2,i1
      cm=1.
      if (i.eq.2.or.i.eq.i1) cm=.5*cm
      if (j.eq.2.or.j.eq.j1) cm=.5*cm
      if(k.eq.2.or.k.eq.k1)cm=.5*cm
      so(i,j,k,kp)=cm*so(i,j,k,kp)+so(i,j,k,kpw)+so(i,j+1,k,kps)
     1+so(i+1,j,k,kpw)+so(i,j,k,kps)+so(i,j,k,kb)+so(i,j,k+1,kb)
      qf(i,j,k)=cm*qf(i,j,k)
  180 continue
      return
      end

