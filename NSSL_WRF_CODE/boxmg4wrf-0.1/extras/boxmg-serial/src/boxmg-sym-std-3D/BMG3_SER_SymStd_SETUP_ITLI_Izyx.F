      SUBROUTINE BMG3_SER_SymStd_SETUP_ITLI_Izyx( 
     &                kgf, kgc, so, soc, sor, ci,
     &                iif, jjf, kkf, iic, jjc, kkc, 
     &                NOG, ifd, NStncl, irelax, zo, yo 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_SETUP_ITLI_Izyx computes an approximation to
C     the true variational coarse-grid operator.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     KGF       Fine grid number.
C     KGC       Coarse grid number.
C     IIF       Number. of grid points in x direction on fine grid,
C               including two fictitious points.
C     JJF       Number. of grid points in y direction on fine grid,
C               including two fictitious points.
C     KKF       Number. of grid points in z direction on fine grid,
C               including two fictitious points.
C     IIC       Number. of grid points in x direction on coarse grid,
C               including two fictitious points.
C     JJC       Number. of grid points in y direction on coarse grid,
C               including two fictitious points.
C     KKC       Number. of grid points in z direction on coarse grid,
C               including two fictitious points.
C
C     NOG       Refer to BMG3_SER_SymStd_SOLVE_boxmg.
C     IFD       Refer to BMG3_SER_SymStd_SOLVE_boxmg.
C     IRELAX    Refer to BMG3_SER_SymStd_SOLVE_boxmg.
C     SO        Refer to BMG3_SER_SymStd_SOLVE_boxmg.
C     SOR       Refer to BMG3_SER_SymStd_SOLVE_boxmg.
C     CI        Refer to BMG3_SER_SymStd_SOLVE_boxmg.
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
C     SOC       SO for coarse grid
C     SORC      SOR for coarse grid
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
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_stencils.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  ifd, iic, iif, irelax, jjc, jjf, kgc, kgf, kkc, kkf,
     &         NOG, NStncl

      REAL*RKIND   ci(iic,jjc,kkc,26),
     &         so(iif,jjf,kkf,NStncl), soc(iic,jjc,kkc,14), 
     &         sor(iif,jjf,kkf,2), 
     &         yo(iif,jjc,kkc+1,14), zo(iif,jjf,kkc,14)

C ----------------------------
C     Local Declarations
C
      INTEGER  ic, i, iicf, iicf1, iic1, iif1, jc, j, jjcf, jjcf1, 
     &         jjc1, jjf1, kc, k, kkcf, kkc1, kkcf1, kkcp1, kkf1, kpz
      REAL*RKIND   a, b, c, eMACH, ep 

C ==========================================================================

      eMACH = 1.d-13

      iic1 = iic-1
      jjc1 = jjc-1
      kkc1 = kkc-1

      iif1 = iif-1
      jjf1 = jjf-1
      kkf1 = kkf-1

      iicf = (iif-2)/2+3
      jjcf = (jjf-2)/2+3
      kkcf = (kkf-2)/2+3

      iicf1 = iicf-1
      jjcf1 = jjcf-1
      kkcf1 = kkcf-1

      kkcp1 = kkc+1



      do kpz=1,14
         do k=1,kkcp1
            do j=1,jjc
               do i=1,iif
                  yo(i,j,k,kpz)= rZERO
               enddo
            enddo
         enddo
      enddo

      do kpz=1,14
         do k=1,kkc
            do j=1,jjf
               do i=1,iif
                  zo(i,j,k,kpz)= rZERO
               enddo
            enddo
         enddo
      enddo


c
c   compute interpolation operator iz, from grid kgfz to grid kgf, where
c   grid kgfz is obtained form grid kgf by coarsening in the z-direction
c   only.
c
      if(kgf.lt.NOG.or.ifd.ne.1) then 
         k=2
         do kc=3,kkcf1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  a=so(i,j+1,k-1,kbnw)+so(i,j+1,k-1,kbn)
     &                 +so(i+1,j+1,k-1,kbne)+so(i,j,k-1,kbw)
     &                 +so(i,j,k-1,kb)+so(i+1,j,k-1,kbe)
     &                 +so(i,j,k-1,kbsw)+so(i,j,k-1,kbs)
     &                 +so(i+1,j,k-1,kbse)
                  b=so(i,j+1,k,kbse)+so(i,j+1,k,kbs)
     &                 +so(i+1,j+1,k,kbsw)+so(i,j,k,kbe)+so(i,j,k,kb)
     &                 +so(i+1,j,k,kbw)+so(i,j,k,kbne)+so(i,j,k,kbn)
     &                 +so(i+1,j,k,kbnw)
                  c=a+b+so(i,j+1,k-1,kpnw)+so(i,j+1,k-1,kps)
     &                 +so(i+1,j+1,k-1,kpsw)+so(i,j,k-1,kpw)
     &                 +so(i+1,j,k-1,kpw)+so(i,j,k-1,kpsw)
     &                 +so(i,j,k-1,kps)+so(i+1,j,k-1,kpnw)
                  ep=MIN(abs(a),abs(b),rONE)
                  c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)
     &                 -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)
     &                 -(rONE+ep)*c)+eMACH)
                  sor(i,j,k-1,mtot)=a/c
                  sor(i,j,k,mtot)=b/c
               enddo
            enddo
         enddo
      else ! if kgf.ge.NOG.and.ifd.eq.1   
         k=2
         do kc=3,kkcf1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  a=so(i,j,k-1,kb)
                  b=so(i,j,k,kb)
                  c=a+b+so(i,j+1,k-1,kps)+so(i,j,k-1,kpw)
     &                 +so(i+1,j,k-1,kpw)+so(i,j,k-1,kps)
                  ep=MIN(abs(a),abs(b),rONE)
                  c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)
     &                 -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)
     &                 -(rONE+ep)*c)+eMACH)
                  sor(i,j,k-1,mtot)=a/c
                  sor(i,j,k,mtot)=b/c
               enddo
            enddo
         enddo
      endif
c
c   compute lz = iz(transpose) l iz, where l is the difference operator
c   on grid kgf.
c
      
      if(kgf.lt.NOG.or.ifd.ne.1) then 
         
         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kpw)=
     &                 so(i,j,k,kpw)
     &                 +so(i,j,k+1,kbe)*sor(i-1,j,k+1,mtot)
     &                 +so(i,j,k,kbw)*sor(i-1,j,k,mtot)
     &                 +sor(i,j,k+1,mtot)
     &                 *( so(i,j,k+1,kpw)*sor(i-1,j,k+1,mtot)
     &                 +so(i,j,k+1,kbw))
     &                 +sor(i,j,k,mtot) 
     &                 *( so(i,j,k-1,kpw)*sor(i-1,j,k,mtot)
     &                 +so(i,j,k,kbe))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kpsw)=so(i,j,k,kpsw)+so(i,j,k+1,kbne)
     &                 *sor(i-1,j-1,k+1,mtot)+so(i,j,k,kbsw)
     &                 *sor(i-1,j-1,k,mtot)+sor(i,j,k,mtot)
     &                 *(so(i,j,k-1,kpsw)*sor(i-1,j-1,k,mtot)
     &                 +so(i,j,k,kbne))+sor(i,j,k+1,mtot)
     &                 *(so(i,j,k+1,kpsw)*sor(i-1,j-1,k+1,mtot)
     &                 +so(i,j,k+1,kbsw))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j+1,kc,kpnw)=so(i,j+1,k,kpnw)+so(i,j+1,k,kbnw)
     &                 *sor(i-1,j+1,k,mtot)+so(i,j+1,k+1,kbse)
     &                 *sor(i-1,j+1,k+1,mtot)+sor(i,j,k,mtot)
     &                 *(so(i,j+1,k-1,kpnw)*sor(i-1,j+1,k,mtot)
     &                 +so(i,j+1,k,kbse))+sor(i,j,k+1,mtot)
     &                 *(so(i,j+1,k+1,kpnw)*sor(i-1,j+1,k+1,mtot)
     &                 +so(i,j+1,k+1,kbnw))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kps)=so(i,j,k,kps)+so(i,j,k+1,kbn)
     &                 *sor(i,j-1,k+1,mtot)+so(i,j,k,kbs)
     &                 *sor(i,j-1,k,mtot)+sor(i,j,k,mtot)
     &                 *(so(i,j,k-1,kps)*sor(i,j-1,k,mtot)
     &                 +so(i,j,k,kbn))+sor(i,j,k+1,mtot)
     &                 *(so(i,j,k+1,kps)*sor(i,j-1,k+1,mtot)
     &                 +so(i,j,k+1,kbs))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kp)=so(i,j,k,kp)-so(i,j,k,kb)
     &                 *sor(i,j,k,mtot)-so(i,j,k+1,kb)
     &                 *sor(i,j,k+1,mtot)-sor(i,j,k,mtot)
     &                 *(-so(i,j,k-1,kp)*sor(i,j,k,mtot)
     &                 +so(i,j,k,kb))-sor(i,j,k+1,mtot)
     &                 *(-so(i,j,k+1,kp)*sor(i,j,k+1,mtot)
     &                 +so(i,j,k+1,kb))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kb)=so(i,j,k,kb)*sor(i,j,k-1,mtot)
     &                 +sor(i,j,k,mtot)*(so(i,j,k-1,kb)
     &                 -so(i,j,k-1,kp)*sor(i,j,k-1,mtot))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kbw)=so(i,j,k,kbw)*sor(i-1,j,k-1,mtot)
     &                 +sor(i,j,k,mtot)*(so(i,j,k-1,kpw)
     &                 *sor(i-1,j,k-1,mtot)+so(i,j,k-1,kbw))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j+1,kc,kbnw)=so(i,j+1,k,kbnw)
     &                 *sor(i-1,j+1,k-1,mtot)+sor(i,j,k,mtot)
     &                 *(so(i,j+1,k-1,kpnw)*sor(i-1,j+1,k-1,mtot)
     &                 +so(i,j+1,k-1,kbnw))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j+1,kc,kbn)=so(i,j+1,k,kbn)
     &                 *sor(i,j+1,k-1,mtot)
     &                 +sor(i,j,k,mtot)*(so(i,j+1,k-1,kps)
     &                 *sor(i,j+1,k-1,mtot)+so(i,j+1,k-1,kbn))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i+1,j+1,kc,kbne)=so(i+1,j+1,k,kbne)
     &                 *sor(i+1,j+1,k-1,mtot)+sor(i,j,k,mtot)
     &                 *(so(i+1,j+1,k-1,kpsw)*sor(i+1,j+1,k-1,mtot)
     &                 +so(i+1,j+1,k-1,kbne))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i+1,j,kc,kbe)=so(i+1,j,k,kbe)
     &                 *sor(i+1,j,k-1,mtot)
     &                 +sor(i,j,k,mtot)*(so(i+1,j,k-1,kpw)
     &                 *sor(i+1,j,k-1,mtot)+so(i+1,j,k-1,kbe))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i+1,j,kc,kbse)=so(i+1,j,k,kbse)
     &                 *sor(i+1,j-1,k-1,mtot)+sor(i,j,k,mtot)
     &                 *(so(i+1,j,k-1,kpnw)*sor(i+1,j-1,k-1,mtot)
     &                 +so(i+1,j,k-1,kbse))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kbs)=so(i,j,k,kbs)*sor(i,j-1,k-1,mtot)
     &                 +sor(i,j,k,mtot)*(so(i,j,k-1,kps)
     &                 *sor(i,j-1,k-1,mtot)+so(i,j,k-1,kbs))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kbsw)=so(i,j,k,kbsw)
     &                 *sor(i-1,j-1,k-1,mtot)
     &                 +sor(i,j,k,mtot)*(so(i,j,k-1,kpsw)
     &                 *sor(i-1,j-1,k-1,mtot)+so(i,j,k-1,kbsw))
               enddo
            enddo
         enddo




      else                      ! if kgf.ge.NOG.and.ifd.eq.1                
         



         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kpw)=so(i,j,k,kpw)+sor(i,j,k+1,mtot)
     &                 *so(i,j,k+1,kpw)*sor(i-1,j,k+1,mtot)
     &                 +sor(i,j,k,mtot)*so(i,j,k-1,kpw)
     &                 *sor(i-1,j,k,mtot)
                  zo(i,j,kc,kpsw)= rZERO
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j+1,kc,kpnw)= rZERO
                  zo(i,j,kc,kps)=so(i,j,k,kps)+sor(i,j,k,mtot)
     &                 *so(i,j,k-1,kps)*sor(i,j-1,k,mtot)
     &                 +sor(i,j,k+1,mtot)*so(i,j,k+1,kps)
     &                 *sor(i,j-1,k+1,mtot)
                  zo(i,j+1,kc,kpnw)= rZERO
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kp)=so(i,j,k,kp)-sor(i,j,k,mtot)
     &                 *(-so(i,j,k-1,kp)*sor(i,j,k,mtot)
     &                 +so(i,j,k,kb))-sor(i,j,k+1,mtot)
     &                 *(-so(i,j,k+1,kp)*sor(i,j,k+1,mtot)
     &                 +so(i,j,k+1,kb))-so(i,j,k,kb)*sor(i,j,k,mtot)
     &                 -so(i,j,k+1,kb)*sor(i,j,k+1,mtot)
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kb)=so(i,j,k,kb)*sor(i,j,k-1,mtot)
     &                 +sor(i,j,k,mtot)*(so(i,j,k-1,kb)
     &                 -so(i,j,k-1,kp)*sor(i,j,k-1,mtot))
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kbw)=sor(i,j,k,mtot)*so(i,j,k-1,kpw)
     &                 *sor(i-1,j,k-1,mtot)
                  zo(i,j+1,kc,kbnw)= rZERO
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j+1,kc,kbn)=sor(i,j,k,mtot)*so(i,j+1,k-1,kps)
     &                 *sor(i,j+1,k-1,mtot)
                  zo(i+1,j+1,kc,kbne)= rZERO
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i+1,j,kc,kbe)=sor(i,j,k,mtot)*so(i+1,j,k-1,kpw)
     &                 *sor(i+1,j,k-1,mtot)
                  zo(i+1,j,kc,kbse)= rZERO
               enddo
            enddo
         enddo

         k=0
         do kc=2,kkc1
            k=k+2
            do j=2,jjf1
               do i=2,iif1
                  zo(i,j,kc,kbs)=sor(i,j,k,mtot)*so(i,j,k-1,kps)
     &                 *sor(i,j-1,k-1,mtot)
                  zo(i,j,kc,kbsw)= rZERO
               enddo
            enddo
         enddo

      endif                     ! of if (kgf.lt.NOG.or.ifd.ne.1)
         



c     compute iy, the interpolation operator form grid kgzy to grid
c     kgz, where grid kgz is obtained from grid kgz by coarsening
c     in the y-direction only.
c     

      do kc=3,kkc
         j=2
         do jc=3,jjcf1
            j=j+2
            do i=2,iif1
               a=zo(i,j-1,kc-1,kbsw)+zo(i,j-1,kc-1,kbs)
     &              +zo(i+1,j-1,kc-1,kbse)+zo(i,j-1,kc-1,kpsw)
     &              +zo(i,j-1,kc-1,kps)+zo(i+1,j-1,kc-1,kpnw)
     &              +zo(i,j-1,kc,kbne)+zo(i,j-1,kc,kbn)
     &              +zo(i+1,j-1,kc,kbnw)
               b=zo(i,j,kc-1,kbnw)+zo(i,j,kc-1,kbn)
     &              +zo(i+1,j,kc-1,kbne)+zo(i,j,kc-1,kpnw)
     &              +zo(i,j,kc-1,kps)+zo(i+1,j,kc-1,kpsw)
     &              +zo(i,j,kc,kbse)+zo(i,j,kc,kbs)
     &              +zo(i+1,j,kc,kbsw)
               c=zo(i,j-1,kc-1,kbw)+zo(i,j-1,kc-1,kb)
     &              +zo(i+1,j-1,kc-1,kbe)+zo(i,j-1,kc-1,kpw)
     &              +a+b+zo(i+1,j-1,kc-1,kpw)
     &              +zo(i,j-1,kc,kbe)+zo(i,j-1,kc,kb)
     &              +zo(i+1,j-1,kc,kbw)
               ep=MIN(abs(a),abs(b),rONE)
               c=a+b+(zo(i,j-1,kc-1,kp)-c)*MAX(zo(i,j-1,kc-1,kp)
     &              -(rONE+ep)*c,rZERO)
     &              /(abs(zo(i,j-1,kc-1,kp)-(rONE+ep)*c)+eMACH)
               sor(i,j-1,kc,msor)=a/c
               sor(i,j,kc,msor)=b/c
            enddo
         enddo
      enddo

c     
c     compute ly = iy(transpose) lz iy.
c     


      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i,jc,kc,kp)=zo(i,j,kc-1,kp)-zo(i,j,kc-1,kps)
     &              *sor(i,j,kc,msor)-zo(i,j+1,kc-1,kps)
     &              *sor(i,j+1,kc,msor)-sor(i,j,kc,msor)
     &              *(-zo(i,j-1,kc-1,kp)*sor(i,j,kc,msor)
     &              +zo(i,j,kc-1,kps))-sor(i,j+1,kc,msor)
     &              *(-zo(i,j+1,kc-1,kp)*sor(i,j+1,kc,msor)
     &              +zo(i,j+1,kc-1,kps))
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
             do i=2,iif1
               yo(i,jc,kc,kpw)=zo(i,j,kc-1,kpw)+zo(i,j,kc-1,kpsw)
     &              *sor(i-1,j,kc,msor)+zo(i,j+1,kc-1,kpnw)
     &              *sor(i-1,j+1,kc,msor)+sor(i,j,kc,msor)
     &              *(zo(i,j-1,kc-1,kpw)*sor(i-1,j,kc,msor)
     &              +zo(i,j,kc-1,kpnw))+sor(i,j+1,kc,msor)
     &              *(zo(i,j+1,kc-1,kpw)*sor(i-1,j+1,kc,msor)
     &              +zo(i,j+1,kc-1,kpsw))
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i,jc,kc,kpnw)=zo(i,j-1,kc-1,kpnw)
     &              *sor(i-1,j,kc,msor)
     &              +sor(i,j-1,kc,msor)*(zo(i,j,kc-1,kpnw)
     &              +zo(i,j-1,kc-1,kpw)*sor(i-1,j,kc,msor))
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i,jc,kc,kps)=zo(i,j,kc-1,kps)*sor(i,j-1,kc,msor)
     &              +sor(i,j,kc,msor)*(zo(i,j-1,kc-1,kps)
     &              -zo(i,j-1,kc-1,kp)*sor(i,j-1,kc,msor))
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i,jc,kc,kpsw)=zo(i,j,kc-1,kpsw)
     &              *sor(i-1,j-1,kc,msor)
     &              +sor(i,j,kc,msor)*(zo(i,j-1,kc-1,kpsw)
     &              +zo(i,j-1,kc-1,kpw)*sor(i-1,j-1,kc,msor))
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i,jc,kc,kb)=zo(i,j,kc-1,kb)+sor(i,j,kc,msor)
     &              *(zo(i,j,kc-1,kbn)+zo(i,j-1,kc-1,kb)
     &              *sor(i,j,kc-1,msor))+sor(i,j+1,kc,msor)
     &              *(zo(i,j+1,kc-1,kbs)+zo(i,j+1,kc-1,kb)
     &              *sor(i,j+1,kc-1,msor))+zo(i,j,kc-1,kbs)
     &              *sor(i,j,kc-1,msor)+zo(i,j+1,kc-1,kbn)
     &              *sor(i,j+1,kc-1,msor)
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i,jc,kc,kbw)=zo(i,j,kc-1,kbw)+zo(i,j,kc-1,kbsw)
     &              *sor(i-1,j,kc-1,msor)+zo(i,j+1,kc-1,kbnw)
     &              *sor(i-1,j+1,kc-1,msor)+sor(i,j,kc,msor)
     &              *(zo(i,j,kc-1,kbnw)+zo(i,j-1,kc-1,kbw)
     &              *sor(i-1,j,kc-1,msor))+sor(i,j+1,kc,msor)
     &              *(zo(i,j+1,kc-1,kbsw)+zo(i,j+1,kc-1,kbw)
     &              *sor(i-1,j+1,kc-1,msor))
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i,jc,kc,kbnw)=zo(i,j-1,kc-1,kbnw)
     &              *sor(i-1,j,kc-1,msor)+sor(i,j-1,kc,msor)
     &              *(zo(i,j-1,kc-1,kbw)*sor(i-1,j,kc-1,msor)
     &              +zo(i,j,kc-1,kbnw))
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i,jc,kc,kbn)=zo(i,j-1,kc-1,kbn)
     &              *sor(i,j,kc-1,msor)
     &              +sor(i,j-1,kc,msor)*(zo(i,j,kc-1,kbn)
     &              +zo(i,j-1,kc-1,kb)*sor(i,j,kc-1,msor))
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i+1,jc,kc,kbne)=zo(i+1,j-1,kc-1,kbne)
     &              *sor(i+1,j,kc-1,msor)+sor(i,j-1,kc,msor)
     &              *(zo(i+1,j,kc-1,kbne)+zo(i+1,j-1,kc-1,kbe)
     &              *sor(i+1,j,kc-1,msor))
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i+1,jc,kc,kbe)=zo(i+1,j,kc-1,kbe)
     &              +zo(i+1,j,kc-1,kbse)
     &              *sor(i+1,j,kc-1,msor)+zo(i+1,j+1,kc-1,kbne)
     &              *sor(i+1,j+1,kc-1,msor)+sor(i,j,kc,msor)
     &              *(zo(i+1,j,kc-1,kbne)+zo(i+1,j-1,kc-1,kbe)
     &              *sor(i+1,j,kc-1,msor))+sor(i,j+1,kc,msor)
     &              *(zo(i+1,j+1,kc-1,kbse)+zo(i+1,j+1,kc-1,kbe)
     &              *sor(i+1,j+1,kc-1,msor))
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i+1,jc,kc,kbse)=zo(i+1,j,kc-1,kbse)
     &              *sor(i+1,j-1,kc-1,msor)+sor(i,j,kc,msor)
     &              *(zo(i+1,j-1,kc-1,kbse)+zo(i+1,j-1,kc-1,kbe)
     &              *sor(i+1,j-1,kc-1,msor))
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i,jc,kc,kbs)=zo(i,j,kc-1,kbs)
     &              *sor(i,j-1,kc-1,msor)
     &              +sor(i,j,kc,msor)*(zo(i,j-1,kc-1,kbs)
     &              +zo(i,j-1,kc-1,kb)*sor(i,j-1,kc-1,msor))
            enddo
         enddo
      enddo

      do kc=3,kkc
         j=0
         do jc=2,jjc1
            j=j+2
            do i=2,iif1
               yo(i,jc,kc,kbsw)=zo(i,j,kc-1,kbsw)
     &              *sor(i-1,j-1,kc-1,msor)+sor(i,j,kc,msor)
     &              *(zo(i,j-1,kc-1,kbsw)+zo(i,j-1,kc-1,kbw)
     &              *sor(i-1,j-1,kc-1,msor))
            enddo
         enddo
      enddo
      
      
c     
c     compute ix, the interpolation operator form grid kgzy to grid kgc
c     

      do kc=4,kkcp1
         do jc=2,jjc1
            
            i=2
            do ic=3,iicf1
               i=i+2
               a=yo(i-1,jc+1,kc-1,kbnw)+yo(i-1,jc,kc-1,kbw)
     &              +yo(i-1,jc,kc-1,kbsw)+yo(i-1,jc+1,kc-1,kpnw)
     &              +yo(i-1,jc,kc-1,kpw)+yo(i-1,jc,kc-1,kpsw)
     &              +yo(i-1,jc+1,kc,kbse)+yo(i-1,jc,kc,kbe)
     &              +yo(i-1,jc,kc,kbne)
               b=yo(i,jc+1,kc-1,kbne)+yo(i,jc,kc-1,kbe)
     &              +yo(i,jc,kc-1,kbse)
     &              +yo(i,jc+1,kc-1,kpsw)+yo(i,jc,kc-1,kpw)
     &              +yo(i,jc,kc-1,kpnw)+yo(i,jc+1,kc,kbsw)
     &              +yo(i,jc,kc,kbw)+yo(i,jc,kc,kbnw)
               c=yo(i-1,jc+1,kc-1,kbn)+yo(i-1,jc,kc-1,kb)
     &              +yo(i-1,jc,kc-1,kbs)+yo(i-1,jc+1,kc-1,kps)+a+b
     &              +yo(i-1,jc,kc-1,kps)+yo(i-1,jc+1,kc,kbs)
     &              +yo(i-1,jc,kc,kb)+yo(i-1,jc,kc,kbn)
               ep=MIN(abs(a),abs(b),rONE)
               c=a+b+(yo(i-1,jc,kc-1,kp)-c)*MAX(yo(i-1,jc,kc-1,kp)
     &              -(rONE+ep)*c,rZERO)
     &              /(abs(yo(i-1,jc,kc-1,kp)-(rONE+ep)*c)+eMACH)
               ci(ic,jc,kc-2,lxyl)=a/c
               ci(ic,jc,kc-2,lxyr)=b/c
            enddo
         enddo
         
      enddo
      
c               
c     compute lc = ix (transpose) ly ix.
c     

      
      do kc=4,kkcp1
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc,kc-2,kp)=yo(i,jc,kc-1,kp)-yo(i,jc,kc-1,kpw)
     &              *ci(ic,jc,kc-2,lxyr)-yo(i+1,jc,kc-1,kpw)
     &              *ci(ic+1,jc,kc-2,lxyl)-ci(ic,jc,kc-2,lxyr)
     &              *(-yo(i-1,jc,kc-1,kp)*ci(ic,jc,kc-2,lxyr)
     &              +yo(i,jc,kc-1,kpw))-ci(ic+1,jc,kc-2,lxyl)
     &              *(-yo(i+1,jc,kc-1,kp)*ci(ic+1,jc,kc-2,lxyl)
     &              +yo(i+1,jc,kc-1,kpw))
            enddo
         enddo
      enddo

      do kc=4,kkcp1
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc,kc-2,kpw)=yo(i,jc,kc-1,kpw)
     &              *ci(ic,jc,kc-2,lxyl)
     &              +ci(ic,jc,kc-2,lxyr)*(yo(i-1,jc,kc-1,kpw)
     &              -yo(i-1,jc,kc-1,kp)*ci(ic,jc,kc-2,lxyl))
            enddo
         enddo
      enddo

      do kc=4,kkcp1         
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc+1,kc-2,kpnw)=yo(i,jc+1,kc-1,kpnw)
     &              *ci(ic,jc+1,kc-2,lxyl)+ci(ic,jc,kc-2,lxyr)
     &              *(yo(i-1,jc+1,kc-1,kpnw)+yo(i-1,jc+1,kc-1,kps)
     &              *ci(ic,jc+1,kc-2,lxyl))
            enddo
         enddo
      enddo

      do kc=4,kkcp1         
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc,kc-2,kpsw)=yo(i,jc,kc-1,kpsw)
     &              *ci(ic,jc-1,kc-2,lxyl)+ci(ic,jc,kc-2,lxyr)
     &              *(yo(i-1,jc,kc-1,kpsw)+yo(i-1,jc,kc-1,kps)
     &              *ci(ic,jc-1,kc-2,lxyl))
            enddo
         enddo
      enddo

      do kc=4,kkcp1         
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc,kc-2,kps)=yo(i,jc,kc-1,kps)+yo(i,jc,kc-1,kpsw)
     &              *ci(ic,jc-1,kc-2,lxyr)+yo(i+1,jc,kc-1,kpnw)
     &              *ci(ic+1,jc-1,kc-2,lxyl)+ci(ic,jc,kc-2,lxyr)
     &              *(yo(i,jc,kc-1,kpnw)+yo(i-1,jc,kc-1,kps)
     &              *ci(ic,jc-1,kc-2,lxyr))+ci(ic+1,jc,kc-2,lxyl)
     &              *(yo(i+1,jc,kc-1,kpsw)+yo(i+1,jc,kc-1,kps)
     &              *ci(ic+1,jc-1,kc-2,lxyl))
            enddo
         enddo
      enddo

      do kc=4,kkcp1         
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc,kc-2,kb)=yo(i,jc,kc-1,kb)
     &              +ci(ic,jc,kc-2,lxyr)
     &              *(yo(i,jc,kc-1,kbe)+yo(i-1,jc,kc-1,kb)
     &              *ci(ic,jc,kc-3,lxyr))+ci(ic+1,jc,kc-2,lxyl)
     &              *(yo(i+1,jc,kc-1,kbw)+yo(i+1,jc,kc-1,kb)
     &              *ci(ic+1,jc,kc-3,lxyl))+yo(i,jc,kc-1,kbw)
     &              *ci(ic,jc,kc-3,lxyr)+yo(i+1,jc,kc-1,kbe)
     &              *ci(ic+1,jc,kc-3,lxyl)
            enddo
         enddo
      enddo

      do kc=4,kkcp1         
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc,kc-2,kbw)=yo(i,jc,kc-1,kbw)
     &              *ci(ic,jc,kc-3,lxyl)
     &              +ci(ic,jc,kc-2,lxyr)*(yo(i-1,jc,kc-1,kbw)
     &              +yo(i-1,jc,kc-1,kb)*ci(ic,jc,kc-3,lxyl))
            enddo
         enddo
      enddo

      do kc=4,kkcp1         
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc+1,kc-2,kbnw)=yo(i,jc+1,kc-1,kbnw)
     &              *ci(ic,jc+1,kc-3,lxyl)+ci(ic,jc,kc-2,lxyr)
     &              *(yo(i-1,jc+1,kc-1,kbnw)+yo(i-1,jc+1,kc-1,kbn)
     &              *ci(ic,jc+1,kc-3,lxyl))
            enddo
         enddo
      enddo

      do kc=4,kkcp1         
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc+1,kc-2,kbn)=yo(i,jc+1,kc-1,kbn)
     &              +yo(i,jc+1,kc-1,kbnw)*ci(ic,jc+1,kc-3,lxyr)
     &              +yo(i+1,jc+1,kc-1,kbne)*ci(ic+1,jc+1,kc-3,lxyl)
     &              +ci(ic,jc,kc-2,lxyr)*(yo(i,jc+1,kc-1,kbne)
     &              +yo(i-1,jc+1,kc-1,kbn)*ci(ic,jc+1,kc-3,lxyr))
     &              +ci(ic+1,jc,kc-2,lxyl)*(yo(i+1,jc+1,kc-1,kbnw)
     &              +yo(i+1,jc+1,kc-1,kbn)*ci(ic+1,jc+1,kc-3,lxyl))
            enddo
         enddo
      enddo

      do kc=4,kkcp1         
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc+1,kc-2,kbne)=yo(i-1,jc+1,kc-1,kbne)
     &              *ci(ic,jc+1,kc-3,lxyr)+ci(ic,jc,kc-2,lxyl)
     &              *(yo(i,jc+1,kc-1,kbne)+yo(i-1,jc+1,kc-1,kbn)
     &              *ci(ic,jc+1,kc-3,lxyr))
            enddo
         enddo
      enddo

      do kc=4,kkcp1         
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc,kc-2,kbe)=yo(i-1,jc,kc-1,kbe)
     &              *ci(ic,jc,kc-3,lxyr)+ci(ic,jc,kc-2,lxyl)
     &              *(yo(i,jc,kc-1,kbe)+yo(i-1,jc,kc-1,kb)
     &              *ci(ic,jc,kc-3,lxyr))
            enddo
         enddo
      enddo

      do kc=4,kkcp1         
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc,kc-2,kbse)=yo(i-1,jc,kc-1,kbse)
     &              *ci(ic,jc-1,kc-3,lxyr)+ci(ic,jc,kc-2,lxyl)
     &              *(yo(i,jc,kc-1,kbse)+yo(i-1,jc,kc-1,kbs)
     &              *ci(ic,jc-1,kc-3,lxyr))
            enddo
         enddo
      enddo

      do kc=4,kkcp1         
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc,kc-2,kbs)=yo(i,jc,kc-1,kbs)+yo(i,jc,kc-1,kbsw)
     &              *ci(ic,jc-1,kc-3,lxyr)+yo(i+1,jc,kc-1,kbse)
     &              *ci(ic+1,jc-1,kc-3,lxyl)+ci(ic,jc,kc-2,lxyr)
     &              *(yo(i,jc,kc-1,kbse)+yo(i-1,jc,kc-1,kbs)
     &              *ci(ic,jc-1,kc-3,lxyr))+ci(ic+1,jc,kc-2,lxyl)
     &              *(yo(i+1,jc,kc-1,kbsw)+yo(i+1,jc,kc-1,kbs)
     &              *ci(ic+1,jc-1,kc-3,lxyl))
            enddo
         enddo
      enddo

      do kc=4,kkcp1         
         do jc=2,jjc1
            i=0
            do ic=2,iic1
               i=i+2
               soc(ic,jc,kc-2,kbsw)=yo(i,jc,kc-1,kbsw)
     &              *ci(ic,jc-1,kc-3,lxyl)+ci(ic,jc,kc-2,lxyr)
     &              *(yo(i-1,jc,kc-1,kbsw)+yo(i-1,jc,kc-1,kbs)
     &              *ci(ic,jc-1,kc-3,lxyl))
            enddo
         enddo
      enddo
      
C ==========================================================================

      return
      end






