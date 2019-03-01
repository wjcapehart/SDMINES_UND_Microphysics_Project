      SUBROUTINE BMG3_SER_SymStd_interp_add(
     &                KCG, KFG, Q ,QC, SO, RES, CI,
     &                IIC, JJC, KKC, IIF, JJF, KKF, NStncl
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   BMG3_SER_SymStd_interp_add.f interpolates Q from the coarse mesh, KCG, to 
C   the fine mesh, KFG, and adds the result to Q on fine mesh.
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
      INCLUDE 'BMG_SER_stencils.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  IIC, IIF, JJC, JJF, KCG, KFG, KKC, KKF, NStncl
      REAL*RKIND   CI(IIC,JJC,KKC,26), Q(IIF,JJF,KKF), QC(IIC,JJC,KKC), 
     &         RES(IIF,JJF,KKF), SO(IIF,JJF,KKF,NStncl)

C ----------------------------
C     Local Declarations
C
      INTEGER  IC, I, IICF, IICF1, IIC1,
     &         JC, J, JJCF, JJCF1, JJC1,
     &         KC, K, KKCF, KKCF1, KKC1
      REAL*RKIND   A, AQ

C ==========================================================================

C -------------------------------------------------
C     Useful index bounds:
C -------------------------------------------------

      iic1=iic-1
      jjc1=jjc-1
      kkc1=kkc-1

      iicf=(iif-2)/2+3
      jjcf=(jjf-2)/2+3
      kkcf=(kkf-2)/2+3

      iicf1=iicf-1
      jjcf1=jjcf-1
      kkcf1=kkcf-1

C --------------------------------------------------
C     NB: division is only possible in the interior
C --------------------------------------------------

      do k=2, kkf-1
         do j=2, jjf-1
            do i=2, iif-1
               RES(i,j,k)=RES(i,j,k)/so(i,j,k,kp)
            enddo
         enddo
      enddo
      
C --------------------------------------------------
C   interpolate answers from coarse to fine mesh 
C   and add to answers on fine mesh.
C --------------------------------------------------

      k=0
      do kc=2,kkc1
         k=k+2


         j=2
         i=2
         q(2,2,k)=q(2,2,k)+qc(2,2,kc)

         do ic=3,iicf1
            i=i+2
            q(i,2,k) = q(i,2,k) + qc(ic,2,kc)
            a = ci(ic,2,kc,lxyr)*qc(ic,2,kc)
     &        + ci(ic,2,kc,lxyl)*qc(ic-1,2,kc)
            q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)
         enddo

         do jc=3,jjcf1
            j=j+2
            i=2

            q(2,j,k) = q(2,j,k) + qc(2,jc,kc)
            aq = ci(2,jc,kc,lxya)*qc(2,jc,kc)
     &         + ci(2,jc,kc,lxyb)*qc(2,jc-1,kc)
            q(2,j-1,k) = q(2,j-1,k) + aq + RES(2,j-1,k)

            do ic=3,iicf1
               i=i+2

               q(i,j,k) = q(i,j,k) + qc(ic,jc,kc)
               a = ci(ic,jc,kc,lxyr)*qc(ic,jc,kc)
     &           + ci(ic,jc,kc,lxyl)*qc(ic-1,jc,kc)
               q(i-1,j,k) = q(i-1,j,k) + a + RES(i-1,j,k)

               aq = ci(ic,jc,kc,lxya)*qc(ic,jc,kc)
     &            + ci(ic,jc,kc,lxyb)*qc(ic,jc-1,kc)
               q(i,j-1,k) = q(i,j-1,k) + aq + RES(i,j-1,k)

               a = ci(ic,jc,kc,lxysw)*qc(ic-1,jc-1,kc)
     &           + ci(ic,jc,kc,lxynw)*qc(ic-1,jc,kc)
     &           + ci(ic,jc,kc,lxyne)*qc(ic,jc,kc)
     &           + ci(ic,jc,kc,lxyse)*qc(ic,jc-1,kc)
               q(i-1,j-1,k) = q(i-1,j-1,k) + a + RES(i-1,j-1,k)

            enddo

         enddo

      enddo



      
      
      k=1
      do kc=3,kkcf1
         k=k+2
         
         j=2
         jc=2
         i=0
         do ic=2,iic1
            i=i+2

            q(i,j,k) = q(i,j,k) + ci(ic,jc,kc,lxza)*qc(ic,jc,kc)
     &           + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1) + RES(i,j,k)
         enddo

         j=2
         do jc=3,jjcf1
            j=j+2
            
            i=0
            do ic=2,iic1
               i=i+2
               q(i,j,k) = q(i,j,k) 
     &                  + ci(ic,jc,kc,lxza)*qc(ic,jc,kc)
     &                  + ci(ic,jc,kc,lxzb)*qc(ic,jc,kc-1)
     &                  + RES(i,j,k)

               q(i,j-1,k) = q(i,j-1,k)
     &                    + ci(ic,jc,kc,lyznw)*qc(ic,jc,kc)
     &                    + ci(ic,jc,kc,lyzne)*qc(ic,jc-1,kc)
     &                    + ci(ic,jc,kc,lyzsw)*qc(ic,jc,kc-1)
     &                    + ci(ic,jc,kc,lyzse)*qc(ic,jc-1,kc-1)
     &                    + RES(i,j-1,k)
            enddo
         enddo
      enddo

      k=1
      do kc=3,kkcf1
         k=k+2
         
         j=2
         jc=2
         i=1
         do ic=3,iicf1
            i=i+2

            q(i,j,k) = q(i,j,k)
     &               + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc)
     &               + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc)
     &               + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1)
     &               + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1)
     &               + RES(i,j,k)

         enddo

         j=2
         do jc=3,jjcf1
            j=j+2

            i=1
            do ic=3,iicf1
               i=i+2

               q(i,j,k) = q(i,j,k) 
     &                  + ci(ic,jc,kc,lxznw)*qc(ic-1,jc,kc)
     &                  + ci(ic,jc,kc,lxzne)*qc(ic,jc,kc)
     &                  + ci(ic,jc,kc,lxzsw)*qc(ic-1,jc,kc-1)
     &                  + ci(ic,jc,kc,lxzse)*qc(ic,jc,kc-1)
     &                  + RES(i,j,k)

               q(i,j-1,k) = q(i,j-1,k) 
     &                    + ci(ic,jc,kc,ltnw)*qc(ic-1,jc,kc)
     &                    + ci(ic,jc,kc,ltne)*qc(ic,jc,kc)
     &                    + ci(ic,jc,kc,ltsw)*qc(ic-1,jc-1,kc)
     &                    + ci(ic,jc,kc,ltse)*qc(ic,jc-1,kc)
     &                    + ci(ic,jc,kc,lbnw)*qc(ic-1,jc,kc-1)
     &                    + ci(ic,jc,kc,lbne)*qc(ic,jc,kc-1)
     &                    + ci(ic,jc,kc,lbsw)*qc(ic-1,jc-1,kc-1)
     &                    + ci(ic,jc,kc,lbse)*qc(ic,jc-1,kc-1)
     &                    + RES(i,j-1,k)
            enddo

         enddo

      enddo

C ==========================================================================

      return
      end
