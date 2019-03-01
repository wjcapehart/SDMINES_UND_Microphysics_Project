      SUBROUTINE BMG3_SymStd_SETUP_ITLI_Izyx( 
     &                KGF, KGC, 
     &                SO, SOC, SOR1, SOR2, CI,
     &                IIF, JJF, KKF, IIC, JJC, KKC, 
     &                NOG, IFD, NStncl, IRELAX, 
     &                zo, yo, NOGm,
     &                IGRD, iWork, NMSGi, pSI_MSG, pMSG, pMSGSO,
     &                BUFFER, NMSGr, NProc,
     &                ProcGrid, MyProcI, MyProcJ, MyProcK, 
     &                NProcI, NProcJ, NProcK, BMG_IOFLAG,
     &                DimX, DimY, DimZ, 
     &                MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_SETUP_ITLI_Izyx computes an approximation to
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
C     NOG       Refer to BMG3_SymStd_SOLVE_boxmg.
C     IFD       Refer to BMG3_SymStd_SOLVE_boxmg.
C     IRELAX    Refer to BMG3_SymStd_SOLVE_boxmg.
C     SO        Refer to BMG3_SymStd_SOLVE_boxmg.
C     SOR       Refer to BMG3_SymStd_SOLVE_boxmg.
C     CI        Refer to BMG3_SymStd_SOLVE_boxmg.
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
      INCLUDE 'mpif.h'
      INCLUDE 'MSG.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER IFD, iic, iif, irelax, jjc, jjf, kgc, KGF, kkc, kkf,
     &        NOG, NOGm, NProc, NProcI, NProcJ, NProcK, NStncl, 
     &        NMSGi, NMSGr

      INTEGER iWork(NMSGi), pMSG(NBMG_pMSG,NOGm), pMSGSO(NBMG_pMSG,NOGm)

      INTEGER IGRD(NOGm,29), ProcGrid(NProcI, NProcJ, NProcK), 
     &        MyProcI, MyProcJ, MyProcK, MPICOMM

      REAL*RKIND  CI(IIC,JJC,KKC,26),
     &        SO(IIF+1,JJF+1,KKF+1,NStncl), 
     &        SOC(IIC+1,JJC+1,KKC+1,14), 
     &        SOR1(IIF,JJF,KKF), SOR2(IIF,JJF,KKC)
      REAL*RKIND  BUFFER(NMSGr)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)
      
      INTEGER DimX(NprocI,NOGm), DimY(NProcJ,NOGm), DimZ(NProcK,NOGm)

C ----------------------------
C     Local Declarations
C
      INTEGER ic, jc, kc, i, j, k, IIF1, JJF1, KKF1, IIC1, JJC1, KKC1, 
     &        KPZ, pSI_MSG, pMSGyo(NBMG_pMSG,2), pMSGzo(NBMG_pMSG,2)

      INTEGER MPI_IERROR, MyProc

      REAL*RKIND  a, b, c, eMACH, ep 
      LOGICAL LEFTPLANE, BOTTOMPLANE, FRONTPLANE,
     &        RIGHTPLANE, TOPPLANE, BACKPLANE
      INTEGER ISTART, JSTART, ICEND, JCEND, ICSTART, JCSTART,
     &        IEND, JEND, KEND,  KSTART, KCSTART, KCEND, 
     &        ISTARTO, JSTARTO, ICENDO, JCENDO, ICSTARTO, JCSTARTO,
     &        KSTARTO, KCSTARTO, KCENDO, 
     &        LXGP, RXGP, LYGP, RYGP, LZGP, RZGP

      REAL*RKIND  YO(IIF,JJC,KKC,14), zo(IIF,JJF,KKC,14)

      INTEGER ierror, ptrn

C ==========================================================================

C ==========================================================================
C      Setup communicators for zo and yo by 
C      calling PtrMSG once for current grid.
C ==========================================================================

      CALL BMG3_SymStd_SETUP_PtrMSG(
     &          IGRD(KGF,idL_BMG_NLx), IGRD(KGF,idL_BMG_NLy), 
     &          IGRD(KGC,idL_BMG_NLz), iZERO, iZERO, iZERO, 
     &          pSI_MSG, NProc, 2, 1, 1, pMSGzo
     &          )

      CALL BMG3_SymStd_SETUP_PtrMSG(
     &          IGRD(KGF,idL_BMG_NLx), IGRD(KGC,idL_BMG_NLy), 
     &          IGRD(KGC,idL_BMG_NLz), iZERO, iZERO, iZERO, 
     &          pSI_MSG, NProc, 2, 1, 1, pMSGyo
     &          )

      eMACH = 1.d-13

      MyProc = ProcGrid(MyProcI,MyProcJ,MyProcK)

      IIC1 = IIC-1
      JJC1 = JJC-1
      KKC1 = KKC-1

      IIF1 = IIF-1
      JJF1 = JJF-1
      KKF1 = KKF-1

      DO KPZ=1,14
         DO k=1,KKC
            DO j=1,JJC
               DO i=1,IIF
                  yo(i,j,k,KPZ)= rZERO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO KPZ=1,14
         DO k=1,KKC
            DO j=1,JJF
               DO i=1,IIF
                  zo(i,j,k,KPZ)= rZERO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO KPZ=1,14
         DO k=1,KKC+1
            DO j=1,JJC+1
               DO i=1,IIC+1
                  soc(i,j,k,KPZ)= rZERO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C     find the location of our node in the processor grid
C

      LEFTPLANE   = ( IGRD(KGF,idL_BMG_ICoord).EQ.1 )
      BOTTOMPLANE = ( IGRD(KGF,idL_BMG_JCoord).EQ.1 )
      FRONTPLANE  = ( IGRD(KGF,idL_BMG_KCoord).EQ.1 )

      !
      ! If ****PLANE = 0, then it is not on ****PLANE
      ! IF ****PLANE = 1, then it is on ****PLANE.
      !

      RIGHTPLANE  = (IGRD(KGF,idL_BMG_ICoord)+IGRD(KGF,idL_BMG_NLx)-1)
     &     .EQ.IGRD(KGF,idL_BMG_NGx)
      TOPPLANE    = (IGRD(KGF,idL_BMG_JCoord)+IGRD(KGF,idL_BMG_NLy)-1)
     &     .EQ.IGRD(KGF,idL_BMG_NGy)
      BACKPLANE   = (IGRD(KGF,idL_BMG_KCoord)+IGRD(KGF,idL_BMG_NLz)-1)
     &     .EQ.IGRD(KGF,idL_BMG_NGz)

c     ------------------------------------

      IF (LEFTPLANE) THEN
         ISTARTO  = 2
         ICSTARTO = 3
      ELSE
         IF (mod(IGRD(KGF,idL_BMG_ICoord),2).EQ.1) THEN
            ISTARTO = 0
         ELSE
            ISTARTO = 1
         ENDIF
         ICSTARTO = 2
      ENDIF
      IF ( (IGRD(KGC,idL_BMG_ICoord)+IIC-3)*2 .EQ.   ! expresses whether or not
     &     (IGRD(KGF,idL_BMG_ICoord)+IIF-2)  ) THEN  ! RT-most interior FG point 
         ICENDO = IIC1                               ! is same as RT-most interior
      ELSE                                           ! CG point
         ICENDO = IIC1+1
      ENDIF
      IF (mod(IGRD(KGF,idL_BMG_ICoord),2).EQ.1) THEN
         ISTART = 0
      ELSE
         ISTART = 1
      ENDIF
      ICSTART = 2
      ICEND = IIC1
      IEND  = IIF1
      
      
c     -----------------------------------
            
      IF (BOTTOMPLANE) THEN
         JSTARTO  = 2
         JCSTARTO = 3
      ELSE
         IF (mod(IGRD(KGF,idL_BMG_JCoord),2).EQ.1) THEN
            JSTARTO = 0
         ELSE
            JSTARTO = 1
         ENDIF
         JCSTARTO = 2
      ENDIF
      IF ( (IGRD(KGC,idL_BMG_JCoord)+JJC-3)*2 .EQ.   ! expresses whether or not
     &     (IGRD(KGF,idL_BMG_JCoord)+JJF-2)  ) THEN  ! RT-most interior FG point
         JCENDO = JJC1                               ! is same as RT-most interior 
      ELSE                                           ! CG point
         JCENDO = JJC1+1
      ENDIF
      IF (mod(IGRD(KGF,idL_BMG_JCoord),2).EQ.1) THEN
         JSTART = 0
      ELSE
         JSTART = 1
      ENDIF
      JCSTART = 2
      JCEND = JJC1
      JEND  = JJF1
      
c     -----------------------------------        
      
      IF (FRONTPLANE) THEN
         KSTARTO  = 2
         KCSTARTO = 3
      ELSE
         IF (mod(IGRD(KGF,idL_BMG_KCoord),2).EQ.1) THEN
            KSTARTO = 0
         ELSE
            KSTARTO = 1
         ENDIF
         KCSTARTO = 2
      ENDIF
      IF ( (IGRD(KGC,idL_BMG_KCoord)+KKC-3)*2 .EQ.   ! expresses whether or not
     &     (IGRD(KGF,idL_BMG_KCoord)+KKF-2)  ) THEN  ! RT-most interior FG point
         KCENDO = KKC1                               ! is same as RT-most interior 
      ELSE                                           ! CG point
         KCENDO = KKC1+1
      ENDIF
      IF (mod(IGRD(KGF,idL_BMG_KCoord),2).EQ.1) THEN
         KSTART = 0
      ELSE
         KSTART = 1
      ENDIF
      KCSTART = 2
      KCEND = KKC1
      KEND  = KKF1

c     -----------------------------------


c
c   compute interpolation operator iz, from grid kgfz to grid kgf, where
c   grid kgfz is obtained form grid kgf by coarsening in the z-direction
c   only.
c
      IF( KGF.LT.NOG .OR. IFD.NE.1 ) THEN 

         k=KSTARTO
         DO kc=KCSTARTO,KCENDO
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
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
                  SOR1(i,j,k-1)=a/c
                  SOR1(i,j,k)=b/c
               ENDDO
            ENDDO
         ENDDO

      ELSE ! if KGF.ge.NOG.and.IFD.eq.1   

         k=KSTARTO
         DO kc=KCSTARTO,KCENDO
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  a=so(i,j,k-1,kb)
                  b=so(i,j,k,kb)
                  c=a+b+so(i,j+1,k-1,kps)+so(i,j,k-1,kpw)
     &                 +so(i+1,j,k-1,kpw)+so(i,j,k-1,kps)
                  ep=MIN(abs(a),abs(b),rONE)
                  c=a+b+(so(i,j,k-1,kp)-c)*MAX(so(i,j,k-1,kp)
     &                 -(rONE+ep)*c,rZERO)/(abs(so(i,j,k-1,kp)
     &                 -(rONE+ep)*c)+eMACH)
                  SOR1(i,j,k-1)=a/c
                  SOR1(i,j,k)=b/c
               ENDDO
            ENDDO
         ENDDO

      ENDIF ! if KGF.ge.NOG.and.IFD.eq.1 
  
c
c     update the ghost points for sor here 
c
c     -------------------------------


      ptrn = 1
      CALL MSG_tbdx_send(SOR1(1,1,1), BUFFER, 
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KGF)),
     &     iWork(pMSG(ipL_MSG_Proc,KGF)),
     &     iWork(pMSG(ipL_MSG_Ipr,KGF)),
     &     iWork(pMSG(ipL_MSG_Index,KGF)),
     &     ptrn, ierror)
            
      CALL MSG_tbdx_receive(SOR1(1,1,1), BUFFER,
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KGF)),
     &     iWork(pMSG(ipL_MSG_Proc,KGF)),
     &     iWork(pMSG(ipL_MSG_Ipr,KGF)),
     &     iWork(pMSG(ipL_MSG_Index,KGF)),
     &     ptrn, ierror)
            
      CALL MSG_tbdx_close(SOR1(1,1,1), BUFFER,
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KGF)),
     &     iWork(pMSG(ipL_MSG_Proc,KGF)),
     &     iWork(pMSG(ipL_MSG_Ipr,KGF)),
     &     iWork(pMSG(ipL_MSG_Index,KGF)),
     &     ptrn, ierror)
 
c     -------------------------------


c
c   compute lz = iz(transpose) l iz, where l is the difference operator
c   on grid kgf.
c

      IF( KGF.LT.NOG .OR. IFD.NE.1 ) THEN 
         
         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kpw)=
     &                 so(i,j,k,kpw)
     &                 +so(i,j,k+1,kbe)*SOR1(i-1,j,k+1)
     &                 +so(i,j,k,kbw)*SOR1(i-1,j,k)
     &                 +SOR1(i,j,k+1)
     &                 *( so(i,j,k+1,kpw)*SOR1(i-1,j,k+1)
     &                 +so(i,j,k+1,kbw))
     &                 +SOR1(i,j,k) 
     &                 *( so(i,j,k-1,kpw)*SOR1(i-1,j,k)
     &                 +so(i,j,k,kbe))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kpsw)=so(i,j,k,kpsw)+so(i,j,k+1,kbne)
     &                 *SOR1(i-1,j-1,k+1)+so(i,j,k,kbsw)
     &                 *SOR1(i-1,j-1,k)+SOR1(i,j,k)
     &                 *(so(i,j,k-1,kpsw)*SOR1(i-1,j-1,k)
     &                 +so(i,j,k,kbne))+SOR1(i,j,k+1)
     &                 *(so(i,j,k+1,kpsw)*SOR1(i-1,j-1,k+1)
     &                 +so(i,j,k+1,kbsw))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j+1,kc,kpnw)=so(i,j+1,k,kpnw)+so(i,j+1,k,kbnw)
     &                 *SOR1(i-1,j+1,k)+so(i,j+1,k+1,kbse)
     &                 *SOR1(i-1,j+1,k+1)+SOR1(i,j,k)
     &                 *(so(i,j+1,k-1,kpnw)*SOR1(i-1,j+1,k)
     &                 +so(i,j+1,k,kbse))+SOR1(i,j,k+1)
     &                 *(so(i,j+1,k+1,kpnw)*SOR1(i-1,j+1,k+1)
     &                 +so(i,j+1,k+1,kbnw))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kps)=so(i,j,k,kps)+so(i,j,k+1,kbn)
     &                 *SOR1(i,j-1,k+1)+so(i,j,k,kbs)
     &                 *SOR1(i,j-1,k)+SOR1(i,j,k)
     &                 *(so(i,j,k-1,kps)*SOR1(i,j-1,k)
     &                 +so(i,j,k,kbn))+SOR1(i,j,k+1)
     &                 *(so(i,j,k+1,kps)*SOR1(i,j-1,k+1)
     &                 +so(i,j,k+1,kbs))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kp)=so(i,j,k,kp)-so(i,j,k,kb)
     &                 *SOR1(i,j,k)-so(i,j,k+1,kb)
     &                 *SOR1(i,j,k+1)-SOR1(i,j,k)
     &                 *(-so(i,j,k-1,kp)*SOR1(i,j,k)
     &                 +so(i,j,k,kb))-SOR1(i,j,k+1)
     &                 *(-so(i,j,k+1,kp)*SOR1(i,j,k+1)
     &                 +so(i,j,k+1,kb))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kb)=so(i,j,k,kb)*SOR1(i,j,k-1)
     &                 +SOR1(i,j,k)*(so(i,j,k-1,kb)
     &                 -so(i,j,k-1,kp)*SOR1(i,j,k-1))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kbw)=so(i,j,k,kbw)*SOR1(i-1,j,k-1)
     &                 +SOR1(i,j,k)*(so(i,j,k-1,kpw)
     &                 *SOR1(i-1,j,k-1)+so(i,j,k-1,kbw))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j+1,kc,kbnw)=so(i,j+1,k,kbnw)
     &                 *SOR1(i-1,j+1,k-1)+SOR1(i,j,k)
     &                 *(so(i,j+1,k-1,kpnw)*SOR1(i-1,j+1,k-1)
     &                 +so(i,j+1,k-1,kbnw))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j+1,kc,kbn)=so(i,j+1,k,kbn)
     &                 *SOR1(i,j+1,k-1)
     &                 +SOR1(i,j,k)*(so(i,j+1,k-1,kps)
     &                 *SOR1(i,j+1,k-1)+so(i,j+1,k-1,kbn))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i+1,j+1,kc,kbne)=so(i+1,j+1,k,kbne)
     &                 *SOR1(i+1,j+1,k-1)+SOR1(i,j,k)
     &                 *(so(i+1,j+1,k-1,kpsw)*SOR1(i+1,j+1,k-1)
     &                 +so(i+1,j+1,k-1,kbne))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i+1,j,kc,kbe)=so(i+1,j,k,kbe)
     &                 *SOR1(i+1,j,k-1)
     &                 +SOR1(i,j,k)*(so(i+1,j,k-1,kpw)
     &                 *SOR1(i+1,j,k-1)+so(i+1,j,k-1,kbe))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i+1,j,kc,kbse)=so(i+1,j,k,kbse)
     &                 *SOR1(i+1,j-1,k-1)+SOR1(i,j,k)
     &                 *(so(i+1,j,k-1,kpnw)*SOR1(i+1,j-1,k-1)
     &                 +so(i+1,j,k-1,kbse))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kbs)=so(i,j,k,kbs)*SOR1(i,j-1,k-1)
     &                 +SOR1(i,j,k)*(so(i,j,k-1,kps)
     &                 *SOR1(i,j-1,k-1)+so(i,j,k-1,kbs))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kbsw)=so(i,j,k,kbsw)
     &                 *SOR1(i-1,j-1,k-1)
     &                 +SOR1(i,j,k)*(so(i,j,k-1,kpsw)
     &                 *SOR1(i-1,j-1,k-1)+so(i,j,k-1,kbsw))
               ENDDO
            ENDDO
         ENDDO

      ELSE ! ===============================================

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kpw)=so(i,j,k,kpw)+SOR1(i,j,k+1)
     &                 *so(i,j,k+1,kpw)*SOR1(i-1,j,k+1)
     &                 +SOR1(i,j,k)*so(i,j,k-1,kpw)
     &                 *SOR1(i-1,j,k)
                  zo(i,j,kc,kpsw)= rZERO
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j+1,kc,kpnw)= rZERO
                  zo(i,j,kc,kps)=so(i,j,k,kps)+SOR1(i,j,k)
     &                 *so(i,j,k-1,kps)*SOR1(i,j-1,k)
     &                 +SOR1(i,j,k+1)*so(i,j,k+1,kps)
     &                 *SOR1(i,j-1,k+1)
                  zo(i,j+1,kc,kpnw)= rZERO
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kp)=so(i,j,k,kp)-SOR1(i,j,k)
     &                 *(-so(i,j,k-1,kp)*SOR1(i,j,k)
     &                 +so(i,j,k,kb))-SOR1(i,j,k+1)
     &                 *(-so(i,j,k+1,kp)*SOR1(i,j,k+1)
     &                 +so(i,j,k+1,kb))-so(i,j,k,kb)*SOR1(i,j,k)
     &                 -so(i,j,k+1,kb)*SOR1(i,j,k+1)
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kb)=so(i,j,k,kb)*SOR1(i,j,k-1)
     &                 +SOR1(i,j,k)*(so(i,j,k-1,kb)
     &                 -so(i,j,k-1,kp)*SOR1(i,j,k-1))
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kbw)=SOR1(i,j,k)*so(i,j,k-1,kpw)
     &                 *SOR1(i-1,j,k-1)
                  zo(i,j+1,kc,kbnw)= rZERO
               ENDDO
            ENDDO
         ENDDO
         
         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j+1,kc,kbn)=SOR1(i,j,k)*so(i,j+1,k-1,kps)
     &                 *SOR1(i,j+1,k-1)
                  zo(i+1,j+1,kc,kbne)= rZERO
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i+1,j,kc,kbe)=SOR1(i,j,k)*so(i+1,j,k-1,kpw)
     &                 *SOR1(i+1,j,k-1)
                  zo(i+1,j,kc,kbse)= rZERO
               ENDDO
            ENDDO
         ENDDO

         k=KSTART
         DO kc=KCSTART,KCEND
            k=k+2
            DO j=2,JJF1
               DO i=2,IIF1
                  zo(i,j,kc,kbs)=SOR1(i,j,k)*so(i,j,k-1,kps)
     &                 *SOR1(i,j-1,k-1)
                  zo(i,j,kc,kbsw)= rZERO
               ENDDO
            ENDDO
         ENDDO

      ENDIF ! of if (KGF.lt.NOG.or.IFD.ne.1)

c
c     Update zo at the ghost boundaries
c

      CALL BMG3_SymStd_update_zo(
     &          zo, IIF, JJF, KKC,
     &          pMSGzo, iWork, NMSGi, IGRD, BUFFER, NMSGr,
     &          KGF, KGC, NOGm, NProc, MyProc, 
     &          DimX, DimY, DimZ, ProcGrid, 
     &          NProcI, NProcJ, NProcK,
     &          MPICOMM
     &          )

c     compute iy, the interpolation operator form grid kgzy to grid
c     kgz, where grid kgz is obtained from grid kgz by coarsening
c     in the y-direction only.
     

      DO kc=KCSTART,KCEND  
         j=JSTARTO         
         DO jc=JCSTARTO,JCENDO 
            j=j+2
            DO i=2,IIF1
               a=zo(i,j-1,kc,kbsw)+zo(i,j-1,kc,kbs)
     &              +zo(i+1,j-1,kc,kbse)+zo(i,j-1,kc,kpsw)
     &              +zo(i,j-1,kc,kps)+zo(i+1,j-1,kc,kpnw)
     &              +zo(i,j-1,kc+1,kbne)+zo(i,j-1,kc+1,kbn)
     &              +zo(i+1,j-1,kc+1,kbnw)
               b=zo(i,j,kc,kbnw)+zo(i,j,kc,kbn)
     &              +zo(i+1,j,kc,kbne)+zo(i,j,kc,kpnw)
     &              +zo(i,j,kc,kps)+zo(i+1,j,kc,kpsw)
     &              +zo(i,j,kc+1,kbse)+zo(i,j,kc+1,kbs)
     &              +zo(i+1,j,kc+1,kbsw)
               c=zo(i,j-1,kc,kbw)+zo(i,j-1,kc,kb)
     &              +zo(i+1,j-1,kc,kbe)+zo(i,j-1,kc,kpw)
     &              +a+b+zo(i+1,j-1,kc,kpw)
     &              +zo(i,j-1,kc+1,kbe)+zo(i,j-1,kc+1,kb)
     &              +zo(i+1,j-1,kc+1,kbw)
                  
               ep=MIN(abs(a),abs(b),rONE)
               c=a+b+(zo(i,j-1,kc,kp)-c)*MAX(zo(i,j-1,kc,kp)
     &              -(rONE+ep)*c,rZERO)
     &              /(abs(zo(i,j-1,kc,kp)-(rONE+ep)*c)+eMACH)
               SOR2(i,j-1,kc)=a/c
               SOR2(i,j,kc)=b/c
            ENDDO
         ENDDO
      ENDDO

c
c     update the ghost points for sor here by using pMSG just like
c     when you update the points for the residual calculation
c

      CALL BMG3_SymStd_update_sor2(
     &          SOR2, IIF, JJF, KKC,
     &          pMSGzo, iWork, NMSGi, IGRD, BUFFER, NMSGr,
     &          KGF, KGC, NOGm, NProc, MyProc, 
     &          DimX, DimY, DimZ, ProcGrid, 
     &          NProcI, NProcJ, NProcK,
     &          MPICOMM
     &          )

c     
c     compute ly = iy(transpose) lz iy.
c     

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i,jc,kc,kp)=zo(i,j,kc,kp)-zo(i,j,kc,kps)
     &              *SOR2(i,j,kc)-zo(i,j+1,kc,kps)
     &              *SOR2(i,j+1,kc)-SOR2(i,j,kc)
     &              *(-zo(i,j-1,kc,kp)*SOR2(i,j,kc)
     &              +zo(i,j,kc,kps))-SOR2(i,j+1,kc)
     &              *(-zo(i,j+1,kc,kp)*SOR2(i,j+1,kc)
     &              +zo(i,j+1,kc,kps))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
             DO i=2,IIF1
               yo(i,jc,kc,kpw)=zo(i,j,kc,kpw)+zo(i,j,kc,kpsw)
     &              *SOR2(i-1,j,kc)+zo(i,j+1,kc,kpnw)
     &              *SOR2(i-1,j+1,kc)+SOR2(i,j,kc)
     &              *(zo(i,j-1,kc,kpw)*SOR2(i-1,j,kc)
     &              +zo(i,j,kc,kpnw))+SOR2(i,j+1,kc)
     &              *(zo(i,j+1,kc,kpw)*SOR2(i-1,j+1,kc)
     &              +zo(i,j+1,kc,kpsw))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i,jc,kc,kpnw) = zo(i,j-1,kc,kpnw)*SOR2(i-1,j,kc) +
     &            SOR2(i,j-1,kc)*zo(i,j-1,kc,kpw )*SOR2(i-1,j,kc)+
     &            SOR2(i,j-1,kc)*zo(i,j,kc,kpnw)
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i,jc,kc,kps)=zo(i,j,kc,kps)*SOR2(i,j-1,kc)
     &              +SOR2(i,j,kc)*(zo(i,j-1,kc,kps)
     &              -zo(i,j-1,kc,kp)*SOR2(i,j-1,kc))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i,jc,kc,kpsw)=
     &                SOR2(i,j,kc)*zo(i,j-1,kc,kpsw)
     &           + SOR2(i,j,kc)*zo(i,j-1,kc,kpw)*SOR2(i-1,j-1,kc)
     &           + zo(i,j,kc,kpsw)*SOR2(i-1,j-1,kc)
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i,jc,kc,kb)=zo(i,j,kc,kb)
     &             + SOR2(i, j ,kc)*zo(i, j ,kc,kbn)
     &             + SOR2(i,j+1,kc)*zo(i,j+1,kc,kbs) 
     &             + SOR2(i, j ,kc)*zo(i,j-1,kc,kb)*SOR2(i, j ,kc-1)
     &             + SOR2(i,j+1,kc)*zo(i,j+1,kc,kb)*SOR2(i,j+1,kc-1) 
     &             + zo(i, j ,kc,kbs)*SOR2(i, j ,kc-1)
     &             + zo(i,j+1,kc,kbn)*SOR2(i,j+1,kc-1)
    
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i,jc,kc,kbw)=zo(i,j,kc,kbw)+zo(i,j,kc,kbsw)
     &              *SOR2(i-1,j,kc-1)+zo(i,j+1,kc,kbnw)
     &              *SOR2(i-1,j+1,kc-1)+SOR2(i,j,kc)
     &              *(zo(i,j,kc,kbnw)+zo(i,j-1,kc,kbw)
     &              *SOR2(i-1,j,kc-1))+SOR2(i,j+1,kc)
     &              *(zo(i,j+1,kc,kbsw)+zo(i,j+1,kc,kbw)
     &              *SOR2(i-1,j+1,kc-1))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i,jc,kc,kbnw)=zo(i,j-1,kc,kbnw)
     &              *SOR2(i-1,j,kc-1)+SOR2(i,j-1,kc)
     &              *(zo(i,j-1,kc,kbw)*SOR2(i-1,j,kc-1)
     &              +zo(i,j,kc,kbnw))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i,jc,kc,kbn)=
     &               zo(i,j-1,kc,kbn)*SOR2(i,j,kc-1)
     &              +SOR2(i,j-1,kc)*(zo(i,j,kc,kbn)
     &              +zo(i,j-1,kc,kb)*SOR2(i,j,kc-1))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i+1,jc,kc,kbne)=zo(i+1,j-1,kc,kbne)
     &              *SOR2(i+1,j,kc-1)+SOR2(i,j-1,kc)
     &              *(zo(i+1,j,kc,kbne)+zo(i+1,j-1,kc,kbe)
     &              *SOR2(i+1,j,kc-1))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i+1,jc,kc,kbe)=zo(i+1,j,kc,kbe)
     &              +zo(i+1,j,kc,kbse)
     &              *SOR2(i+1,j,kc-1)+zo(i+1,j+1,kc,kbne)
     &              *SOR2(i+1,j+1,kc-1)+SOR2(i,j,kc)
     &              *(zo(i+1,j,kc,kbne)+zo(i+1,j-1,kc,kbe)
     &              *SOR2(i+1,j,kc-1))+SOR2(i,j+1,kc)
     &              *(zo(i+1,j+1,kc,kbse)+zo(i+1,j+1,kc,kbe)
     &              *SOR2(i+1,j+1,kc-1))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i+1,jc,kc,kbse)=zo(i+1,j,kc,kbse)
     &              *SOR2(i+1,j-1,kc-1)+SOR2(i,j,kc)
     &              *(zo(i+1,j-1,kc,kbse)+zo(i+1,j-1,kc,kbe)
     &              *SOR2(i+1,j-1,kc-1))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i,jc,kc,kbs)=zo(i,j,kc,kbs)
     &              *SOR2(i,j-1,kc-1)
     &              +SOR2(i,j,kc)*(zo(i,j-1,kc,kbs)
     &              +zo(i,j-1,kc,kb)*SOR2(i,j-1,kc-1))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         j=JSTART
         DO jc=JCSTART,JCEND
            j=j+2
            DO i=2,IIF1
               yo(i,jc,kc,kbsw)=zo(i,j,kc,kbsw)
     &              *SOR2(i-1,j-1,kc-1)+SOR2(i,j,kc)
     &              *(zo(i,j-1,kc,kbsw)+zo(i,j-1,kc,kbw)
     &              *SOR2(i-1,j-1,kc-1))
            ENDDO
         ENDDO
      ENDDO
      
c
c     Update yo at the ghost boundaries
c

      CALL BMG3_SymStd_update_yo(
     &         yo, IIF, JJC, KKC,
     &         pMSGyo, iWork, NMSGi, IGRD, BUFFER, NMSGr,
     &         KGF, KGC, NOGm, NProc, MyProc, 
     &         DimX, DimY, DimZ, ProcGrid, 
     &         NProcI, NProcJ, NProcK,
     &         MPICOMM
     &         )

c     
c     compute ix, the interpolation operator form grid kgzy to grid kgc
c     

      DO kc=KCSTART,KCEND ! 4,KKCP1 
         DO jc=JCSTART,JCEND
            i=ISTARTO
            DO ic=ICSTARTO,ICENDO
               i=i+2
               a=yo(i-1,jc+1,kc,kbnw)+yo(i-1,jc,kc,kbw)
     &              +yo(i-1,jc,kc,kbsw)+yo(i-1,jc+1,kc,kpnw)
     &              +yo(i-1,jc,kc,kpw)+yo(i-1,jc,kc,kpsw)
     &              +yo(i-1,jc+1,kc+1,kbse)+yo(i-1,jc,kc+1,kbe)
     &              +yo(i-1,jc,kc+1,kbne)
               b=yo(i,jc+1,kc,kbne)+yo(i,jc,kc,kbe)
     &              +yo(i,jc,kc,kbse)
     &              +yo(i,jc+1,kc,kpsw)+yo(i,jc,kc,kpw)
     &              +yo(i,jc,kc,kpnw)+yo(i,jc+1,kc+1,kbsw)
     &              +yo(i,jc,kc+1,kbw)+yo(i,jc,kc+1,kbnw)
               c=yo(i-1,jc+1,kc,kbn)+yo(i-1,jc,kc,kb)
     &              +yo(i-1,jc,kc,kbs)+yo(i-1,jc+1,kc,kps)+a+b
     &              +yo(i-1,jc,kc,kps)+yo(i-1,jc+1,kc+1,kbs)
     &              +yo(i-1,jc,kc+1,kb)+yo(i-1,jc,kc+1,kbn)
               ep=MIN(abs(a),abs(b),rONE)
               c=a+b+(yo(i-1,jc,kc,kp)-c)*MAX(yo(i-1,jc,kc,kp)
     &              -(rONE+ep)*c,rZERO)
     &              /(abs(yo(i-1,jc,kc,kp)-(rONE+ep)*c)+eMACH)
               CI(ic,jc,kc,lxyl)=a/c
               CI(ic,jc,kc,lxyr)=b/c
            ENDDO
         ENDDO
      ENDDO   
  

c
c     Update CI at the ghost boundaries
c
      DO kpz=1,2

         ptrn = 1
         CALL MSG_tbdx_send(CI(1,1,1,kpz), BUFFER, 
     &        iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),
     &        iWork(pMSG(ipL_MSG_Proc,KGC)),
     &        iWork(pMSG(ipL_MSG_Ipr,KGC)),
     &        iWork(pMSG(ipL_MSG_Index,KGC)),
     &        ptrn, ierror)
            
         CALL MSG_tbdx_receive(CI(1,1,1,kpz), BUFFER,
     &        iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),
     &        iWork(pMSG(ipL_MSG_Proc,KGC)),
     &        iWork(pMSG(ipL_MSG_Ipr,KGC)),
     &        iWork(pMSG(ipL_MSG_Index,KGC)),
     &        ptrn, ierror)
            
         CALL MSG_tbdx_close(CI(1,1,1,kpz), BUFFER,
     &        iWork(pMSG(ipL_MSG_NumAdjProc,KGC)),
     &        iWork(pMSG(ipL_MSG_Proc,KGC)),
     &        iWork(pMSG(ipL_MSG_Ipr,KGC)),
     &        iWork(pMSG(ipL_MSG_Index,KGC)),
     &        ptrn, ierror)      
      enddo

c               
c     compute lc = ix (transpose) ly ix.
c     
      DO kc=KCSTART,KCEND
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc,kc,kp)=yo(i,jc,kc,kp)-yo(i,jc,kc,kpw)
     &              *CI(ic,jc,kc,lxyr)-yo(i+1,jc,kc,kpw)
     &              *CI(ic+1,jc,kc,lxyl)-CI(ic,jc,kc,lxyr)
     &              *(-yo(i-1,jc,kc,kp)*CI(ic,jc,kc,lxyr)
     &              +yo(i,jc,kc,kpw))-CI(ic+1,jc,kc,lxyl)
     &              *(-yo(i+1,jc,kc,kp)*CI(ic+1,jc,kc,lxyl)
     &              +yo(i+1,jc,kc,kpw))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc,kc,kpw)=yo(i,jc,kc,kpw)
     &              *CI(ic,jc,kc,lxyl)
     &              +CI(ic,jc,kc,lxyr)*(yo(i-1,jc,kc,kpw)
     &              -yo(i-1,jc,kc,kp)*CI(ic,jc,kc,lxyl))
            ENDDO
         ENDDO
      ENDDO

C
C     ** Needs different message scheme **
C
      DO kc=KCSTART,KCEND         
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc+1,kc,kpnw)=yo(i,jc+1,kc,kpnw)
     &              *CI(ic,jc+1,kc,lxyl)+CI(ic,jc,kc,lxyr)
     &              *(yo(i-1,jc+1,kc,kpnw)+yo(i-1,jc+1,kc,kps)
     &              *CI(ic,jc+1,kc,lxyl))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND         
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc,kc,kpsw)=yo(i,jc,kc,kpsw)
     &              *CI(ic,jc-1,kc,lxyl)+CI(ic,jc,kc,lxyr)
     &              *(yo(i-1,jc,kc,kpsw)+yo(i-1,jc,kc,kps)
     &              *CI(ic,jc-1,kc,lxyl))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND         
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc,kc,kps)=yo(i,jc,kc,kps)+yo(i,jc,kc,kpsw)
     &              *CI(ic,jc-1,kc,lxyr)+yo(i+1,jc,kc,kpnw)
     &              *CI(ic+1,jc-1,kc,lxyl)+CI(ic,jc,kc,lxyr)
     &              *(yo(i,jc,kc,kpnw)+yo(i-1,jc,kc,kps)
     &              *CI(ic,jc-1,kc,lxyr))+CI(ic+1,jc,kc,lxyl)
     &              *(yo(i+1,jc,kc,kpsw)+yo(i+1,jc,kc,kps)
     &              *CI(ic+1,jc-1,kc,lxyl))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND         
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc,kc,kb)=yo(i,jc,kc,kb)
     &              +CI(ic,jc,kc,lxyr)
     &              *(yo(i,jc,kc,kbe)+yo(i-1,jc,kc,kb)
     &              *CI(ic,jc,kc-1,lxyr))+CI(ic+1,jc,kc,lxyl)
     &              *(yo(i+1,jc,kc,kbw)+yo(i+1,jc,kc,kb)
     &              *CI(ic+1,jc,kc-1,lxyl))+yo(i,jc,kc,kbw)
     &              *CI(ic,jc,kc-1,lxyr)+yo(i+1,jc,kc,kbe)
     &              *CI(ic+1,jc,kc-1,lxyl)
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND         
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc,kc,kbw)=yo(i,jc,kc,kbw)
     &              *CI(ic,jc,kc-1,lxyl)
     &              +CI(ic,jc,kc,lxyr)*(yo(i-1,jc,kc,kbw)
     &              +yo(i-1,jc,kc,kb)*CI(ic,jc,kc-1,lxyl))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND         
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc+1,kc,kbnw)=yo(i,jc+1,kc,kbnw)
     &              *CI(ic,jc+1,kc-1,lxyl)+CI(ic,jc,kc,lxyr)
     &              *(yo(i-1,jc+1,kc,kbnw)+yo(i-1,jc+1,kc,kbn)
     &              *CI(ic,jc+1,kc-1,lxyl))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND         
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc+1,kc,kbn)=yo(i,jc+1,kc,kbn)
     &              +yo(i,jc+1,kc,kbnw)*CI(ic,jc+1,kc-1,lxyr)
     &              +yo(i+1,jc+1,kc,kbne)*CI(ic+1,jc+1,kc-1,lxyl)
     &              +CI(ic,jc,kc,lxyr)*(yo(i,jc+1,kc,kbne)
     &              +yo(i-1,jc+1,kc,kbn)*CI(ic,jc+1,kc-1,lxyr))
     &              +CI(ic+1,jc,kc,lxyl)*(yo(i+1,jc+1,kc,kbnw)
     &              +yo(i+1,jc+1,kc,kbn)*CI(ic+1,jc+1,kc-1,lxyl))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND         
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc+1,kc,kbne)=yo(i-1,jc+1,kc,kbne)
     &              *CI(ic,jc+1,kc-1,lxyr)+CI(ic,jc,kc,lxyl)
     &              *(yo(i,jc+1,kc,kbne)+yo(i-1,jc+1,kc,kbn)
     &              *CI(ic,jc+1,kc-1,lxyr))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND         
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc,kc,kbe)=yo(i-1,jc,kc,kbe)
     &              *CI(ic,jc,kc-1,lxyr)+CI(ic,jc,kc,lxyl)
     &              *(yo(i,jc,kc,kbe)+yo(i-1,jc,kc,kb)
     &              *CI(ic,jc,kc-1,lxyr))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND         
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc,kc,kbse)=yo(i-1,jc,kc,kbse)
     &              *CI(ic,jc-1,kc-1,lxyr)+CI(ic,jc,kc,lxyl)
     &              *(yo(i,jc,kc,kbse)+yo(i-1,jc,kc,kbs)
     &              *CI(ic,jc-1,kc-1,lxyr))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND         
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc,kc,kbs)=yo(i,jc,kc,kbs)+yo(i,jc,kc,kbsw)
     &              *CI(ic,jc-1,kc-1,lxyr)+yo(i+1,jc,kc,kbse)
     &              *CI(ic+1,jc-1,kc-1,lxyl)+CI(ic,jc,kc,lxyr)
     &              *(yo(i,jc,kc,kbse)+yo(i-1,jc,kc,kbs)
     &              *CI(ic,jc-1,kc-1,lxyr))+CI(ic+1,jc,kc,lxyl)
     &              *(yo(i+1,jc,kc,kbsw)+yo(i+1,jc,kc,kbs)
     &              *CI(ic+1,jc-1,kc-1,lxyl))
            ENDDO
         ENDDO
      ENDDO

      DO kc=KCSTART,KCEND         
         DO jc=JCSTART,JCEND
            i=ISTART
            DO ic=ICSTART,ICEND
               i=i+2
               SOC(ic,jc,kc,kbsw)=yo(i,jc,kc,kbsw)
     &              *CI(ic,jc-1,kc-1,lxyl)+CI(ic,jc,kc,lxyr)
     &              *(yo(i-1,jc,kc,kbsw)+yo(i-1,jc,kc,kbs)
     &              *CI(ic,jc-1,kc-1,lxyl))
            ENDDO
         ENDDO
      ENDDO

      DO kpz=1,14

         IF ( kpz.NE.kpnw .AND. kpz.NE.kbnw .AND. 
     &	      kpz.NE.kbn  .AND. kpz.NE.kbne ) THEN

            ptrn = 1 
            call MSG_tbdx_send(SOC(1,1,1,kpz), buffer, 
     &           iWork(pMSGSO(ipL_MSG_NumAdjProc,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Proc,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Ipr,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Index,KGC)),
     &           ptrn, ierror)
         
            call MSG_tbdx_receive(SOC(1,1,1,kpz), buffer,
     &           iWork(pMSGSO(ipL_MSG_NumAdjProc,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Proc,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Ipr,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Index,KGC)),
     &           ptrn, ierror)
         
            call MSG_tbdx_close(SOC(1,1,1,kpz), buffer,
     &           iWork(pMSGSO(ipL_MSG_NumAdjProc,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Proc,KGC)),
     &           iWork(pMSGSO(ipL_MSG_Ipr,KGC)), 
     &           iWork(pMSGSO(ipL_MSG_Index,KGC)),
     &           ptrn, ierror)

         ENDIF

      ENDDO 

      LXGP = 1
      RXGP = 2
      LYGP = 2
      RYGP = 1
      LZGP = 1
      RZGP = 2

      CALL  BMG3_SymStd_SETUP_nMSGGrid(
     &           IGRD(KGC,idL_BMG_NGx), IGRD(KGC,idL_BMG_NGy),
     &           IGRD(KGC,idL_BMG_NGz),
     &           iONE, iONE, iONE,
     &           LXGP, RXGP, LYGP, RYGP, LZGP, RZGP,
     &           iWork(pMSGzo(ipL_MSG_LocalArraySize,1)),
     &           iWork(pMSGzo(ipL_MSG_GlobalCoordLocalData,1)),
     &           iWork(pMSGzo(ipL_MSG_GlobalCoordActData,1)),
     &           iWork(pMSGzo(ipL_MSG_ActDataStart,1)),
     &           DimX, DimY, DimZ, ProcGrid,
     &           NProc, NProcI, NProcJ, NProcK, NOGm, KGF
     &           )


      iWork(pMSGzo(ipL_MSG_NumAdjProc,1)) = 0

      CALL MSG_tp_setup (
     &        iWork(pMSGzo(ipL_MSG_LocalArraySize,1)),
     &        iWork(pMSGzo(ipL_MSG_ActDataStart,1)),
     &        iWork(pMSGzo(ipL_MSG_GlobalCoordLocalData,1)),
     &        iWork(pMSGzo(ipL_MSG_GlobalCoordActData,1)),
     &        NProc, MyProc,
     &        iWork(pMSGzo(ipL_MSG_NumAdjProc,1)),
     &        iWork(pMSGzo(ipL_MSG_Proc,1)),
     &        iWork(pMSGzo(ipL_MSG_Ipr,1)),
     &        iWork(pMSGzo(ipL_MSG_Index,1)),
     &        1,1,ierror)

      ptrn = 1 
      kpz = kpnw
      call MSG_tbdx_send(SOC(1,1,1,kpz), buffer,
     &     iWork(pMSGzo(ipL_MSG_NumAdjProc,1)),
     &     iWork(pMSGzo(ipL_MSG_Proc,1)),
     &     iWork(pMSGzo(ipL_MSG_Ipr,1)),
     &     iWork(pMSGzo(ipL_MSG_Index,1)),
     &     ptrn, ierror)

      call MSG_tbdx_receive(SOC(1,1,1,kpz), buffer,
     &     iWork(pMSGzo(ipL_MSG_NumAdjProc,1)),
     &     iWork(pMSGzo(ipL_MSG_Proc,1)),
     &     iWork(pMSGzo(ipL_MSG_Ipr,1)),
     &     iWork(pMSGzo(ipL_MSG_Index,1)),
     &     ptrn, ierror)

      call MSG_tbdx_close(SOC(1,1,1,kpz), buffer,
     &     iWork(pMSGzo(ipL_MSG_NumAdjProc,1)),
     &     iWork(pMSGzo(ipL_MSG_Proc,1)),
     &     iWork(pMSGzo(ipL_MSG_Ipr,1)),
     &     iWork(pMSGzo(ipL_MSG_Index,1)),
     &     ptrn, ierror)


      DO kpz=8,10

         ptrn = 1 
         call MSG_tbdx_send(SOC(1,1,1,kpz), buffer,
     &        iWork(pMSGzo(ipL_MSG_NumAdjProc,1)),
     &        iWork(pMSGzo(ipL_MSG_Proc,1)),
     &        iWork(pMSGzo(ipL_MSG_Ipr,1)),
     &        iWork(pMSGzo(ipL_MSG_Index,1)),
     &        ptrn, ierror)

         call MSG_tbdx_receive(SOC(1,1,1,kpz), buffer,
     &        iWork(pMSGzo(ipL_MSG_NumAdjProc,1)),
     &        iWork(pMSGzo(ipL_MSG_Proc,1)),
     &        iWork(pMSGzo(ipL_MSG_Ipr,1)),
     &        iWork(pMSGzo(ipL_MSG_Index,1)),
     &        ptrn, ierror)

         call MSG_tbdx_close(SOC(1,1,1,kpz), buffer,
     &        iWork(pMSGzo(ipL_MSG_NumAdjProc,1)),
     &        iWork(pMSGzo(ipL_MSG_Proc,1)),
     &        iWork(pMSGzo(ipL_MSG_Ipr,1)),
     &        iWork(pMSGzo(ipL_MSG_Index,1)),
     &        ptrn, ierror)

      ENDDO

C ==========================================================================
      
      RETURN
      END





