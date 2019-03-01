      SUBROUTINE BMG3_SymStd_relax_GS( 
     &                KG, SO, QF, Q, RES, SOR,
     &                NLx, NLy, NLz, NGx, NGy, NGz, 
     &                RES_L2, BMG_IOFLAG, NOG, NOGm,
     &                IFD, NStncl, NSORv, IRELAX, IRELAX_SYM, UPDOWN,
     &                iGs, jGs, kGs,
     &                MyProcI, MyProcJ, MyProcK, MyProc,
     &                ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &                iWork, NMSGi, pMSG, MSG_Buffer, NMSGr, MPICOMM
     &               ,ISEP)

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_relax_GS performs one sweep of Gauss Seidel (with the
C     correct ordering depending on whether we are on the way down, or
C     up in the symmetric cycling case)
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

C ------------------------------------------------
C     Includes
C
      INCLUDE 'mpif.h'
      INCLUDE 'MSG.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

C ------------------------------------------------
C     Argument Declarations
C
      INTEGER NGx, NGy, NGz, NLx, NLy, NLz

      INTEGER IFD, IRELAX, IRELAX_SYM, KG, NOG, NOGm,
     &        NSORv, NStncl, UPDOWN, NMSGi, NMSGr

      INTEGER iGs, jGs, kGs
      INTEGER MyProc, MyProcI, MyProcJ, MyProcK,
     &        NProc, NProcI, NProcJ, NProcK
      INTEGER ProcGrid(NProcI,NProcJ,NProcK)

      INTEGER MPICOMM, iWork(NMSGi), pMSG(NBMG_pMSG,NOGm)
      REAL*RKIND  Q(NLx,NLy,NLz), QF(NLx,NLy,NLz), 
     &        RES(NLx,NLy,NLz),
     &        SO(NLx+1,NLy+1,NLz+1,NStncl), SOR(NLx,NLy,NLz,NSORv)
      
      REAL*8 RES_L2
      
      REAL*RKIND  MSG_Buffer(NMSGr)

      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

C ------------------------------------------------
C     Local Declarations
C
      INTEGER i, i1, j ,j1, k, k1
      
      INTEGER ie, iw, js, jn, kbot, kt

      INTEGER pts, ibeg, iend, ptsred, ptsblk
      INTEGER ptstart, ptend, ptstride
      INTEGER ptrn, ierror
      INTEGER relaxtest
      LOGICAL comhalf
      
      REAL*RKIND SOx(NLx+1),SOy(NLy+1),SOz(NLz+1)
      
      INTEGER ISEP

C ==========================================================================

      j1 = NLy-1
      i1 = NLx-1
      k1 = NLz-1
      
c      write(0,*) 'mproc, mpi,mpj,mpk, igs,jgs,kgs = ',myproc,myproci,
c     :   myprocj,myprock,igs,jgs,kgs
c      write(0,*) 'mproc, isep, BMG_SEPARABLE = ',myproc,isep,
c     :   BMG_SEPARABLE, BMG_NONSEPARABLE
      

c      write(0,*) 'RELAX_GS: kg,nog,ifd = ',KG,NOG,ifd

      comhalf = .false.

      
      IF ( KG.LT.NOG .OR. ifd.NE.1 ) THEN 

         !
         !   27-point relaxation (8-color)
         !

C         IF ((UPDOWN.eq.BMG_UP).or.(IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) 
C     &              THEN
C            ptstart = 1
C            ptend   = 8
C            ptstride = 1
C         ELSE
C            ptstart = 8
C            ptend = 1
C            ptstride = -1
C         ENDIF         
C        
C         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<
C
C            DO k=2+mod(mod((pts-1)/4,2)+mod(kGs+1,2),2),K1,2
C               !
C               DO j=2+mod(mod((pts-1)/2,2)+mod(jGs+1,2),2),J1,2
C                  !
C                  DO i=2+mod(mod(pts-1,2)+mod(iGs+1,2),2),I1,2
                     !

         IF ( (UPDOWN.eq.BMG_UP)
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF         
        
        
        relaxtest = 2
        
        IF ( relaxtest .eq. 0 ) THEN ! default
        
         DO pts = ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

c            DO k=2+mod(mod((pts-1)/4,2)+mod(kGs+1,2),2),K1,2
c               !
c               DO j=2+mod(mod((pts-1)/2,2)+mod(jGs+1,2),2),J1,2
c                  !
c                  DO i=2+mod(mod(pts-1,2)+mod(iGs+1,2),2),I1,2
             DO k = 2,k1
              
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                     !
                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            !
            
            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
            !
            
            ENDIF
          ENDDO

        ELSEIF ( relaxtest .eq. 1 ) THEN ! non-redblack GS

             DO k = 2,k1
              
               DO j=2, j1
                  !
                  !
                  DO i = 2,i1
                     !
                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO

        ELSEIF ( relaxtest .eq. 2 ) THEN ! layered redblack GS

c         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

c do 'red' points on first plane
            ptsred = ptstart
            ptsblk = ptend
            DO  k=2, 2
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsred,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO

c do 'red' points on k plane and 'black' points on k-1 plane

            DO  k=3, k1

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsred,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO

            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
            !
            
            ENDIF

c black points
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-2)+ptsblk,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k-1) = ( QF(i,j,k-1) 
     &                    + SO(i,j,k-1,kpw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k-1,kpnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k-1,kps)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k-1,kpsw)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k-1,kpw)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k-1,kpnw)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k-1,kps)*Q(i,j-1,k-1)
     &                    + SO(i,j,k-1,kpsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k-1,kb)*Q(i,j,k-1-1)
     &                    + SO(i,j,k-1,kbw)*Q(i-1,j,k-1-1)
     &                    + SO(i,j+1,k-1,kbnw)*Q(i-1,j+1,k-1-1)
     &                    + SO(i,j+1,k-1,kbn)*Q(i,j+1,k-1-1)
     &                    + SO(i+1,j+1,k-1,kbne)*Q(i+1,j+1,k-1-1)
     &                    + SO(i+1,j,k-1,kbe)*Q(i+1,j,k-1-1)
     &                    + SO(i+1,j,k-1,kbse)*Q(i+1,j-1,k-1-1)
     &                    + SO(i,j,k-1,kbs)*Q(i,j-1,k-1-1)
     &                    + SO(i,j,k-1,kbsw)*Q(i-1,j-1,k-1-1)
     &                    + SO(i,j,k-1+1,kb)*Q(i,j,k-1+1)
     &                    + SO(i,j,k-1+1,kbe)*Q(i-1,j,k-1+1)
     &                    + SO(i,j+1,k-1+1,kbse)*Q(i-1,j+1,k-1+1)
     &                    + SO(i,j+1,k-1+1,kbs)*Q(i,j+1,k-1+1)
     &                    + SO(i+1,j+1,k-1+1,kbsw)*Q(i+1,j+1,k-1+1)
     &                    + SO(i+1,j,k-1+1,kbw)*Q(i+1,j,k-1+1)
     &                    + SO(i+1,j,k-1+1,kbnw)*Q(i+1,j-1,k-1+1)
     &                    + SO(i,j,k-1+1,kbn)*Q(i,j-1,k-1+1)
     &                    + SO(i,j,k-1+1,kbne)*Q(i-1,j-1,k-1+1)
     &                    )*SOR(i,j,k-1,msor)
                     !
                  ENDDO
                  !
               ENDDO

               !
            ENDDO

c do 'black' points on last plane
            DO  k=k1, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+ptsblk,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !

                     Q(i,j,k) = ( QF(i,j,k) 
     &                    + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                    + SO(i,j+1,k,kpnw)*Q(i-1,j+1,k)
     &                    + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                    + SO(i+1,j+1,k,kpsw)*Q(i+1,j+1,k)
     &                    + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                    + SO(i+1,j,k,kpnw)*Q(i+1,j-1,k)
     &                    + SO(i,j,k,kps)*Q(i,j-1,k)
     &                    + SO(i,j,k,kpsw)*Q(i-1,j-1,k)
     &                    + SO(i,j,k,kb)*Q(i,j,k-1)
     &                    + SO(i,j,k,kbw)*Q(i-1,j,k-1)
     &                    + SO(i,j+1,k,kbnw)*Q(i-1,j+1,k-1)
     &                    + SO(i,j+1,k,kbn)*Q(i,j+1,k-1)
     &                    + SO(i+1,j+1,k,kbne)*Q(i+1,j+1,k-1)
     &                    + SO(i+1,j,k,kbe)*Q(i+1,j,k-1)
     &                    + SO(i+1,j,k,kbse)*Q(i+1,j-1,k-1)
     &                    + SO(i,j,k,kbs)*Q(i,j-1,k-1)
     &                    + SO(i,j,k,kbsw)*Q(i-1,j-1,k-1)
     &                    + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                    + SO(i,j,k+1,kbe)*Q(i-1,j,k+1)
     &                    + SO(i,j+1,k+1,kbse)*Q(i-1,j+1,k+1)
     &                    + SO(i,j+1,k+1,kbs)*Q(i,j+1,k+1)
     &                    + SO(i+1,j+1,k+1,kbsw)*Q(i+1,j+1,k+1)
     &                    + SO(i+1,j,k+1,kbw)*Q(i+1,j,k+1)
     &                    + SO(i+1,j,k+1,kbnw)*Q(i+1,j-1,k+1)
     &                    + SO(i,j,k+1,kbn)*Q(i,j-1,k+1)
     &                    + SO(i,j,k+1,kbne)*Q(i-1,j-1,k+1)
     &                    )*SOR(i,j,k,msor)
                     !
                  ENDDO
                  !
               ENDDO
               !
            ENDDO
            
       ENDIF ! relaxtest
       !
c         ENDDO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<

        IF ( relaxtest .ne. 0 ) THEN
 
            ptrn = 1
            !
            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
     
        ENDIF
            !
!         ENDDO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<


         !
      ELSE
         !  
         !  7 point relaxation (2-color)
         !         
         IF ( (UPDOWN.eq.BMG_UP)
     &        .OR. (IRELAX_SYM.EQ.BMG_RELAX_NONSYM)) THEN
            ptstart  = 0
            ptend    = 1
            ptstride = 1
         ELSE
            ptstart  = 1
            ptend    = 0
            ptstride = -1
         ENDIF         

      IF ( ISEP .eq. BMG_SEPARABLE ) THEN
      
      DO  k=1, k1+1
        SOz(k) = SO(3,3,k,kb)
c       IF ( myproc == 1 ) write(0,*) 'k, SOz = ',k,SOz(k)
      ENDDO
c        SOz(2) = SOz(3)
c        SOz(k1+1) = SOz(k1)
        
      
      DO j=1, j1+1
        SOy(j) = SO(3,j,3,kps)
c        IF ( myproc == 1 ) write(0,*) 'j, SOy = ',j,SOy(j)
c        write(0,*) 'j, SOy, myproc = ',j,SOy(j),myproc
      ENDDO
c        SOy(2) = SOy(3)
c        SOy(j1+1) = SOy(j1)
      
      
      DO i=1,i1+1
        SOx(i) = SO(i,3,3,kpw)
c        IF ( myproc == 1 ) write(0,*) 'i, SOx = ',i,SOx(i)
      ENDDO
c        SOx(2) = SOx(3)
c        SOx(i1+1) = SOx(i1)
      
      ENDIF

         DO pts=ptstart, ptend, ptstride ! >>> BEGIN: loop over colors <<<

          IF ( ISEP .eq. BMG_SEPARABLE ) THEN

          IF ( .true. ) THEN
          
          
            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SOx(i)*Q(i-1,j,k)
     &                        + SOy(j+1)*Q(i,j+1,k)
     &                        + SOx(i+1)*Q(i+1,j,k)
     &                        + SOy(j)*Q(i,j-1,k)
     &                        + SOz(k)*Q(i,j,k-1)
     &                        + SOz(k+1)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)
                  END DO
                  !
               END DO
               !
            END DO
            
          ELSE ! old try

          kbot = 0
          kt = 0
          js = 0
          jn = 0
          iw = 0
          ie = 0
          
          IF (  MyProcK .eq. 1 )     kbot = 1
          IF (  MyProcK .eq. NProcK) kt = 1
          IF (  MyProcJ .eq. 1)      js = 1
          IF (  MyProcJ .eq. NProcJ) jn = 1
          IF (  MyProcI .eq. 1)      iw = 1
          IF (  MyProcI .eq. NProcI) ie = 1
          
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  
                  IF ( ibeg .eq. 2 ) ibeg = ibeg + 2*iw
                  IF ( iend .eq. i1 ) iend = iend - 2*ie
                  
                  !
                  DO i=ibeg,iend,2
                  !
c                  IF ( .not. (( (k.eq.2      .and. MyProcK.eq.1 ) .or.
c     &                        ( k .eq. k1    .and. MyProcK .eq. NProcK) .or.
c     &                        ( j .eq. 2     .and. MyProcJ .eq. 1)  .or.
c     &                        ( j .eq. j1    .and. MyProcJ .eq. NProcJ) .or.
c     &                        ( i .eq. ibeg  .and. MyProcI .eq. 1)  .or. 
c     &                        ( i .eq. iend  .and. MyProcI .eq. NProcI) ) ) THEN
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SOx(i)*Q(i-1,j,k)
     &                        + SOy(j+1)*Q(i,j+1,k)
     &                        + SOx(i+1)*Q(i+1,j,k)
     &                        + SOy(j)*Q(i,j-1,k)
     &                        + SOz(k)*Q(i,j,k-1)
     &                        + SOz(k+1)*Q(i,j,k+1)
     &                        )/(SOx(i) + SOx(i+1) + SOy(j) + 
     &                           SOy(j+1) + SOz(k) + SOz(k+1))

                  END DO
                  !
               END DO
               !
            END DO

c edges:

           IF ( kbot .eq. 1 ) THEN
            k=2

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg 
                  !
                  DO i=ibeg,iend,2

                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)


                 END DO
              END DO
           
           ENDIF

           IF ( kt .eq. 1 ) THEN
            k=k1

               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg 
                  !
                  DO i=ibeg,iend,2

                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)


                 END DO
              END DO
           
           ENDIF
           
           IF ( js .eq. 1 ) THEN
                      
            DO  k=2+kbot, k1-kt
              j = 2
              ibeg = mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend = 2*((i1-ibeg)/2)+ibeg 

               DO i=ibeg,iend,2
    
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)

                END DO
             END DO
          ENDIF

           
           IF ( jn .eq. 1 ) THEN
                      
            DO  k=2+kbot, k1-kt
              j = j1
              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg 
                  DO i=ibeg,iend,2

    
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)

                END DO
             END DO
          ENDIF

           
           IF ( ie .eq. 1 ) THEN
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn

              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg 
              i = ibeg
               IF ( ibeg .eq. 2 ) THEN
    
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)

                ENDIF
                END DO
             END DO
          ENDIF

           IF ( iw .eq. 1 ) THEN
            DO  k=2+kbot, k1-kt
               DO j=2+js, j1-jn

              ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
              iend=2*((i1-ibeg)/2)+ibeg 
              
              i = iend
              
               IF ( i .eq. i1 ) THEN
    
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)

                ENDIF
                END DO
             END DO
          ENDIF
          
         ENDIF
         

         ELSE  ! Non separable

            DO  k=2, k1
               DO j=2, j1
                  !
                  ibeg=mod((iGs-1)+(j+jGs-1)+(k+kGs-1)+pts,2) + 2
                  iend=2*((i1-ibeg)/2)+ibeg
                  !
                  DO i=ibeg,iend,2
                  !
                     Q(i,j,k) = ( QF(i,j,k)
     &                        + SO(i,j,k,kpw)*Q(i-1,j,k)
     &                        + SO(i,j+1,k,kps)*Q(i,j+1,k)
     &                        + SO(i+1,j,k,kpw)*Q(i+1,j,k)
     &                        + SO(i,j,k,kps)*Q(i,j-1,k)
     &                        + SO(i,j,k,kb)*Q(i,j,k-1)
     &                        + SO(i,j,k+1,kb)*Q(i,j,k+1)
     &                        )*SOR(i,j,k,msor)
                  END DO
                  !
               END DO
               !
            END DO

          ENDIF
            !

            IF ( comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
            !
            
            ENDIF
            
            
            !
         END DO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<
            !
            IF ( .not. comhalf ) THEN

            ptrn = 1

            CALL MSG_tbdx_send(Q, MSG_Buffer, 
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)

            CALL MSG_tbdx_receive(Q, MSG_Buffer,
     &           iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &           iWork(pMSG(ipL_MSG_Proc,KG)),
     &           iWork(pMSG(ipL_MSG_Ipr,KG)),
     &           iWork(pMSG(ipL_MSG_Index,KG)),
     &           ptrn, ierror)
            !
            
            ENDIF
            !
C         END DO   ! >>> END: loop over colors <<<<<<<<<<<<<<<<<<<<<<<<<<

      ENDIF

      ptrn = 1

      CALL MSG_tbdx_close(Q, MSG_Buffer,
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &     iWork(pMSG(ipL_MSG_Proc,KG)),
     &     iWork(pMSG(ipL_MSG_Ipr,KG)),
     &     iWork(pMSG(ipL_MSG_Index,KG)),
     &     ptrn, ierror)

      IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) ) THEN
         !
         CALL BMG3_SymStd_residual( 
     &             KG, NOG, IFD,
     &             Q, QF, SO, RES, NLx, NLy, NLz, NStncl, 
     &             iWork, NMSGi, pMSG, 
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )
         CALL BMG3_SymStd_UTILS_norm_l2( 
     &             RES, NLx, NLy, NLz, RES_L2, MPICOMM
     &             )
         !
      ENDIF

C ==========================================================================

      RETURN
      END

