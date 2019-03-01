      SUBROUTINE BMG3_SymStd_residual(
     &                KG, NOG, IFD,
     &                Q, QF, SO, RES, II, JJ, KK, NStncl,
     &                iWork, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_residual computes the residual on level kg.
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
      INCLUDE 'mpif.h'
      INCLUDE 'MSG.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  IFD, II, JJ, KG, KK, MyProc, NOG, NMSGi, NMSGr
      INTEGER  NStncl, MPICOMM, iWork(NMSGi), pMSG(NBMG_pMSG,NOG)
      REAL*RKIND   q(II,JJ,KK), qf(II,JJ,KK), RES(II,JJ,KK)
      REAL*RKIND   so(II+1,JJ+1,KK+1,NStncl)
      REAL*RKIND   MSG_Buffer(NMSGr)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, i1, j, j1, k, k1
      INTEGER  ierror, ptrn, MPI_IERROR

      REAL*RKIND SOx(II+1),SOy(JJ+1),SOz(KK+1)
      
      INTEGER ISEP

C ==========================================================================

      ISEP = BMG_NONSEPARABLE


      i1=II-1
      j1=JJ-1
      k1=KK-1
      
      IF( KG.LT.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN

         DO k=2,k1
            DO j=2,j1
               DO i=2,i1
                  RES(i,j,k) = qf(i,j,k)
     &                       + so(i,j,k,kpw)*q(i-1,j,k)
     &                       + so(i,j+1,k,kpnw)*q(i-1,j+1,k)
     &                       + so(i,j+1,k,kps)*q(i,j+1,k)
     &                       + so(i+1,j+1,k,kpsw) *q(i+1,j+1,k)
     &                       + so(i+1,j,k,kpw)*q(i+1,j,k)
     &                       + so(i+1,j,k,kpnw)*q(i+1,j-1,k)
     &                       + so(i,j,k,kps)*q(i,j-1,k)
     &                       + so(i,j,k,kpsw)*q(i-1,j-1,k)
     &                       + so(i,j,k,kb)*q(i,j,k-1)
     &                       + so(i,j,k,kbw)*q(i-1,j,k-1)
     &                       + so(i,j+1,k,kbnw)*q(i-1,j+1,k-1)
     &                       + so(i,j+1,k,kbn)*q(i,j+1,k-1)
     &                       + so(i+1,j+1,k,kbne)*q(i+1,j+1,k-1)
     &                       + so(i+1,j,k,kbe)*q(i+1,j,k-1)
     &                       + so(i+1,j,k,kbse)*q(i+1,j-1,k-1)
     &                       + so(i,j,k,kbs)*q(i,j-1,k-1)
     &                       + so(i,j,k,kbsw)*q(i-1,j-1,k-1)
     &                       + so(i,j,k+1,kb)*q(i,j,k+1)
     &                       + so(i,j,k+1,kbe)*q(i-1,j,k+1)
     &                       + so(i,j+1,k+1,kbse)*q(i-1,j+1,k+1)
     &                       + so(i,j+1,k+1,kbs)*q(i,j+1,k+1)
     &                       + so(i+1,j+1,k+1,kbsw)*q(i+1,j+1,k+1)
     &                       + so(i+1,j,k+1,kbw)*q(i+1,j,k+1)
     &                       + so(i+1,j,k+1,kbnw)*q(i+1,j-1,k+1)
     &                       + so(i,j,k+1,kbn)*q(i,j-1,k+1)
     &                       + so(i,j,k+1,kbne)*q(i-1,j-1,k+1)
     &                       - so(i,j,k,kp)*q(i,j,k)
               ENDDO
            ENDDO
         ENDDO

      ELSE
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
      


         DO k=2,k1
            DO j=2,j1
               DO i=2,i1
                  RES(i,j,k) = qf(i,j,k)
     &                        + SOx(i)*Q(i-1,j,k)
     &                        + SOy(j+1)*Q(i,j+1,k)
     &                        + SOx(i+1)*Q(i+1,j,k)
     &                        + SOy(j)*Q(i,j-1,k)
     &                        + SOz(k)*Q(i,j,k-1)
     &                        + SOz(k+1)*Q(i,j,k+1)
     &                       - so(i,j,k,kp)*q(i,j,k)
               ENDDO
            ENDDO
         ENDDO

       ELSE

         DO k=2,k1
            DO j=2,j1
               DO i=2,i1
                  RES(i,j,k) = qf(i,j,k)
     &                       + so(i,j,k,kpw)*q(i-1,j,k)
     &                       + so(i,j+1,k,kps)*q(i,j+1,k)
     &                       + so(i+1,j,k,kpw)*q(i+1,j,k)
     &                       + so(i,j,k,kps)*q(i,j-1,k)
     &                       + so(i,j,k,kb)*q(i,j,k-1)
     &                       + so(i,j,k+1,kb)*q(i,j,k+1)
     &                       - so(i,j,k,kp)*q(i,j,k)
               ENDDO
            ENDDO
         ENDDO
         
      ENDIF
      
      ENDIF


      ptrn = 1

      CALL MSG_tbdx_send(RES, MSG_Buffer, 
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &     iWork(pMSG(ipL_MSG_Proc,KG)),
     &     iWork(pMSG(ipL_MSG_Ipr,KG)),
     &     iWork(pMSG(ipL_MSG_Index,KG)),
     &     ptrn, ierror)

      CALL MSG_tbdx_receive(RES, MSG_Buffer,
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &     iWork(pMSG(ipL_MSG_Proc,KG)),
     &     iWork(pMSG(ipL_MSG_Ipr,KG)),
     &     iWork(pMSG(ipL_MSG_Index,KG)),
     &     ptrn, ierror)

      CALL MSG_tbdx_close(RES, MSG_Buffer,
     &     iWork(pMSG(ipL_MSG_NumAdjProc,KG)),
     &     iWork(pMSG(ipL_MSG_Proc,KG)),
     &     iWork(pMSG(ipL_MSG_Ipr,KG)),
     &     iWork(pMSG(ipL_MSG_Index,KG)),
     &     ptrn, ierror)


C ==========================================================================

      RETURN
      END
