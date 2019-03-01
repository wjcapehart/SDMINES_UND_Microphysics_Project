      SUBROUTINE BMG2_SymStd_SETUP_ITLI_bl( 
     &                KF, KC, SO, SOC, CI, 
     &                IIF, JJF, IIC, JJC, iGs, jGs, 
     &                NOG, IFD, NStncl, BMG_iPARMS, BMG_IOFLAG,
     &                iWork, NMSGi, pMSGSO, 
     &                BUFFER, NMSGr, MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SETUP_cg_ITLI constructs the Galerkin (variational)
C     coarse-grid operator on the coasre grid, KC, given the fine-grid
C     stencil, SO, and interpolation stencil, CI, on the fine-grid, KF.
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

#include      "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'
      
C ---------------------------
C    Argument Declarations:
C
      INTEGER   iGs, IIC, IIF, IFD, jGs, JJC, JJF, KC, KF,
     &          NMSGi, NMSGr, NOG, NStncl, MPICOMM

      INTEGER   iWork(NMSGi), pMSGSO(NBMG_pMSG,NOG)
      REAL*RKIND    CI(IIC,JJC,8), SO(IIF+1,JJF+1,NStncl), 
     &          SOC(IIC+1,JJC+1,5), 
     &          BUFFER(NMSGr)

      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER   BMG_iPARMS(NBMG_IPARMS)

C --------------------------
C     Externals
C      
      REAL*RKIND    BMG2_SymStd_UTILS_reswt
      EXTERNAL  BMG2_SymStd_UTILS_reswt

C --------------------------
C     Local Declarations:
C
      INTEGER   IC, I, IIC1, IIF1, IL, ILG, ILGP, INDA, INDB, INDC,
     &          JC, J, JJC1, JJF1, JL, JLG, JLGP,
     &          ptrn, ierror
      REAL*RKIND    CONW, CENW, CSENW, CFNW, CSNW, Q(7,7,2), R(7,7,2)

      INTEGER  MyProc, ISTART, JSTART

C ==========================================================================

C ----------------------------------
C     Sanity Check:
C ----------------------------------
      
      IF (KF-1.NE.KC ) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,*) 'ERROR: BMG2_SymStd_SETUP_cg_ITLI   .... '
            WRITE(*,*) '*****  KC = ', KC
            WRITE(*,*) '*****  KF = ', KF
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,16)
         RETURN

      ENDIF

      MyProc = MSG_MyProc(MPICOMM)

C -----------------------------------
C     Useful indexing bounds:
C -----------------------------------

      IIC1=IIC-1
      JJC1=JJC-1

      IIF1=IIF-1
      JJF1=JJF-1

      if (mod(iGs,2).eq.1) THEN
         ISTART = 0
      else
         ISTART = 1
      endif
      
      if (mod(jGs,2).eq.1) then
         JSTART = 0
      else
         JSTART = 1
      endif
      

      IF ( IFD.NE.1 .OR. KF.LT.NOG ) THEN

C******************************
C     begin computation of grid kc difference operator when kf difference
C     operator is nine point unless kc. ge. NOG
C     
         J=JSTART

         DO JC=2,JJC1

            J=J+2

            INDA = 1
            INDB = 2

            I  = ISTART
            IC = 1

            CALL BMG2_SymStd_UTILS_bbody9(
     &                IIF, JJF, IIC, JJC, IC, JC, I, J, NStncl,
     &                INDA, CI, SO, Q, R)           

            DO IC=2,IIC1

               I=I+2

               CALL BMG2_SymStd_UTILS_bbody9(
     &              IIF, JJF, IIC, JJC, IC, JC, I,J, NStncl,
     &              INDB, CI, SO, Q, R)

               !
               ! Apply residual weighting
               !
               SOC(IC,JC,KO) 
     &         =   BMG2_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, IC, JC, 4, 4, INDB 
     &                  )
               SOC(IC,JC,KW) 
     &         = - BMG2_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, IC-1, JC, 2, 4, INDB 
     &                  )
               SOC(IC,JC,KS) 
     &         = - BMG2_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, IC, JC-1, 4, 2, INDB 
     &                  )
               SOC(IC,JC,KSW) 
     &         = - BMG2_SymStd_UTILS_reswt(
     &                  CI, R, IIC, JJC, IC-1, JC-1, 2, 2, INDB
     &                  )
               SOC(IC,JC,KNW) 
     &              = - BMG2_SymStd_UTILS_reswt(
     &              CI, R, IIC, JJC, IC, JC-1, 6, 2, INDA
     &              )
            !

            INDC = INDB
            INDB = INDA
            INDA = INDC

            ENDDO
         ENDDO

C     end of computation of kc difference operator when kf difference
C     operator is nine point
C******************************

      ELSE

C******************************
C   begin computation of kc difference operator when kf difference
C   operator is five point unless kc.ge.NOG
C
         J=JSTART

         DO JC=2,JJC1

            J=J+2

            INDA = 1
            INDB = 2

            I = ISTART
            IC = 1

            CALL BMG2_SymStd_UTILS_bbody5(
     &                IIF, JJF, IIC, JJC, IC, JC, I, J, NStncl,
     &                INDA, CI, SO, Q, R)           

            DO IC=2,IIC1

               I=I+2

               CALL BMG2_SymStd_UTILS_bbody5(
     &              IIF, JJF, IIC, JJC, IC, JC, I, J, NStncl,
     &              INDB, CI, SO, Q, R)

              !
              !  Apply residual weighting
              !
              SOC(IC,JC,KO) 
     &        =   BMG2_SymStd_UTILS_reswt(
     &                 CI, R, IIC, JJC, IC, JC, 4, 4, INDB 
     &                 )
              SOC(IC,JC,KW) 
     &        = - BMG2_SymStd_UTILS_reswt(
     &                 CI, R, IIC, JJC, IC-1, JC, 2, 4, INDB 
     &                 )
              SOC(IC,JC,KS) 
     &        = - BMG2_SymStd_UTILS_reswt(
     &                 CI, R, IIC, JJC, IC, JC-1, 4, 2, INDB
     &                 )
              SOC(IC,JC,KSW) 
     &        = - BMG2_SymStd_UTILS_reswt(
     &                 CI, R, IIC, JJC, IC-1, JC-1, 2, 2, INDB
     &             )
              SOC(IC,JC,KNW)
     &        = - BMG2_SymStd_UTILS_reswt(
     &                 CI, R, IIC, JJC, IC, JC-1, 6, 2, INDA
     &            )

              INDC = INDB
              INDB = INDA
              INDA = INDC

            ENDDO
          ENDDO

        ENDIF  


      !
      ! update the ghost boundaries of SOC
      ! note that the coarse stencil is always 
      ! nine point
      !
      DO I=1,5
         ptrn = 6
         
         call MSG_tbdx_send(SOC(1,1,I), buffer, 
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,KC)),
     &        iWork(pMSGSO(ipL_MSG_Proc,KC)),
     &        iWork(pMSGSO(ipL_MSG_Ipr,KC)),
     &        iWork(pMSGSO(ipL_MSG_Index,KC)),
     &        ptrn, ierror)
         
         call MSG_tbdx_receive(SOC(1,1,I), buffer,
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,KC)),
     &        iWork(pMSGSO(ipL_MSG_Proc,KC)),
     &        iWork(pMSGSO(ipL_MSG_Ipr,KC)),
     &        iWork(pMSGSO(ipL_MSG_Index,KC)),
     &        ptrn, ierror)
         
         call MSG_tbdx_close(SOC(1,1,I), buffer,
     &        iWork(pMSGSO(ipL_MSG_NumAdjProc,KC)),
     &        iWork(pMSGSO(ipL_MSG_Proc,KC)),
     &        iWork(pMSGSO(ipL_MSG_Ipr,KC)),
     &        iWork(pMSGSO(ipL_MSG_Index,KC)),
     &        ptrn, ierror)
         
      ENDDO


C   end of computation of grid kc difference operator, when kf
C   difference operator is five point
C******************************

      RETURN
      END
