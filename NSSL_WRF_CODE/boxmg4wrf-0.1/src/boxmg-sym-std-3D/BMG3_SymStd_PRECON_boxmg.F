      SUBROUTINE BMG3_SymStd_PRECON_boxmg( 
     &                NLx, NLy, NLz, NGx, NGy, NGz, 
     &                iGs, jGs, kGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, NFm, 
     &                SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                BMG_rWORK_PL, NBMG_rWORK_PLm,
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                NProc, MSG_Buffer, NMSGrm, MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_PRECON_boxmg is used to wrap the call to the BoxMG
C     solver, BMG3_SymStd_SOLVE_boxmg, for preconditioning in
C     BMG3_SymStd_SOLVE_pcg.
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
#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  NBMG_iWORKm, NBMG_rWORKm, NOGm, NFm, NSOm,
     &         iGs, jGs, kGs, NGx, NGy, NGz, NLx, NLy, NLz,
     &         NBMG_iWORK_CSm, NBMG_rWORK_CSm,
     &         NBMG_iWORK_PLm, NBMG_rWORK_PLm

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_iWORK(NBMG_iWORKm),
     &         BMG_pWORK(NBMG_pWORK)

      INTEGER  BMG_iWORK_CS(NBMG_iWORK_CSm)
      REAL*RKIND   BMG_rWORK_CS(NBMG_rWORK_CSm)

      INTEGER  BMG_iWORK_PL(NBMG_iWORK_PLm)
      REAL*RKIND   BMG_rWORK_PL(NBMG_rWORK_PLm)

      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS), BMG_rWORK(NBMG_rWORKm),
     &         Q(NFm), QF(NFm), SO(NSOm)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      !
      ! MPI-Related Argument Declarations
      !
      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  MPICOMM, NProc, NMSGim, NMSGrm
      INTEGER  iWorkMSG(NMSGim), pMSG(NBMG_pMSG,NOGm), 
     &         pMSGSO(NBMG_pMSG,NOGm), pLS(NBMG_pLS*NOGm)
      REAL*RKIND   MSG_Buffer(NMSGrm) 

C ----------------------------
C     Local Declarations
C
      INTEGER  NOG, NF, NC, NSO, NCI, NSOR, NCBW, NCU,
     &         p_CI, p_CSO, p_CU, p_iG, p_RES, p_SOR

C ==========================================================================

      !
      !  BoxMG dimensional parameters
      !
      NOG  = BMG_iPARMS(id_BMG3_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG3_DIM_NF)
      NC   = BMG_iPARMS(id_BMG3_DIM_NC)
      NSO  = BMG_iPARMS(id_BMG3_DIM_NSO)
      NCI  = BMG_iPARMS(id_BMG3_DIM_NCI)
      NSOR = BMG_iPARMS(id_BMG3_DIM_NSOR)

      !
      !  BoxMG pointers
      !
      p_RES = BMG_pWORK(ip_RES)
      p_SOR = BMG_pWORK(ip_SOR)
      p_CI  = BMG_pWORK(ip_CI)
      p_CSO = BMG_pWORK(ip_CSO)
      p_CU  = BMG_pWORK(ip_CU)
      p_iG  = BMG_pWORK(ip_iG)


      CALL BMG3_SymStd_SOLVE_boxmg(
     &                 NLx, NLy, NLz, NGx, NGy, NGz, 
     &                 iGs, jGs, kGs,
     &                 BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                 Q, QF, BMG_rWORK(p_RES), NF, NC, 
     &                 SO, NSO, BMG_rWORK(p_SOR), NSOR, 
     &                 BMG_rWORK(p_CI), NCI, 
     &                 BMG_iWORK(p_IG), NOGm, NOG, 
     &                 BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                 BMG_rWORK_PL, NBMG_rWORK_PLm,
     &                 BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                 BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                 iWorkMSG, NMSGim, pMSG, pMSGSO,
     &                 BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                 NProc, MSG_Buffer, NMSGrm, MPICOMM
     &                 )
      IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
         RETURN
      END IF

C ==========================================================================

      RETURN
      END

