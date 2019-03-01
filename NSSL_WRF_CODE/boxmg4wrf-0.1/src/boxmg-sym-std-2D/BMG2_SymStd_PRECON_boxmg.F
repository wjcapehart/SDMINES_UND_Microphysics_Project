      SUBROUTINE BMG2_SymStd_PRECON_boxmg( 
     &                NLx, NLy, NGx, NGy, iGs, jGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, NFm, SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &                MSG_Buffer, NMSGrm, MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_PRECON_boxmg is used to wrap the call to the BoxMG
C     solver, BMG2_SymStd_SOLVE_boxmg, for preconditioning in
C     BMG2_SymStd_SOLVE_pcg.
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

C ---------------------------------------------
C     Includes
C ---------------------------------------------

#include      "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ---------------------------------------------
C     Argument Declarations
C ---------------------------------------------

      INTEGER  iGs, jGs, NGx, NGy, NLx, NLy

C ---------------------------------------------
C     BOXMG Specific Argument Declarations
C ---------------------------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER   BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND    BMG_rPARMS(NBMG_rPARMS)
      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG)

      !
      ! Workspace: Generic
      !
      INTEGER   NBMG_iWORKm, NBMG_rWORKm
      INTEGER   BMG_iWORK(NBMG_iWORKm)
      REAL*RKIND    BMG_rWORK(NBMG_rWORKm)

      !
      ! Workspace: Coarse-grid Solve
      !
      INTEGER   NBMG_iWORK_CSm, NBMG_rWORK_CSm
      INTEGER   BMG_iWORK_CS(NBMG_iWORK_CSm)
      REAL*RKIND    BMG_rWORK_CS(NBMG_rWORK_CSm)

      !
      ! Workspace: pointers
      !
      INTEGER   BMG_pWORK(NBMG_pWORK)

      !
      ! Solution, RHS, and Stencil
      !
      INTEGER   NFm, NOGm, NSOm
      REAL*RKIND    Q(NFm), QF(NFm), SO(NSOm)

C ---------------------------------------------
C     MPI-Related Argument Declarations
C ---------------------------------------------

      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  MPICOMM, NProc, NMSGim, NMSGrm
      INTEGER  iWorkMSG(NMSGim), pMSG(NBMG_pMSG,NOGm),
     &         pMSGSO(NBMG_pMSG,NOGm), pLS(NBMG_pLS*NOGm)
      REAL*RKIND   MSG_Buffer(NMSGrm) 

C ---------------------------------------------
C     Local Declarations
C ---------------------------------------------

      INTEGER  NOG, NF, NC, NSO, NCI, NSOR, NCBW, NCU,
     &         p_CI, p_CSO, p_CU, p_iG, p_RES, p_SOR

C ==========================================================================

      !
      !  BoxMG dimensional parameters
      !
      NOG  = BMG_iPARMS(id_BMG2_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG2_DIM_NF)
      NC   = BMG_iPARMS(id_BMG2_DIM_NC)
      NSO  = BMG_iPARMS(id_BMG2_DIM_NSO)
      NCI  = BMG_iPARMS(id_BMG2_DIM_NCI)
      NSOR = BMG_iPARMS(id_BMG2_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG2_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG2_DIM_NCU)

      !
      !  BoxMG pointers
      !
      p_RES = BMG_pWORK(ip_RES)
      p_SOR = BMG_pWORK(ip_SOR)
      p_CI  = BMG_pWORK(ip_CI)
      p_CSO = BMG_pWORK(ip_CSO)
      p_CU  = BMG_pWORK(ip_CU)
      p_iG  = BMG_pWORK(ip_iG)

      CALL BMG2_SymStd_SOLVE_boxmg(
     &          NLx, NLy, NGx, NGy, iGs, jGs,
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          Q, QF, BMG_rWORK(p_RES), NF, NC, 
     &          SO, NSO, BMG_rWORK(p_SOR), NSOR, 
     &          BMG_rWORK(p_CI), NCI, 
     &          BMG_iWORK(p_IG), NOGm, NOG, 
     &          BMG_iWORK_CS, NBMG_iWORK_CSm,
     &          BMG_rWORK_CS, NBMG_rWORK_CSm,
     &          BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &          iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &          MSG_Buffer, NMSGrm, MPICOMM
     &          )

C ==========================================================================

      RETURN
      END

