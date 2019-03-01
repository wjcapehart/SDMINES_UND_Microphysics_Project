      SUBROUTINE BMG3_SER_SymStd_PRECON_boxmg( 
     &                       Nx, Ny, Nz,
     &                       BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                       Q, QF, NFm, 
     &                       SO, NSOm, NOGm,
     &                       BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                       BMG_rWORK, NBMG_rWORKm,
     &                       BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                       BMG_rWORK_PL, NBMG_rWORK_PLm
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_PRECON_boxmg is used to wrap the call to the BoxMG
C     solver, BMG3_SER_SymStd_SOLVE_boxmg, for preconditioning in
C     BMG3_SER_SymStd_SOLVE_pcg.
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
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_stencils.h'
      INCLUDE 'BMG_SER_workspace.h'
      INCLUDE 'BMG_SER_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  NBMG_iWORKm, NBMG_rWORKm, NOGm, NFm, NSOm, Nx, Ny, Nz,
     &         NBMG_iWORK_PLm, NBMG_rWORK_PLm

      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS), BMG_iWORK(NBMG_iWORKm),
     &         BMG_pWORK(NBMG_SER_pWORK), BMG_iWORK_PL(NBMG_iWORK_PLm)
      REAL*RKIND   BMG_rPARMS(NBMG_SER_rPARMS), BMG_rWORK(NBMG_rWORKm),
     &         Q(NFm), QF(NFm), SO(NSOm), BMG_rWORK_PL(NBMG_rWORK_PLm)
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  NOG, NF, NC, NSO, NCI, NSOR, NCBW, NCU,
     &         p_CI, p_CSO, p_CU, p_iG, p_RES, p_SOR

C ==========================================================================

      !
      !  BoxMG dimensional parameters
      !
      NOG  = BMG_iPARMS(id_BMG3_SER_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG3_SER_DIM_NF)
      NC   = BMG_iPARMS(id_BMG3_SER_DIM_NC)
      NSO  = BMG_iPARMS(id_BMG3_SER_DIM_NSO)
      NCI  = BMG_iPARMS(id_BMG3_SER_DIM_NCI)
      NSOR = BMG_iPARMS(id_BMG3_SER_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG3_SER_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG3_SER_DIM_NCU)

      !
      !  BoxMG pointers
      !
      p_RES = BMG_pWORK(ip_RES)
      p_SOR = BMG_pWORK(ip_SOR)
      p_CI  = BMG_pWORK(ip_CI)
      p_CSO = BMG_pWORK(ip_CSO)
      p_CU  = BMG_pWORK(ip_CU)
      p_iG  = BMG_pWORK(ip_iG)
      

      CALL BMG3_SER_SymStd_SOLVE_boxmg( 
     &                 Nx, Ny, Nz,
     &                 BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                 Q, QF, BMG_rWORK(p_RES), NF, NC, 
     &                 SO, NSO, BMG_rWORK(p_SOR), NSOR,
     &                 BMG_rWORK(p_CI), NCI,
     &                 BMG_rWORK(p_CSO), BMG_rWORK(p_CU), NCBW, NCU,
     &                 BMG_iWORK(p_iG), NOGm, NOG,
     &                 BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                 BMG_rWORK_PL, NBMG_rWORK_PLm
     &                 )
      IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
         RETURN
      END IF

C ==========================================================================

      RETURN
      END

