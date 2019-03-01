      SUBROUTINE BMG2_SymStd_SETUP_PtrLS( 
     &                       pLS, pSI_MSG, NProcI, NProcJ, NOG
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Setup pointers into the line-solve workspace.
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

      IMPLICIT   NONE

C -----------------------------
C     Includes
C
      INCLUDE    'BMG_workspace.h'
#include      "BMG_constants.h"

C ---------------------------
C     Argument Declarations:
C
      INTEGER  NOG, NProcI, NProcJ, pSI_MSG
      INTEGER  pLS(NBMG_pLS,NOG)

C ----------------------------
C     Local Declarations:
C
      INTEGER kg, ptr

C ==========================================================================

      ptr = pSI_MSG

      DO kg=NOG, 1, -1

         pLS(ipL_LS_XDataDist,kg) = ptr
         ptr = ptr+2*NProcI
         
         pLS(ipL_LS_YDataDist,kg) = ptr
         ptr = ptr+2*NProcJ
         
      ENDDO
      
      pSI_MSG = ptr

C ========================================================================

      RETURN
      END


