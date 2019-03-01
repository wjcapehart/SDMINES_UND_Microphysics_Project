      REAL*RKIND FUNCTION BMG2_SymStd_UTILS_reswt(
     &                     CI, R, IIC, JJC, ICL, JCL, IL, JL, IND 
     &                     )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_UTILS_reswt.f local residual weighting routing used
C     by BMG2_SymStd_cg_ITLI_bl.f.
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
C   OUTPUT:
C  --------------------
C
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

C ----------------------------
C     Includes
C

      INCLUDE 'mpif.h'
      INCLUDE 'MSG.h'

#include      "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'

C ----------------------------
C     Argument Declarations
C
      INTEGER ICL, IIC, IL, IND, JCL, JJC, JL
      REAL*RKIND  CI(IIC,JJC,8), R(7,7,2)

C ----------------------------
C     Local Variables
C     
      INTEGER ICP, JCP

C ==========================================================================

      ICP = MIN(ICL+1,IIC)
      JCP = MIN(JCL+1,JJC)
  
      BMG2_SymStd_UTILS_reswt = CI(ICL,JCL,LNE)*R(IL-1,JL-1,IND)
     &                        + CI(ICL,JCL,LA)*R(IL,JL-1,IND)
     &                        + CI(ICP,JCL,LNW)*R(IL+1,JL-1,IND)
     &                        + CI(ICL,JCL,LR)*R(IL-1,JL,IND)
     &                        + R(IL,JL,IND)
     &                        + CI(ICP,JCL,LL)*R(IL+1,JL,IND)
     &                        + CI(ICL,JCP,LSE)*R(IL-1,JL+1,IND)
     &                        + CI(ICL,JCP,LB)*R(IL,JL+1,IND)
     &                        + CI(ICP,JCP,LSW)*R(IL+1,JL+1,IND)

C ==========================================================================

      RETURN
      END
