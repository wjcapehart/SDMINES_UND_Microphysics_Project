      REAL*RKIND FUNCTION BMG3_SymStd_UTILS_reswt(
     &                     CI, R, IIC, JJC, KKC, ICL, JCL, KCL, 
     &                     IL, JL, KL 
     &                     )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_UTILS_reswt local residual weighting routing used
C     by BMG3_SymStd_ITLI_bl.
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
      
#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'

C ----------------------------
C     Argument Declarations
C
      INTEGER ICL, IIC, IL, JCL, JJC, JL, KCL, KKC, KL
      REAL*RKIND  CI(IIC,JJC,KKC,26), R(7,7,7)

C ----------------------------
C     Local Variables
C     
      INTEGER ICP, JCP, KCP

C ==========================================================================

      ICP = MIN(ICL+1,IIC)
      JCP = MIN(JCL+1,JJC)
      KCP = MIN(KCL+1,KKC)
  
      BMG3_SymStd_UTILS_reswt = 
     &       CI(ICL,JCL,KCL,LXYNE)*R(IL-1,JL-1,KL)
     &     + CI(ICL,JCL,KCL,LXYA)*R(IL,JL-1,KL)
     &     + CI(ICP,JCL,KCL,LXYNW)*R(IL+1,JL-1,KL)
     &     + CI(ICL,JCL,KCL,LXYR)*R(IL-1,JL,KL)
     &     + R(IL,JL,KL)
     &     + CI(ICP,JCL,KCL,LXYL)*R(IL+1,JL,KL)
     &     + CI(ICL,JCP,KCL,LXYSE)*R(IL-1,JL+1,KL)
     &     + CI(ICL,JCP,KCL,LXYB)*R(IL,JL+1,KL)
     &     + CI(ICP,JCP,KCL,LXYSW)*R(IL+1,JL+1,KL)
     &     + CI(ICL,JCL,KCL,LTNE)*R(IL-1,JL-1,KL-1)
     &     + CI(ICL,JCL,KCL,LYZNW)*R(IL,JL-1,KL-1)
     &     + CI(ICP,JCL,KCL,LTNW)*R(IL+1,JL-1,KL-1)
     &     + CI(ICL,JCL,KCL,LXZNE)*R(IL-1,JL,KL-1)
     &     + CI(ICL,JCL,KCL,LXZA)*R(IL,JL,KL-1)
     &     + CI(ICP,JCL,KCL,LXZNW)*R(IL+1,JL,KL-1)
     &     + CI(ICL,JCP,KCL,LTSE)*R(IL-1,JL+1,KL-1)
     &     + CI(ICL,JCP,KCL,LYZNE)*R(IL,JL+1,KL-1)
     &     + CI(ICP,JCP,KCL,LTSW)*R(IL+1,JL+1,KL-1)
     &     + CI(ICL,JCL,KCP,LBNE)*R(IL-1,JL-1,KL+1)
     &     + CI(ICL,JCL,KCP,LYZSW)*R(IL,JL-1,KL+1)
     &     + CI(ICP,JCL,KCP,LBNW)*R(IL+1,JL-1,KL+1)
     &     + CI(ICL,JCL,KCP,LXZSE)*R(IL-1,JL,KL+1)
     &     + CI(ICL,JCL,KCP,LXZB)*R(IL,JL,KL+1)
     &     + CI(ICP,JCL,KCP,LXZSW)*R(IL+1,JL,KL+1)
     &     + CI(ICL,JCP,KCP,LBSE)*R(IL-1,JL+1,KL+1)
     &     + CI(ICL,JCP,KCP,LYZSE)*R(IL,JL+1,KL+1)
     &     + CI(ICP,JCP,KCP,LBSW)*R(IL+1,JL+1,KL+1)

C ==========================================================================

      RETURN
      END
