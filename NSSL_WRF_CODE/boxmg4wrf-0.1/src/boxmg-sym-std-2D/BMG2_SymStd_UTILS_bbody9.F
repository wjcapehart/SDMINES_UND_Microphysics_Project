      SUBROUTINE BMG2_SymStd_UTILS_bbody9(
     &                IIF, JJF, IIC, JJC, IC, JC, I, J, 
     &                NStncl, IND, CI, SO, Q, R
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Used by SETUP_ITLI_bl to preform the body of the local 
C     computations in the case of 9-point fine-grid stencil.
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
      INTEGER I, IC, IIC, IIF, IND, J, JC, JJC, JJF, NStncl
      REAL*RKIND CI(IIC,JJC,8), SO(IIF+1,JJF+1,Nstncl),Q(7,7,2),R(7,7,2)

C ----------------------------
C     Local Variables
C
      INTEGER IL, ILG, ILGP, JL, JLG, JLGP

C ==========================================================================

      DO JL = 1,7
         DO IL = 1,7
            Q(IL,JL,IND) = rZERO
            R(IL,JL,IND) = rZERO
         ENDDO
      ENDDO 
      
      Q(4,4,IND) = rONE
      Q(3,4,IND) = CI(IC,JC,LR)
      Q(5,4,IND) = CI(IC+1,JC,LL)
      Q(4,3,IND) = CI(IC,JC,LA)
      Q(4,5,IND) = CI(IC,JC+1,LB)
      Q(3,3,IND) = CI(IC,JC,LNE)
      Q(3,5,IND) = CI(IC,JC+1,LSE)
      Q(5,3,IND) = CI(IC+1,JC,LNW)
      Q(5,5,IND) = CI(IC+1,JC+1,LSW)
      
      !
      ! Apply fine-grid operator
      !
      DO JL = 2,6
         DO IL = 2,6

            ILG = I + IL - 4
            JLG = J + JL - 4
            
            ILGP = MIN(ILG+1,IIF+1)
            JLGP = MIN(JLG+1,JJF+1)
            
            ILGP = MAX(ILGP,1)
            JLGP = MAX(JLGP,1)
            
            ILG = MIN(ILG,IIF+1)
            JLG = MIN(JLG,JJF+1)
            
            ILG = MAX(ILG,1)
            JLG = MAX(JLG,1)
            
            R(IL,JL,IND) = SO(ILG,JLG,KO)*Q(IL,JL,IND)
     &                   - SO(ILG,JLG,KW)*Q(IL-1,JL,IND)
     &                   - SO(ILG,JLGP,KS)*Q(IL,JL+1,IND)
     &                   - SO(ILGP,JLG,KW)*Q(IL+1,JL,IND)
     &                   - SO(ILG,JLG,KS)*Q(IL,JL-1,IND)
     &                   - SO(ILG,JLGP,KNW)*Q(IL-1,JL+1,IND)
     &                   - SO(ILGP,JLGP,KSW)*Q(IL+1,JL+1,IND)
     &                   - SO(ILGP,JLG,KNW)*Q(IL+1,JL-1,IND)
     &                   - SO(ILG,JLG,KSW)*Q(IL-1,JL-1,IND)


         ENDDO   
      ENDDO   

C ==========================================================================

      RETURN
      END
