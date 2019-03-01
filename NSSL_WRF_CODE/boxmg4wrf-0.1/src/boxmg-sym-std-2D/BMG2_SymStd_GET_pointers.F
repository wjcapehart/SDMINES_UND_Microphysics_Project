      SUBROUTINE BMG2_SymStd_GET_pointers( 
     &                       k, IGRD, NOG,
     &                       p_U, p_SO, p_SOR, p_CI, NLx, NLy, NGx, NGy 
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_GET_pointers retrieves pointers to data at grid
C     level k within the internal pointer array IGRD.
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

C ---------------------------
C    Includes:
C
      INCLUDE 'BMG_workspace.h'

C ---------------------------
C    Argument Declarations:
C     
      INTEGER k, NLx, NLy, NGx, NGy, NOG
      INTEGER IGRD(NOG,NBMG_pIGRD), p_CI, p_SO, p_SOR, p_U

C ==========================================================================

      NLx = IGRD(k,idL_BMG_NLx)
      NLy = IGRD(k,idL_BMG_NLy)

      NGx = IGRD(k,idL_BMG_NGx)
      NGy = IGRD(k,idL_BMG_NGy)
      

      p_U   = IGRD(k,ipL_BMG_U)
      p_SO  = IGRD(k,ipL_BMG_SO)
      p_SOR = IGRD(k,ipL_BMG_SOR)
      p_CI  = IGRD(k,ipL_BMG_CI)

C ==========================================================================

      RETURN
      END
