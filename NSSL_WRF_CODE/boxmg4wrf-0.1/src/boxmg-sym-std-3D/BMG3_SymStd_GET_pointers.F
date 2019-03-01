      SUBROUTINE BMG3_SymStd_GET_pointers( 
     &                       k, IGRD, NOGm,
     &                       p_U, p_SO, p_SOR, p_CI, 
     &                       NLx, NLy, NLz, NGx, NGy, NGz 
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_GET_pointers retrieves pointers to data at grid
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
      INTEGER k, NOGm, NGx, NGy, NGz,  NLx, NLy, NLz
      INTEGER IGRD(NOGm,NBMG_pIGRD), p_CI, p_SO, p_SOR, p_U

C ==========================================================================

      NLx = igrd(k,idL_BMG_NLx)
      NLy = igrd(k,idL_BMG_NLy)
      NLz = igrd(k,idL_BMG_NLz)

      NGx = igrd(k,idL_BMG_NGx)
      NGy = igrd(k,idL_BMG_NGy)
      NGz = igrd(k,idL_BMG_NGz)

      p_U   = IGRD(k,ipL_BMG_U)
      p_SO  = IGRD(k,ipL_BMG_SO)
      p_SOR = IGRD(k,ipL_BMG_SOR)
      p_CI  = IGRD(k,ipL_BMG_CI)

C ==========================================================================

      RETURN
      END

