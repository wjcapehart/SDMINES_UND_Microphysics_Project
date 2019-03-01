      SUBROUTINE BMG2_SER_SymStd_GET_pointers(
     &                k, IGRD, NOG, Nx, Ny, p_U, p_SO, p_SOR, p_CI
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_SymStd_GET_pointers.f retrieves pointers to data at grid
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
      INCLUDE 'BMG_SER_workspace.h'

C ---------------------------
C    Argument Declarations:
C     
      INTEGER k, Nx, Ny, NOG

      INTEGER IGRD(NOG,9), p_CI, p_SO, p_SOR, p_U

C ===========================================================================

      Nx = IGRD(k,idL_BMG_SER_Nx)
      Ny = IGRD(k,idL_BMG_SER_Ny)

      p_U   = IGRD(k,ipL_BMG_SER_U)
      p_SO  = IGRD(k,ipL_BMG_SER_SO)
      p_SOR = IGRD(k,ipL_BMG_SER_SOR)
      p_CI  = IGRD(k,ipL_BMG_SER_CI)

C ===========================================================================

      RETURN
      END
