      SUBROUTINE BMG3_SER_SymStd_GET_pointers( 
     &                k, IGRD, NOG,
     &                p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_GET_pointers is a subroutine provided to access the
C     data structures used in BMG3D.  It returns pointers for grid
C     K, as well as the dimensions of grid K.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     k       Grid index.
C
C     NOG     The number of grids and leading dimension of the
C             workspace array IGRD.
C
C     IGRD    An integer array containing the pointers for different 
C             grids.
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
C     Nx      The number of grid points in the x direction for grid k, 
C             including two fictitious points.
C     Ny      The number of grid points in the y directiom for grid K, 
C             including two fictitious points.
C     Nz      The number of grid points in the z direction for grid K,
C             including two fictitious points
C
C     p_CI    Pointer to the interpolation on grid k, in CI.
C     p_SO    Pointer to the stencil on grid k, in SO.
C     p_SOR   Pointer to the stencil workspace on grid k, in SOR.
C     p_U     Pointer to vectors on grid k, in Q and QF
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ===========================================================================

      IMPLICIT NONE

C ---------------------------
C    Includes:
C
      INCLUDE 'BMG_SER_workspace.h'

C ---------------------------
C    Argument Declarations:
C     
      INTEGER k, NOG, Nx, Ny, Nz

      INTEGER IGRD(NOG,NBMG_SER_pIGRD), p_CI, p_SO, p_SOR, p_U

C ===========================================================================

      Nx = igrd(k,idL_BMG_SER_Nx)
      Ny = igrd(k,idL_BMG_SER_Ny)
      Nz = igrd(k,idL_BMG_SER_Nz)

      p_U   = IGRD(k,ipL_BMG_SER_U)
      p_SO  = IGRD(k,ipL_BMG_SER_SO)
      p_SOR = IGRD(k,ipL_BMG_SER_SOR)
      p_CI  = IGRD(k,ipL_BMG_SER_CI)

C ===========================================================================

      RETURN
      END

