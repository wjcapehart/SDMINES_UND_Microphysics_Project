C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Parameter declarations for multigrid/workspace memory allocation.
C     
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   VARIABLES:
C  --------------------
C
C
C =======================================================================
 
      INTEGER   Kmax, Lmax, Mmax, NXYCmax, NXYCmin
      PARAMETER ( Kmax=30, Lmax=30, Mmax=60, NXYCmax=3, NXYCmin=3 )

      INTEGER   NFm, NOGm, NSOm
      PARAMETER ( NFm=60452, NSOm=339798, NOGm=5 )

      INTEGER   NBMG_iWORK, NBMG_rWORK
      PARAMETER (  NBMG_iWORK=40, NBMG_rWORK=436535 )

      INTEGER   NBMG_iWORK_PL, NBMG_rWORK_PL
      PARAMETER ( NBMG_iWORK_PL=15799, NBMG_rWORK_PL=3270324 )

C =======================================================================

