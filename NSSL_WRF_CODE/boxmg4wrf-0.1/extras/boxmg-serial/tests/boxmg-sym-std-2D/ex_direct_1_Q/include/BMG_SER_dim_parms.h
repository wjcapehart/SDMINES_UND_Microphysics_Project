C ------------------------------------------------
C     Multigrid/Workspace Memory Allocation: 
C ------------------------------------------------
 
      INTEGER   Lmax, Mmax, NXYCmax, NXYCmin
      PARAMETER ( Lmax=30, Mmax=60, NXYCmax=3, NXYCmin=3 )

      INTEGER   NFm, NOGm, NSOm
      PARAMETER ( NFm=51533, NSOm=292712, NOGm=5 )

      INTEGER   NBMG_iWORK, NBMG_rWORK
      PARAMETER (  NBMG_iWORK=58900, NBMG_rWORK=1616441 )

