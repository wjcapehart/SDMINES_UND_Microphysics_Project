C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This file provides a single resource for defining commonly used
C     parameters in the pcg code used with the BOXMG family of codes.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   VARIABLES:
C  --------------------
C
C
C ==========================================================================

C ==========================================================================
C ----------------------------------------------------
C     Parameter Indexing
C ---------------------------------

C ---------------------------------
C     INTEGER Parameters
C ---------------------------------

      INTEGER   id_BMG_SER_PCG_NMG_CYCLES,
     &          id_BMG_SER_PCG_STOP_TEST,
     &          id_BMG_SER_PCG_PRECON,
     &          id_BMG_SER_BMG_iPARMS0_SETUP,
     &          id_BMG_SER_PCG_MAX_ITERS

      PARAMETER ( id_BMG_SER_PCG_NMG_CYCLES   = 1,
     &            id_BMG_SER_PCG_STOP_TEST    = 2,
     &            id_BMG_SER_PCG_PRECON       = 3,
     &            id_BMG_SER_BMG_iPARMS0_SETUP    = 4,
     &            id_BMG_SER_PCG_MAX_ITERS    = 5  )

      INTEGER   NBMG_SER_SER_PCG_iPARMS 
      PARAMETER ( NBMG_SER_SER_PCG_iPARMS = 5 )


C -------------------------------
C     REAL Parameters
C -------------------------------

      INTEGER   id_BMG_SER_PCG_STOP_TOL
      PARAMETER ( id_BMG_SER_PCG_STOP_TOL = 1 )

      INTEGER   NBMG_SER_SER_PCG_rPARMS 
      PARAMETER ( NBMG_SER_SER_PCG_rPARMS = 2 )

C --------------------------------
C     Stopping Critieria
C --------------------------------

      INTEGER BMG_SER_PCG_STOP_ABS_RES_L2,
     &        BMG_SER_PCG_STOP_REL_RES_L2,
     &        BMG_SER_PCG_STOP_ABS_RES_M2,
     &        BMG_SER_PCG_STOP_REL_RES_M2

      PARAMETER( BMG_SER_PCG_STOP_ABS_RES_L2 = 0,
     &           BMG_SER_PCG_STOP_REL_RES_L2 = 1,
     &           BMG_SER_PCG_STOP_ABS_RES_M2 = 2,
     &           BMG_SER_PCG_STOP_REL_RES_M2 = 3   ) 

C --------------------------------
C     Preconditioners
C --------------------------------

      INTEGER BMG_SER_PCG_PRECON_NONE,
     &        BMG_SER_PCG_PRECON_DIAG,
     &        BMG_SER_PCG_PRECON_BMG

      PARAMETER( BMG_SER_PCG_PRECON_NONE  = 1,
     &           BMG_SER_PCG_PRECON_DIAG  = 2,
     &           BMG_SER_PCG_PRECON_BMG   = 3  ) 

C -------------------------------
C     BMG Setup options:
C -------------------------------

      INTEGER BMG_SER_BMG_iPARMS0_SETUP_all,
     &        BMG_SER_BMG_iPARMS0_SETUP_none

      PARAMETER ( BMG_SER_BMG_iPARMS0_SETUP_all = 0,
     &            BMG_SER_BMG_iPARMS0_SETUP_none = 1 )

C ==========================================================================

