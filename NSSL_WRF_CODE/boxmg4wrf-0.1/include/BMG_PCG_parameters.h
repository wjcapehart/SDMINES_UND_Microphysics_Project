! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     This file provides a single resource for defining commonly used
!     parameters in the pcg code used with the BOXMG family of codes.
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   VARIABLES:
!  --------------------
!
!
! ==========================================================================
! ----------------------------------------------------
!     Parameter Indexing
! ---------------------------------

! ---------------------------------
!     INTEGER Parameters
! ---------------------------------

      INTEGER   id_BMG_PCG_NMG_CYCLES,                                  &
     &          id_BMG_PCG_STOP_TEST,                                   &
     &          id_BMG_PCG_PRECON,                                      &
     &          id_BMG_PCG_BMG_SETUP,                                   &
     &          id_BMG_PCG_MAX_ITERS,                                   &
     &          id_BMG_PCG_OUTPUT,                                      &
     &          id_BMG_PCG_OUT_FREQ
     

      PARAMETER ( id_BMG_PCG_NMG_CYCLES   = 1,                          &
     &            id_BMG_PCG_STOP_TEST    = 2,                          &
     &            id_BMG_PCG_PRECON       = 3,                          &
     &            id_BMG_PCG_BMG_SETUP    = 4,                          &
     &            id_BMG_PCG_MAX_ITERS    = 5,                          &
     &            id_BMG_PCG_OUTPUT       = 6,                          &
     &            id_BMG_PCG_OUT_FREQ     = 7  )

      INTEGER   NBMG_PCG_iPARMS 
      PARAMETER ( NBMG_PCG_iPARMS = 7 )


! -------------------------------
!     REAL Parameters
! -------------------------------

      INTEGER   id_BMG_PCG_STOP_TOL
      PARAMETER ( id_BMG_PCG_STOP_TOL = 1 )

      INTEGER   NBMG_PCG_rPARMS 
      PARAMETER ( NBMG_PCG_rPARMS = 2 )

! --------------------------------
!     Stopping Critieria
! --------------------------------

      INTEGER BMG_PCG_STOP_ABS_RES_L2,                                  &
     &        BMG_PCG_STOP_REL_RES_L2,                                  &
     &        BMG_PCG_STOP_ABS_RES_M2,                                  &
     &        BMG_PCG_STOP_REL_RES_M2

      PARAMETER( BMG_PCG_STOP_ABS_RES_L2 = 0,                           &
     &           BMG_PCG_STOP_REL_RES_L2 = 1,                           &
     &           BMG_PCG_STOP_ABS_RES_M2 = 2,                           &
     &           BMG_PCG_STOP_REL_RES_M2 = 3   ) 

! --------------------------------
!     Preconditioners
! --------------------------------

      INTEGER BMG_PCG_PRECON_NONE,                                      &
     &        BMG_PCG_PRECON_DIAG,                                      &
     &        BMG_PCG_PRECON_BMG

      PARAMETER( BMG_PCG_PRECON_NONE  = 1,                              &
     &           BMG_PCG_PRECON_DIAG  = 2,                              &
     &           BMG_PCG_PRECON_BMG   = 3  ) 

! -------------------------------
!     BMG Setup options:
! -------------------------------

      INTEGER BMG_PCG_BMG_SETUP_all,                                    &
     &        BMG_PCG_BMG_SETUP_none

      PARAMETER ( BMG_PCG_BMG_SETUP_all = 0,                            &
     &            BMG_PCG_BMG_SETUP_none = 1 )

! =======================================================================
! ----------------------------------------------------
!     IOFLAG Indexing
! ---------------------------------

      INTEGER  BMG_PCG_OUT_NONE,                                        &
     &         BMG_PCG_OUT_INIT_RES,                                    &
     &         BMG_PCG_OUT_FIN_RES,                                     &
     &         BMG_PCG_OUT_INIT_FIN_RES,                                &
     &         BMG_PCG_OUT_ITS,                                         &
     &         BMG_PCG_OUT_ALL
  
      PARAMETER ( BMG_PCG_OUT_NONE         = 1,                         &
     &            BMG_PCG_OUT_INIT_RES     = 2,                         &
     &            BMG_PCG_OUT_FIN_RES      = 3,                         &
     &            BMG_PCG_OUT_INIT_FIN_RES = 4,                         &
     &            BMG_PCG_OUT_ITS          = 5,                         &
     &            BMG_PCG_OUT_ALL          = 6 )





