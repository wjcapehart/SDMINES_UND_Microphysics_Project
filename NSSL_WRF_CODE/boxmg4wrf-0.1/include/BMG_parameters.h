! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     This include file provides a single resource for defining commonly
!     used parameters in the BOXMG family of codes.
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

      INTEGER   id_BMG2_DIM_NOG,                                        &
     &          id_BMG2_DIM_NF,                                         &
     &          id_BMG2_DIM_NC,                                         &
     &          id_BMG2_DIM_NSO,                                        &
     &          id_BMG2_DIM_NSOR,                                       &
     &          id_BMG2_DIM_NCI,                                        &
     &          id_BMG2_DIM_NCBW,                                       &
     &          id_BMG2_DIM_NCU,                                        & 
     &          id_BMG2_DIM_NMSGi,                                      & 
     &          id_BMG2_DIM_NMSGr,                                      & 
     &          id_BMG2_POINTERS,                                       &
     &          id_BMG2_STENCIL,                                        &
     &          id_BMG2_BC,                                             &
     &          id_BMG2_SETUP,                                          &
     &          id_BMG2_RELAX,                                          &
     &          id_BMG2_RELAX_SYM,                                      &
     &          id_BMG2_NRELAX_DOWN,                                    &
     &          id_BMG2_NRELAX_UP,                                      &
     &          id_BMG2_NRELAX_FG,                                      &
     &          id_BMG2_CYCLE_CLASS,                                    &
     &          id_BMG2_NCYCLE_TYPE,                                    &
     &          id_BMG2_FMG_NNCYCLE,                                    &
     &          id_BMG2_MAX_ITERS,                                      &
     &          id_BMG2_STOP_TEST,                                      &
     &          id_BMG2_MIN_NOG,                                        &
     &          id_BMG2_CG_MIN_DIM,                                     &
     &          id_BMG2_CG_TYPE,                                        &
     &          id_BMG2_CG_CONSTRUCT,                                   &
     &          id_BMG2_CG_COMM,                                        &
     &          id_BMG2_CG_SOLVER,                                      &
     &          id_BMG2_Err_Code,                                       &
     &          id_BMG2_Ext_Err_Code,                                   &
     &          id_BMG2_SYNC_INITIAL_GUESS,                             &
     &          id_BMG2_LINE_SOLVE_COMM_TYPE,                           &
     &          id_BMG3_DIM_NOG,                                        &
     &          id_BMG3_DIM_NF,                                         &
     &          id_BMG3_DIM_NC,                                         &
     &          id_BMG3_DIM_NSO,                                        &
     &          id_BMG3_DIM_NSOR,                                       &
     &          id_BMG3_DIM_NCI,                                        &
     &          id_BMG3_DIM_NCBW,                                       &
     &          id_BMG3_DIM_NCU,                                        &
     &          id_BMG3_DIM_NMSGi,                                      & 
     &          id_BMG3_DIM_NMSGr,                                      & 
     &          id_BMG3_DIM_NSO3,                                       & 
     &          id_BMG3_DIM_NSR3,                                       & 
     &          id_BMG3_DIM_NCI3,                                       & 
     &          id_BMG3_DIM_NCPL,                                       & 
     &          id_BMG3_POINTERS,                                       &
     &          id_BMG3_STENCIL,                                        &
     &          id_BMG3_BC,                                             &
     &          id_BMG3_SETUP,                                          &
     &          id_BMG3_RELAX,                                          &
     &          id_BMG3_RELAX_SYM,                                      &
     &          id_BMG3_NRELAX_DOWN,                                    &
     &          id_BMG3_NRELAX_UP,                                      &
     &          id_BMG3_NRELAX_FG,                                      &
     &          id_BMG3_CYCLE_CLASS,                                    &
     &          id_BMG3_NCYCLE_TYPE,                                    &
     &          id_BMG3_FMG_NNCYCLE,                                    &
     &          id_BMG3_MAX_ITERS,                                      &
     &          id_BMG3_STOP_TEST,                                      &
     &          id_BMG3_MIN_NOG,                                        &
     &          id_BMG3_CG_MIN_DIM,                                     &
     &          id_BMG3_CG_TYPE,                                        &
     &          id_BMG3_CG_CONSTRUCT,                                   &
     &          id_BMG3_CG_COMM,                                        &
     &          id_BMG3_CG_SOLVER,                                      &
     &          id_BMG3_Err_Code,                                       &
     &          id_BMG3_Ext_Err_Code,                                   &
     &          id_BMG3_SYNC_INITIAL_GUESS,                             &
     &          NBMG2_iPARMS, NBMG3_iPARMS, NBMG_iPARMS,                &
     &          id_BMG3_SEPARABLE


      PARAMETER ( id_BMG2_DIM_NOG      =  1,                            &
     &            id_BMG2_DIM_NF       =  2,                            &
     &            id_BMG2_DIM_NC       =  3,                            &
     &            id_BMG2_DIM_NSO      =  4,                            &
     &            id_BMG2_DIM_NSOR     =  5,                            &
     &            id_BMG2_DIM_NCI      =  6,                            &
     &            id_BMG2_DIM_NCBW     =  7,                            &
     &            id_BMG2_DIM_NCU      =  8,                            & 
     &            id_BMG2_DIM_NMSGi    =  9,                            & 
     &            id_BMG2_DIM_NMSGr    =  10,                           & 
     &            id_BMG2_POINTERS     =  11,                           &
     &            id_BMG2_STENCIL      =  12,                           &
     &            id_BMG2_BC           =  13,                           & 
     &            id_BMG2_SETUP        =  14,                           &
     &            id_BMG2_RELAX        =  15,                           &
     &            id_BMG2_RELAX_SYM    =  16,                           &
     &            id_BMG2_NRELAX_DOWN  =  17,                           &
     &            id_BMG2_NRELAX_UP    =  18,                           &
     &            id_BMG2_NRELAX_FG    =  19,                           &
     &            id_BMG2_CYCLE_CLASS  =  20,                           &
     &            id_BMG2_NCYCLE_TYPE  =  21,                           &
     &            id_BMG2_FMG_NNCYCLE  =  22,                           &
     &            id_BMG2_MAX_ITERS    =  23,                           &
     &            id_BMG2_STOP_TEST    =  24,                           &    
     &            id_BMG2_MIN_NOG      =  25,                           &
     &            id_BMG2_CG_MIN_DIM   =  26,                           &
     &            id_BMG2_CG_TYPE      =  27,                           &
     &            id_BMG2_CG_CONSTRUCT =  28,                           &
     &            id_BMG2_CG_COMM      =  29,                           &
     &            id_BMG2_CG_SOLVER    =  30,                           &
     &            id_BMG2_Err_Code     =  31,                           &
     &            id_BMG2_Ext_Err_Code =  32,                           &
     &            id_BMG2_SYNC_INITIAL_GUESS = 33,                      &
     &            id_BMG2_LINE_SOLVE_COMM_TYPE = 34,                    &
     &            id_BMG3_DIM_NOG      =  35,                           &
     &            id_BMG3_DIM_NF       =  36,                           &
     &            id_BMG3_DIM_NC       =  37,                           &
     &            id_BMG3_DIM_NSO      =  38,                           &
     &            id_BMG3_DIM_NSOR     =  39,                           &
     &            id_BMG3_DIM_NCI      =  40,                           &
     &            id_BMG3_DIM_NCBW     =  41,                           &
     &            id_BMG3_DIM_NCU      =  42,                           & 
     &            id_BMG3_DIM_NMSGi    =  43,                           & 
     &            id_BMG3_DIM_NMSGr    =  44,                           & 
     &            id_BMG3_DIM_NSO3     =  45,                           & 
     &            id_BMG3_DIM_NSR3     =  46,                           & 
     &            id_BMG3_DIM_NCI3     =  47,                           & 
     &            id_BMG3_DIM_NCPL     =  48,                           & 
     &            id_BMG3_POINTERS     =  49,                           &
     &            id_BMG3_STENCIL      =  50,                           &
     &            id_BMG3_BC           =  51,                           &
     &            id_BMG3_SETUP        =  52,                           &
     &            id_BMG3_RELAX        =  53,                           &
     &            id_BMG3_RELAX_SYM    =  54,                           &
     &            id_BMG3_NRELAX_DOWN  =  55,                           &
     &            id_BMG3_NRELAX_UP    =  56,                           &
     &            id_BMG3_NRELAX_FG    =  57,                           &
     &            id_BMG3_CYCLE_CLASS  =  58,                           &
     &            id_BMG3_NCYCLE_TYPE  =  59,                           &
     &            id_BMG3_FMG_NNCYCLE  =  60,                           &
     &            id_BMG3_MAX_ITERS    =  61,                           &
     &            id_BMG3_STOP_TEST    =  62,                           &
     &            id_BMG3_MIN_NOG      =  63,                           &
     &            id_BMG3_CG_MIN_DIM   =  64,                           &
     &            id_BMG3_CG_TYPE      =  65,                           &
     &            id_BMG3_CG_CONSTRUCT =  66,                           &
     &            id_BMG3_CG_COMM      =  67,                           &
     &            id_BMG3_CG_SOLVER    =  68,                           &
     &            id_BMG3_Err_Code     =  69,                           &
     &            id_BMG3_Ext_Err_Code =  70,                           &
     &            id_BMG3_SYNC_INITIAL_GUESS = 71,                      &
     &            NBMG2_iPARMS = 34,                                    &
!     &            NBMG3_iPARMS = 37,                                    &
!     &            NBMG_iPARMS  = 71                                     &
!     &            )
     &            NBMG3_iPARMS = 38,                                    &
     &            id_BMG3_SEPARABLE = 72,                               &
     &            NBMG_iPARMS  = 72                                     &
     &            )

! -------------------------------
!     REAL Parameters
! -------------------------------

      INTEGER   id_BMG2_STOP_TOL,                                       &
     &          id_BMG2_TIME_SETUP_FINE_STENCIL,                        &
     &          id_BMG2_TIME_SETUP_CG_ITLI,                             &
     &          id_BMG2_TIME_SETUP_INTERP_OI,                           &
     &          id_BMG2_TIME_SETUP_RELAX,                               &
     &          id_BMG2_TIME_SETUP_CG_LU,                               &
     &          id_BMG2_TIME_SOLVE_CG,                                  &
     &          id_BMG2_TIME_relax,                                     &
     &          id_BMG2_TIME_restrict,                                  &
     &          id_BMG2_TIME_interp_add,                                &
     &          id_BMG2_TIME_SETUP_MSG,                                 &
     &          id_BMG2_TIME_SETUP_LS      ,                            &
     &          id_BMG2_TIME_SETUP_PTR_GRID,                            &
     &          id_BMG2_TIME_SETUP_PARTS,                               &
     &          id_BMG2_TIME_SETUP_TOTAL,                               &
     &          id_BMG2_TIME_SOLVE_TOTAL,                               &
     &          id_BMG2_TIME_PCG_TOTAL,                                 &
     &          id_BMG2_TIME_PCG_PRECON,                                &  
     &          id_BMG3_STOP_TOL,                                       &
     &          id_BMG3_TIME_SETUP_FINE_STENCIL,                        &
     &          id_BMG3_TIME_SETUP_CG_ITLI,                             &
     &          id_BMG3_TIME_SETUP_INTERP_OI,                           &
     &          id_BMG3_TIME_SETUP_RELAX,                               &
     &          id_BMG3_TIME_SETUP_CG_LU,                               &
     &          id_BMG3_TIME_SOLVE_CG,                                  &
     &          id_BMG3_TIME_relax,                                     &
     &          id_BMG3_TIME_restrict,                                  &
     &          id_BMG3_TIME_interp_add,                                &
     &          id_BMG3_TIME_SETUP_MSG,                                 &
     &          id_BMG3_TIME_SETUP_PTR_GRID,                            &
     &          id_BMG3_TIME_SETUP_PARTS,                               &
     &          id_BMG3_TIME_SETUP_TOTAL,                               &
     &          id_BMG3_TIME_SOLVE_TOTAL,                               &
     &          id_BMG3_TIME_PCG_TOTAL,                                 &
     &          id_BMG3_TIME_PCG_PRECON,                                &
     &          NBMG2_rPARMS,                                           &
     &          NBMG3_rPARMS,                                           &
     &          NBMG_rPARMS 

      PARAMETER ( id_BMG2_STOP_TOL                = 1,                  &
     &            id_BMG2_TIME_SETUP_FINE_STENCIL = 2,                  &
     &            id_BMG2_TIME_SETUP_CG_ITLI      = 3,                  &
     &            id_BMG2_TIME_SETUP_INTERP_OI    = 4,                  &
     &            id_BMG2_TIME_SETUP_RELAX        = 5,                  &
     &            id_BMG2_TIME_SETUP_CG_LU        = 6,                  &
     &            id_BMG2_TIME_SOLVE_CG           = 7,                  &
     &            id_BMG2_TIME_relax              = 8,                  &
     &            id_BMG2_TIME_restrict           = 9,                  &
     &            id_BMG2_TIME_interp_add         = 10,                 &
     &            id_BMG2_TIME_SETUP_MSG          = 11,                 &
     &            id_BMG2_TIME_SETUP_LS           = 12,                 &
     &            id_BMG2_TIME_SETUP_PTR_GRID     = 13,                 &
     &            id_BMG2_TIME_SETUP_PARTS        = 14,                 &
     &            id_BMG2_TIME_SETUP_TOTAL        = 15,                 &
     &            id_BMG2_TIME_SOLVE_TOTAL        = 16,                 &
     &            id_BMG2_TIME_PCG_TOTAL          = 17,                 &
     &            id_BMG2_TIME_PCG_PRECON         = 18,                 &  
     &            id_BMG3_STOP_TOL                = 19,                 &
     &            id_BMG3_TIME_SETUP_FINE_STENCIL = 20,                 &
     &            id_BMG3_TIME_SETUP_CG_ITLI      = 21,                 &
     &            id_BMG3_TIME_SETUP_INTERP_OI    = 22,                 &
     &            id_BMG3_TIME_SETUP_RELAX        = 23,                 &
     &            id_BMG3_TIME_SETUP_CG_LU        = 24,                 &
     &            id_BMG3_TIME_SOLVE_CG           = 25,                 &
     &            id_BMG3_TIME_relax              = 26,                 &
     &            id_BMG3_TIME_restrict           = 27,                 &
     &            id_BMG3_TIME_interp_add         = 28,                 &
     &            id_BMG3_TIME_SETUP_MSG          = 29,                 &
     &            id_BMG3_TIME_SETUP_PTR_GRID     = 30,                 &
     &            id_BMG3_TIME_SETUP_PARTS        = 31,                 &
     &            id_BMG3_TIME_SETUP_TOTAL        = 32,                 &
     &            id_BMG3_TIME_SOLVE_TOTAL        = 33,                 &
     &            id_BMG3_TIME_PCG_TOTAL          = 34,                 &
     &            id_BMG3_TIME_PCG_PRECON         = 35,                 &
     &            NBMG2_rPARMS = 18,                                    &
     &            NBMG3_rPARMS = 17,                                    &
     &            NBMG_rPARMS = 35       )

! ----------------------------------------------------
! ==========================================================================
! ----------------------------------------------------
!     Parameter Values
! -------------------------------

! -------------------------------
!     Stencil Type
! -------------------------------

      INTEGER BMG_STENCIL_5pt,                                          &
     &        BMG_STENCIL_9pt,                                          &
     &        BMG_STENCIL_7pt,                                          &
     &        BMG_STENCIL_27pt

      PARAMETER ( BMG_STENCIL_5pt  = 1,                                 &
     &            BMG_STENCIL_9pt  = 2,                                 &
     &            BMG_STENCIL_7pt  = 1,                                 &
     &            BMG_STENCIL_27pt = 2  )

! -------------------------------
!     Periodicity
! -------------------------------

      INTEGER  BMG_BCs_definite,                                        &
     &         BMG_BCs_def_per_x,                                       &
     &         BMG_BCs_def_per_y,                                       &
     &         BMG_BCs_def_per_xy,                                      &
     &         BMG_BCs_indef_per_x,                                     &
     &         BMG_BCs_indef_per_y,                                     &
     &         BMG_BCs_indef_per_xy,                                    &
     &         BMG_BCs_indef_nonper

      PARAMETER ( BMG_BCs_definite     =  0,                            &
     &            BMG_BCs_def_per_x    =  2,                            &
     &            BMG_BCs_def_per_y    =  1,                            &
     &            BMG_BCs_def_per_xy   =  3,                            &
     &            BMG_BCs_indef_per_x  = -2,                            &
     &            BMG_BCs_indef_per_y  = -1,                            &
     &            BMG_BCs_indef_per_xy = -3,                            &
     &            BMG_BCs_indef_nonper = -4  )

! -------------------------------
!     Setup options:
! -------------------------------

      INTEGER BMG_SETUP_only,                                           &
     &        BMG_SETUP_none,                                           &
     &        BMG_SETUP_opers,                                          &
     &        BMG_SETUP_ptrs,                                           &
     &        BMG_SETUP_ptrs_opers 

      PARAMETER ( BMG_SETUP_only  = 3,                                  &
     &            BMG_SETUP_none  = 2,                                  &
     &            BMG_SETUP_ptrs  = 4,                                  &
     &            BMG_SETUP_opers = 1,                                  &
     &            BMG_SETUP_ptrs_opers = 0 )

! -------------------------------
!     Memory Allocation:
! -------------------------------

      INTEGER BMG_USE_pointers,                                         &
     &        BMG_NO_pointers

      PARAMETER ( BMG_USE_pointers = 1,                                 &
     &            BMG_NO_pointers = 0   )

! -------------------------------
!     Relaxation:
! -------------------------------
      
      INTEGER  BMG_GS_RB_point,                                         &
     &         BMG_GS_RB_x_lines,                                       &
     &         BMG_GS_RB_y_lines,                                       &
     &         BMG_GS_RB_x_y_lines,                                     &
     &         BMG_GS_RB_planes_xy_yz_xz

      PARAMETER( BMG_GS_RB_point     = 1,                               &
     &           BMG_GS_RB_x_lines   = 2,                               &
     &           BMG_GS_RB_y_lines   = 3,                               &
     &           BMG_GS_RB_x_y_lines = 4,                               &
     &           BMG_GS_RB_planes_xy_yz_xz = 5 )


! --------------------------------
!     Symmetry of the MG n-cycle:
! --------------------------------

      INTEGER  BMG_RELAX_NONSYM,                                        &
     &         BMG_RELAX_SYM 
      PARAMETER( BMG_RELAX_NONSYM = 0,                                  &
     &           BMG_RELAX_SYM = 1    )

      INTEGER  BMG_DOWN,                                                &
     &         BMG_UP
      PARAMETER(  BMG_DOWN = 0,                                         &
     &            BMG_UP = 1   )

      INTEGER  BMG_NONSEPARABLE,                                        &
     &         BMG_SEPARABLE

      PARAMETER(  BMG_NONSEPARABLE = 0,                                 &
     &            BMG_SEPARABLE = 1  )

! --------------------------------
!     Stopping Critieria
! --------------------------------

      INTEGER BMG_STOP_ABS_RES_L2,                                      &
     &        BMG_STOP_REL_RES_L2

      PARAMETER( BMG_STOP_ABS_RES_L2 = 0,                               &
     &           BMG_STOP_REL_RES_L2 = 1 ) 

! --------------------------------
!     Cycle Class and Type
! --------------------------------
      
      INTEGER   BMG_FMG_CYCLE,                                          &
     &          BMG_N_CYCLE,                                            &
     &          BMG_V_CYCLE,                                            &
     &          BMG_W_CYCLE

      PARAMETER ( BMG_FMG_CYCLE = 0,                                    &
     &            BMG_N_CYCLE   = 1,                                    &
     &            BMG_V_CYCLE   = 1,                                    &
     &            BMG_W_CYCLE   = 2 )

! --------------------------------
!     Coarse-Grid Operator Type
! --------------------------------

      INTEGER   BMG_CG_ITLI_IzIyIx,                                     &
     &          BMG_CG_ITLI,                                            &
     &          BMG_CG_USER

      PARAMETER ( BMG_CG_ITLI_IzIyIx = 1,                               &
     &            BMG_CG_ITLI        = 2,                               &
     &            BMG_CG_USER        = 3 )

! --------------------------------
!     I^{T} L I Construction
! --------------------------------

      INTEGER   BMG_CG_CONS_explicit,                                   &
     &          BMG_CG_CONS_block

      PARAMETER ( BMG_CG_CONS_explicit = 1,                             &
     &            BMG_CG_CONS_block    = 2  )

! --------------------------------
!     Coarse Grid Commm Pattern
! --------------------------------

      INTEGER BMG_CG_ALLGATHER,                                         &
     &        BMG_CG_GATHER_SCATTER

      PARAMETER ( BMG_CG_ALLGATHER = 0,                                 &
     &            BMG_CG_GATHER_SCATTER = 1 )
     
! --------------------------------
!     Coarse Grid Solve Type
! --------------------------------

      INTEGER BMG_CG_SOLVE_LU,                                          &
     &        BMG_CG_SOLVE_BOXMG

      PARAMETER ( BMG_CG_SOLVE_LU = 0,                                  &
     &            BMG_CG_SOLVE_BOXMG = 1 )

! --------------------------------
!     Syncronize initial guess
! --------------------------------

      INTEGER BMG_SYNC_INITIAL_GUESS,                                   &
     &        BMG_SYNC_INITIAL_GUESS_none

      PARAMETER ( BMG_SYNC_INITIAL_GUESS = 0,                           &
     &            BMG_SYNC_INITIAL_GUESS_none = 1 )

! --------------------------------
!     Coarse Grid Solve Type
! --------------------------------

      INTEGER BMG_LINE_SOLVE_COMM_TRADITIONAL,                          &
     &        BMG_LINE_SOLVE_COMM_TUNED

      PARAMETER ( BMG_LINE_SOLVE_COMM_TRADITIONAL = 0,                  &
     &            BMG_LINE_SOLVE_COMM_TUNED = 1 )
! ----------------------------------------------------
! ==========================================================================
! ----------------------------------------------------
!     IOFLAG Indexing
! ---------------------------------

      INTEGER  iBMG2_BUG_STENCIL_FG,                                    &
     &         iBMG2_BUG_STENCIL_CG,                                    &
     &         iBMG2_BUG_STENCIL_CG1,                                   &
     &         iBMG2_BUG_RESTRICT,                                      &
     &         iBMG2_BUG_INTERP,                                        &
     &         iBMG2_BUG_RES_CG_SOLVE,                                  &
     &         iBMG2_BUG_RES_INTERP,                                    &
     &         iBMG2_BUG_RES_RELAX,                                     &
     &         iBMG2_BUG_RES_RESTRICT,                                  &
     &         iBMG2_BUG_PARAMETERS,                                    &
     &         iBMG2_OUT_ITERATIONS,                                    &
     &         iBMG2_OUT_STENCIL_TTY,                                   &
     &         iBMG2_OUT_RESTRICT_TTY,                                  &
     &         iBMG2_OUT_INTERP_TTY,                                    &
     &         iBMG2_OUT_TIME_CYCLING,                                  &
     &         iBMG2_OUT_TIME_SETUP,                                    &
     &         iBMG2_OUT_TIME_TOTAL,                                    &
     &         iBMG2_OUT_WSPACE_SIZE,                                   &
     &         iBMG2_OUT_WSPACE_POINT,                                  &
     &         iBMG2_WARN_ZERO_RESIDUAL,                                &
     &         iBMG2_OUT_STOP_ERROR,                                    &
     &         iBMG3_BUG_STENCIL_FG,                                    &
     &         iBMG3_BUG_STENCIL_CG,                                    &
     &         iBMG3_BUG_STENCIL_CG1,                                   &
     &         iBMG3_BUG_RESTRICT,                                      &
     &         iBMG3_BUG_INTERP,                                        &
     &         iBMG3_BUG_RES_CG_SOLVE,                                  &
     &         iBMG3_BUG_RES_INTERP,                                    &
     &         iBMG3_BUG_RES_RELAX,                                     &
     &         iBMG3_BUG_RES_RESTRICT,                                  &
     &         iBMG3_BUG_PARAMETERS,                                    &
     &         iBMG3_OUT_ITERATIONS,                                    &
     &         iBMG3_OUT_STENCIL_TTY,                                   &
     &         iBMG3_OUT_RESTRICT_TTY,                                  &
     &         iBMG3_OUT_INTERP_TTY,                                    &
     &         iBMG3_OUT_TIME_CYCLING,                                  &
     &         iBMG3_OUT_TIME_SETUP,                                    &
     &         iBMG3_OUT_TIME_TOTAL,                                    &
     &         iBMG3_OUT_WSPACE_SIZE,                                   &
     &         iBMG3_OUT_WSPACE_POINT,                                  &
     &         iBMG3_WARN_ZERO_RESIDUAL,                                &
     &         iBMG3_OUT_STOP_ERROR,                                    &
     &         iBMG_OUT_RHS,                                            &
     &         iBMG_OUT_SOLUTION,                                       &
     &         NBMG2_IOFLAG, NBMG3_IOFLAG, NBMG_IOFLAG

      PARAMETER( iBMG2_BUG_STENCIL_FG     = 1,                          &
     &           iBMG2_BUG_STENCIL_CG     = 2,                          & 
     &           iBMG2_BUG_STENCIL_CG1    = 3,                          &
     &           iBMG2_BUG_RESTRICT       = 4,                          &
     &           iBMG2_BUG_INTERP         = 5,                          &
     &           iBMG2_BUG_RES_INTERP     = 6,                          &
     &           iBMG2_BUG_RES_RESTRICT   = 7,                          &
     &           iBMG2_BUG_RES_RELAX      = 8,                          &
     &           iBMG2_BUG_RES_CG_SOLVE   = 9,                          &
     &           iBMG2_BUG_PARAMETERS     = 10,                         &
     &           iBMG2_OUT_WSPACE_SIZE    = 11,                         &
     &           iBMG2_OUT_WSPACE_POINT   = 12,                         &
     &           iBMG2_OUT_TIME_SETUP     = 13,                         &
     &           iBMG2_OUT_TIME_CYCLING   = 14,                         &
     &           iBMG2_OUT_TIME_TOTAL     = 15,                         &
     &           iBMG2_OUT_ITERATIONS     = 16,                         &
     &           iBMG2_OUT_STENCIL_TTY    = 17,                         &
     &           iBMG2_OUT_RESTRICT_TTY   = 18,                         &
     &           iBMG2_OUT_INTERP_TTY     = 19,                         &
     &           iBMG2_WARN_ZERO_RESIDUAL = 20,                         &
     &           iBMG2_OUT_STOP_ERROR     = 21,                         &
     &           iBMG3_BUG_STENCIL_FG     = 22,                         &
     &           iBMG3_BUG_STENCIL_CG     = 23,                         & 
     &           iBMG3_BUG_STENCIL_CG1    = 24,                         &
     &           iBMG3_BUG_RESTRICT       = 25,                         &
     &           iBMG3_BUG_INTERP         = 26,                         &
     &           iBMG3_BUG_RES_INTERP     = 27,                         &
     &           iBMG3_BUG_RES_RESTRICT   = 28,                         &
     &           iBMG3_BUG_RES_RELAX      = 29,                         &
     &           iBMG3_BUG_RES_CG_SOLVE   = 30,                         &
     &           iBMG3_BUG_PARAMETERS     = 31,                         &
     &           iBMG3_OUT_WSPACE_SIZE    = 32,                         &
     &           iBMG3_OUT_WSPACE_POINT   = 33,                         &
     &           iBMG3_OUT_TIME_SETUP     = 34,                         &
     &           iBMG3_OUT_TIME_CYCLING   = 35,                         &
     &           iBMG3_OUT_TIME_TOTAL     = 36,                         &
     &           iBMG3_OUT_ITERATIONS     = 37,                         & 
     &           iBMG3_OUT_STENCIL_TTY    = 38,                         &
     &           iBMG3_OUT_RESTRICT_TTY   = 39,                         &
     &           iBMG3_OUT_INTERP_TTY     = 40,                         &
     &           iBMG3_WARN_ZERO_RESIDUAL = 41,                         &
     &           iBMG3_OUT_STOP_ERROR     = 42,                         &
     &           iBMG_OUT_SOLUTION        = 43,                         & 
     &           iBMG_OUT_RHS             = 44,                         &
     &           NBMG2_IOFLAG = 22,                                     &
     &           NBMG3_IOFLAG = 22,                                     & 
     &           NBMG_IOFLAG  = 44                                      &
     &           )

! ----------------------------------------------------
! ==========================================================================




