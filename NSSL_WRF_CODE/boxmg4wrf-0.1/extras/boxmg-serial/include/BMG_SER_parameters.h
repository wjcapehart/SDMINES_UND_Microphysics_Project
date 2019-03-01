C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This include file provides a single resource for defining commonly
C     used parameters in the BOXMG family of codes.
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
C ----------------------------------------------------
C     Parameter Indexing
C ---------------------------------

C ---------------------------------
C     INTEGER Parameters
C ---------------------------------

      INTEGER   id_BMG2_SER_DIM_NOG,
     &          id_BMG2_SER_DIM_NF,
     &          id_BMG2_SER_DIM_NC,
     &          id_BMG2_SER_DIM_NSO,
     &          id_BMG2_SER_DIM_NSOR,
     &          id_BMG2_SER_DIM_NCI,
     &          id_BMG2_SER_DIM_NCBW,
     &          id_BMG2_SER_DIM_NCU, 
     &          id_BMG2_SER_POINTERS,
     &          id_BMG2_SER_STENCIL,
     &          id_BMG2_SER_BC,
     &          id_BMG2_SER_SETUP,
     &          id_BMG2_SER_RELAX,
     &          id_BMG2_SER_RELAX_SYM,
     &          id_BMG2_SER_NRELAX_DOWN,
     &          id_BMG2_SER_NRELAX_UP,
     &          id_BMG2_SER_NRELAX_FG,
     &          id_BMG2_SER_CYCLE_CLASS,
     &          id_BMG2_SER_NCYCLE_TYPE,
     &          id_BMG2_SER_FMG_NNCYCLE,
     &          id_BMG2_SER_MAX_ITERS,
     &          id_BMG2_SER_STOP_TEST,
     &          id_BMG2_SER_MIN_NOG,
     &          id_BMG2_SER_CG_MIN_DIM,
     &          id_BMG2_SER_CG_TYPE,
     &          id_BMG2_SER_CG_CONSTRUCT,
     &          id_BMG2_SER_Err_Code,
     &          id_BMG2_SER_Ext_Err_Code,
     &          id_BMG3_SER_DIM_NOG,
     &          id_BMG3_SER_DIM_NF,
     &          id_BMG3_SER_DIM_NC,
     &          id_BMG3_SER_DIM_NSO,
     &          id_BMG3_SER_DIM_NSOR,
     &          id_BMG3_SER_DIM_NCI,
     &          id_BMG3_SER_DIM_NCBW,
     &          id_BMG3_SER_DIM_NCU, 
     &          id_BMG3_SER_DIM_NSO3, 
     &          id_BMG3_SER_DIM_NSR3, 
     &          id_BMG3_SER_DIM_NCI3, 
     &          id_BMG3_SER_DIM_NCPL, 
     &          id_BMG3_SER_POINTERS,
     &          id_BMG3_SER_STENCIL,
     &          id_BMG3_SER_BC,
     &          id_BMG3_SER_SETUP,
     &          id_BMG3_SER_RELAX,
     &          id_BMG3_SER_RELAX_SYM,
     &          id_BMG3_SER_NRELAX_DOWN,
     &          id_BMG3_SER_NRELAX_UP,
     &          id_BMG3_SER_NRELAX_FG,
     &          id_BMG3_SER_CYCLE_CLASS,
     &          id_BMG3_SER_NCYCLE_TYPE,
     &          id_BMG3_SER_FMG_NNCYCLE,
     &          id_BMG3_SER_MAX_ITERS,
     &          id_BMG3_SER_STOP_TEST,
     &          id_BMG3_SER_MIN_NOG,
     &          id_BMG3_SER_CG_MIN_DIM,
     &          id_BMG3_SER_CG_TYPE,
     &          id_BMG3_SER_CG_CONSTRUCT,
     &          id_BMG3_SER_Err_Code,
     &          id_BMG3_SER_Ext_Err_Code,
     &          NBMG2_SER_iPARMS, NBMG3_SER_iPARMS, NBMG_SER_iPARMS 

      PARAMETER ( id_BMG2_SER_DIM_NOG      =  1,
     &            id_BMG2_SER_DIM_NF       =  2,
     &            id_BMG2_SER_DIM_NC       =  3,
     &            id_BMG2_SER_DIM_NSO      =  4,
     &            id_BMG2_SER_DIM_NSOR     =  5,
     &            id_BMG2_SER_DIM_NCI      =  6,
     &            id_BMG2_SER_DIM_NCBW     =  7,
     &            id_BMG2_SER_DIM_NCU      =  8, 
     &            id_BMG2_SER_POINTERS     =  9,
     &            id_BMG2_SER_STENCIL      =  10,
     &            id_BMG2_SER_BC           =  11, 
     &            id_BMG2_SER_SETUP        =  12,
     &            id_BMG2_SER_RELAX        =  13,
     &            id_BMG2_SER_RELAX_SYM    =  14,
     &            id_BMG2_SER_NRELAX_DOWN  =  15,
     &            id_BMG2_SER_NRELAX_UP    =  16,
     &            id_BMG2_SER_NRELAX_FG    =  17,
     &            id_BMG2_SER_CYCLE_CLASS  =  18,
     &            id_BMG2_SER_NCYCLE_TYPE  =  19,
     &            id_BMG2_SER_FMG_NNCYCLE  =  20,
     &            id_BMG2_SER_MAX_ITERS    =  21,
     &            id_BMG2_SER_STOP_TEST    =  22,    
     &            id_BMG2_SER_MIN_NOG      =  23,
     &            id_BMG2_SER_CG_MIN_DIM   =  24,
     &            id_BMG2_SER_CG_CONSTRUCT =  25,
     &            id_BMG2_SER_CG_TYPE      =  26,
     &            id_BMG2_SER_Err_Code     =  27,
     &            id_BMG2_SER_Ext_Err_Code =  28,
     &            id_BMG3_SER_DIM_NOG      =  29,
     &            id_BMG3_SER_DIM_NF       =  30,
     &            id_BMG3_SER_DIM_NC       =  31,
     &            id_BMG3_SER_DIM_NSO      =  32,
     &            id_BMG3_SER_DIM_NSOR     =  33,
     &            id_BMG3_SER_DIM_NCI      =  34,
     &            id_BMG3_SER_DIM_NCBW     =  35,
     &            id_BMG3_SER_DIM_NCU      =  36, 
     &            id_BMG3_SER_DIM_NSO3     =  37, 
     &            id_BMG3_SER_DIM_NSR3     =  38, 
     &            id_BMG3_SER_DIM_NCI3     =  39, 
     &            id_BMG3_SER_DIM_NCPL     =  40, 
     &            id_BMG3_SER_POINTERS     =  41,
     &            id_BMG3_SER_STENCIL      =  42,
     &            id_BMG3_SER_BC           =  43,
     &            id_BMG3_SER_SETUP        =  44,
     &            id_BMG3_SER_RELAX        =  45,
     &            id_BMG3_SER_RELAX_SYM    =  46,
     &            id_BMG3_SER_NRELAX_DOWN  =  47,
     &            id_BMG3_SER_NRELAX_UP    =  48,
     &            id_BMG3_SER_NRELAX_FG    =  49,
     &            id_BMG3_SER_CYCLE_CLASS  =  50,
     &            id_BMG3_SER_NCYCLE_TYPE  =  51,
     &            id_BMG3_SER_FMG_NNCYCLE  =  52,
     &            id_BMG3_SER_MAX_ITERS    =  53,
     &            id_BMG3_SER_STOP_TEST    =  54,
     &            id_BMG3_SER_MIN_NOG      =  55,
     &            id_BMG3_SER_CG_MIN_DIM   =  56,
     &            id_BMG3_SER_CG_TYPE      =  57,
     &            id_BMG3_SER_CG_CONSTRUCT =  58,
     &            id_BMG3_SER_Err_Code     =  59,
     &            id_BMG3_SER_Ext_Err_Code =  60,
     &            NBMG2_SER_iPARMS = 28,      ! Number of 2D parameters
     &            NBMG3_SER_iPARMS = 32,      ! Number of 3D parameters
     &            NBMG_SER_iPARMS  = 60       ! Dimension of BMG_SER_iPARM
     &            )

C -------------------------------
C     REAL Parameters
C -------------------------------

      INTEGER   id_BMG2_SER_STOP_TOL,
     &          id_BMG3_SER_STOP_TOL,
     &          NBMG_SER_rPARMS 

      PARAMETER ( id_BMG2_SER_STOP_TOL  = 1,
     &            id_BMG3_SER_STOP_TOL  = 2,
     &            NBMG_SER_rPARMS = 2       )


C ----------------------------------------------------
C ==========================================================================
C ----------------------------------------------------
C     Parameter Values
C -------------------------------

C -------------------------------
C     Stencil Type
C -------------------------------

      INTEGER BMG_SER_STENCIL_5pt,
     &        BMG_SER_STENCIL_9pt,
     &        BMG_SER_STENCIL_7pt,
     &        BMG_SER_STENCIL_27pt

      PARAMETER ( BMG_SER_STENCIL_5pt  = 1,
     &            BMG_SER_STENCIL_9pt  = 2,
     &            BMG_SER_STENCIL_7pt  = 1,
     &            BMG_SER_STENCIL_27pt = 2  )


C -------------------------------
C     Periodicity
C -------------------------------

      INTEGER  BMG_SER_BCs_definite,
     &         BMG_SER_BCs_def_per_x,
     &         BMG_SER_BCs_def_per_y,
     &         BMG_SER_BCs_def_per_xy,
     &         BMG_SER_BCs_indef_per_x,
     &         BMG_SER_BCs_indef_per_y,
     &         BMG_SER_BCs_indef_per_xy,
     &         BMG_SER_BCs_indef_nonper

      PARAMETER ( BMG_SER_BCs_definite     =  0,
     &            BMG_SER_BCs_def_per_x    =  2,
     &            BMG_SER_BCs_def_per_y    =  1,
     &            BMG_SER_BCs_def_per_xy   =  3,
     &            BMG_SER_BCs_indef_per_x  = -2,
     &            BMG_SER_BCs_indef_per_y  = -1,
     &            BMG_SER_BCs_indef_per_xy = -3,
     &            BMG_SER_BCs_indef_nonper = -4  )

C -------------------------------
C     Setup options:
C -------------------------------

      INTEGER BMG_SER_SETUP_only,
     &        BMG_SER_SETUP_none,
     &        BMG_SER_SETUP_opers, 
     &        BMG_SER_SETUP_ptrs_opers 

      PARAMETER ( BMG_SER_SETUP_only  = 3,
     &            BMG_SER_SETUP_none  = 2,
     &            BMG_SER_SETUP_opers = 1,
     &            BMG_SER_SETUP_ptrs_opers = 0 )

C -------------------------------
C     Memory Allocation:
C -------------------------------

      INTEGER BMG_SER_USE_pointers,
     &        BMG_SER_NO_pointers

      PARAMETER ( BMG_SER_USE_pointers = 1,
     &            BMG_SER_NO_pointers = 0   )

C -------------------------------
C     Relaxation:
C -------------------------------
      
      INTEGER  BMG_SER_GS_RB_point, 
     &         BMG_SER_GS_RB_x_lines,
     &         BMG_SER_GS_RB_y_lines,
     &         BMG_SER_GS_RB_x_y_lines,
     &         BMG_SER_GS_RB_planes_xy_yz_xz 

      PARAMETER( BMG_SER_GS_RB_point     = 1,
     &           BMG_SER_GS_RB_x_lines   = 2,
     &           BMG_SER_GS_RB_y_lines   = 3,
     &           BMG_SER_GS_RB_x_y_lines = 4,
     &           BMG_SER_GS_RB_planes_xy_yz_xz = 5 )


C --------------------------------
C     Symmetry of the MG n-cycle:
C --------------------------------

      INTEGER  BMG_SER_RELAX_NONSYM,
     &         BMG_SER_RELAX_SYM 
      PARAMETER( BMG_SER_RELAX_NONSYM = 0, 
     &           BMG_SER_RELAX_SYM = 1    )

      INTEGER  BMG_SER_DOWN, 
     &         BMG_SER_UP
      PARAMETER(  BMG_SER_DOWN = 0, 
     &            BMG_SER_UP = 1   )

C --------------------------------
C     Stopping Critieria
C --------------------------------

      INTEGER BMG_SER_STOP_ABS_RES_L2,
     &        BMG_SER_STOP_REL_RES_L2

      PARAMETER( BMG_SER_STOP_ABS_RES_L2 = 0, 
     &           BMG_SER_STOP_REL_RES_L2 = 1 ) 

C --------------------------------
C     Cycle Class and Type
C --------------------------------
      
      INTEGER   BMG_SER_FMG_CYCLE, 
     &          BMG_SER_N_CYCLE, 
     &          BMG_SER_V_CYCLE, 
     &          BMG_SER_W_CYCLE

      PARAMETER ( BMG_SER_FMG_CYCLE = 0,
     &            BMG_SER_N_CYCLE   = 1,
     &            BMG_SER_V_CYCLE   = 1,
     &            BMG_SER_W_CYCLE   = 2 )


C --------------------------------
C     Coarse-Grid Operator Type
C --------------------------------

      INTEGER   BMG_SER_CG_ITLI_IzIyIx,
     &          BMG_SER_CG_ITLI,
     &          BMG_SER_CG_USER

      PARAMETER ( BMG_SER_CG_ITLI_IzIyIx = 1,
     &            BMG_SER_CG_ITLI        = 2,
     &            BMG_SER_CG_USER        = 3 )


C --------------------------------
C     I^{T} L I Construction
C --------------------------------

      INTEGER   BMG_SER_CG_CONS_explicit,
     &          BMG_SER_CG_CONS_block

      PARAMETER ( BMG_SER_CG_CONS_explicit = 1,
     &            BMG_SER_CG_CONS_block    = 2  )
                 

C ----------------------------------------------------
C ==========================================================================
C ----------------------------------------------------
C     IOFLAG Indexing
C ---------------------------------

      INTEGER  iBMG2_SER_BUG_STENCIL_FG,
     &         iBMG2_SER_BUG_STENCIL_CG,
     &         iBMG2_SER_BUG_STENCIL_CG1,
     &         iBMG2_SER_BUG_RESTRICT,
     &         iBMG2_SER_BUG_INTERP,
     &         iBMG2_SER_BUG_RES_CG_SOLVE,         
     &         iBMG2_SER_BUG_RES_INTERP,
     &         iBMG2_SER_BUG_RES_RELAX,
     &         iBMG2_SER_BUG_RES_RESTRICT,
     &         iBMG2_SER_BUG_PARAMETERS,
     &         iBMG2_SER_OUT_ITERATIONS,
     &         iBMG2_SER_OUT_STENCIL_TTY,
     &         iBMG2_SER_OUT_RESTRICT_TTY,
     &         iBMG2_SER_OUT_INTERP_TTY,
     &         iBMG2_SER_OUT_TIME_CYCLING,
     &         iBMG2_SER_OUT_TIME_SETUP,
     &         iBMG2_SER_OUT_TIME_TOTAL,
     &         iBMG2_SER_OUT_WSPACE_SIZE,
     &         iBMG2_SER_OUT_WSPACE_POINT,
     &         iBMG2_SER_WARN_ZERO_RESIDUAL,
     &         iBMG2_SER_OUT_STOP_ERROR,
     &         iBMG3_SER_BUG_STENCIL_FG,
     &         iBMG3_SER_BUG_STENCIL_CG,
     &         iBMG3_SER_BUG_STENCIL_CG1,
     &         iBMG3_SER_BUG_RESTRICT,
     &         iBMG3_SER_BUG_INTERP,
     &         iBMG3_SER_BUG_RES_CG_SOLVE,         
     &         iBMG3_SER_BUG_RES_INTERP,
     &         iBMG3_SER_BUG_RES_RELAX,         
     &         iBMG3_SER_BUG_RES_RESTRICT,
     &         iBMG3_SER_BUG_PARAMETERS,
     &         iBMG3_SER_OUT_ITERATIONS,
     &         iBMG3_SER_OUT_STENCIL_TTY,
     &         iBMG3_SER_OUT_RESTRICT_TTY,
     &         iBMG3_SER_OUT_INTERP_TTY,
     &         iBMG3_SER_OUT_TIME_CYCLING,
     &         iBMG3_SER_OUT_TIME_SETUP,
     &         iBMG3_SER_OUT_TIME_TOTAL,
     &         iBMG3_SER_OUT_WSPACE_SIZE,
     &         iBMG3_SER_OUT_WSPACE_POINT,
     &         iBMG3_SER_WARN_ZERO_RESIDUAL,
     &         iBMG3_SER_OUT_STOP_ERROR,
     &         iBMG_SER_OUT_RHS,
     &         iBMG_SER_OUT_SOLUTION,
     &         NBMG2_SER_IOFLAG, NBMG3_SER_IOFLAG, NBMG_SER_IOFLAG

      PARAMETER( iBMG2_SER_BUG_STENCIL_FG     = 1,
     &           iBMG2_SER_BUG_STENCIL_CG     = 2, 
     &           iBMG2_SER_BUG_STENCIL_CG1    = 3,
     &           iBMG2_SER_BUG_RESTRICT       = 4,
     &           iBMG2_SER_BUG_INTERP         = 5,
     &           iBMG2_SER_BUG_RES_INTERP     = 6,
     &           iBMG2_SER_BUG_RES_RESTRICT   = 7,
     &           iBMG2_SER_BUG_RES_RELAX      = 8,
     &           iBMG2_SER_BUG_RES_CG_SOLVE   = 9,
     &           iBMG2_SER_BUG_PARAMETERS     = 10,
     &           iBMG2_SER_OUT_WSPACE_SIZE    = 11,
     &           iBMG2_SER_OUT_WSPACE_POINT   = 12,
     &           iBMG2_SER_OUT_TIME_SETUP     = 13,
     &           iBMG2_SER_OUT_TIME_CYCLING   = 14,
     &           iBMG2_SER_OUT_TIME_TOTAL     = 15,
     &           iBMG2_SER_OUT_ITERATIONS     = 16,
     &           iBMG2_SER_OUT_STENCIL_TTY    = 17,
     &           iBMG2_SER_OUT_RESTRICT_TTY   = 18,
     &           iBMG2_SER_OUT_INTERP_TTY     = 19,
     &           iBMG2_SER_WARN_ZERO_RESIDUAL = 20,
     &           iBMG2_SER_OUT_STOP_ERROR     = 21,
     &           iBMG3_SER_BUG_STENCIL_FG     = 22,
     &           iBMG3_SER_BUG_STENCIL_CG     = 23, 
     &           iBMG3_SER_BUG_STENCIL_CG1    = 24,
     &           iBMG3_SER_BUG_RESTRICT       = 25,
     &           iBMG3_SER_BUG_INTERP         = 26,
     &           iBMG3_SER_BUG_RES_INTERP     = 27,
     &           iBMG3_SER_BUG_RES_RESTRICT   = 28,
     &           iBMG3_SER_BUG_RES_RELAX      = 29,
     &           iBMG3_SER_BUG_RES_CG_SOLVE   = 30,
     &           iBMG3_SER_BUG_PARAMETERS     = 31,
     &           iBMG3_SER_OUT_WSPACE_SIZE    = 32,
     &           iBMG3_SER_OUT_WSPACE_POINT   = 33,
     &           iBMG3_SER_OUT_TIME_SETUP     = 34,
     &           iBMG3_SER_OUT_TIME_CYCLING   = 35,
     &           iBMG3_SER_OUT_TIME_TOTAL     = 36,
     &           iBMG3_SER_OUT_ITERATIONS     = 37,
     &           iBMG3_SER_OUT_STENCIL_TTY    = 38,
     &           iBMG3_SER_OUT_RESTRICT_TTY   = 39,
     &           iBMG3_SER_OUT_INTERP_TTY     = 40,
     &           iBMG3_SER_WARN_ZERO_RESIDUAL = 41,
     &           iBMG3_SER_OUT_STOP_ERROR     = 42,
     &           iBMG_SER_OUT_SOLUTION        = 43, 
     &           iBMG_SER_OUT_RHS             = 44,
     &           NBMG2_SER_IOFLAG = 22,       ! Number of 2D I/O FLAGS
     &           NBMG3_SER_IOFLAG = 22,       ! Number of 3D I/O FLAGS
     &           NBMG_SER_IOFLAG  = 44        ! Dimension of BMG_IOFLAG
     &           )

C ----------------------------------------------------
C ==========================================================================



