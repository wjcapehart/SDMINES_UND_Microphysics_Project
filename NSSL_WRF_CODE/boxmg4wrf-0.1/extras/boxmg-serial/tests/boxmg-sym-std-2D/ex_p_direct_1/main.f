      PROGRAM MAIN

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This solves Example 1 (periodic case) given in Victor Bandy's BBMG
C     interface guide.  But it uses a direct call to the symmetric
C     two-dimensional BoxMG solver (BMG2_SER_SymStd_SOLVE_boxmg), and it
C     uses a W-cycle instead of F-cycles.  The discretization (putf),
C     and the pre/post processing (normal/checker) are the subroutines
C     from the BBMG package, and consquently they fall under the BBMG
C     public-domain license and copyright (see extras/bbmg/README).
C
C =======================================================================
C $license_flag$
C =======================================================================
C ==========================================================================

      IMPLICIT  NONE

      INCLUDE   'BMG_SER_workspace.h'
      INCLUDE   'BMG_SER_parameters.h'

C --------------------------------------------
C     Critical Parameters:
C
      INTEGER   Lm, Mm, NXYCm
      PARAMETER ( Lm=10, Mm=10, NXYCm=3 )

C --------------------------------------------
C     Multigrid/Workspace Memory Allocation:
C
      INTEGER   NFm, NOGm, NSOm
      PARAMETER ( NFm=218, NSOm=802, NOGm=4 )

      INTEGER   NBMG_iWORK, NBMG_rWORK
      PARAMETER ( NBMG_iWORK=27, NBMG_rWORK=10534 )

      INTEGER   pSI, pSR, BMG_pWORK(NBMG_SER_pWORK)
      LOGICAL   BMG_InWORK(NBMG_SER_InWORK)

C --------------------------------------------
C     Variable Declarations:
C      
      INTEGER   BMG_iPARMS(NBMG_SER_iPARMS),
     &          BMG_iWORK(NBMG_iWORK),
     &          i, KG, L, M, NC, NCI, NCU,
     &          NF, NOG, NSO, NCBW, NSOR
      REAL*8    BMG_rPARMS(NBMG_SER_rPARMS),
     &          hx, hy, Pi, RHS(NFm),
     &          BMG_rWORK(NBMG_rWORK),
     &          SO(NSOm), u(NFm), utemp(NFm)
      LOGICAL   BMG_IOFLAG(NBMG_SER_IOFLAG)

C ===========================================================================

C ---------------------------------------------------------
C     Initialize the arguments: 
C ----------------------------------
C     Domain and grid:
C ----------------------------------

      L=10                     ! # of grid points in x
      M=10                     ! # of grid points in y

      Pi=4D0*ATAN(1D0)         ! Pi=3.1415926....

      hx=2*PI/L                ! grid spacing in x
      hy=PI/M                  ! grid spacing in y

      WRITE(*,200) '*** Solving Victor''s periodic example ***'

C ========================================================================
C     >>>>>>>>>>>>  Setup the first case: FMG with Line Relaxation
C ========================================================================

C ------------------------------
C     Multigird parameters:
C  -----------------------------
C     NB: parameters affect
C         space requirements!!
C ------------------------------

      CALL BMG2_SER_SymStd_SETUP_parms( 
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG
     &          )

      BMG_iPARMS(id_BMG2_SER_STENCIL)  = BMG_SER_STENCIL_5pt
      BMG_iPARMS(id_BMG2_SER_BC)       = BMG_SER_BCs_indef_per_x
      BMG_rPARMS(id_BMG2_SER_STOP_TOL) = 1D-6

      BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS) = BMG_SER_FMG_CYCLE
      BMG_iPARMS(id_BMG2_SER_NCYCLE_TYPE) = BMG_SER_V_CYCLE
      BMG_iPARMS(id_BMG2_SER_RELAX_SYM )  = BMG_SER_RELAX_NONSYM

C -------------------------------------
C     I/O Parameters
C -------------------------------------

      BMG_IOFLAG(iBMG2_SER_OUT_WSPACE_SIZE)  = .TRUE.
      BMG_IOFLAG(iBMG2_SER_OUT_ITERATIONS)   = .TRUE.

      BMG_IOFLAG(iBMG2_SER_OUT_TIME_SETUP)   = .TRUE.
      BMG_IOFLAG(iBMG2_SER_OUT_TIME_CYCLING) = .TRUE.
      BMG_IOFLAG(iBMG2_SER_OUT_TIME_TOTAL)   = .TRUE.

      BMG_IOFLAG(iBMG2_SER_BUG_RES_RELAX)    = .TRUE.
      BMG_IOFLAG(iBMG2_SER_BUG_RES_CG_SOLVE) = .TRUE.

C -------------------------------------
C     Space requirements and pointers:
C -------------------------------------

      pSR=1
      pSI=1


      DO i=1, NBMG_SER_InWORK
         BMG_InWORK(i) = .FALSE.
      ENDDO

      BMG_InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      BMG_InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      BMG_InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      BMG_InWORK(i_InWORK_RES) = .TRUE.     ! we don't have RES yet!!!

      CALL BMG2_SER_SymStd_SETUP_PtrWork( 
     &                 L, M, BMG_iPARMS, 
     &                 NOGm, NFm, NSOm, NBMG_iWORK, NBMG_rWORK,
     &                 BMG_pWORK, BMG_InWORK, pSR, pSI 
     &                 ) 

      NOG  = BMG_iPARMS(id_BMG2_SER_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG2_SER_DIM_NF)
      NC   = BMG_iPARMS(id_BMG2_SER_DIM_NC)
      NSO  = BMG_iPARMS(id_BMG2_SER_DIM_NSO)
      NCI  = BMG_iPARMS(id_BMG2_SER_DIM_NCI)
      NSOR = BMG_iPARMS(id_BMG2_SER_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG2_SER_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG2_SER_DIM_NCU)

C --------------------------------------
C     Zero everyting: 
C       of particular importance is 
C       the rhs, as we are doing FMG.
C ---------------------------------------

      DO 10 i=1, NSO
         SO(i)=0
 10   CONTINUE

      DO 20 i=1, NF
         u(i)=0
         rhs(i)=0
 20   CONTINUE

C --------------------------------------
C     Define the coefficient matrix
C --------------------------------------

      CALL PUTF(SO, rhs, L, M, hx, hy, 0, BMG_iPARMS(id_BMG2_SER_BC) )

C --------------------------------------
C     Solve the system
C --------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_boxmg(
     &          L, M, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          u, rhs, BMG_rWORK(BMG_pWORK(ip_RES)), NF, NC,
     &          SO, NSO, BMG_rWORK(BMG_pWORK(ip_SOR)), NSOR,
     &          BMG_rWORK(BMG_pWORK(ip_CI)), NCI, 
     &          BMG_rWORK(BMG_pWORK(ip_CSO)),
     &          BMG_rWORK(BMG_pWORK(ip_CU)), NCBW, NCU,
     &          BMG_iWORK(BMG_pWORK(ip_iG)), NOG, KG                         
     &          )

C --------------------------------------
C     Normailize the computed solution
C --------------------------------------

      CALL NORMAL(u, hx, hy, L, M)

C --------------------------------------
C     Compare with the exact solution
C --------------------------------------

      CALL CHECKER(u, utemp, L, M, hx, hy, BMG_iPARMS(id_BMG2_SER_BC) )

C ========================================================================
C     >>>>>>>>>>>>  Do it all again but N-cycle with Point Relaxation:
C ========================================================================

C ------------------------------
C     Multigird parameters:
C ------------------------------

      BMG_iPARMS(id_BMG2_SER_RELAX)       = BMG_SER_GS_RB_point
      BMG_iPARMS(id_BMG2_SER_RELAX_SYM )  = BMG_SER_RELAX_NONSYM
      BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS) = BMG_SER_N_CYCLE
      BMG_iPARMS(id_BMG2_SER_NCYCLE_TYPE) = BMG_SER_V_CYCLE

      BMG_iPARMS(id_BMG2_SER_MAX_ITERS) = 20
      BMG_rPARMS(id_BMG2_SER_STOP_TOL)  = 1D-7

C -------------------------------------
C     Space requirements and pointers:
C -------------------------------------

      pSR=1
      pSI=1

      BMG_InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      BMG_InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      BMG_InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      BMG_InWORK(i_InWORK_RES) = .TRUE.    ! we don't have RES yet!!!

      CALL BMG2_SER_SymStd_SETUP_PtrWork( 
     &                 L, M, BMG_iPARMS, 
     &                 NOGm, NFm, NSOm, NBMG_iWORK, NBMG_rWORK,
     &                 BMG_pWORK, BMG_InWORK, pSR, pSI 
     &                 ) 

      NOG  = BMG_iPARMS(id_BMG2_SER_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG2_SER_DIM_NF)
      NC   = BMG_iPARMS(id_BMG2_SER_DIM_NC)
      NSO  = BMG_iPARMS(id_BMG2_SER_DIM_NSO)
      NCI  = BMG_iPARMS(id_BMG2_SER_DIM_NCI)
      NSOR = BMG_iPARMS(id_BMG2_SER_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG2_SER_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG2_SER_DIM_NCU)

C --------------------------------------
C     Zero everyting: 
C       of particular importance is 
C       the rhs, as we are doing FMG.
C ---------------------------------------

      DO i=1, NSO
         SO(i)=0
      ENDDO

      DO i=1, NF
         u(i)=0
         rhs(i)=0
      ENDDO
      
C --------------------------------------
C     Define the coefficient matrix
C --------------------------------------

      CALL PUTF(SO, rhs, L, M, hx, hy, 0, BMG_iPARMS(id_BMG2_SER_BC) )

C --------------------------------------
C     Solve the system
C --------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_boxmg( 
     &          L, M, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          u, rhs, BMG_rWORK(BMG_pWORK(ip_RES)), NF, NC,
     &          SO, NSO, BMG_rWORK(BMG_pWORK(ip_SOR)), NSOR,
     &          BMG_rWORK(BMG_pWORK(ip_CI)), NCI, 
     &          BMG_rWORK(BMG_pWORK(ip_CSO)),
     &          BMG_rWORK(BMG_pWORK(ip_CU)), NCBW, NCU,
     &          BMG_iWORK(BMG_pWORK(ip_iG)), NOG, KG
     &          )

C --------------------------------------
C     Normailize the computed solution
C --------------------------------------

      CALL NORMAL(u, hx, hy, L, M)

C --------------------------------------
C     Compare with the exact solution
C --------------------------------------

      CALL CHECKER(u, utemp, L, M, hx, hy, BMG_iPARMS(id_BMG2_SER_BC) )


C ========================================================================
C     >>>>>>>>>>>>  Do it all again but skip the setup:
C ========================================================================

C ------------------------------
C     Multigird parameters:
C ------------------------------

      BMG_iPARMS(id_BMG2_SER_SETUP)    = BMG_SER_SETUP_none
      BMG_rPARMS(id_BMG2_SER_STOP_TOL) = 1D-7

C --------------------------------------
C     Zero the Initial Guess:
C ---------------------------------------

      DO i=1, NF
         u(i)=0
      ENDDO

C --------------------------------------
C     Solve the system
C --------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_boxmg(
     &          L, M, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          u, rhs, BMG_rWORK(BMG_pWORK(ip_RES)), NF, NC,
     &          SO, NSO, BMG_rWORK(BMG_pWORK(ip_SOR)), NSOR,
     &          BMG_rWORK(BMG_pWORK(ip_CI)), NCI, 
     &          BMG_rWORK(BMG_pWORK(ip_CSO)),
     &          BMG_rWORK(BMG_pWORK(ip_CU)), NCBW, NCU,
     &          BMG_iWORK(BMG_pWORK(ip_iG)), NOG, KG
     &          )

C --------------------------------------
C     Normailize the computed solution
C --------------------------------------

      CALL NORMAL(u, hx, hy, L, M)

C --------------------------------------
C     Compare with the exact solution
C --------------------------------------

      CALL CHECKER(u, utemp, L, M, hx, hy, BMG_iPARMS(id_BMG2_SER_BC) )


C ==========================================================================

C ========================================================================
C     >>>>>>>>>>>>  Do it all again but for the degenerate case NOG=1!
C ========================================================================

      WRITE(*,*) ' Testing degenerate case of NOG = 1'

C ------------------------------
C     Multigird parameters:
C ------------------------------

      BMG_iPARMS(id_BMG2_SER_SETUP)      = BMG_SER_SETUP_ptrs_opers
      BMG_iPARMS(id_BMG2_SER_BC)         = BMG_SER_BCs_indef_per_x
      BMG_iPARMS(id_BMG2_SER_STENCIL)    = BMG_SER_STENCIL_5pt
      BMG_iPARMS(id_BMG2_SER_MIN_NOG)    = 1
      BMG_iPARMS(id_BMG2_SER_CG_MIN_DIM) = L
      BMG_rPARMS(id_BMG2_SER_STOP_TOL)   = 1D-6

C -------------------------------------
C     Space requirements and pointers:
C -------------------------------------

      pSR=1
      pSI=1

      BMG_InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      BMG_InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      BMG_InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      BMG_InWORK(i_InWORK_RES) = .TRUE.     ! we don't have RES yet!!!

      CALL BMG2_SER_SymStd_SETUP_PtrWork( 
     &                 L, M, BMG_iPARMS, 
     &                 NOGm, NFm, NSOm, NBMG_iWORK, NBMG_rWORK,
     &                 BMG_pWORK, BMG_InWORK, pSR, pSI 
     &                 ) 

      NOG  = BMG_iPARMS(id_BMG2_SER_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG2_SER_DIM_NF)
      NC   = BMG_iPARMS(id_BMG2_SER_DIM_NC)
      NSO  = BMG_iPARMS(id_BMG2_SER_DIM_NSO)
      NCI  = BMG_iPARMS(id_BMG2_SER_DIM_NCI)
      NSOR = BMG_iPARMS(id_BMG2_SER_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG2_SER_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG2_SER_DIM_NCU)

C --------------------------------------
C     Zero everyting: 
C       of particular importance is 
C       the rhs, as we are doing FMG.
C ---------------------------------------

      DO i=1, NSO
         SO(i)=0
      ENDDO

      DO i=1, NF
         u(i)=0
         rhs(i)=0
      ENDDO
      
C --------------------------------------
C     Define the coefficient matrix
C --------------------------------------

      CALL PUTF(SO, rhs, L, M, hx, hy, 0, BMG_iPARMS(id_BMG2_SER_BC) )

C --------------------------------------
C     Solve the system
C --------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_boxmg( 
     &          L, M, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          u, rhs, BMG_rWORK(BMG_pWORK(ip_RES)), NF, NC,
     &          SO, NSO, BMG_rWORK(BMG_pWORK(ip_SOR)), NSOR,
     &          BMG_rWORK(BMG_pWORK(ip_CI)), NCI, 
     &          BMG_rWORK(BMG_pWORK(ip_CSO)),
     &          BMG_rWORK(BMG_pWORK(ip_CU)), NCBW, NCU,
     &          BMG_iWORK(BMG_pWORK(ip_iG)), NOG, KG
     &          )

C --------------------------------------
C     Normailize the computed solution
C --------------------------------------

      CALL NORMAL(u, hx, hy, L, M)

C --------------------------------------
C     Compare with the exact solution
C --------------------------------------

      CALL CHECKER(u, utemp, L, M, hx, hy, BMG_iPARMS(id_BMG2_SER_BC) )

C ==========================================================================

 200  FORMAT (/,2X,A,2X,/)

C ======================================

      END

      



