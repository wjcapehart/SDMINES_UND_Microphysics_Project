      PROGRAM MAIN

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This solves Example 1 given in Victor Bandy's BBMG interface guide
C     using a direct call to the Preconditioned Conjugate Gradients
C     (PCG) solver with (BMG2_SER_SymStd_SOLVE_pcg). It uses BoxMG as a
C     preconditioner.  The discretization (putf), and the post
C     processing (diff) are the subroutines from the BBMG package, and
C     consquently they fall under the BBMG public-domain license and
C     copyright (see extras/bbmg).
C
C =======================================================================
C $license_flag$
C =======================================================================
C ==========================================================================

      IMPLICIT  NONE
      
      INCLUDE   'BMG_SER_workspace.h'
      INCLUDE   'BMG_SER_parameters.h'
      INCLUDE   'BMG_SER_PCG_parameters.h' 

C ---------------------------------------------
C     Dimensional Parameters:
C
      INTEGER   Lm, Mm, NXYCm
      PARAMETER ( Lm=30, Mm=60, NXYCm=3 )

C --------------------------------------------
C     Multigrid/Workspace Memory Allocation:
C
      INTEGER   NFm, NOGm, NSOm
      PARAMETER ( NFm=2758, NSOm=9822, NOGm=4 )

      INTEGER   NIWORK, NRWORK
      PARAMETER ( NIWORK=36, NRWORK=28480 )

      INTEGER   pSI, pSR, pWORK(NBMG_SER_pWORK)
      LOGICAL   InWORK(NBMG_SER_InWORK)

      INTEGER   MAX_PCG_ITERSm
      PARAMETER ( MAX_PCG_ITERSm = 1000 )

C --------------------------------------------
C     Variable Declarations:
C      
      INTEGER   BMG_iPARMS(NBMG_SER_iPARMS),
     &          i, IWORK(NIWORK), L, M,
     &          MAX_PCG_ITERS, NF, NOG, NSO
      REAL*8    BMG_rPARMS(NBMG_SER_rPARMS),
     &          hx, hy, MAXnorm, RHS(NFm),
     &          RWORK(NRWORK), SO(NSOm), u(NFm)
      LOGICAL   BMG_IOFLAG(NBMG_SER_IOFLAG)

      INTEGER   BMG_SER_PCG_iPARMS(NBMG_SER_SER_PCG_iPARMS)
      REAL*8    BMG_SER_PCG_rPARMS(NBMG_SER_SER_PCG_rPARMS),
     &          BMG_SER_PCG_RES(MAX_PCG_ITERSm)

C ==========================================================================


C ---------------------------------------------------------
C     Initialize the arguments
C ---------------------------------- 
C     Domain and grid:
C ----------------------------------

      WRITE(*,200) 'Enter the number of points in x ... '
      READ(*,*) L
      WRITE(*,200) 'Enter the number of points in y ... '
      READ(*,*) M

      IF (L.GT.Lm) THEN
         WRITE(*,*) '***Error: memory allocation insufficient in x ***'
         STOP
      ELSE IF (M.GT.Mm) THEN
         WRITE(*,*) '***Error: memory allocation insufficient in y ***'
         STOP
      ENDIF

      hx=1D0/L                  ! grid spacing in x
      hy=2D0/M                  ! grid spacing in y


      WRITE(*,*) 
      WRITE(*,*) 
      WRITE(*,*) ' TESTING:  Preconditioned Conjugate Gradients ... '
      WRITE(*,*) 
      WRITE(*,*) ' **** PCG with no preconditioning ... '
      WRITE(*,*) 

C ------------------------------
C     PCG parameters:
C ------------------------------

      BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_PRECON) 
     & = BMG_SER_PCG_PRECON_NONE

      BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_NMG_CYCLES) = 1
      BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_STOP_TEST)  = 1
      BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_MAX_ITERS)  = 1000
      BMG_SER_PCG_rPARMS(id_BMG_SER_PCG_STOP_TOL)   = 1D-5

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

C -----------------------------------
C     Multigrid parameters necessary
C     to work as preconditioner
C -----------------------------------

      BMG_iPARMS(id_BMG2_SER_RELAX)       = BMG_SER_GS_RB_x_y_lines
      BMG_iPARMS(id_BMG2_SER_RELAX_SYM)   = BMG_SER_RELAX_SYM
      BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS) = BMG_SER_N_CYCLE
      BMG_iPARMS(id_BMG2_SER_NCYCLE_TYPE) = BMG_SER_V_CYCLE
      BMG_iPARMS(id_BMG2_SER_MAX_ITERS)   = 1
      
C -------------------------------------
C     I/O Parameters
C -------------------------------------

      BMG_IOFLAG(iBMG2_SER_OUT_WSPACE_SIZE)  = .FALSE.
      BMG_IOFLAG(iBMG2_SER_OUT_ITERATIONS)   = .FALSE.

      BMG_IOFLAG(iBMG2_SER_OUT_TIME_SETUP)   = .FALSE.
      BMG_IOFLAG(iBMG2_SER_OUT_TIME_CYCLING) = .FALSE.
      BMG_IOFLAG(iBMG2_SER_OUT_TIME_TOTAL)   = .FALSE.

      BMG_IOFLAG(iBMG2_SER_BUG_RES_RELAX)    = .FALSE.
      BMG_IOFLAG(iBMG2_SER_BUG_RES_CG_SOLVE) = .FALSE.
      
      BMG_IOFLAG(iBMG2_SER_OUT_STOP_ERROR)   = .TRUE.

C -------------------------------------
C     Space requirements and pointers:
C -------------------------------------

      pSR=1
      pSI=1      

      InWORK(i_InWORK_SO)    = .FALSE.    ! use a separate array for SO
      InWORK(i_InWORK_U)     = .FALSE.    ! use a separate array for Q
      InWORK(i_InWORK_Q)     = .FALSE.    ! use a separate array for QF
      InWORK(i_InWORK_RES)   = .TRUE.     ! store RES in work array
      InWORK(i_InWORK_PCG_P) = .TRUE.     ! store PCG vector P in work array
      InWORK(i_InWORK_PCG_R) = .TRUE.     ! store PCG vector R in work array
      InWORK(i_InWORK_PCG_Z) = .TRUE.     ! store PCG vector Z in work array

      CALL BMG2_SER_SymStd_SETUP_PtrWork( 
     &                 L, M, BMG_iPARMS,
     &                 NOGm, NFm, NSOm, NIWORK, NRWORK,
     &                 pWORK, InWORK, pSR, pSI, BMG_IOFLAG 
     &                 ) 

      NOG  = BMG_iPARMS(id_BMG2_SER_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG2_SER_DIM_NF)
      NSO  = BMG_iPARMS(id_BMG2_SER_DIM_NSO)

      MAX_PCG_ITERS = BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_MAX_ITERS)

C --------------------------------------
C     Zero everyting: 
C     ( of particular importance is the 
C       rhs, if one does FMG! )
C --------------------------------------

      DO 10 i=1, NSO
         SO(i)=0
 10   CONTINUE

      DO 20 i=1, NF
         u(i)=0
         rhs(i)=0
 20   CONTINUE

      DO 30 i=1, NRWORK
         RWORK(i)=0
 30   CONTINUE

C -------------------------------------
C     Define the coefficient matrix
C -------------------------------------

      CALL PUTF(SO, rhs, L, M, hx, hy, 0)

C -------------------------------------
C     Solve the system with pcg
C -------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_pcg(
     &          L, M, 
     &          BMG_SER_PCG_iPARMS, BMG_SER_PCG_rPARMS,
     &          BMG_SER_PCG_RES, MAX_PCG_ITERS,
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &          u, RHS, NF, SO, NSO, NOG,
     &          pWORK, IWORK, NIWORK, RWORK, NRWORK
     &          )

C -------------------------------------
C     Compute the errors
C -------------------------------------

      CALL DIFF(L, M, hx, hy, u, MAXnorm, .FALSE.)
        
      WRITE(*,100) 'The max-norm of the difference between the exact ',
     &             'and the computed solution is ', MAXnorm


C ------------------------------------------------------------------------
C     Do it all again but: 
C       - skip the stencil definition
C       - use diagonal preconditionig
C ------------------------------------------------------------------------

      WRITE(*,*) 
      WRITE(*,*) ' **** PCG with diagonal preconditioning ... '
      WRITE(*,*) 

C ------------------------------
C     PCG parameters:
C ------------------------------

      BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_PRECON) 
     & = BMG_SER_PCG_PRECON_DIAG

C -------------------------------------
C     Zero the initial guess.
C -------------------------------------

      DO 40 i=1, NF
         u(i)=0
 40   CONTINUE

C -------------------------------------
C     Solve the system with pcg
C -------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_pcg(
     &          L, M, 
     &          BMG_SER_PCG_iPARMS, BMG_SER_PCG_rPARMS,
     &          BMG_SER_PCG_RES, MAX_PCG_ITERS,
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &          u, RHS, NF, SO, NSO, NOG,
     &          pWORK, IWORK, NIWORK, RWORK, NRWORK
     &          )

C -------------------------------------
C     Compute the errors
C -------------------------------------

      CALL DIFF(L, M, hx, hy, u, MAXnorm, .FALSE.)
        
      WRITE(*,100) 'The max-norm of the difference between the exact ',
     &             'and the computed solution is ', MAXnorm


C ------------------------------------------------------------------------
C     Do it all again but: 
C       - skip the stencil definition
C       - use boxmg n-cycle for preconditionig
C ------------------------------------------------------------------------

      WRITE(*,*) 
      WRITE(*,*) ' ***** PCG with BOXMG n-cycle preconditioning ... '
      WRITE(*,*) '       ( includes setup of boxmg operators ) '
      WRITE(*,*) 

C ------------------------------
C     PCG parameters:
C ------------------------------

      BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_PRECON)
     &  = BMG_SER_PCG_PRECON_BMG
      BMG_SER_PCG_iPARMS(id_BMG_SER_BMG_iPARMS0_SETUP)
     &  = BMG_SER_BMG_iPARMS0_SETUP_all

C -------------------------------------
C     Zero the initial guess.
C -------------------------------------

      DO 50 i=1, NF
         u(i)=0
 50   CONTINUE

C -------------------------------------
C     Solve the system with pcg
C -------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_pcg(
     &          L, M, 
     &          BMG_SER_PCG_iPARMS, BMG_SER_PCG_rPARMS,
     &          BMG_SER_PCG_RES, MAX_PCG_ITERS,
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &          u, RHS, NF, SO, NSO, NOG,
     &          pWORK, IWORK, NIWORK, RWORK, NRWORK
     &          )

C -------------------------------------
C     Compute the errors
C -------------------------------------

      CALL DIFF(L, M, hx, hy, u, MAXnorm, .FALSE.)
        
      WRITE(*,100) 'The max-norm of the difference between the exact ',
     &             'and the computed solution is ', MAXnorm


C ------------------------------------------------------------------------
C     Do it all again but: 
C       - skip the stencil definition
C       - use boxmg n-cycle for preconditionig
C       - skip boxmg setup
C ------------------------------------------------------------------------

      WRITE(*,*) 
      WRITE(*,*) ' ***** PCG with BOXMG n-cycle preconditioning ... '
      WRITE(*,*) '       ( skips setup of boxmg operators ) '
      WRITE(*,*) 

C ------------------------------
C     PCG parameters:
C ------------------------------

      BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_PRECON)   
     & = BMG_SER_PCG_PRECON_BMG
      BMG_SER_PCG_iPARMS(id_BMG_SER_BMG_iPARMS0_SETUP)
     & = BMG_SER_BMG_iPARMS0_SETUP_none

C -------------------------------------
C     Zero the initial guess.
C -------------------------------------

      DO 60 i=1, NF
         u(i)=0
 60   CONTINUE

C -------------------------------------
C     Solve the system with pcg
C -------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_pcg(
     &          L, M, 
     &          BMG_SER_PCG_iPARMS, BMG_SER_PCG_rPARMS,
     &          BMG_SER_PCG_RES, MAX_PCG_ITERS,
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &          u, RHS, NF, SO, NSO, NOG,
     &          pWORK, IWORK, NIWORK, RWORK, NRWORK
     &          )

C -------------------------------------
C     Compute the errors
C -------------------------------------

      CALL DIFF(L, M, hx, hy, u, MAXnorm, .FALSE.)
        
      WRITE(*,100) 'The max-norm of the difference between the exact ',
     &             'and the computed solution is ', MAXnorm


C =========================================================================

 100  FORMAT (1X,A,/,1X,A,E15.7,/)
 200  FORMAT (2X,A,$)

C =======================================

      END

      






