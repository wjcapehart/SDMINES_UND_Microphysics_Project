      PROGRAM MAIN

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This solves Example 1 given in Victor Bandy's BBMG interface
C     guide.  But it uses a direct call to the symmetric two-dimensional
C     BoxMG solver (BMG2_SER_SymStd_SOLVE_boxmg), and it uses a W-cycle
C     instead of F-cycles.  The discretization (putf), and the post
C     processing (diff) are the subroutines from the BBMG package, and
C     consquently they fall under the BBMG public-domain license and
C     copyright (see extras/bbmg/README).
C
C =======================================================================
C $license_flag$
C =======================================================================
C ==========================================================================

      IMPLICIT  NONE

      INCLUDE   'BMG_SER_workspace.h'
      INCLUDE   'BMG_SER_parameters.h'

C ---------------------------------------------
C     Dimensional Parameters:
C
      INTEGER     Lmax, Mmax, NXYCmax, NXYCmin
      PARAMETER ( Lmax=1900, Mmax=1900, NXYCmax=3, NXYCmin=3 )

C --------------------------------------------
C     Multigrid/Workspace Memory Allocation:
C
      INTEGER   NFmax, NOGmax, NOGmin, NSOmax
      PARAMETER ( NFmax=2758, NSOmax=9822, NOGmin=4, NOGmax=4 )

      INTEGER   NIWORK, NRWORK
      PARAMETER ( NIWORK=36, NRWORK=20206 )

      INTEGER   pSI, pSR, pWORK(NBMG_SER_pWORK)
      LOGICAL   InWORK(NBMG_SER_InWORK)

C --------------------------------------------
C     Variable Declarations:
C      
      INTEGER   BMG_iPARMS(NBMG_SER_iPARMS),
     &          i, IWORK(NIWORK), KG, L, M,
     &          NC, NCI, NCU, NF, NOG,
     &          NSO, NCBW, NSOR
      REAL*8    BMG_rPARMS(NBMG_SER_rPARMS),
     &          hx, hy, MAXnorm, RHS(NFmax),
     &          RWORK(NRWORK), SO(NSOmax), u(NFmax)
      LOGICAL   BMG_IOFLAG(NBMG_SER_IOFLAG)

C ===========================================================================


C ---------------------------------------------------------
C     Initialize the arguments
C ---------------------------------- 
C     Domain and grid:
C ----------------------------------

      WRITE(*,200) 'Enter the number of points in x ... '
      READ(*,*) L
      WRITE(*,200) 'Enter the number of points in y ... '
      READ(*,*) M

      IF (L.GT.Lmax) THEN
         WRITE(*,*) '***Error: memory allocation insufficient in x ***'
         STOP
      ELSE IF (M.GT.Mmax) THEN
         WRITE(*,*) '***Error: memory allocation insufficient in y ***'
         STOP
      ENDIF

      hx=1D0/L                  ! grid spacing in x
      hy=2D0/M                  ! grid spacing in y

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
      BMG_rPARMS(id_BMG2_SER_STOP_TOL) = 1D-6

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

      BMG_IOFLAG(iBMG2_SER_OUT_STOP_ERROR)   = .TRUE.
      
      !
      ! Verbose operator debugging
      !
      BMG_IOFLAG(iBMG2_SER_BUG_STENCIL_FG)  = .FALSE.
      BMG_IOFLAG(iBMG2_SER_BUG_STENCIL_CG)  = .FALSE.
      BMG_IOFLAG(iBMG2_SER_BUG_STENCIL_CG1) = .FALSE.
      BMG_IOFLAG(iBMG2_SER_OUT_STENCIL_TTY) = .FALSE.

      BMG_IOFLAG(iBMG2_SER_BUG_RESTRICT)     = .FALSE.
      BMG_IOFLAG(iBMG2_SER_OUT_RESTRICT_TTY) = .FALSE.

C -------------------------------------
C     Space requirements and pointers:
C -------------------------------------

      pSR=1
      pSI=1      

      !
      ! Start with nothing stored in the workspace array
      !
      DO i=1, NBMG_SER_InWORK
         InWORK(i) = .FALSE.
      ENDDO

      !
      ! Explicitly mark arrays to be in the workspace
      !
      InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      InWORK(i_InWORK_RES) = .TRUE.     ! store RES in work array

      CALL BMG2_SER_SymStd_SETUP_PtrWork( 
     &                 L, M, BMG_iPARMS,
     &                 NOGmax, NFmax, NSOmax, NIWORK, NRWORK,
     &                 pWORK, InWORK, pSR, pSI, BMG_IOFLAG
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

C -------------------------------------
C     Define the coefficient matrix
C -------------------------------------

      CALL PUTF(SO, rhs, L, M, hx, hy, 0)

C -------------------------------------
C     Solve the system
C -------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_boxmg( 
     &          L, M, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          u, rhs, RWORK(pWORK(ip_RES)), NF, NC, SO, NSO, 
     &          RWORK(pWORK(ip_SOR)), NSOR, RWORK(pWORK(ip_CI)), NCI, 
     &          RWORK(pWORK(ip_CSO)), RWORK(pWORK(ip_CU)), NCBW, NCU,
     &          IWORK(pWORK(ip_iG)), NOG, KG                         
     &          )

      IF (BMG_iPARMS(id_BMG2_SER_Err_Code) .ne. 0) THEN
         WRITE(*,*) "ERROR CODE = ", BMG_iPARMS(id_BMG2_SER_Err_Code)
         STOP
      END IF

C -------------------------------------
C     Compute the errors
C -------------------------------------

      CALL DIFF(L, M, hx, hy, u, MAXnorm, .FALSE.)
        
      WRITE(*,100) 'The max-norm of the difference between the exact ',
     &             'and the computed solution is ', MAXnorm


C ------------------------------------------------------------------------
C     Do it all again but skip the stencil definition:
C ------------------------------------------------------------------------

C ------------------------------
C     Multigird parameters:
C ------------------------------

      BMG_iPARMS(id_BMG2_SER_SETUP)    = BMG_SER_SETUP_none
      BMG_rPARMS(id_BMG2_SER_STOP_TOL) = 1D-6

C -------------------------------------
C     Zero the initial guess.
C -------------------------------------

      DO 40 i=1, NF
         u(i)=0
 40   CONTINUE

C -------------------------------------
C     Solve the system
C -------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_boxmg(
     &          L, M, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          u, rhs, RWORK(pWORK(ip_RES)), NF, NC, SO, NSO, 
     &          RWORK(pWORK(ip_SOR)), NSOR, RWORK(pWORK(ip_CI)), NCI, 
     &          RWORK(pWORK(ip_CSO)), RWORK(pWORK(ip_CU)), NCBW, NCU,
     &          IWORK(pWORK(ip_iG)), NOG, KG                         
     &          )

C -------------------------------------
C     Compute the errors
C -------------------------------------

      CALL DIFF(L, M, hx, hy, u, MAXnorm, .FALSE.)
        
      WRITE(*,100) 'The max-norm of the difference between the exact ',
     &             'and the computed solution is ', MAXnorm
 

C ------------------------------------------------------------------------
C     Do it all again but skip recompute the stencil and 
C     and change # of relaxations and do V-cycles
C ------------------------------------------------------------------------
C ------------------------------
C     Multigird parameters:
C ------------------------------

      BMG_iPARMS(id_BMG2_SER_SETUP)       = BMG_SER_SETUP_opers
      BMG_iPARMS(id_BMG2_SER_NCYCLE_TYPE) = BMG_SER_V_CYCLE
      BMG_rPARMS(id_BMG2_SER_STOP_TOL)    = 1D-6

C --------------------------------------
C     Zero everyting: 
C     ( of particular importance is the 
C       rhs, if one does FMG! )
C --------------------------------------

      DO 50 i=1, NSO
         SO(i)=0
 50   CONTINUE

      DO 60 i=1, NF
         u(i)=0
         rhs(i)=0
 60   CONTINUE

C -------------------------------------
C     Define the coefficient matrix
C -------------------------------------

      CALL PUTF(SO, rhs, L, M, hx, hy, 0)

C -------------------------------------
C     Solve the system
C -------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_boxmg( 
     &          L, M, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          u, rhs, RWORK(ip_RES), NF, NC, SO, NSO, 
     &          RWORK(pWORK(ip_SOR)), NSOR, RWORK(pWORK(ip_CI)), NCI, 
     &          RWORK(pWORK(ip_CSO)), RWORK(pWORK(ip_CU)), NCBW, NCU,
     &          IWORK(pWORK(ip_iG)), NOG, KG                         
     &          )

C -------------------------------------
C     Compute the errors
C -------------------------------------

      CALL DIFF(L, M, hx, hy, u, MAXnorm, .TRUE.)
        
      WRITE(*,100) 'The max-norm of the difference between the exact ',
     &             'and the computed solution is ', MAXnorm

C ------------------------------------------------------------------------
C     Do it all again but reset dimensions to test the
C     degenerate case of NOG=1.
C ------------------------------------------------------------------------
      
      WRITE(*,*) ' Testing degenerate case of NOG = 1'

      L=8
      M=16

      hx=1D0/L                  ! grid spacing in x
      hy=2D0/M                  ! grid spacing in y

C ------------------------------
C     Multigird parameters:
C  -----------------------------
C     NB: parameters affect
C         space requirements!!
C ------------------------------

      CALL BMG2_SER_SymStd_SETUP_parms( 
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG
     &          )

      BMG_iPARMS(id_BMG2_SER_STENCIL)     = BMG_SER_STENCIL_5pt
      BMG_iPARMS(id_BMG2_SER_MIN_NOG)     = 1
      BMG_iPARMS(id_BMG2_SER_CG_MIN_DIM)  = L
      BMG_rPARMS(id_BMG2_SER_STOP_TOL)    = 1D-6

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

      InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      InWORK(i_InWORK_RES) = .TRUE.     ! store RES in work array

      CALL BMG2_SER_SymStd_SETUP_PtrWork( 
     &                 L, M, BMG_iPARMS,
     &                 NOGmax, NFmax, NSOmax, NIWORK, NRWORK,
     &                 pWORK, InWORK, pSR, pSI, BMG_IOFLAG
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
C     ( of particular importance is the 
C       rhs, if one does FMG! )
C --------------------------------------

      DO 70 i=1, NSO
         SO(i)=0
 70   CONTINUE

      DO 80 i=1, NF
         u(i)=0
         rhs(i)=0
 80   CONTINUE

C -------------------------------------
C     Define the coefficient matrix
C -------------------------------------

      CALL PUTF(SO, rhs, L, M, hx, hy, 0)

C -------------------------------------
C     Solve the system
C -------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_boxmg(
     &          L, M, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          u, rhs, RWORK(pWORK(ip_RES)), NF, NC, SO, NSO, 
     &          RWORK(pWORK(ip_SOR)), NSOR, RWORK(pWORK(ip_CI)), NCI, 
     &          RWORK(pWORK(ip_CSO)), RWORK(pWORK(ip_CU)), NCBW, NCU,
     &          IWORK(pWORK(ip_iG)), NOG, KG                         )

C -------------------------------------
C     Compute the errors
C -------------------------------------

      CALL DIFF(L, M, hx, hy, u, MAXnorm, .TRUE.)
        
      WRITE(*,100) 'The max-norm of the difference between the exact ',
     &             'and the computed solution is ', MAXnorm


C ------------------------------------------------------------------------
C     Do it all again but reset dimensions to test the
C     special case of NOG=2.
C ------------------------------------------------------------------------
      
      WRITE(*,*) ' Testing special case of NOG = 2'

      L=8
      M=16

      hx=1D0/L                  ! grid spacing in x
      hy=2D0/M                  ! grid spacing in y

C ------------------------------
C     Multigird parameters:
C  -----------------------------
C     NB: parameters affect
C         space requirements!!
C ------------------------------

      CALL BMG2_SER_SymStd_SETUP_parms( 
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG
     &          )

      BMG_iPARMS(id_BMG2_SER_STENCIL)     = BMG_SER_STENCIL_5pt
      BMG_iPARMS(id_BMG2_SER_MIN_NOG)     = 1
      BMG_iPARMS(id_BMG2_SER_CG_MIN_DIM)  = L/2
      BMG_rPARMS(id_BMG2_SER_STOP_TOL)    = 1D-6

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

      InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      InWORK(i_InWORK_RES) = .TRUE.     ! store RES in work array

      CALL BMG2_SER_SymStd_SETUP_PtrWork( 
     &                 L, M, BMG_iPARMS,
     &                 NOGmax, NFmax, NSOmax, NIWORK, NRWORK,
     &                 pWORK, InWORK, pSR, pSI, BMG_IOFLAG
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
C     ( of particular importance is the 
C       rhs, if one does FMG! )
C --------------------------------------

      DO i=1, NSO
         SO(i)=0
      END DO

      DO i=1, NF
         u(i)=0
         rhs(i)=0
      END DO

C -------------------------------------
C     Define the coefficient matrix
C -------------------------------------

      CALL PUTF(SO, rhs, L, M, hx, hy, 0)

C -------------------------------------
C     Solve the system
C -------------------------------------

      CALL BMG2_SER_SymStd_SOLVE_boxmg(
     &          L, M, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          u, rhs, RWORK(pWORK(ip_RES)), NF, NC, SO, NSO, 
     &          RWORK(pWORK(ip_SOR)), NSOR, RWORK(pWORK(ip_CI)), NCI, 
     &          RWORK(pWORK(ip_CSO)), RWORK(pWORK(ip_CU)), NCBW, NCU,
     &          IWORK(pWORK(ip_iG)), NOG, KG                         )

C -------------------------------------
C     Compute the errors
C -------------------------------------

      CALL DIFF(L, M, hx, hy, u, MAXnorm, .TRUE.)
        
      WRITE(*,100) 'The max-norm of the difference between the exact ',
     &             'and the computed solution is ', MAXnorm

C ==========================================================================

 100  FORMAT (1X,A,/,1X,A,E15.7,/)
 200  FORMAT (2X,A,$)


C ===================================

      END

      






