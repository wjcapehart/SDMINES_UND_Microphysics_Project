      PROGRAM MAIN

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This solves a simple Poisson problem with the source (RHS)
C     computed from the discrete approximation to the solution.  So this
C     is not a real example, but is handy for debugging.  It uses a
C     direct direct call to the symmetric two-dimensional BoxMG solver
C     (BMG2_SymStd_SOLVE_boxmg).
C
C =======================================================================
C $license_flag$
C =======================================================================
C ==========================================================================

      IMPLICIT  NONE

      INCLUDE   'mpif.h'
      INCLUDE   'MSG.h'

      INCLUDE   'BMG_workspace.h'
      INCLUDE   'BMG_parameters.h'

C ------------------------------------------------
C     PDE Grid Parameters:
C ------------------------------------------------

      INTEGER   Lm, Mm, NXYCm
      PARAMETER ( Lm=8192, Mm=8192, NXYCm=3 )

      REAL*8    x(0:Lm+1), y(0:Mm+1), u_ex((Lm+2)*(Mm+2))

C ------------------------------------------------
C     Multigrid/Workspace Memory Allocation: 
C ------------------------------------------------

      INTEGER   NFm, NOGm, NSOm
      PARAMETER ( NFm=1341482, NSOm=4723702, NOGm=9 )

      INTEGER   NBMG_iWORK, NBMG_rWORK
      PARAMETER ( NBMG_iWORK=725638, NBMG_rWORK=6743430 )

      INTEGER   NBMG_iWORK_CS, NBMG_rWORK_CS
      PARAMETER ( NBMG_iWORK_CS=18, NBMG_rWORK_CS=692393 )

      INTEGER    NBMG_MSG_iGRIDm
      PARAMETER ( NBMG_MSG_iGRIDm=303 )

      !
      ! Workspace pointers and logicals
      !
      INTEGER   BMG_pWORK(NBMG_pWORK)
      LOGICAL   BMG_InWORK(NBMG_InWORK)

      !
      ! Workspace pointer shift variables
      !
      INTEGER   pSI, pSR 

C -------------------------------------------------
C     Multigrid:  Variables
C -------------------------------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER   BMG_iPARMS(NBMG_iPARMS)
      REAL*8    BMG_rPARMS(NBMG_rPARMS)
      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG)

      !
      ! Workspace: Coarse-grid Solve
      !
      INTEGER   BMG_iWORK_CS(NBMG_iWORK_CS)
      REAL*8    BMG_rWORK_CS(NBMG_rWORK_CS)

      !
      ! Workspace: Generic
      !
      INTEGER   BMG_iWORK(NBMG_iWORK)
      REAL*8    BMG_rWORK(NBMG_rWORK)

      !
      ! Solution, Source/RHS, and Stencil
      !
      REAL*8    u(NFm), rhs(NFm), SO(NSOm)

      !
      ! Miscellaneous
      !
      INTEGER   NC, NCBW, NCI, NCU, NF, NOG, NOGc, NSO, NSOR
      REAL*8    hx, hy, MaxNorm

      INTEGER   iGs, jGs, NGx, NGy, NLx, NLy

      INTEGER   pMSG(NBMG_pMSG*NOGm), pMSGSO(NBMG_pMSG*NOGm),
     &          pLS(NBMG_pLS*NOGm)

C -------------------------------------------------
C     MPI/MSG:  Variables
C -------------------------------------------------

      INTEGER   NBMG_MSG_iGRID
      
      INTEGER   BMG_MSG_iGRID(NBMG_MSG_iGRIDm),
     &          BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER   NMSGi, NMSGr
      INTEGER   MPI_IERROR, NProc, MyProc
      INTEGER   NProcI, NProcJ

C --------------------------------------------------
C     Local Variables:
C --------------------------------------------------

      INTEGER      i, j
      CHARACTER*80 GRIDFILEi

C ==========================================================================

      CALL MPI_INIT (MPI_IERROR)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, NProc, MPI_IERROR )      
      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, MyProc, MPI_IERROR )

C ==========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

C -----------------------------------
C     Cycle Parameters:
C -----------------------------------

      CALL BMG2_SymStd_SETUP_parms( 
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG
     &          )

      BMG_iPARMS(id_BMG2_STENCIL)   = BMG_STENCIL_5pt
      BMG_iPARMS(id_BMG2_RELAX  )   = BMG_GS_RB_x_y_lines
      BMG_iPARMS(id_BMG2_MAX_ITERS) = 10
      BMG_rPARMS(id_BMG2_STOP_TOL)  = 1D-10

C -------------------------------------
C     I/O Parameters
C -------------------------------------

      BMG_IOFLAG(iBMG2_OUT_WSPACE_SIZE)  = .TRUE.
      BMG_IOFLAG(iBMG2_OUT_ITERATIONS)   = .TRUE.

      BMG_IOFLAG(iBMG2_OUT_TIME_SETUP)   = .TRUE.
      BMG_IOFLAG(iBMG2_OUT_TIME_CYCLING) = .TRUE.
      BMG_IOFLAG(iBMG2_OUT_TIME_TOTAL)   = .TRUE.

      BMG_IOFLAG(iBMG2_BUG_RES_RELAX)    = .TRUE.
      BMG_IOFLAG(iBMG2_BUG_RES_CG_SOLVE) = .TRUE. !.TRUE.

C ---------------------------------------------------------
C     Initialize the arguments
C ---------------------------------- 
C     Domain and grid:
C ----------------------------------

C ---------------------------------- 
C     Read File Names:
C ----------------------------------

      IF ( MyProc.EQ.0 ) THEN 
         READ(*,'(A)') GRIDFILEi      ! Processor grid data
      ENDIF

C ---------------------------------- 
C     Logical Grid:
C ----------------------------------

      CALL BMG2_SymStd_SETUP_ProcGrid_file( 
     &                 GRIDFILEi, BMG_MSG_iGRID, BMG_MSG_pGRID,
     &                 BMG_iPARMS, NBMG_MSG_iGRIDm, NBMG_MSG_iGRID,
     &                 MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR 
     &                 )

      NGx = BMG_MSG_iGRID(id_BMG_MSG_NGx)
      NGy = BMG_MSG_iGRID(id_BMG_MSG_NGy)

      NLx = BMG_MSG_iGRID(id_BMG_MSG_NLx)
      NLy = BMG_MSG_iGRID(id_BMG_MSG_NLy)

      iGs = BMG_MSG_iGRID(id_BMG_MSG_iGs)
      JGs = BMG_MSG_iGRID(id_BMG_MSG_jGs)

      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)

C ==========================================================================
C     >>>>>>>>>>>>>>>>     END: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

C ==========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

C ------------------------------
C     Workspace Allocation:
C ------------------------------

      pSR=1
      pSI=1      

      !
      ! Start with nothing stored in the workspace array
      !
      DO i=1, NBMG_InWORK
         BMG_InWORK(i) = .FALSE.
      ENDDO

      BMG_InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      BMG_InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      BMG_InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      BMG_InWORK(i_InWORK_RES) = .TRUE.     ! store RES in work array

      CALL BMG2_SymStd_SETUP_PtrWork(
     &          NLx, NLy, NGx, NGy, iGs, jGs, 
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &          MyProc, NProc, NProcI, NProcJ, 
     &          MPI_COMM_WORLD,
     &          NOGm, NFm, NSOm,
     &          NBMG_iWORK, NBMG_rWORK,
     &          NBMG_iWORK_CS, NBMG_rWORK_CS,
     &          BMG_pWORK, BMG_InWORK, pSR, pSI
     &          )
      
      NOG  = BMG_iPARMS(id_BMG2_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG2_DIM_NF)
      NC   = BMG_iPARMS(id_BMG2_DIM_NC)
      NCI  = BMG_iPARMS(id_BMG2_DIM_NCI)
      NSO  = BMG_iPARMS(id_BMG2_DIM_NSO)
      NSOR = BMG_iPARMS(id_BMG2_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG2_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG2_DIM_NCU)

      NMSGi = BMG_iPARMS(id_BMG2_DIM_NMSGi)
      NMSGr = BMG_iPARMS(id_BMG2_DIM_NMSGr)

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

      hx=1D0/(NGx+1)                  ! grid spacing in x
      hy=1D0/(NGy+1)                  ! grid spacing in y

      DO i = 0, NLx+1
        x(i) = hx*(iGs-1+i)
      ENDDO

      DO j = 0, NLy+1
        y(j) = hy*(jGs-1+j)
      ENDDO

      CALL ExactSol( u_ex, NLx, NLy, iGs, jGs, x, y, hx, hy )

      CALL PUTF( SO, rhs, u_ex,
     &           NLx, NLy, iGs, jGs, NGx, NGy, x, y, hx, hy,
     &           MPI_COMM_WORLD )

C -------------------------------------
C     Solve the system
C -------------------------------------

      CALL BMG2_SymStd_SOLVE_boxmg(
     &          NLx, NLy, NGx, NGy, iGs, jGs, 
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &          u, rhs, BMG_rWORK(BMG_pWORK(ip_RES)), NF, NC, 
     &          SO, NSO, BMG_rWORK(BMG_pWORK(ip_SOR)), NSOR, 
     &          BMG_rWORK(BMG_pWORK(ip_CI)), NCI,
     &          BMG_iWORK(BMG_pWORK(ip_iG)), NOG, NOGc,
     &          BMG_iWORK_CS, NBMG_iWORK_CS,
     &          BMG_rWORK_CS, NBMG_rWORK_CS,
     &          BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &          BMG_iWORK(BMG_pWORK(ip_MSG)), NMSGi, pMSG, pMSGSO, pLS,
     &          BMG_rWORK(BMG_pWORK(ip_MSG_BUF)), NMSGr,
     &          MPI_COMM_WORLD
     &          )

C -------------------------------------
C     Compute the errors
C -------------------------------------

      CALL DIFF( 
     &           NLx, NLy, u, u_ex, MaxNorm, MPI_COMM_WORLD
     &         )
      
      IF (MyProc.eq.1) THEN
         WRITE(*,100)'The max-norm of the difference between the exact '
     &        ,'and the computed solution is ', MAXnorm
      ENDIF

C ---------------------------------
C     Close MPI
C ---------------------------------

      CALL MPI_FINALIZE(MPI_IERROR)

C ==========================================================================
C     >>>>>>>>>>>>>>>>     END: SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

 100  FORMAT (1X,A,/,1X,A,E15.7,/)
 200  FORMAT (2X,A,$)

C ===================================

      END

      






