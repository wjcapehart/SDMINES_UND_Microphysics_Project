      PROGRAM MAIN

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This example demonstrates a direct call to the Preconditioned 
C     Conjugate Gradients (PCG) solver, BMG3_SymStd_SOLVE_pcg, with
C     BoxMG as the preconditioner.  It demonstrates a modular
C     setup of a discretized PDE and the BoxMG solve.  It uses the 
C     subroutine BMG3_SymStd_SETUP_PtrWork to compute the dimensions
C     of the required workspace and the necessary pointers into that 
C     workspace. In addition, it uses BMG3_SymStd_SETUP_ProcGrid_file
C     to read data from a file and define the processor grid.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   VARIABLES:
C  --------------------
C 
C     -- boundaries of region "i" for a tensor product grid
C
C     x1(i)    - minimum in x
C     x2(i)    - maximum in x
C     y1(i)    - minimum in y
C     y2(i)    - maximum in y
C     z1(i)    - minimum in z
C     z2(i)    - maximum in z
C
C     -- entries of the diagonal diffusion tensor in region "i"
C
C     dix(i)   - entry D_xx
C     diy(i)   - entry D_yy
C     diz(i)   - entry D_zz
C
C     -- other parameters (assumed constant) for region "i"
C
C     si(i)    - zeroth order term (zero if there isn't one)
C     fi(i)    - source term
C
C     -- boundary condition indeces for the global domain
C
C     ibl      - index for boundary in x; left
C     ibr      - index for boundary in x; right
C     ibyb     - index for boundary in y; bottom
C     ibyt     - index for boundary in y; top
C     ibzb     - index for boundary in z; bottom
C     ibzt     - index for boundary in z; top
C
C     -- discretization dimensions
C
C     Nx       - number of points in x (excluding ghost points)
C     Ny       - number of points in y (excluding ghost points)
C     Nz       - number of points in z (excluding ghost points)
C
C     hx       - grid spacing in x (assumed constant)
C     hy       - grid spacing in y (assumed constant)
C     hz       - grid spacing in z (assumed constant)
C
C ==========================================================================

      IMPLICIT NONE

C ------------------------------------------------
C     Includes
C ------------------------------------------------

      INCLUDE   'mpif.h'
      INCLUDE   'MSG.h'

      INCLUDE   'BMG_constants.h'
      INCLUDE   'BMG_workspace.h'
      INCLUDE   'BMG_parameters.h'
      INCLUDE   'BMG_PCG_parameters.h'  ! new include for pcg

      INCLUDE   'common2.h'

C ------------------------------------------------
C     Multigrid/Workspace Memory Allocation: 
C ------------------------------------------------
  
      INCLUDE   'BMG_dim_parms.h'

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
      ! Workspace: Plane Relaxation
      !
      INTEGER   BMG_iWORK_PL(NBMG_iWORK_PL)
      REAL*8    BMG_rWORK_PL(NBMG_rWORK_PL)

      !
      ! Workspace: Plane Relaxation
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
      INTEGER   NC, NCBW, NCI, NCU, NF, NOG, NSO, NSOR
      REAL*8    hx, hy, hz, TOL

      INTEGER   iGs, jGs, kGs, NGx, NGy, NGz, NLx, NLy, NLz

      INTEGER   pMSG(NBMG_pMSG*NOGm), pMSGSO(NBMG_pMSG*NOGm),
     &          pLS(NBMG_pLS*NOGm)

C -------------------------------------------------
C     PCG:  Variables
C -------------------------------------------------

      INTEGER   MAX_PCG_ITERSm, MAX_PCG_ITERS
      PARAMETER ( MAX_PCG_ITERSm = 1000 )

      INTEGER   BMG_PCG_iPARMS(NBMG_PCG_iPARMS)
      REAL*8    BMG_PCG_rPARMS(NBMG_PCG_rPARMS)
      REAL*8    BMG_PCG_RES(MAX_PCG_ITERSm)
      
C -------------------------------------------------
C     MPI/MSG:  Variables
C -------------------------------------------------

      INTEGER   NBMG_MSG_iGRID
      
      INTEGER   BMG_MSG_iGRID(NBMG_MSG_iGRIDm),
     &          BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER   NMSGim, NMSGrm, 
     &          MPI_IERROR, NProc, MPI_MyProc, MyProc,
     &          MyProcI, MyProcJ, MyProcK, NProcI, NProcJ, NProcK

C --------------------------------------------------
C     Local Variables:
C
      INTEGER      i, p_ProcGrid
      CHARACTER*40 GRIDFILEi, CYCLEFILEi, PDEFILEi

C ===========================================================================

      CALL MPI_INIT( MPI_IERROR )

      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NProc, MPI_IERROR )      
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, MPI_MyProc, MPI_IERROR )

      MyProc = MPI_MyProc

C ===========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
C ===========================================================================

C ---------------------------------- 
C     Read File Names:
C ----------------------------------

      IF ( MyProc.EQ.0 ) THEN 
         READ(*,'(A)') GRIDFILEi      ! Processor grid data
         READ(*,'(A)') CYCLEFILEi     ! BMG cycling parameters
         READ(*,'(A)') PDEFILEi       ! PDE parameters
      ENDIF

C -----------------------------------
C     PCG Parameters:
C -----------------------------------

      BMG_PCG_iPARMS(id_BMG_PCG_NMG_CYCLES) = 1
      BMG_PCG_iPARMS(id_BMG_PCG_STOP_TEST)  = BMG_PCG_STOP_REL_RES_L2
      BMG_PCG_iPARMS(id_BMG_PCG_PRECON)     = BMG_PCG_PRECON_BMG
      BMG_PCG_iPARMS(id_BMG_PCG_MAX_ITERS)  = 10
      BMG_PCG_iPARMS(id_BMG_PCG_BMG_SETUP)  = BMG_PCG_BMG_SETUP_all

      BMG_PCG_iPARMS(id_BMG_PCG_OUTPUT)     = BMG_PCG_OUT_ALL
      BMG_PCG_iPARMS(id_BMG_PCG_OUT_FREQ)   = 1
      BMG_PCG_rPARMS(id_BMG_PCG_STOP_TOL)   = 1D-5

C -----------------------------------
C     BoxMG Cycle Parameters:
C -----------------------------------

      CALL EX_SETUP_BMG_parms( 
     &              CYCLEFILEi, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &              MPI_MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR
     &              )

      !
      !  Multigrid parameters necessary
      !  to work as preconditioner
      !
      BMG_iPARMS(id_BMG3_RELAX)       = BMG_GS_RB_point
      BMG_iPARMS(id_BMG3_RELAX_SYM)   = BMG_RELAX_SYM
      BMG_iPARMS(id_BMG3_CYCLE_CLASS) = BMG_N_CYCLE
      BMG_iPARMS(id_BMG3_NCYCLE_TYPE) = BMG_V_CYCLE
      BMG_iPARMS(id_BMG3_MAX_ITERS)   = 1

C -----------------------------------
C     PDE Parameters:
C -----------------------------------

      PDEFILEi='input/in.pde-poisson'
      CALL EX_SETUP_PDE_parms(
     &              PDEFILEi,
     &              MPI_MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR
     &              )

C ---------------------------------- 
C     Logical Grid:
C ----------------------------------

      CALL BMG3_SymStd_SETUP_ProcGrid_file( 
     &                 GRIDFILEi, BMG_MSG_iGRID, BMG_MSG_pGRID,
     &                 BMG_iPARMS, NBMG_MSG_iGRIDm, NBMG_MSG_iGRID,
     &                 MPI_MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR 
     &                 )

      NGx = BMG_MSG_iGRID(id_BMG_MSG_NGx)
      NGy = BMG_MSG_iGRID(id_BMG_MSG_NGy)
      NGz = BMG_MSG_iGRID(id_BMG_MSG_NGz)

      NLx = BMG_MSG_iGRID(id_BMG_MSG_NLx)
      NLy = BMG_MSG_iGRID(id_BMG_MSG_NLy)
      NLz = BMG_MSG_iGRID(id_BMG_MSG_NLz)

      iGs = BMG_MSG_iGRID(id_BMG_MSG_iGs)
      jGs = BMG_MSG_iGRID(id_BMG_MSG_jGs)
      kGs = BMG_MSG_iGRID(id_BMG_MSG_kGs)

      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)

C ===========================================================================
C     >>>>>>>>>>>>>>>>     END: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
C ===========================================================================

C ===========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
C ===========================================================================

C ------------------------------
C     Workspace Allocation:
C ------------------------------

      pSI=1
      pSR=1

      !
      ! Start with nothing stored in the workspace array
      !
      DO i=1, NBMG_InWORK
         BMG_InWORK(i) = .FALSE.
      ENDDO

      BMG_InWORK(i_InWORK_SO)    = .FALSE. ! use a separate array for SO
      BMG_InWORK(i_InWORK_U)     = .FALSE. ! use a separate array for Q
      BMG_InWORK(i_InWORK_Q)     = .FALSE. ! use a separate array for QF
      BMG_InWORK(i_InWORK_RES)   = .TRUE.  ! store RES in work array
      BMG_InWORK(i_InWORK_PCG_P) = .TRUE.  ! store PCG vector P in work array
      BMG_InWORK(i_InWORK_PCG_R) = .TRUE.  ! store PCG vector R in work array
      BMG_InWORK(i_InWORK_PCG_Z) = .TRUE.  ! store PCG vector Z in work array

      CALL BMG3_SymStd_SETUP_PtrWork(
     &          NLx, NLy, NLz, NGx, NGy, NGz, iGs, jGs, kGs,
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          MyProc, NProc, NProcI, NProcJ, NProcK,
     &          MPI_COMM_WORLD,
     &          NOGm, NFm, NSOm,
     &          NBMG_iWORK, NBMG_rWORK,
     &          NBMG_iWORK_PL, NBMG_rWORK_PL,
     &          NBMG_iWORK_CS, NBMG_rWORK_CS,
     &          BMG_pWORK, BMG_InWORK, pSR, pSI
     &          )

      NOG  = BMG_iPARMS(id_BMG3_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG3_DIM_NF)
      NC   = BMG_iPARMS(id_BMG3_DIM_NC)
      NCI  = BMG_iPARMS(id_BMG3_DIM_NCI)
      NSO  = BMG_iPARMS(id_BMG3_DIM_NSO)
      NSOR = BMG_iPARMS(id_BMG3_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG3_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG3_DIM_NCU)

      MAX_PCG_ITERS = BMG_PCG_iPARMS(id_BMG_PCG_MAX_ITERS)

      NMSGim = BMG_iPARMS(id_BMG3_DIM_NMSGi)
      NMSGrm = BMG_iPARMS(id_BMG3_DIM_NMSGr)

C ------------------------------
C     Zero the Stencil:
C ------------------------------
      
      DO i=1, NSO
         SO(i)=rZERO
      END DO

      DO i=1, NF
         u(i)=rZERO
         rhs(i)=rZERO
      END DO

      DO i=1, NBMG_rWORK
         BMG_rWORK(i)=rZERO
      END DO

      DO i=1, NBMG_iWORK
         BMG_iWORK(i)=iZERO
      END DO

C ------------------------------
C     Compute the Stencil:
C ------------------------------

      hx = (xGf-xGs)/(NGx-1)
      hy = (yGf-yGs)/(NGy-1)
      hz = (zGf-zGs)/(NGz-1)

      CALL PUTF( SO, rhs, u, 
     &           NLx, NLy, NLz, NGx, NGy, NGz, 
     &           iGs, jGs, kGs, hx, hy, hz, 0 
     &          )

C ------------------------------
C     Solve the System:    
C ------------------------------

      CALL BMG3_SymStd_SOLVE_pcg(
     &          NLx, NLy, NLz, NGx, NGy, NGz, 
     &          iGs, jGs, kGs,
     &          BMG_PCG_iPARMS, BMG_PCG_rPARMS,
     &          BMG_PCG_RES, MAX_PCG_ITERSm,
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &          u, rhs, NF, SO, NSO, NOG,
     &          BMG_pWORK, BMG_iWORK, NBMG_iWORK,
     &          BMG_rWORK, NBMG_rWORK,
     &          BMG_iWORK_PL, NBMG_iWORK_PL, 
     &          BMG_rWORK_PL, NBMG_rWORK_PL,
     &          BMG_iWORK_CS, NBMG_iWORK_CS, 
     &          BMG_rWORK_CS, NBMG_rWORK_CS,
     &          BMG_iWORK(BMG_pWORK(ip_MSG)), NMSGim,
     &          pMSG, pMSGSO, pLS,
     &          BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &          NProc, BMG_rWORK(BMG_pWORK(ip_MSG_BUF)), 
     &          NMSGrm, MPI_COMM_WORLD
     &          )

      TOL = BMG_rPARMS(id_BMG3_STOP_TOL)

C ===========================================================================
C     >>>>>>>>>>>>>>>>     END: SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
C ===========================================================================
      
      CALL MPI_FINALIZE(MPI_IERROR)
      STOP

C --------------------------------
C     Region Parameters
C
 1000 FORMAT(/,2X,A,I3)
 1005 FORMAT(/,2X,A,$)
 1010 FORMAT(/,2X,A,I3)
 1020 FORMAT(2X,A,6(F5.2,A))
 1030 FORMAT(2X,A,T35,1P,3(E14.7))
 1040 FORMAT(2X,A,T35,1P,E14.7)

C --------------------------------
C     Boundary Condition Indeces
C
 1050 FORMAT(/,2X,30("="),/,4X,A,/,2X,30("="),/)
 1060 FORMAT(6X,A,T20,I2)

C --------------------------------
C     Grid Parameters
C
 1070 FORMAT(/,2X,A,3(1X,I3))
 1080 FORMAT(2X,A,1P,3(1X,E12.5))

C ---------------------------------
C     Cycle Parameters
C
 1090 FORMAT(/,2X,40("="),/,4X,A,/,2X,40("="),/)
 1100 FORMAT(4X,A,T26,1P,E12.5)
 1110 FORMAT(4X,A,T26,1P,I5)

 1180 FORMAT(/,1X,A,1P,E14.7,/)

C ===========================================================================

      END


