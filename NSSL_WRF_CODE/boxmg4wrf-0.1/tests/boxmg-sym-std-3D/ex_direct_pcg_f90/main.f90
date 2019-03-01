      PROGRAM MAIN

! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     This example demonstrates the use of a allocatable arrays in
!     f90/f95 to allocate memory as required for a direct call to 
!     the Preconditioned Conjugate Gradients (PCG) solver,
!     BMG3_SymStd_SOLVE_pcg.  The problem and setup are the same
!     as ex_direct_pcg, and just as in that case it uses the
!     subroutine BMG3_SymStd_SETUP_PtrWork to compute the size of the
!     required workspace, and the pointers into that workspace.  Then a
!     memory usage table is output, and these workspace arrays are
!     allocated. In addition, it uses BMG3_SymStd_SETUP_ProcGrid_file to
!     read data from a file and define the processor grid.
!
!     ** NOTE:
!
!     One important difference between this example and ex_direct_pcg
!     is that the source is set to zero here, and a random initial guess
!     is used. This facilitates accurate estimates of the asymptotic 
!     convergence rate and more accurate timing estimates (i.e., with
!     larger number of iterations).
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   VARIABLES:
!  --------------------
! 
!     -- boundaries of region "i" for a tensor product grid
!
!     x1(i)    - minimum in x
!     x2(i)    - maximum in x
!     y1(i)    - minimum in y
!     y2(i)    - maximum in y
!     z1(i)    - minimum in z
!     z2(i)    - maximum in z
!
!     -- entries of the diagonal diffusion tensor in region "i"
!
!     dix(i)   - entry D_xx
!     diy(i)   - entry D_yy
!     diz(i)   - entry D_zz
!
!     -- other parameters (assumed constant) for region "i"
!
!     si(i)    - zeroth order term (zero if there isn't one)
!     fi(i)    - source term
!
!     -- boundary condition indeces for the global domain
!
!     ibl      - index for boundary in x; left
!     ibr      - index for boundary in x; right
!     ibyb     - index for boundary in y; bottom
!     ibyt     - index for boundary in y; top
!     ibzb     - index for boundary in z; bottom
!     ibzt     - index for boundary in z; top
!
!     -- discretization dimensions
!
!     Nx       - number of points in x (excluding ghost points)
!     Ny       - number of points in y (excluding ghost points)
!     Nz       - number of points in z (excluding ghost points)
!
!     hx       - grid spacing in x (assumed constant)
!     hy       - grid spacing in y (assumed constant)
!     hz       - grid spacing in z (assumed constant)
!
! ===========================================================================

      IMPLICIT NONE

! ------------------------------------------------
!     Includes
! ------------------------------------------------

      INCLUDE   'mpif.h'

      INCLUDE   'BMG_constants.h'
      INCLUDE   'BMG_workspace.h'
      INCLUDE   'BMG_parameters.h'
      INCLUDE   'BMG_PCG_parameters.h'  ! new include for pcg

      INCLUDE   'common2.h'

! ------------------------------------------------
!     Multigrid/Workspace Memory Allocation: 
! ------------------------------------------------
 
      INTEGER   NFm, NOGm, NSOm
      INTEGER   NIWORKm, NRWORKm
      INTEGER   NIWORK_PLm, NRWORK_PLm
      INTEGER   pSI, pSR, pWORK(NBMG_pWORK)
      LOGICAL   InWORK(NBMG_InWORK)

! -------------------------------------------------
!     Multigrid:  Variables
! -------------------------------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER   BMG_iPARMS(NBMG_iPARMS)
      REAL*8    BMG_rPARMS(NBMG_rPARMS)
      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG)

      !
      ! Workspace: Plane Relaxation
      !
      INTEGER, allocatable, dimension(:) :: IWORKPL
      REAL*8, allocatable, dimension(:) ::  RWORK_PL

      INTEGER, allocatable, dimension(:) :: BMG_iWORK_CS
      REAL*8, allocatable, dimension(:) ::  BMG_rWORK_CS

      !
      ! Workspace: Generic
                            
      INTEGER,allocatable, dimension(:) ::  IWORKG
      REAL*8,allocatable, dimension(:) :: RWORK

      !
      ! Solution, Source/RHS, and Stencil
      !
      REAL*8 , allocatable, dimension(:):: u, rhs, SO



       INTEGER  NBMG_iWORK_CS, NBMG_rWORK_CS


      !
      ! Miscellaneous
      !
      INTEGER   NC, NCBW, NCI, NCU, NF, NOG, NOGc, NSO, NSOR
      REAL*8    hx, hy, hz, TOL, TOL_SAVE

      INTEGER   iGs, jGs, kGs, NGx, NGy, NGz, NLx, NLy, NLz

      INTEGER, allocatable, dimension(:) :: pMSG, pMSGSO, pLS

! -------------------------------------------------
!     PCG:  Variables
! -------------------------------------------------

      INTEGER   MAX_PCG_ITERSm, MAX_PCG_ITERS
      PARAMETER ( MAX_PCG_ITERSm = 10000 )

      INTEGER   BMG_PCG_iPARMS(NBMG_PCG_iPARMS)
      REAL*8    BMG_PCG_rPARMS(NBMG_PCG_rPARMS)
      REAL*8    BMG_PCG_RES(MAX_PCG_ITERSm)
      
      
! -------------------------------------------------
!     MPI/MSG:  Variables
! -------------------------------------------------

      INTEGER   NBMG_MSG_iGRIDm, NBMG_MSG_iGRID
      
      INTEGER, allocatable, dimension(:) :: BMG_MSG_iGRID
      INTEGER  BMG_MSG_pGRID(NBMG_MSG_pGRID)



      INTEGER   NMSGim, NMSGrm
      INTEGER   MPI_IERROR, NProc, MPI_MyProc, MyProc
      INTEGER   MyProcI, MyProcJ, MyProcK, NProcI, NProcJ, NProcK

! --------------------------------------------------
!     Local Variables:
!
      INTEGER       i, isl, p_ProcGrid
      CHARACTER*40  GRIDFILEi, CYCLEFILEi, PDEFILEi, PCGFILEi

      REAL*8    T2, T1, T_SOLVE_LOOP

! ===========================================================================

      CALL MPI_INIT( MPI_IERROR )

      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NProc, MPI_IERROR )      
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, MPI_MyProc, MPI_IERROR )

      MyProc = MPI_MyProc

! ===========================================================================
!     >>>>>>>>>>>>>>>>     BEGIN: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
! ===========================================================================

      IF ( MyProc.EQ.0 ) THEN
         READ(*,'(A)') GRIDFILEi      ! Processor grid data
         READ(*,'(A)') CYCLEFILEi     ! BMG cycling parameters
         READ(*,'(A)') PDEFILEi       ! PDE parameters
	 READ(*,'(A)') PCGFILEi       ! PCG parameters
      ENDIF



! ---------------------------------------------------------
!     Initialize the arguments
! ---------------------------------- 
!     Domain and grid:
! ----------------------------------

      BMG_iPARMS(id_BMG3_POINTERS) = BMG_USE_pointers
      CALL BMG3_SymStd_SETUP_ProcGrid_file(                        &
                      GRIDFILEi, BMG_MSG_iGRID, BMG_MSG_pGRID,     &
                      BMG_iPARMS, NBMG_MSG_iGRIDm, NBMG_MSG_iGRID, &
                      MPI_MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR &
                      )
      allocate(BMG_MSG_iGRID(1:NBMG_MSG_iGRID))
      BMG_iPARMS(id_BMG3_POINTERS) = BMG_NO_pointers
      NBMG_MSG_iGRIDm = NBMG_MSG_iGRID
      CALL BMG3_SymStd_SETUP_ProcGrid_file(                        &
                      GRIDFILEi, BMG_MSG_iGRID, BMG_MSG_pGRID,     &
                      BMG_iPARMS, NBMG_MSG_iGRIDm, NBMG_MSG_iGRID, &
                      MPI_MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR &
                      )
  

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

! -----------------------------------
!     Cycle Parameters:
! -----------------------------------

      CALL EX_SETUP_BMG_parms(  &
                   CYCLEFILEi, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, &
                   MPI_MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR &
                   )

!
!     To test NOG = 1 
!     
!     BMG_iPARMS(id_BMG3_MIN_NOG)     = 1
!     BMG_iPARMS(id_BMG3_CG_MIN_DIM)  = NLx

!
!     Tired of watching to many iterations ...
!
!     BMG_iPARMS(id_BMG3_MAX_ITERS)   = 1

! -----------------------------------
!     PCG Parameters:
! -----------------------------------

      CALL EX_SETUP_PCG_parms(  &
                   PCGFILEi, BMG_PCG_iPARMS, BMG_PCG_rPARMS, &
                   MPI_MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR &
                   )
      
      MAX_PCG_ITERS = BMG_PCG_iPARMS(id_BMG_PCG_MAX_ITERS)

! -----------------------------------
!     Multigrid parameters necessary
!     to work as preconditioner
! -----------------------------------

      BMG_iPARMS(id_BMG3_RELAX)       = BMG_GS_RB_point
      BMG_iPARMS(id_BMG3_RELAX_SYM)   = BMG_RELAX_SYM
      BMG_iPARMS(id_BMG3_CYCLE_CLASS) = BMG_N_CYCLE
      BMG_iPARMS(id_BMG3_NCYCLE_TYPE) = BMG_V_CYCLE
      BMG_iPARMS(id_BMG3_MAX_ITERS)   = 1

! -----------------------------------
!     PDE Parameters:
! -----------------------------------

      CALL EX_SETUP_PDE_parms( &
                   PDEFILEi, &
                   MPI_MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR &
                    )

      hx = (xGf-xGs)/(NGx-1)
      hy = (yGf-yGs)/(NGy-1)
      hz = (zGf-zGs)/(NGz-1)

! ---------------------------------------
!     I/O Parameters:  Override defaults
! ---------------------------------------

      BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE)  = .FALSE. !.TRUE.
      BMG_IOFLAG(iBMG3_OUT_ITERATIONS)   = .FALSE.

      BMG_IOFLAG(iBMG3_OUT_TIME_SETUP)   = .FALSE. !.TRUE.
      BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING) = .FALSE.
      BMG_IOFLAG(iBMG3_OUT_TIME_TOTAL)   = .FALSE. !TRUE.

      BMG_IOFLAG(iBMG3_BUG_RES_RELAX)    = .FALSE.
      BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) = .FALSE.

      BMG_IOFLAG(iBMG3_BUG_STENCIL_FG)   = .FALSE.
      BMG_IOFLAG(iBMG3_BUG_STENCIL_CG)   = .FALSE.
      BMG_IOFLAG(iBMG3_BUG_STENCIL_CG1)  = .FALSE.
      BMG_IOFLAG(iBMG3_BUG_RESTRICT)     = .FALSE.
  
      BMG_PCG_iPARMS(id_BMG_PCG_OUTPUT)   = BMG_PCG_OUT_ITS
      BMG_PCG_iPARMS(id_BMG_PCG_OUT_FREQ) = 1

! ===========================================================================
!     >>>>>>>>>>>>>>>>     END: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
! ===========================================================================

!C ===========================================================================
!C     >>>>>>>>>>>>>>>>     BEGIN: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<<
!C ===========================================================================

!C ------------------------------
!C     Workspace Allocation:
!C ------------------------------

      pSI=1
      pSR=1

      !
      ! Start with nothing stored in the workspace array
      !
      DO i=1, NBMG_InWORK
         InWORK(i) = .FALSE.
      ENDDO

      InWORK(i_InWORK_SO)    = .FALSE.    ! use a separate array for SO
      InWORK(i_InWORK_U)     = .FALSE.    ! use a separate array for Q
      InWORK(i_InWORK_Q)     = .FALSE.    ! use a separate array for QF
      InWORK(i_InWORK_RES)   = .TRUE.     ! store RES in work array
      InWORK(i_InWORK_PCG_P) = .TRUE.     ! store PCG vector P in work array
      InWORK(i_InWORK_PCG_R) = .TRUE.     ! store PCG vector R in work array
      InWORK(i_InWORK_PCG_Z) = .TRUE.     ! store PCG vector Z in work array

      BMG_iPARMS(id_BMG3_POINTERS) = BMG_USE_pointers

      CALL BMG3_SymStd_SETUP_PtrWork(                   &
                NLx, NLy, NLz, NGx, NGy, NGz,           &
                iGs, jGs, kGs,                          &
                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,     &
                MyProc, NProc, NProcI, NProcJ, NProcK,  &
                MPI_COMM_WORLD,                         &
                NOGm, NFm, NSOm,                        &
                NIWORKm, NRWORKm,                       &
                NIWORK_PLm, NRWORK_PLm,                 &
                NBMG_iWORK_CS, NBMG_rWORK_CS,           &
                pWORK, InWORK, pSR, pSI                 &
                )

      NOG  = BMG_iPARMS(id_BMG3_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG3_DIM_NF)
      NC   = BMG_iPARMS(id_BMG3_DIM_NC)
      NCI  = BMG_iPARMS(id_BMG3_DIM_NCI)
      NSO  = BMG_iPARMS(id_BMG3_DIM_NSO)
      NSOR = BMG_iPARMS(id_BMG3_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG3_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG3_DIM_NCU)

      NMSGim = BMG_iPARMS(id_BMG3_DIM_NMSGi)
      NMSGrm = BMG_iPARMS(id_BMG3_DIM_NMSGr)

!C ------------------------------
!C    Report space requirements
!C ------------------------------

      IF (MyProc .eq. 0) THEN
         WRITE(*,*)
         WRITE(*,*) 'Memory allocation in MB'
         WRITE(*,*) '======================='
         WRITE(*,1210) 'SO           : ',NSO*8./(1024*1024)
         WRITE(*,1210) 'Q            : ',NF*8./(1024*1024)
         WRITE(*,1210) 'QF           : ',NF*8./(1024*1024)
         WRITE(*,1210) 'RWORK        : ',NRWORKm*8./(1024*1024)
         WRITE(*,1210) 'IWORK        : ',NIWORKm*8./(1024*1024)
         WRITE(*,1210) 'pMSG         : ',NBMG_pMSG*NOG*8./(1024*1024)
         WRITE(*,1210) 'pMSGSO       : ',NBMG_pMSG*NOG*8./(1024*1024)
         WRITE(*,1210) 'pLS          : ',NBMG_pLS*NOG*8./(1024*1024)
         WRITE(*,1210) 'NIWORK_PL    : ',NIWORK_PLm*8./(1024*1024)
         WRITE(*,1210) 'RWORK_PL     : ',NRWORK_PLm*8./(1024*1024)
         WRITE(*,1210) 'BMG_iWORK_CS : ',NBMG_iWORK_CS*8./(1024*1024)
         WRITE(*,1210) 'BMG_rWORK_CS : ',NBMG_rWORK_CS*8./(1024*1024)
         WRITE(*,*)
      ENDIF

!C ------------------------------
!C    Allocate the space
!C ------------------------------

      allocate(SO(1:NSO))
      allocate(u(1:NF))
      allocate(rhs(1:NF))
      allocate(RWORK(1:NRWORKm))
      allocate(IWORKG(1:NIWORKm))
      allocate(IWORKPL(1:NIWORK_PLm))
      allocate(RWORK_PL(1:NRWORK_PLm))
      allocate(pMSG(1:NBMG_pMSG*NOG))
      allocate(pMSGSO(1:NBMG_pMSG*NOG))
      allocate(pLS(1:NBMG_pLS*NOG))
      allocate(BMG_iWORK_CS(1:NBMG_iWORK_CS))
      allocate(BMG_rWORK_CS(1:NBMG_rWORK_CS))


!C ===========================================================================
!C     >>>>>>>>>>>>>>>>     END: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<<<<
!C ===========================================================================
      
!C ===========================================================================
!C     >>>>>>>>>>>>>>>>     BEGIN: SOLVE   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!C ===========================================================================

      !
      !  Loop to ensure enough timing data is generated.
      !

      T_SOLVE_LOOP = rZERO

      DO isl=1, 3   

         ! ------------------------------
         !   Zero the Stencil:
         ! ------------------------------
         
         DO i=1, NSO
            SO(i)=rZERO
         END DO

         DO i=1, NF
            u(i)=rZERO
            rhs(i)=rZERO
         END DO

         DO i=1, NRWORKm
            RWORK(i)=rZERO
         END DO
         
         DO i=1, NIWORKm
            IWORKG(i)=iZERO
         END DO

         ! ------------------------------
         !     Compute the Stencil:
         ! ------------------------------

         CALL PUTF( SO, rhs, u,                                         &
                    NLx, NLy, NLz, NGx, NGy, NGz,                       &
                    iGs, jGs, kGs, hx, hy, hz, 0                        &
                  )

         ! ------------------------------
         !     Solve the System:    
         ! ------------------------------

         CALL BMG3_SymStd_UTILS_zero_times(BMG_rPARMS)

         CALL MPI_Barrier(MPI_COMM_WORLD, MPI_ierror)
         T1 = MPI_Wtime()

         CALL BMG3_SymStd_SOLVE_pcg(                                    &
                   NLx, NLy, NLz, NGx, NGy, NGz,                        &   
                   iGs, jGs, kGs,                                       &
                   BMG_PCG_iPARMS, BMG_PCG_rPARMS,                      &
                   BMG_PCG_RES, MAX_PCG_ITERS,                          &
                   BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,                  & 
                   u, rhs, NF, SO, NSO, NOG,                            &
                   pWORK, IWORKG, NIWORKm, RWORK, NRWORKm,              &
                   IWORKPL, NIWORK_PLm,                                 &
                   RWORK_PL, NRWORK_PLm,                                &
                   BMG_iWORK_CS, NBMG_iWORK_CS,                         &
                   BMG_rWORK_CS, NBMG_rWORK_CS,                         &
                   IWORKG(pWORK(ip_MSG)), NMSGim,                       &
                   pMSG, pMSGSO, pLS,                                   &
                   BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,        &
                   NProc, RWORK(pWORK(ip_MSG_BUF)),                     &
                   NMSGrm, MPI_COMM_WORLD                               &
                   )

         CALL MPI_Barrier(MPI_COMM_WORLD, MPI_ierror)
         T2 = MPI_Wtime()

         IF (MyProc .eq. 0) THEN
            WRITE(*,1200) 'Total Time = ', T2 - T1
         ENDIF
         
         CALL BMG3_SymStd_OUTPUT_times(                                  &
                   BMG_rPARMS, MyProc, NProc, MPI_COMM_WORLD             &
                   )

         !
         !  Update total time
         !  
         T_SOLVE_LOOP = T_SOLVE_LOOP + T2 - T1       

         !
         !  Run for a minimum time of 50 seconds
         !
         IF ( T_SOLVE_LOOP .GT. 50 ) EXIT

      ENDDO

! ===========================================================================
!     >>>>>>>>>>>>>>>>     END: SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
! ===========================================================================

      CALL MPI_FINALIZE(MPI_IERROR)

! ===========================================================================
      
! --------------------------------
!     Region Parameters
!
 1000 FORMAT(/,2X,A,I3)
 1005 FORMAT(/,2X,A,$)
 1010 FORMAT(/,2X,A,I3)
 1020 FORMAT(2X,A,6(F5.2,A))
 1030 FORMAT(2X,A,T35,1P,3(E14.7))
 1040 FORMAT(2X,A,T35,1P,E14.7)

! --------------------------------
!     Boundary Condition Indeces
!
 1050 FORMAT(/,2X,30("="),/,4X,A,/,2X,30("="),/)
 1060 FORMAT(6X,A,T20,I2)

! --------------------------------
!     Grid Parameters
!
 1070 FORMAT(/,2X,A,3(1X,I3))
 1080 FORMAT(2X,A,1P,3(1X,E12.5))

! ---------------------------------
!     Cycle Parameters
!
 1090 FORMAT(/,2X,40("="),/,4X,A,/,2X,40("="),/)
 1100 FORMAT(4X,A,T26,1P,E12.5)
 1110 FORMAT(4X,A,T26,1P,I5)

 1180 FORMAT(/,1X,A,1P,E14.7,/)

 1200 FORMAT(/,/,1X,32("-"),/,2X,A,1X,F16.4,/,1X,32("-"),/)

 1210 FORMAT(4X,A,1X,F18.6,1X,'MBytes')

! ===========================================================================

      END PROGRAM MAIN


