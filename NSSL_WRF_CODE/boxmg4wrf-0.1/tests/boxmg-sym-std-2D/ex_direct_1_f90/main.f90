      PROGRAM MAIN

! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     This example demonstrates the use of a allocatable arrays in
!     f90/f95 to allocate memory as required.  The problem and setup are
!     the same as ex_direct_1_Q, and just as in that case it uses the
!     subroutine BMG2_SymStd_SETUP_PtrWork to compute the size of the
!     required workspace, and the pointers into that workspace.  Then a
!     memory usage table is output, and these workspace arrays are
!     allocated. In addition, it uses BMG2_SymStd_SETUP_ProcGrid_file to
!     read data from a file and define the processor grid.
!
!     ** NOTE:
!
!     One important difference between this example and ex_direct_1_Q
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
!
!     -- entries of the diagonal diffusion tensor in region "i"
!
!     dix(i)   - entry D_xx
!     diy(i)   - entry D_yy
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
!     ibb      - index for boundary in y; bottom
!     ibt      - index for boundary in y; top
!
!     -- discretization dimensions
!
!     Nx       - number of points in x (excluding ghost points)
!     Ny       - number of points in y (excluding ghost points)
!
!     hx       - grid spacing in x (assumed constant)
!     hy       - grid spacing in y (assumed constant)
!
! ==========================================================================

      IMPLICIT NONE

! ------------------------------------------------
!     Includes
! ------------------------------------------------

      INCLUDE   'mpif.h'

      INCLUDE   'BMG_constants.h'
      INCLUDE   'BMG_workspace.h'
      INCLUDE   'BMG_parameters.h'

      INCLUDE   'common2.h'

! ------------------------------------------------
!     Multigrid/Workspace Memory Allocation: 
! ------------------------------------------------
 
      !
      ! Workspace pointers and logicals
      !
      INTEGER   BMG_pWORK(NBMG_pWORK)
      LOGICAL   BMG_InWORK(NBMG_InWORK)

      !
      ! Workspace pointer shift variables
      !
      INTEGER   pSI, pSR 

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
      ! Workspace: Generic
      !
      INTEGER,  allocatable, dimension(:) ::  BMG_iWORK
      REAL*8,   allocatable, dimension(:) ::  BMG_rWORK

      !
      ! Workspace: Generic
      !
      INTEGER,  allocatable, dimension(:) ::  BMG_iWORK_CS
      REAL*8,   allocatable, dimension(:) ::  BMG_rWORK_CS

      !
      ! Solution, Source/RHS, and Stencil
      !
      REAL*8,   allocatable, dimension(:) ::   Q, QF, SO

      !
      ! Miscellaneous
      !
      INTEGER   NC, NCBW, NCI, NCU, NF, NOG, NOGc, NSO, NSOR
      REAL*8    hx, hy, TOL, TOL_SAVE

      INTEGER   iGs, jGs, NGx, NGy, NLx, NLy

      INTEGER,  allocatable, dimension(:) ::  pMSG, pMSGSO, pLS

! -------------------------------------------------
!     MPI/MSG:  Variables
! -------------------------------------------------

      INTEGER   NBMG_MSG_iGRID, NBMG_MSG_iGRIDm
      
      INTEGER,  allocatable, dimension(:) ::   BMG_MSG_iGRID
      INTEGER   BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER   NMSGi, NMSGr
      INTEGER   MPI_IERROR, NProc, MPI_MyProc, MyProc
      INTEGER   MyProcI, MyProcJ, NProcI, NProcJ

! -------------------------------------------------
!     Local Variables:
! -------------------------------------------------
      INTEGER      i, isl, p_ProcGrid
      CHARACTER*80 GRIDFILEi, CYCLEFILEi, PDEFILEi

      INTEGER   NFm, NOGm, NSOm
      INTEGER   NBMG_iWORK, NBMG_rWORK
      INTEGER   NBMG_iWORK_CS, NBMG_rWORK_CS

      REAL*8    T1, T2, T_SOLVE_LOOP

! ==========================================================================

      CALL MPI_INIT( MPI_IERROR )

      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NProc, MPI_IERROR )      
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, MPI_MyProc, MPI_IERROR )

      !
      !  Note:
      !
      !  Internally BoxMG uses the MSG numbering of processors:
      !
      !        MyProc = MPI_MyProc + 1.
      !
      !  However, this is hidden from the user in BMG_MSG_iGRID,
      !  and hence, the MyProc argument in user level BoxMG 
      !  routines is expected to be MPI_MyProc.
      !
      MyProc = MPI_MyProc

! ==========================================================================
!     >>>>>>>>>>>>>>>>     BEGIN: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

! ---------------------------------- 
!     Read File Names:
! ----------------------------------

      IF ( MyProc.EQ.0 ) THEN 
         READ(*,'(A)') GRIDFILEi      ! Processor grid data
         READ(*,'(A)') CYCLEFILEi     ! BMG cycling parameters
         READ(*,'(A)') PDEFILEi       ! PDE parameters
      ENDIF

! ---------------------------------- 
!     Logical Grid:
! ----------------------------------

      BMG_iPARMS(ID_BMG2_POINTERS) = BMG_USE_pointers

      CALL BMG2_SymStd_SETUP_ProcGrid_file(                            &
                       GRIDFILEi, BMG_MSG_iGRID, BMG_MSG_pGRID,        &
                       BMG_iPARMS, NBMG_MSG_iGRIDm, NBMG_MSG_iGRID,    &
                       MPI_MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR   &
                       )

      ! allocate space for BMG_MSG_iGRID
      allocate(BMG_MSG_iGRID(1:NBMG_MSG_iGRID))

      BMG_iPARMS(id_BMG2_POINTERS) = BMG_NO_pointers
      NBMG_MSG_iGRIDm = NBMG_MSG_iGRID

      CALL BMG2_SymStd_SETUP_ProcGrid_file(                            &
                       GRIDFILEi, BMG_MSG_iGRID, BMG_MSG_pGRID,        &
                       BMG_iPARMS, NBMG_MSG_iGRIDm, NBMG_MSG_iGRID,    &
                       MPI_MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR   &
                       )

      NGx = BMG_MSG_iGRID(id_BMG_MSG_NGx)
      NGy = BMG_MSG_iGRID(id_BMG_MSG_NGy)

      NLx = BMG_MSG_iGRID(id_BMG_MSG_NLx)
      NLy = BMG_MSG_iGRID(id_BMG_MSG_NLy)

      iGs = BMG_MSG_iGRID(id_BMG_MSG_iGs)
      jGs = BMG_MSG_iGRID(id_BMG_MSG_jGs)

      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)

      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)

! -----------------------------------
!     Cycle Parameters:
! -----------------------------------

      CALL EX_SETUP_BMG_parms(                                          &
              CYCLEFILEi, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,           &
              MPI_MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR             &
              )

      BMG_IOFLAG(iBMG2_BUG_RES_RELAX)    =  .FALSE.
      BMG_IOFLAG(iBMG2_BUG_RES_CG_SOLVE) =  .FALSE.
      BMG_IOFLAG(iBMG2_OUT_TIME_SETUP)   =  .FALSE.
      BMG_IOFLAG(iBMG2_OUT_TIME_CYCLING) =  .FALSE.
      BMG_IOFLAG(iBMG2_OUT_TIME_TOTAL)   =  .FALSE.

!     To test NOG = 1 
!     
!     BMG_iPARMS(id_BMG2_MIN_NOG)     = 1
!     BMG_iPARMS(id_BMG2_CG_MIN_DIM)  = NLx

!     Tired of watching too many iterations ...
!
!     BMG_iPARMS(id_BMG2_MAX_ITERS)   = 1

! -----------------------------------
!     PDE Parameters:
! -----------------------------------

      CALL EX_SETUP_PDE_parms(                                          &
              PDEFILEi,                                                 &
              MPI_MyProc, NProc, MPI_COMM_WORLD, MPI_IERROR             &
              )


! ==========================================================================
!     >>>>>>>>>>>>>>>>     END: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

! ==========================================================================
!     >>>>>>>>>>>>>>>>     BEGIN: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

! ------------------------------
!     Workspace Allocation:
! ------------------------------

      pSI=1
      pSR=1

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

      BMG_iPARMS(id_BMG2_POINTERS) = BMG_USE_pointers

      CALL BMG2_SymStd_SETUP_PtrWork(                                   &
                NLx, NLy, NGx, NGy, iGs, jGs,                           &
                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,                     &
                MPI_MyProc, NProc, NProcI, NProcJ,                      &
                MPI_COMM_WORLD,                                         &
                NOGm, NFm, NSOm,                                        &
                NBMG_iWORK, NBMG_rWORK,                                 &
                NBMG_iWORK_CS, NBMG_rWORK_CS,                           &
                BMG_pWORK, BMG_InWORK, pSR, pSI                         &
                )

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

! ------------------------------
!    Report space requirements
! ------------------------------

      IF ( MyProc.EQ.0 ) THEN
         WRITE(*,*)
         WRITE(*,*) 'Memory allocation in MB'
         WRITE(*,*) '======================='
         WRITE(*,1210) 'SO           : ',NSO*8./(1024*1024)
         WRITE(*,1210) 'Q            : ',NF*8./(1024*1024)
         WRITE(*,1210) 'QF           : ',NF*8./(1024*1024)
         WRITE(*,1210) 'BMG_rWORK    : ',NBMG_rWORK*8./(1024*1024)
         WRITE(*,1210) 'BMG_iWORK    : ',NBMG_iWORK*8./(1024*1024)
         WRITE(*,1210) 'BMG_rWORK_CS : ',NBMG_rWORK*8./(1024*1024)
         WRITE(*,1210) 'BMG_iWORK_CS : ',NBMG_iWORK*8./(1024*1024)
         WRITE(*,1210) 'pMSG         : ',NBMG_pMSG*NOG*8./(1024*1024)
         WRITE(*,1210) 'pMSGSO       : ',NBMG_pMSG*NOG*8./(1024*1024)
         WRITE(*,1210) 'pLS          : ',NBMG_pLS*NOG*8./(1024*1024)
         WRITE(*,1210)
      ENDIF

! ------------------------------
!    Allocate the space
! ------------------------------

      allocate(SO(1:NSO))
      allocate(Q(1:NF))
      allocate(QF(1:NF))

      allocate(BMG_rWORK(1:NBMG_rWORK))
      allocate(BMG_iWORK(1:NBMG_iWORK))

      allocate(BMG_rWORK_CS(1:NBMG_rWORK_CS))
      allocate(BMG_iWORK_CS(1:NBMG_iWORK_CS))

      allocate(pMSG(1:NBMG_pMSG*NOG))
      allocate(pMSGSO(1:NBMG_pMSG*NOG))
      allocate(pLS(1:NBMG_pLS*NOG))

! ==========================================================================
!     >>>>>>>>>>>>>>>>     END: WORKSPACE SETUP   <<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================
      
! ==========================================================================
!     >>>>>>>>>>>>>>>>     BEGIN: SOLVE   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================


      !
      !  Loop to ensure enough timing data is generated.
      !

      T_SOLVE_LOOP = rZERO

      DO isl=1, 1   

         ! ------------------------------
         !     Zero the Stencil:
         ! ------------------------------
      
         DO i=1, NSO
            SO(i)=rZERO
         END DO

         DO i=1, NF
            Q(i)  =rZERO
            QF(i)=rZERO
         END DO
      
         ! ------------------------------
         !     Zero the workspace:
         ! ------------------------------

         DO i=1, NBMG_rWORK
            BMG_rWORK(i)=rZERO
         END DO

         DO i=1, NBMG_iWORK
            BMG_iWORK(i)=iZERO
         END DO

         DO i=1, NBMG_rWORK_CS
            BMG_rWORK_CS(i)=rZERO
         END DO

         DO i=1, NBMG_iWORK_CS
            BMG_iWORK_CS(i)=iZERO
         END DO

         ! ------------------------------
         !     Compute the Stencil:
         ! ------------------------------

         !
         !  Effective uniform grid spacing
         !
         hx = (xGf-xGs)/(NGx-1)
         hy = (yGf-yGs)/(NGy-1)
      
         CALL PUTF( MyProc, SO, QF, Q,                                  &
                    NLx, NLy, NGx, NGy,                                 &  
                    iGs, jGs, hx, hy, 0                                 &
                  )

         ! ------------------------------
         !     Solve the System:    
         ! ------------------------------

         !
         !  Override defaults (.FALSE.) for debugging:
         !
         BMG_IOFLAG(iBMG2_BUG_STENCIL_FG)  = .FALSE.
         BMG_IOFLAG(iBMG2_BUG_STENCIL_CG)  = .FALSE.
         BMG_IOFLAG(iBMG2_BUG_STENCIL_CG1) = .FALSE.
         BMG_IOFLAG(iBMG2_BUG_RESTRICT)    = .FALSE.

         TOL_SAVE = BMG_rPARMS(id_BMG2_STOP_TOL)

         CALL BMG2_SymStd_UTILS_zero_times(BMG_rPARMS)

         
         CALL MPI_Barrier(MPI_COMM_WORLD, MPI_ierror)
         T1 = MPI_Wtime()

         CALL BMG2_SymStd_SOLVE_boxmg(                                  &
                   NLx, NLy, NGx, NGy, iGs, jGs,                        &
                   BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,                  &
                   Q, QF, BMG_rWORK(BMG_pWORK(ip_RES)), NF, NC,         &
                   SO, NSO, BMG_rWORK(BMG_pWORK(ip_SOR)), NSOR,         &
                   BMG_rWORK(BMG_pWORK(ip_CI)), NCI,                    &
                   BMG_iWORK(BMG_pWORK(ip_iG)), NOG, NOGc,              &
                   BMG_iWORK_CS, NBMG_iWORK_CS,                         &
                   BMG_rWORK_CS, NBMG_rWORK_CS,                         &
                   BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,        &
                   BMG_iWORK(BMG_pWORK(ip_MSG)), NMSGi,                 &
                   pMSG, pMSGSO, pLS,                                   &
                   BMG_rWORK(BMG_pWORK(ip_MSG_BUF)), NMSGr,             &
                   MPI_COMM_WORLD                                       &
                   )

         CALL MPI_Barrier(MPI_COMM_WORLD, MPI_ierror)
         T2 = MPI_Wtime()

         IF ( MPI_MyProc.EQ.0 ) THEN
            WRITE(*,1200) 'Total time = ', T2 - T1
         ENDIF

         CALL BMG2_SymStd_OUTPUT_times(                                 &
                   BMG_rPARMS, MyProc, NProc, MPI_COMM_WORLD            &
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

      CALL MPI_FINALIZE(MPI_IERROR)

! ==========================================================================
!     >>>>>>>>>>>>>>>>     END: SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! ==========================================================================

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

! ==========================================================================

      END PROGRAM MAIN


