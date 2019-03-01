      PROGRAM MAIN

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This example uses a direct call to the symmetric three-dimensional
C     BoxMG solver (BMG3_SER_SymStd_SOLVE_boxmg).  It uses subroutines
C     EX_SETUP_BMG_SER_parms, and EX_SETUP_PDE_parms to demonstrate the
C     typical setup of the BoxMG and PDE discretization parameters.  In
C     addition, it uses the subroutine BMG3_SER_SymStd_SETUP_PtrWork to
C     compute the dimensions of the required workspace and the necessary
C     pointers into that workspace.
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
C     -- multigrid parameters
C 
C     TOL      - iteration tolerance 
C
C     iFD      - iFD=1 => 7-points, iFD>1 => 27-points
C     iU       - number of relaxation sweeps after interpolation
C     iD       - number of relaxation sweeps after restriction
C     iM       - number of fine grid relaxation sweeps
C     iSTRT    - ISTRT>0 => FMG, ISTRT<0 => .NOT.FMG
C     iRELAX   - 1=>RB-GS, 5=>Alternating planes
C     iTAU     - ITAU=0 => compute an estimate of truncation error
C     iVW      - IVW=1=>V-cycles, IVW=2=>W-cycles
C     MCYCL    - NOT USED
C     iSKIP    - ISKIP.NE.0 initial setup phase is skipped
C     iSTRT2   - istrt for boxmg (usually -1)
C     iVW2     - ivw for boxmg  (usually 1)
C     MCYCL2   - mcycl for boxmg (usually 1)
C
C ==========================================================================

      IMPLICIT NONE

C ------------------------------------------------
C     Includes
C ------------------------------------------------

      INCLUDE 'BMG_SER_constants.h'
      INCLUDE 'BMG_SER_workspace.h'
      INCLUDE 'BMG_SER_parameters.h'

      INCLUDE 'common2.h'

C ------------------------------------------------
C     Multigrid/Workspace Memory Allocation: 
C ------------------------------------------------

      INCLUDE   'BMG_SER_dim_parms.h'
      !
      ! Workspace pointers and logicals
      !
      INTEGER   BMG_pWORK(NBMG_SER_pWORK)
      LOGICAL   BMG_InWORK(NBMG_SER_InWORK)

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
      INTEGER   BMG_iPARMS(NBMG_SER_iPARMS)
      REAL*8    BMG_rPARMS(NBMG_SER_rPARMS)
      LOGICAL   BMG_IOFLAG(NBMG_SER_IOFLAG)

      !
      ! Workspace: Plane Relaxation
      !
      INTEGER   BMG_iWORK_PL(NBMG_iWORK_PL)
      REAL*8    BMG_rWORK_PL(NBMG_rWORK_PL)

      !
      ! Workspace: Generic
      !
      INTEGER   BMG_iWORK(NBMG_iWORK)
      REAL*8    BMG_rWORK(NBMG_rWORK)

      !
      ! Solution, Source/RHS, and Stencil
      !
      REAL*8    Q(NFm), QF(NFm), SO(NSOm)

      !
      ! Miscellaneous
      !
      INTEGER   NC, NCBW, NCI, NCU, NF, NOG, NOGc, NSO, NSOR
      REAL*8    hx, hy, hz, TOL, TOL_SAVE

      INTEGER   iSTRT, iSTRT2, Nx, Ny, Nz

C --------------------------------------------------
C     Local Variables:
C --------------------------------------------------

      INTEGER i
      CHARACTER*80 GRIDFILEi, CYCLEFILEi, PDEFILEi

C ==========================================================================

C ==========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

C ---------------------------------- 
C     Read File Names:
C ----------------------------------

      READ(*,'(A)') GRIDFILEi   ! Processor grid data
      READ(*,'(A)') CYCLEFILEi  ! BMG cycling parameters
      READ(*,'(A)') PDEFILEi    ! PDE parameters

C -----------------------------------
C     Cycle Parameters:
C -----------------------------------

      CALL EX_SETUP_BMG_SER_parms( 
     &              CYCLEFILEi, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG
     &              )

C -----------------------------------
C     PDE Parameters:
C -----------------------------------

      CALL EX_SETUP_PDE_parms( PDEFILEi )

C ---------------------------------- 
C     Problem Size:
C ----------------------------------

      CALL EX_SETUP_GRID_size ( GRIDFILEi, Nx, Ny, Nz )

C ==========================================================================
C     >>>>>>>>>>>>>>>>     END: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

C ==========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

C ------------------------------
C     Workspace Allocation:
C ------------------------------

      pSI=1
      pSR=1

      DO i=1, NBMG_SER_InWORK
         BMG_InWORK(i)=.FALSE.
      ENDDO

      BMG_InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      BMG_InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      BMG_InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      BMG_InWORK(i_InWORK_RES) = .TRUE.     ! store RES in work array

      CALL BMG3_SER_SymStd_SETUP_PtrWork(
     &                 Nx, Ny, Nz, BMG_iPARMS,
     &                 NOGm, NFm, NSOm, NBMG_iWORK, NBMG_rWORK,
     &                 NBMG_iWORK_PL, NBMG_rWORK_PL,
     &                 BMG_pWORK, BMG_InWORK, pSR, pSI, BMG_IOFLAG
     &                 )


      NOG  = BMG_iPARMS(id_BMG3_SER_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG3_SER_DIM_NF)
      NC   = BMG_iPARMS(id_BMG3_SER_DIM_NC)
      NCI  = BMG_iPARMS(id_BMG3_SER_DIM_NCI)
      NSO  = BMG_iPARMS(id_BMG3_SER_DIM_NSO)
      NSOR = BMG_iPARMS(id_BMG3_SER_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG3_SER_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG3_SER_DIM_NCU)

C ------------------------------
C     Zero the Stencil:
C ------------------------------
      
      DO i=1, NSO
         SO(i)=rZERO
      END DO

      DO i=1, NF
         Q(i)=rZERO
         QF(i)=rZERO
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

      hx = (xGf-xGs)/(Nx-1)
      hy = (yGf-yGs)/(Ny-1)
      hz = (zGf-zGs)/(Nz-1)

      CALL PUTF( SO, QF, Q, 
     &           Nx+2, Ny+2, Nz+2, 
     &           hx, hy, hz, 0
     &         )

C ------------------------------
C     Solve the System:    
C ------------------------------

      !
      !  Override defaults (.FALSE.) for debugging:
      !
      BMG_IOFLAG(iBMG3_SER_BUG_STENCIL_FG)  = .TRUE.
      BMG_IOFLAG(iBMG2_SER_BUG_STENCIL_CG)  = .FALSE.
      BMG_IOFLAG(iBMG2_SER_BUG_STENCIL_CG1) = .FALSE.
      BMG_IOFLAG(iBMG2_SER_OUT_STENCIL_TTY) = .FALSE.

      !BMG_IOFLAG(iBMG2_SER_OUT_RESTRICT_TTY)  = .FALSE.

      TOL_SAVE = BMG_rPARMS(id_BMG3_SER_STOP_TOL)

      CALL BMG3_SER_SymStd_SOLVE_boxmg(
     &          Nx, Ny, Nz, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &          Q, QF, BMG_rWORK(BMG_pWORK(ip_RES)), NF, NC, SO, NSO,
     &          BMG_rWORK(BMG_pWORK(ip_SOR)), NSOR, 
     &          BMG_rWORK(BMG_pWORK(ip_CI)), NCI, 
     &          BMG_rWORK(BMG_pWORK(ip_CSO)),
     &          BMG_rWORK(BMG_pWORK(ip_CU)), NCBW, NCU,
     &          BMG_iWORK(BMG_pWORK(ip_iG)), NOG, NOGc,
     &          BMG_iWORK_PL, NBMG_iWORK_PL, 
     &          BMG_rWORK_PL, NBMG_rWORK_PL
     &          )

      TOL = BMG_rPARMS(id_BMG3_SER_STOP_TOL)

      IF ( BMG_iPARMS(id_BMG3_SER_STOP_TEST)
     &    .EQ.BMG_SER_STOP_REL_RES_L2 ) THEN
         !
         WRITE(*,1180) 'Final Relative Residual (l2-norm) = ', TOL
         !
      ELSEIF ( BMG_iPARMS(id_BMG3_SER_STOP_TEST)
     &        .EQ.BMG_SER_STOP_ABS_RES_L2 ) THEN
         !
         WRITE(*,1180) 'Final Absolute Residual (l2-norm) = ', TOL
         !
      ELSE
         WRITE(*,*) 'WARNING: Unknown stopping criteria: TOL = ', TOL
      ENDIF

C ==========================================================================
C     >>>>>>>>>>>>>>>>     END: SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

C ==========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: ISKIP-SKIP   <<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

      BMG_iPARMS(id_BMG3_SER_SETUP)    = BMG_SER_SETUP_none
      BMG_rPARMS(id_BMG3_SER_STOP_TOL) = TOL_SAVE

C ----------------------
C     Initial Guess
C ----------------------

      DO i=1, NF
         Q(i)=rZERO
      END DO

C ------------------------------
C     Solve the System:    
C ------------------------------

      CALL BMG3_SER_SymStd_SOLVE_boxmg(
     &          Nx, Ny, Nz, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &          Q, QF, BMG_rWORK(BMG_pWORK(ip_RES)), NF, NC, SO, NSO,
     &          BMG_rWORK(BMG_pWORK(ip_SOR)), NSOR, 
     &          BMG_rWORK(BMG_pWORK(ip_CI)), NCI, 
     &          BMG_rWORK(BMG_pWORK(ip_CSO)),
     &          BMG_rWORK(BMG_pWORK(ip_CU)), NCBW, NCU,
     &          BMG_iWORK(BMG_pWORK(ip_iG)), NOG, NOGc,
     &          BMG_iWORK_PL, NBMG_iWORK_PL, 
     &          BMG_rWORK_PL, NBMG_rWORK_PL
     &          )

      TOL = BMG_rPARMS(id_BMG3_SER_STOP_TOL)

      IF ( BMG_iPARMS(id_BMG3_SER_STOP_TEST)
     &    .EQ.BMG_SER_STOP_REL_RES_L2 ) THEN 
         !
         WRITE(*,1180) 'Final Relative Residual (l2-norm) = ', TOL
         !
      ELSEIF ( BMG_iPARMS(id_BMG3_SER_STOP_TEST)
     &        .EQ.BMG_SER_STOP_ABS_RES_L2 ) THEN
         !
         WRITE(*,1180) 'Final Absolute Residual (l2-norm) = ', TOL
         !
      ELSE
         WRITE(*,*) 'WARNING: Unknown stopping criteria: TOL = ', TOL
      ENDIF

C ==========================================================================
C     >>>>>>>>>>>>>>>>     END: ISKIP-SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

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

C ==========================================================================

      END


