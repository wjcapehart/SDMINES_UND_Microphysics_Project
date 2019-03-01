      PROGRAM MAIN

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This example uses a direct call to the symmetric three-dimensional
C     BoxMG solver (BMG3_SER_SymStd_SOLVE_boxmg).  It uses the subroutine
C     BMG3_SER_SymStd_SETUP_PtrWork to compute the dimensions of the
C     required workspace and the necessary pointers into that workspace.
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

#include      'BMG_SER_constants.h'
      INCLUDE 'BMG_SER_workspace.h'
      INCLUDE 'BMG_SER_parameters.h'
#include      'common2.h'

C ------------------------------------------------
C     Multigrid/Workspace Memory Allocation: 
C ------------------------------------------------

      INTEGER   NFm, NOGm, NSOm
      PARAMETER ( NFm=60452, NSOm=339798, NOGm=5 )

      INTEGER   NBMG_iWORK, NBMG_rWORK
      PARAMETER (  NBMG_iWORK=40, NBMG_rWORK=436535 )

      INTEGER   NBMG_iWORK_PL, NBMG_rWORK_PL
      PARAMETER ( NBMG_iWORK_PL=15799, NBMG_rWORK_PL=3270324 )

      INTEGER   pSI, pSR, pWORK(NBMG_SER_pWORK)
      LOGICAL   InWORK(NBMG_SER_InWORK)

C -------------------------------------------------
C     Multigrid:  Variables
C -------------------------------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER   BMG_iPARMS(NBMG_SER_iPARMS)
      REAL*RKIND  BMG_rPARMS(NBMG_SER_rPARMS)
      LOGICAL   BMG_IOFLAG(NBMG_SER_IOFLAG)

      !
      ! Workspace: Plane Relaxation
      !
      INTEGER   BMG_iWORK_PL(NBMG_iWORK_PL)
      REAL*RKIND    BMG_rWORK_PL(NBMG_rWORK_PL)

      !
      ! Workspace: Generic
      !
      INTEGER   BMG_iWORK(NBMG_iWORK)
      REAL*RKIND    BMG_rWORK(NBMG_rWORK)

      !
      ! Solution, Source/RHS, and Stencil
      !
      REAL*RKIND    Q(NFm), QF(NFm), SO(NSOm)

      !
      ! Miscellaneous
      !
      INTEGER   NC, NCBW, NCI, NCU, NF, NOG, NOGc, NSO, NSOR
      REAL*RKIND    hx, hy, hz, TOL, TOL_SAVE

      INTEGER   IPR, iSTRT, iSTRT2, iVERTEX, Nx, Ny, Nz

C --------------------------------------------------
C     Local Variables:
C --------------------------------------------------

      INTEGER   i, IO_in, IO_keyboard
      CHARACTER FILE_in*20

C =========================================================================

C =========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
C =========================================================================

C -----------------------------------
C     Default cycle parameters:
C -----------------------------------

      CALL BMG3_SER_SymStd_SETUP_parms( 
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG
     &          )

C -------------------------------------
C     Enter the name of the input file:
C -------------------------------------

      IO_keyboard=5
      IO_in=5

      IF (IO_in.NE.IO_keyboard) THEN

         WRITE(*,1005) 'Enter the name of the input file ... '
         READ(IO_keyboard,*) FILE_in
         
         OPEN( IO_in, FILE=FILE_in, STATUS='OLD' ) 

      ENDIF

C -------------------------------------
C     The number of regions:
C -------------------------------------

      READ(IO_in,*) iREG
      WRITE(*,1000) 'The number of regions is ', iREG

C  ----------------------------------
C     Read data for each region:
C  ----------------------------------
      
      DO 10 i=1, iREG

         READ(IO_in,*) x1(i),x2(i),y1(i),y2(i),z1(i),z2(i)
         READ(IO_in,*) dix(i),diy(i),diz(i),si(i),fi(i)
         WRITE(*,1010) 'Parameters for region ', i
         WRITE(*,1020) 'Physical Domain:  [',x1(i),',',x2(i),'] x [',
     &                  y1(i),',',y2(i),'] x [',
     &                  z1(i),',',z2(i),']' 
         WRITE(*,1030) '(Diagonal) Diffusion Tensor:',
     &                 Dix(i), Diy(i), Diz(i)
         WRITE(*,1040) 'Absorption Coefficient:', si(i)
         WRITE(*,1040) 'Source Strength:', fi(i)

   10 CONTINUE

      READ(IO_in,*) ibl,ibr,ibyb,ibyt,ibzb,ibzt
      WRITE(*,1050) 'Boundary Condition Indeces'
      WRITE(*,1060) 'x=Xs', ibl
      WRITE(*,1060) 'x=Xf', ibr
      WRITE(*,1060) 'y=Ys', ibyb
      WRITE(*,1060) 'y=Yf', ibyt
      WRITE(*,1060) 'z=Zs', ibzb
      WRITE(*,1060) 'z=Zf', ibzt

C ----------------------------------
C     Discretization Parameters:
C ----------------------------------

      READ(IO_in,*) Nx, Ny, Nz, iVERTEX,
     &              BMG_rPARMS(id_BMG3_SER_STOP_TOL)
      IF (iREG.EQ.1) THEN
         hx=(x2(1)-x1(1))/(Nx-iVERTEX)
         hy=(y2(1)-y1(1))/(Ny-iVERTEX)
         hz=(z2(1)-z1(1))/(Nz-iVERTEX)
      ELSE
         WRITE(*,*) 'ERROR: Code cannot handle this case!'
         STOP
      ENDIF
      
      WRITE(*,1070) 'Global Grid Dimensions:', Nx, Ny, Nz
      WRITE(*,1080) 'Global h-parameters: ', hx, hy, hz

C -----------------------------------
C     Cycle Parameters:
C -----------------------------------
      
      READ(IO_in,*) BMG_iPARMS(id_BMG3_SER_STENCIL),
     &              BMG_iPARMS(id_BMG3_SER_NRELAX_DOWN),
     &              BMG_iPARMS(id_BMG3_SER_NRELAX_UP),
     &              BMG_iPARMS(id_BMG3_SER_NRELAX_FG),
     &              iSTRT,
     &              BMG_iPARMS(id_BMG3_SER_RELAX)
      READ(IO_in,*) BMG_iPARMS(id_BMG3_SER_NCYCLE_TYPE),
     &              BMG_iPARMS(id_BMG3_SER_FMG_NNCYCLE)
      READ(IO_in,*) iPR, BMG_iPARMS(id_BMG3_SER_SETUP)
      READ(IO_in,*) iSTRT2,
     &              BMG_iPARMS(id_BMG2_SER_NCYCLE_TYPE),
     &              BMG_iPARMS(id_BMG2_SER_FMG_NNCYCLE)

C ------------------------------------
C     Coarse-Grid Operators:
C ------------------------------------

      READ(IO_in,*) BMG_iPARMS(id_BMG3_SER_CG_TYPE),
     &              BMG_iPARMS(id_BMG3_SER_CG_CONSTRUCT)

      IF (IO_in.NE.IO_keyboard) THEN
         CLOSE(IO_in)
      ENDIF

C -----------------------------------
C     Override Defaults:
C -----------------------------------

      BMG_iPARMS(id_BMG3_SER_STOP_TEST) = BMG_SER_STOP_REL_RES_L2
      BMG_iPARMS(id_BMG3_SER_MAX_ITERS) = ABS(ISTRT)
      !
      IF ( ISTRT.GT.0 ) THEN
         BMG_iPARMS(id_BMG3_SER_CYCLE_CLASS) = BMG_SER_FMG_CYCLE
      ELSE IF ( ISTRT.LT.0 ) THEN 
         BMG_iPARMS(id_BMG3_SER_CYCLE_CLASS) = BMG_SER_N_CYCLE
      ELSE 
         WRITE(*,*) 'ERROR: ISTRT out of range ... '
         STOP
      ENDIF
      !
      IF ( ISTRT2.GT.0 ) THEN
         BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS) = BMG_SER_FMG_CYCLE
      ELSE IF ( ISTRT2.LT.0 ) THEN 
         BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS) = BMG_SER_N_CYCLE
      ELSE 
         WRITE(*,*) 'ERROR: ISTRT out of range ... '
         STOP
      ENDIF
      !
      BMG_iPARMS(id_BMG3_SER_RELAX_SYM ) = BMG_SER_RELAX_NONSYM
      BMG_iPARMS(id_BMG2_SER_RELAX_SYM ) = BMG_SER_RELAX_NONSYM

      !
      ! Dump the cycling parameters
      !
      WRITE(*,1090) 'Multigrid Parameters'
      WRITE(*,1100) 'Convergence Criteria',
     &              BMG_rPARMS(id_BMG3_SER_STOP_TOL)
      WRITE(*,1110) 'Stencil Index', 
     &              BMG_iPARMS(id_BMG3_SER_STENCIL)
      WRITE(*,1110) 'Up relaxations', 
     &              BMG_iPARMS(id_BMG3_SER_NRELAX_UP)
      WRITE(*,1110) 'Down relaxations', 
     &              BMG_iPARMS(id_BMG3_SER_NRELAX_DOWN)
      WRITE(*,1110) 'Fine Grid relaxations',
     &              BMG_iPARMS(id_BMG3_SER_NRELAX_FG)
      WRITE(*,1110) 'iSTRT', iSTRT
      WRITE(*,1110) 'Relaxation Index',
     &              BMG_iPARMS(id_BMG3_SER_RELAX)
      WRITE(*,1110) 'Cycle Recursion', 
     &              BMG_iPARMS(id_BMG3_SER_NCYCLE_TYPE)
      WRITE(*,1110) 'Mcycl', 
     &              BMG_iPARMS(id_BMG3_SER_FMG_NNCYCLE)
      WRITE(*,1110) 'Output Level', iPR
      WRITE(*,1110) 'CG construction', 
     &              BMG_iPARMS(id_BMG3_SER_SETUP)
      WRITE(*,1110) 'iSTRT2', iSTRT2
      WRITE(*,1110) 'Cycle Recursion',
     &              BMG_iPARMS(id_BMG2_SER_NCYCLE_TYPE)
      WRITE(*,1110) 'Mcycl2',
     &              BMG_iPARMS(id_BMG2_SER_FMG_NNCYCLE)
      WRITE(*,*)

C
C     Uncomment to test NOG=1 case, only
C     practical for small problems 
C
C     BMG_iPARMS(id_BMG3_SER_MIN_NOG)     = 1
C     BMG_iPARMS(id_BMG3_SER_CG_MIN_DIM)  = Nx

C
C     Tired of watching to many iterations ...
C     
C     BMG_iPARMS(id_BMG3_SER_MAX_ITERS)   = 1


C -------------------------------------
C     I/O Parameters
C -------------------------------------

      BMG_IOFLAG(iBMG3_SER_OUT_WSPACE_SIZE)  = .TRUE.
      BMG_IOFLAG(iBMG3_SER_OUT_ITERATIONS)   = .TRUE.

      BMG_IOFLAG(iBMG3_SER_OUT_TIME_SETUP)   = .TRUE.
      BMG_IOFLAG(iBMG3_SER_OUT_TIME_CYCLING) = .TRUE.
      BMG_IOFLAG(iBMG3_SER_OUT_TIME_TOTAL)   = .TRUE.

      BMG_IOFLAG(iBMG3_SER_BUG_RES_RELAX)    = .TRUE.
      BMG_IOFLAG(iBMG3_SER_BUG_RES_CG_SOLVE) = .TRUE.

      BMG_IOFLAG(iBMG3_SER_BUG_STENCIL_FG)   = .FALSE.
      BMG_IOFLAG(iBMG3_SER_BUG_STENCIL_CG)   = .FALSE.
      BMG_IOFLAG(iBMG3_SER_BUG_STENCIL_CG1)  = .FALSE.
      BMG_IOFLAG(iBMG3_SER_BUG_RESTRICT)     = .FALSE.

      BMG_IOFLAG(iBMG3_SER_OUT_STOP_ERROR)   = .TRUE.

C =========================================================================
C     >>>>>>>>>>>>>>>>     END: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
C =========================================================================

C =========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
C =========================================================================

C ------------------------------
C     Workspace Allocation:
C ------------------------------

      pSI=1
      pSR=1

      DO i=1, NBMG_SER_InWORK
         InWORK(i)=.FALSE.
      ENDDO

      InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      InWORK(i_InWORK_RES) = .TRUE.     ! store RES in work array

      CALL BMG3_SER_SymStd_SETUP_PtrWork(
     &                 Nx, Ny, Nz, BMG_iPARMS,
     &                 NOGm, NFm, NSOm, NBMG_iWORK, NBMG_rWORK,
     &                 NBMG_iWORK_PL, NBMG_rWORK_PL,
     &                 pWORK, InWork, pSR, pSI, BMG_IOFLAG
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

      CALL PUTF( SO, QF, Q, Nx+2, Ny+2, Nz+2, hx, hy, hz, 0)

C ------------------------------
C     Solve the System:    
C ------------------------------

      TOL_SAVE = BMG_rPARMS(id_BMG3_SER_STOP_TOL)

      CALL BMG3_SER_SymStd_SOLVE_boxmg(
     &          Nx, Ny, Nz, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &          Q, QF, BMG_rWORK(pWORK(ip_RES)), NF, NC, SO, NSO,
     &          BMG_rWORK(pWORK(ip_SOR)), NSOR, 
     &          BMG_rWORK(pWORK(ip_CI)), NCI, 
     &          BMG_rWORK(pWORK(ip_CSO)),
     &          BMG_rWORK(pWORK(ip_CU)), NCBW, NCU,
     &          BMG_iWORK(pWORK(ip_iG)), NOG, NOGc,
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

C =========================================================================
C     >>>>>>>>>>>>>>>>     END: SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
C =========================================================================

C =========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: ISKIP-SKIP   <<<<<<<<<<<<<<<<<<<<<<<<<
C =========================================================================

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
     &          Q, QF, BMG_rWORK(pWORK(ip_RES)), NF, NC, SO, NSO,
     &          BMG_rWORK(pWORK(ip_SOR)), NSOR, 
     &          BMG_rWORK(pWORK(ip_CI)), NCI, 
     &          BMG_rWORK(pWORK(ip_CSO)),
     &          BMG_rWORK(pWORK(ip_CU)), NCBW, NCU,
     &          BMG_iWORK(pWORK(ip_iG)), NOG, NOGc,
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

C =========================================================================
C     >>>>>>>>>>>>>>>>     END: ISKIP-SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
C =========================================================================



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


