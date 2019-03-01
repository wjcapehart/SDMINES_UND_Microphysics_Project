      PROGRAM MAIN

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This is an example that tests the ability of the BoxMG solver to
C     perform V-cycles in a symmetric fashion.  This is necessary for
C     using BoxMG as a preconditioner for a symmetric Krylov method such
C     as Conjugate Gradients (CG).  It uses the subroutine
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
C     MCYCL    - ??? set to 1
C     ALPHl    - truncation crap
C     ALPHm    - truncation crap
C     iPR      - output is generated if ipr.eq.0
C     iSKIP    - ISKIP.NE.0 initial setup phase is skipped
C     iSTRT2   - istrt for boxmg (usually -1)
C     iVW2     - ivw for boxmg  (usually 1)
C     MCYCL2   - mcycl for boxmg (usually 1)
C
C ===========================================================================

      IMPLICIT NONE

C ------------------------------------------------
C     Includes
C
      INCLUDE 'BMG_SER_constants.h'
      INCLUDE 'BMG_SER_workspace.h'
      INCLUDE 'BMG_SER_parameters.h'
      INCLUDE 'common2.h'

C ------------------------------------------------
C     Multigrid/Workspace Memory Allocation: 
C 
      INTEGER   NFm, NOGm, NSOm
      PARAMETER ( NFm=1799, NSOm=11876, NOGm=3 )

      INTEGER   NIWORK, NRWORK
      PARAMETER (  NIWORK=2153, NRWORK=24766 )

      INTEGER   NBMG_iWORK_PL, NBMG_rWORK_PL
      PARAMETER ( NBMG_iWORK_PL=16, NBMG_rWORK_PL=23408 )

      INTEGER   pSI, pSR, pWORK(NBMG_SER_pWORK)
      LOGICAL   InWORK(NBMG_SER_InWORK)

C -------------------------------------------------
C     Multigrid:  Variables
C
      INTEGER   BMG_iPARMS(NBMG_SER_iPARMS),
     &          BMG_iWORK_PL(NBMG_iWORK_PL),
     &          iPR, iSTRT, iSTRT2, iVERTEX,
     &          IWORK(NIWORK), NC, NCBW, NCI,
     &          NCU, NF, NOG, NOGc, NSO, NSOR,
     &          Nx, Ny, Nz
      REAL*8    BMG_rPARMS(NBMG_SER_rPARMS),
     &          BMG_rWORK_PL(NBMG_rWORK_PL),
     &          hx, hy, hz, prod1, prod2,
     &          Q(NFm), QF(NFm), RWORK(NRWORK),
     &          SO(NSOm), TOL_SAVE, 
     &          u(NFm), v(NFm), bu(NFm), bv(NFm)
      LOGICAL   BMG_IOFLAG(NBMG_SER_IOFLAG)

C --------------------------------------------------
C     Local Variables:
C
      INTEGER i

C ==========================================================================

C ==========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

C -----------------------------------
C     Default Cycling:
C -----------------------------------

      CALL BMG3_SER_SymStd_SETUP_parms( 
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG
     &          )

C -------------------------------------
C     The number of regions:
C -------------------------------------

      READ(*,*) iREG
C      WRITE(*,1000) 'The number of regions is ', iREG

C  ----------------------------------
C     Read data for each region:
C  ----------------------------------
      
      DO 10 i=1, iREG

         READ(*,*) x1(i),x2(i),y1(i),y2(i),z1(i),z2(i)
         READ(*,*) dix(i),diy(i),diz(i),si(i),fi(i)
C         WRITE(*,1010) 'Parameters for region ', i
C         WRITE(*,1020) 'Physical Domain:  [',x1(i),',',x2(i),'] x [',
C     &                  y1(i),',',y2(i),'] x [',
C     &                  z1(i),',',z2(i),']' 
C         WRITE(*,1030) '(Diagonal) Diffusion Tensor:',
C     &                 Dix(i), Diy(i), Diz(i)
C         WRITE(*,1040) 'Absorption Coefficient:', si(i)
C         WRITE(*,1040) 'Source Strength:', fi(i)

   10 CONTINUE

      READ(*,*) ibl,ibr,ibyb,ibyt,ibzb,ibzt
C      WRITE(*,1050) 'Boundary Condition Indeces'
C      WRITE(*,1060) 'x=Xs', ibl
C      WRITE(*,1060) 'x=Xf', ibr
C      WRITE(*,1060) 'y=Ys', ibyb
C      WRITE(*,1060) 'y=Yf', ibyt
C      WRITE(*,1060) 'z=Zs', ibzb
C      WRITE(*,1060) 'z=Zf', ibzt

C ----------------------------------
C     Discretization Parameters:
C ----------------------------------

      READ(*,*) Nx, Ny, Nz, iVERTEX, BMG_rPARMS(id_BMG3_SER_STOP_TOL)
      IF (iREG.EQ.1) THEN
         hx=(x2(1)-x1(1))/(Nx-iVERTEX)
         hy=(y2(1)-y1(1))/(Ny-iVERTEX)
         hz=(z2(1)-z1(1))/(Nz-iVERTEX)
      ELSE
         WRITE(*,*) 'ERROR: Code cannot handle this case!'
         STOP
      ENDIF
      
C      WRITE(*,1070) 'Global Grid Dimensions:', Nx, Ny, Nz
C      WRITE(*,1080) 'Global h-parameters: ', hx, hy, hz

C -----------------------------------
C     Cycle Parameters:
C -----------------------------------
      
      READ(*,*) BMG_iPARMS(id_BMG3_SER_STENCIL),
     &          BMG_iPARMS(id_BMG3_SER_NRELAX_DOWN),
     &          BMG_iPARMS(id_BMG3_SER_NRELAX_UP),
     &          BMG_iPARMS(id_BMG3_SER_NRELAX_FG),
     &          iSTRT,
     &          BMG_iPARMS(id_BMG3_SER_RELAX)
      READ(*,*) BMG_iPARMS(id_BMG3_SER_NCYCLE_TYPE),
     &          BMG_iPARMS(id_BMG3_SER_FMG_NNCYCLE)
      READ(*,*) ipr, BMG_iPARMS(id_BMG3_SER_SETUP)
      READ(*,*) iSTRT2,
     &          BMG_iPARMS(id_BMG2_SER_NCYCLE_TYPE),
     &          BMG_iPARMS(id_BMG2_SER_FMG_NNCYCLE)

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
      BMG_iPARMS(id_BMG3_SER_RELAX_SYM ) = BMG_SER_RELAX_SYM
      BMG_iPARMS(id_BMG2_SER_RELAX_SYM ) = BMG_SER_RELAX_SYM

C -------------------------------------
C     I/O Parameters
C -------------------------------------

      DO i=1, NBMG_SER_IOFLAG
         BMG_IOFLAG(i)=.FALSE.
      ENDDO

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

      InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      InWORK(i_InWORK_RES) = .TRUE.     ! store RES in work array

      CALL BMG3_SER_SymStd_SETUP_PtrWork(
     &                 Nx, Ny, Nz, BMG_iPARMS,
     &                 NOGm, NFm, NSOm, NIWORK, NRWORK,
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

      DO i=1, NRWORK
         RWORK(i)=rZERO
      END DO

      DO i=1, NIWORK
         IWORK(i)=iZERO
      END DO

C ------------------------------
C     Compute the Stencil:
C ------------------------------

      CALL PUTF( SO, QF, Q, Nx+2, Ny+2, Nz+2, hx, hy, hz, 0)

      DO 20 i=1, NF
         u(i)=0
         v(i)=0
         bu(i)=0
         bv(i)=0
         QF(i)=0
 20   CONTINUE

C -------------------------------------
C     initialize the vectors
C -------------------------------------

      CALL INITVECS (u, v, Nx+2, Ny+2, Nz+2)

C -------------------------------------
C     call boxmg 
C -------------------------------------

      TOL_SAVE = BMG_rPARMS(id_BMG3_SER_STOP_TOL)

      CALL BMG3_SER_SymStd_SOLVE_boxmg(
     &          Nx, Ny, Nz, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &          bu, u, RWORK(pWORK(ip_RES)), NF, NC, SO, NSO,
     &          RWORK(pWORK(ip_SOR)), NSOR, 
     &          RWORK(pWORK(ip_CI)), NCI, 
     &          RWORK(pWORK(ip_CSO)), RWORK(pWORK(ip_CU)), NCBW, NCU,
     &          IWORK(pWORK(ip_iG)), NOG, NOGc,
     &          BMG_iWORK_PL, NBMG_iWORK_PL, 
     &          BMG_rWORK_PL, NBMG_rWORK_PL
     &          )

C ==========================================================================
C     >>>>>>>>>>>>>>>>     END: SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

C ==========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: ISKIP SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================

      BMG_iPARMS(id_BMG3_SER_SETUP)    = BMG_SER_SETUP_none
      BMG_rPARMS(id_BMG3_SER_STOP_TOL) = TOL_SAVE

C ------------------------------
C     Solve the System:    
C ------------------------------

      CALL BMG3_SER_SymStd_SOLVE_boxmg(
     &          Nx, Ny, Nz, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &          bv, v, RWORK(pWORK(ip_RES)), NF, NC, SO, NSO,
     &          RWORK(pWORK(ip_SOR)), NSOR, 
     &          RWORK(pWORK(ip_CI)), NCI, 
     &          RWORK(pWORK(ip_CSO)), RWORK(pWORK(ip_CU)), NCBW, NCU,
     &          IWORK(pWORK(ip_iG)), NOG, NOGc,
     &          BMG_iWORK_PL, NBMG_iWORK_PL, 
     &          BMG_rWORK_PL, NBMG_rWORK_PL
     &          )

C ==========================================================================
C     >>>>>>>>>>>>>>>>     END: ISKIP-SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<<<<
C ==========================================================================


C ---------------------------------
C     calculate the inner products
C ---------------------------------

      CALL VECPROD (prod1, bu, v, Nx+2, Ny+2, Nz+2)
      CALL VECPROD (prod2, u, bv, Nx+2, Ny+2, Nz+2)

      WRITE(*,300)'The two inner products:',prod1,prod2
      WRITE(*,310)'Their difference:      ',prod1-prod2


C ==========================================================================

 300  FORMAT(/,2X,A,1X,F22.14,4X,F22.14,/)
 310  FORMAT(2X,A,1X,1P,E22.14,/)

C --------------------------------
C     Region Parameters
C
 1000 FORMAT(/,2X,A,I3)
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


