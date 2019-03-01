      PROGRAM MAIN

C ============================================================================
C
C     Example 2:  A direct call to BMG3D 
C
C   ---------------
C   DESCRIPTION:
C   ---------------
C
C   This example uses a direct call to the symmetric three-dimensional
C   black box multigrid code.  It utilizes INTEGER and REAL work arrays
C   with pointer indexing based on the include file of "boxspace3.h".
C   It checks the symmetry of the BOXMG solver, when symmetric cycling
C   is switched on.   
C
C   ---------------------
C   HISTORY:  
C   ---------------------
C
C   Written 2000/28/02 M.Berndt, the code was mostly taken from 
C                      ex_direct_1 and modified a bit
C 
C   --------------
C   VARIABLES:
C   --------------
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
C     iCOEF    - ICOEF=0 to get operator-induced variational coarsening
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
      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'common2.h'

C ------------------------------------------------
C     Multigrid/Workspace Memory Allocation: 
C 
      INTEGER   NFm, NOGm, NSOm
      PARAMETER ( NFm=1799, NSOm=11876, NOGm=3 )

      INTEGER   NIWORK, NRWORK
      PARAMETER (  NIWORK=2153, NRWORK=24766 )

      INTEGER   NInWORK, NpWORK
      PARAMETER ( NInWORK=4, NpWORK=24 )

      INTEGER   pSI, pSR, pWORK(NpWORK)
      LOGICAL   InWORK(NInWORK)


C -------------------------------------------------
C     Multigrid:  Variables
C
      INTEGER   iCOEF, iD, iFD, iM, iPR, iRELAX, IRELAX_SYM, iSKIP,
     &          ISTOP, iSTRT, iSTRT2, iU, iVERTEX, iVW, iVW2,
     &          IWORK(NIWORK), MCYCL, MCYCL2, NC, NCBW, NCBW2,
     &          NCI, NCI3, NCPL2, NCU, NCU2, NF, NF2, NOG, NOGc,
     &          NSO, NSO3, NSOR, NSR3, Nx, Ny, Nz
      REAL*8    hx, hy, hz, prod1, prod2, Q(NFm), QF(NFm),
     &          RWORK(NRWORK), SO(NSOm), TOL, TOL_SAVE, 
     &          u(NFm), v(NFm), bu(NFm), bv(NFm)
      LOGICAL   IOFLAG(NIOFLAG)

C --------------------------------------------------
C     Local Variables:
C
      INTEGER i

C ---------------------------------------------------
C     Common Blocks:
C
      INTEGER  kdebug, linp, lout, ltty
      COMMON   /iocomm/ kdebug, linp, lout, ltty

C ===========================================================================

C ===========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: DATA INPUT    <<<<<<<<<<<<<<<<<<<<<<<<<
C ===========================================================================

      DATA kdebug /2/, linp /5/, lout /6/, ltty /6/

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

      READ(*,*) Nx, Ny, Nz, iVERTEX, TOL
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
      
      READ(*,*) iFD, iU, iD, iM, iSTRT, iRELAX, iCOEF
      READ(*,*) iVW, Mcycl
      READ(*,*) iPR,iSKIP
      READ(*,*) iSTRT2, iVW2, Mcycl2

      IRELAX_SYM = BMG_RELAX_SYM
      ISTOP = BMG_STOP_REL_RES_L2

C      WRITE(*,1090) 'Multigrid Parameters'
C      WRITE(*,1100) 'Convergence Criteria', TOL
C      WRITE(*,1110) 'Stencil Index', iFD
C      WRITE(*,1110) 'Up relaxations', iU
C      WRITE(*,1110) 'Down relaxations', iD
C      WRITE(*,1110) 'Fine Grid relaxations', iM
C      WRITE(*,1110) 'iSTRT', iSTRT
C      WRITE(*,1110) 'Relaxation Index', iRELAX
C      WRITE(*,1110) 'Truncation Error', iTAU
C      WRITE(*,1110) 'iCOEF', iCOEF
C      WRITE(*,1110) 'Cycle Recursion', iVW
C      WRITE(*,1110) 'Mcycl', Mcycl
C      WRITE(*,1100) 'ALPHl', ALPHl
C      WRITE(*,1100) 'ALPHm', ALPHm
C      WRITE(*,1110) 'Output Level', iPR
C      WRITE(*,1110) 'CG construction', iSKIP
C      WRITE(*,1110) 'iSTRT2', iSTRT2
C      WRITE(*,1110) 'Cycle Recursion', iVW2
C      WRITE(*,1110) 'Mcycl2', Mcycl2
      
C -------------------------------------
C     I/O Parameters
C -------------------------------------

      DO i=1, NIOFLAG
         IOFLAG(i)=.FALSE.
      ENDDO

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

      InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      InWORK(i_InWORK_RES) = .TRUE.     ! store RES in work array

      CALL BOXspace3( Nx, Ny, Nz, iFD, iRELAX, IWORK,
     &                NOGm, NFm, NSOm, NOG, NF, NC, NSO, NSOR,
     &                NCI, NCBW, NCU, NF2, NSO3, NSR3, NCI3, 
     &                NCBW2, NCU2, NCPL2, NIWORK, NRWORK,
     &                pWORK, NpWORK, InWORK, NInWORK, pSR, pSI )

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

      TOL_SAVE = TOL

      CALL BMG3D( Nx, Ny, Nz,
     &            TOL, ISTOP, IFD, iU, iD, iM,
     &            ISTRT, IRELAX, IRELAX_SYM,
     &            ICOEF, IVW, MCYCL, IOFLAG, ISKIP,
     &            bu, u, RWORK(pWORK(ip_RES)), NF, NC, SO, NSO,
     &            RWORK(pWORK(ip_SOR)), NSOR, 
     &            RWORK(pWORK(ip_CI)), NCI, 
     &            RWORK(pWORK(ip_CSO)), RWORK(pWORK(ip_CU)), NCBW, NCU,
     &            IWORK(pWORK(ip_iG)), NOG,
     &            iSTRT2, iVW2, MCYCL2,
     &            RWORK(pWORK(ip_PL_U)), RWORK(pWORK(ip_PL_Q)),
     &            RWORK(pWORK(ip_PL_RES)), NF2, 
     &            RWORK(pWORK(ip_PL_SO)), NSO3, 
     &            RWORK(pWORK(ip_PL_SOR)), NSR3, 
     &            RWORK(pWORK(ip_PL_CI)), NCI3,
     &            RWORK(pWORK(ip_PL_CSO)), 
     &            RWORK(pWORK(ip_PL_CU)), NCBW2, NCU2, NCPL2, 
     &            IWORK(pWORK(ip_PL_iG)), NOGc )

C ===========================================================================
C     >>>>>>>>>>>>>>>>     END: SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C ===========================================================================

C ===========================================================================
C     >>>>>>>>>>>>>>>>     BEGIN: ISKIP SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<
C ===========================================================================

      ISKIP=2
      TOL = TOL_SAVE

C ------------------------------
C     Solve the System:    
C ------------------------------

      CALL BMG3D( Nx, Ny, Nz,
     &            TOL, ISTOP, IFD, iU, iD, iM,
     &            ISTRT, IRELAX, IRELAX_SYM,
     &            ICOEF, IVW, MCYCL, IOFLAG, ISKIP,
     &            bv, v, RWORK(pWORK(ip_RES)), NF, NC, SO, NSO,
     &            RWORK(pWORK(ip_SOR)), NSOR, 
     &            RWORK(pWORK(ip_CI)), NCI, 
     &            RWORK(pWORK(ip_CSO)), RWORK(pWORK(ip_CU)), NCBW, NCU,
     &            IWORK(pWORK(ip_iG)), NOG,
     &            iSTRT2, iVW2, MCYCL2,
     &            RWORK(pWORK(ip_PL_U)), RWORK(pWORK(ip_PL_Q)),
     &            RWORK(pWORK(ip_RES)), NF2, 
     &            RWORK(pWORK(ip_PL_SO)), NSO3, 
     &            RWORK(pWORK(ip_PL_SOR)), NSR3, 
     &            RWORK(pWORK(ip_PL_CI)), NCI3,
     &            RWORK(pWORK(ip_PL_CSO)), 
     &            RWORK(pWORK(ip_PL_CU)), NCBW2, NCU2, NCPL2, 
     &            IWORK(pWORK(ip_PL_iG)), NOGc )


C ===========================================================================
C     >>>>>>>>>>>>>>>>     END: ISKIP-SOLVE    <<<<<<<<<<<<<<<<<<<<<<<<<<<<
C ===========================================================================


C ---------------------------------
C     calculate the inner products
C ---------------------------------

      CALL VECPROD (prod1, bu, v, Nx+2, Ny+2, Nz+2)
      CALL VECPROD (prod2, u, bv, Nx+2, Ny+2, Nz+2)

      WRITE(*,300)'The two inner products:',prod1,prod2
      WRITE(*,310)'Their difference:      ',prod1-prod2


C ===========================================================================

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

C ===========================================================================

      END


