      PROGRAM MAIN

C ==========================================================================
C
C     Author: J. David Moulton
C     Date:   November 19, 1993.
C      
C     This solves Example 1 given in Victor Bandy's BBMG interface guide.
C     But it uses a W-cycle instead of F-cycles.
C
c------------------------------------------------------------------------
c***Author  Bandy, Victor A.
c             Computational Mathematics Group
c             University of Colorado at Denver
c
c   This program solves the PDE given below using the BBMG interface to
c   the BOXMG, black box multigrid, solver where the system of equations
c   is obtained by using central finite difference approximations.
c
c      d [ P du ]   d [ P du ]            2    2
c    - --    --   - --    --   + u = -7 (x  + y  ) - 4
c      dx    dx     dy    dy
c
c     where  P = x**2 + y**2 + 1
c     and divergence flux boundary conditions on all four boundaries,
c
c       normal dot [ P grad u ] = g
c
c     and the exact solution of     u = x**2 + y**2
c     with domain X=[0,1], Y=[1,3]
c
c--variables
c
c  l,m   the number of unknowns in the x and y directions respectively.
c
c  deltax,deltay  the x and y direction grid spacings respectively.
c
c  For the other variables see the documentation header in BBMG.
c
c
c----------------------------------------------------------------------
C
C   However, I use a direct call to BOXMG with the memory allocation
C   parameters computed by the program space.f for max dimesions. For
C   the actual call to BOXMG the parameters are computed by the subroutine
C   BOXspace based on the actaul dimensions of the problem and the desired
C   number of grids.
C
C ==========================================================================

      IMPLICIT  NONE

      INCLUDE   'BMG_SER_constants.h'
      INCLUDE   'BMG_SER_workspace.h'
      INCLUDE   'BMG_SER_parameters.h'

C ---------------------------------------------
C     Critical Parameters:
C
      INTEGER   Lmax, Mmax, NXYCmax, NXYCmin
      PARAMETER ( Lmax=30, Mmax=60, NXYCmax=3, NXYCmin=3 )

C --------------------------------------------
C     Multigrid/Workspace Memory Allocation:
C
      INTEGER   NFmax, NOGmax, NOGmin, NSOmax
      PARAMETER ( NFmax=2758, NSOmax=9822, NOGmin=4, NOGmax=4 )

      INTEGER   NBMG_iWORK, NBMG_rWORK
      PARAMETER ( NBMG_iWORK=36, NBMG_rWORK=20206 )

      INTEGER   pSI, pSR, BMG_pWORK(NBMG_SER_pWORK)
      LOGICAL   BMG_InWORK(NBMG_SER_InWORK)

C --------------------------------------------
C     Variable Declarations:
C      
      INTEGER   BMG_iPARMS(NBMG_SER_iPARMS), i, BMG_iWORK(NBMG_iWORK),
     &          KG, L, M, NC, NCI, NCU, NF, NOG, NSO, NCBW, NSOR
      REAL*8    BMG_rPARMS(NBMG_SER_rPARMS), bu(NFmax), bv(NFmax), 
     &          hx, hy, prod1, prod2, RES(NFmax),
     &          RHS(NFmax), BMG_rWORK(NBMG_rWORK), SO(NSOmax), TOL_SAVE,
     &          u(NFmax), v(NFmax)
      LOGICAL   BMG_IOFLAG(NBMG_SER_IOFLAG)

C ===========================================================================

C ---------------------------------------------------------
C     Initialize the arguments
C ---------------------------------- 
C     Domain and grid:
C ----------------------------------

      WRITE(*,*) 'Enter the number of points in x ... '
      READ(*,*) L
      WRITE(*,*) 'Enter the number of points in y ... '
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

      BMG_iPARMS(id_BMG2_SER_STENCIL)   = BMG_SER_STENCIL_5pt
      BMG_iPARMS(id_BMG2_SER_MAX_ITERS) = 1
      BMG_rPARMS(id_BMG2_SER_STOP_TOL)  = 1D-6

C -------------------------------------
C     Space requirements and pointers:
C -------------------------------------

      pSR=1
      pSI=1      

      DO i=1, NBMG_SER_InWORK
         BMG_InWORK(i) = .FALSE.
      ENDDO

      BMG_InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      BMG_InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      BMG_InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      BMG_InWORK(i_InWORK_RES) = .TRUE.     ! store RES in work array

      CALL BMG2_SER_SymStd_SETUP_PtrWork(
     &                 L, M, BMG_iPARMS, 
     &                 NOGmax, NFmax, NSOmax, NBMG_iWORK, NBMG_rWORK,
     &                 BMG_pWORK, BMG_InWORK, pSR, pSI 
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
         SO(i)=rZERO
 10   CONTINUE

      DO 20 i=1, NF
         u(i)   = rZERO
         v(i)   = rZERO
         bu(i)  = rZERO
         bv(i)  = rZERO
         rhs(i) = rZERO
         RES(i) = rZERO
 20   CONTINUE

C -------------------------------------
C     Define the coefficient matrix
C -------------------------------------

      CALL PUTF( SO, rhs, L, M, hx, hy, 0 )

C -------------------------------------
C     initialize the vectors
C -------------------------------------

      CALL INITVECS ( u, v, L, M )

      do i=1,NF
         bu(i) = rZERO
         bv(i) = rZERO
      enddo

C -------------------------------------
C     call boxmg 
C -------------------------------------

      TOL_SAVE = BMG_rPARMS(id_BMG2_SER_STOP_TOL)

      CALL BMG2_SER_SymStd_SOLVE_boxmg(
     &          L, M, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          bv, v, RES, NF, NC, SO, NSO, 
     &          BMG_rWORK(BMG_pWORK(ip_SOR)), NSOR,
     &          BMG_rWORK(BMG_pWORK(ip_CI)), NCI, 
     &          BMG_rWORK(BMG_pWORK(ip_CSO)),
     &          BMG_rWORK(BMG_pWORK(ip_CU)), NCBW, NCU,
     &          BMG_iWORK(BMG_pWORK(ip_iG)), NOG, KG
     &          )
      
      BMG_iPARMS(id_BMG2_SER_SETUP)    = BMG_SER_SETUP_none   ! skip the setup
      BMG_rPARMS(id_BMG2_SER_STOP_TOL) = TOL_SAVE         ! reset criteria

      CALL BMG2_SER_SymStd_SOLVE_boxmg(
     &          L, M, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          bu, u, RES, NF, NC, SO, NSO, 
     &          BMG_rWORK(BMG_pWORK(ip_SOR)), NSOR,
     &          BMG_rWORK(BMG_pWORK(ip_CI)), NCI, 
     &          BMG_rWORK(BMG_pWORK(ip_CSO)),
     &          BMG_rWORK(BMG_pWORK(ip_CU)), NCBW, NCU,
     &          BMG_iWORK(BMG_pWORK(ip_iG)), NOG, KG
     &          )

C ---------------------------------
C     calculate the inner products
C ---------------------------------

      CALL VECPROD (prod1, bu, v, L, M)
      CALL VECPROD (prod2, u, bv, L, M)

      WRITE(*,300)'The two inner products:',prod1,prod2
      WRITE(*,310)'Their difference:      ',prod1-prod2

 300  FORMAT(/,2X,A,1X,F22.14,4X,F22.14,/)
 310  FORMAT(2X,A,1X,1P,E22.14,/)

      END

      






