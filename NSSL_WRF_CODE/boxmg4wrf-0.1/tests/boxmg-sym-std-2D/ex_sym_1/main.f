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
C   However, I use a direct call to BOXMGP with the memory allocation
C   parameters computed by the program space.f for max dimesions. For
C   the actual call to BOXMG the parameters are computed by the subroutine
C   BOXspace based on the actaul dimensions of the problem and the desired
C   number of grids.
C
C ==========================================================================

      IMPLICIT  NONE

      INCLUDE   'BMG_constants.h'
      INCLUDE   'BMG_workspace.h'
      INCLUDE   'BMG_parameters.h'

C ---------------------------------------------
C     Critical Parameters:
C
      INTEGER   IBC, IRELAX, IStncl, Lmax, Mmax, NXYCmax, NXYCmin
      PARAMETER (IBC=0, IRELAX=4, IStncl=1, Lmax=30, Mmax=60,
     &           NXYCmax=3, NXYCmin=3)

C --------------------------------------------
C     Multigrid/Workspace Memory Allocation:
C
      INTEGER   NFmax, NOGmax, NOGmin, NSOmax
      PARAMETER ( NFmax=2758, NSOmax=9822, NOGmin=4, NOGmax=4 )

      INTEGER   NIWORK, NRWORK
      PARAMETER ( NIWORK=36, NRWORK=20206 )

      INTEGER   NInWORK, NpWORK
      PARAMETER ( NinWORK=4, NpWORK=19 )

      INTEGER   pSI, pSR, pWORK(NBMG_pWORK)
      LOGICAL   InWORK(NBMG_InWORK)

C --------------------------------------------
C     Variable Declarations:
C      
      INTEGER   i, ICG, IFMG, IFMGk, IRELAX_SYM, ISKIP, ISTOP, IVW,
     &          IWORK(NIWORK), KG, L, M, MGCmax, NC, NCI, NCU, NF,
     &          NOG, NSO, NCBW, NSOR, NUd, NUf, NUu, NXYC
      REAL*8    bu(NFmax), bv(NFmax), hx, hy, prod1, prod2, RES(NFmax),
     &          RHS(NFmax), RWORK(NRWORK), SO(NSOmax), TOL, TOL_SAVE,
     &          u(NFmax), v(NFmax)
      LOGICAL   IOFLAG(NBMG_IOFLAG)

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
C ------------------------------

      NUd=1                    ! # of relaxations before moving to grid k-1
      NUu=1                    ! # of relaxations before moving to grid k+1
      NUf=0                    ! # of relaxations on the finest grid
      IFMG=-1                  ! +1 -> FMG
                               ! -1 -> n-cycling from finest grid
      IFMGk=1                  ! # of times grid k is visited before moving
                               ! up to grid k+1
      IVW=2                    ! n-cycles: 1 -> V cycles 2 -> W cycles etc.
      MGCmax=1                 ! max. # of MG cycles
      
      TOL=rZERO                ! while only standard MG cycles are performed

      ISKIP=0                  ! skip nothing

      NXYC=3                   ! min. # of grid points in x or y

      IRELAX_SYM = BMG_RELAX_SYM     ! Symmetric n-cycles

      ISTOP = BMG_STOP_REL_RES_L2    !  Convergence Measure
      TOL   = 1D-6                   !  Convergence Criteria

C -------------------------------------
C     I/O Parameters
C -------------------------------------

      DO i=1, NBMG_IOFLAG
         IOFLAG(i)=.FALSE.
      ENDDO

C -------------------------------------
C     Space requirements and pointers:
C -------------------------------------

      pSR=1
      pSI=1      

      InWORK(i_InWORK_SO)  = .FALSE.    ! use a separate array for SO
      InWORK(i_InWORK_U)   = .FALSE.    ! use a separate array for Q
      InWORK(i_InWORK_Q)   = .FALSE.    ! use a separate array for QF
      InWORK(i_InWORK_RES) = .TRUE.     ! store RES in work array

      CALL BMG2_SymStd_SETUP_PtrWork(
     &                 L, M, NXYc, IStncl, IBC, IRELAX, 
     &                 NOGmax, NFmax, NSOmax,
     &                 NOG, NF, NC, NCI, NSO, NSOR,
     &                 NCBW, NCU, NIWORK, NRWORK,
     &                 pWORK, NpWORK, InWORK, NInWORK, pSR, pSI 
     &                ) 

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

      TOL_SAVE=TOL

      CALL BOXMG( L, M,
     &            TOL, ISTOP, IStncl, NUu, NUd, NUf, IFMG*MGCmax,
     &            IRELAX, IRELAX_SYM, IVW, IFMGk, IOFLAG, ISKIP,
     &            bv, v, RES, NF, NC, SO, NSO, 
     &            RWORK(pWORK(ip_SOR)), NSOR, RWORK(pWORK(ip_CI)), NCI, 
     &            RWORK(pWORK(ip_CSO)), RWORK(pWORK(ip_CU)), NCBW, NCU,
     &            IWORK(pWORK(ip_iG)), NOG, NXYC, KG                  )
        
      ISKIP=2                  ! skip everything
      TOL = TOL_SAVE

      CALL BOXMG( L, M,
     &            TOL, ISTOP, IStncl, NUu, NUd, NUf, IFMG*MGCmax,
     &            IRELAX, IRELAX_SYM, IVW, IFMGk, IOFLAG, ISKIP,
     &            bu, u, RES, NF, NC, SO, NSO, 
     &            RWORK(pWORK(ip_SOR)), NSOR, RWORK(pWORK(ip_CI)), NCI, 
     &            RWORK(pWORK(ip_CSO)), RWORK(pWORK(ip_CU)), NCBW, NCU,
     &            IWORK(pWORK(ip_iG)), NOG, NXYC, KG                  )

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

      






