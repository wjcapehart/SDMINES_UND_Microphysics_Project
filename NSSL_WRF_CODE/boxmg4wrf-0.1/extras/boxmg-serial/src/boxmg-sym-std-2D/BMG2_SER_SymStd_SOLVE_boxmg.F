      SUBROUTINE BMG2_SER_SymStd_SOLVE_boxmg(
     &                Nxm, Nym, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm,
     &                SO, NSOm, SOR, NSORm, CI, NCIm,
     &                ABD, BBD, NCBWm, NCUm,
     &                IGRD, NOGm, NOG 
     &                )


C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This subroutine is a black box multigrid solver (BoxMG) for
C     discretizations of second order elliptic partial differential
C     equations that generate, at most, a 9-point stencil on a logically
C     rectangular grid.  It may be applied to other, similarly
C     structured, problems.
C
C     BMG2_SER_SymStd_SOLVE_boxmg takes a matrix defined by the user from a
C     given (fine) grid and constructs coarser grids and their
C     associated coefficient matrices. BOXMG then performs the multigrid
C     algorithm and returns the solution vector. The matrix is set by a
C     user written subroutine. The difference stencil at (i,j) is
C     assumed to be of the form
C
C                      xnw xn xne
C                      xw  xo xe
C                      xsw xs xse
C
C     where 
C
C         xo  = SO(i,j,ko) 
C         xw  = - SO(i,j,kw)
C         xs  = - SO(i,j,ks)
C         xsw = - SO(i,j,ksw)
C         xse = - SO(i+1,j,knw)
C
C     and where xe, xn, xnw, xn and xne, are specified by symmetry.  The
C     indeces k* are defined in the include file 'BMG_SER_stencils.h'.  In
C     the case of a five point operator only xo, xw, and xs need to be
C     set. The difference scheme is assumed to be positive definite.
C     
C     It is important to note that the fictitious points:
C
C        (1,j), j=1,jj, (ii,j),j=1,jj, (i,1), i=1,ii, (i,jj),i=1,ii
C
C     are used for ease of programming.  In particular, it is assumed
C     that the user has set any difference coefficient referring  to 
C     these points to zero. For example, if the above difference
C     template is centered at (2,2), then xnw, xw, xs, and xse
C     would be assumed zero.
C
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     Nx      x-dimension of the grid, excluding fictitious points
C     Ny      y-dimension of the grid, excluding fictitious points
C
C     TOL      Convergence tolerance.
C
C     ISTOP    Indicates the stopping criteria. At present there
C              is a relative residual test
C
C                ISTOP = BMG_SER_STOP_REL_RES_L2 
C
C                => iteration has converged if the l2-norm of the
C                   current residual, divided by the l2-norm of the 
C                   initial residual is less than TOL.
C
C              and an absolute residual test
C  
C                ISTOP = BMG_SER_STOP_REL_ABS_L2
C
C                => iteration has converged if the l2-norm of the 
C                  current residual is less than TOL.
C
C              NOTE: if the selected test fails and the number of
C              iterations is greater than ABS(ISTRT) then the iteration
C              has failed => set TOL=-RES_L2 and return.
C
C
C     IFD      Indicator for difference scheme. IFD=1 means a five
C              point operator. I.e., in the template given above, xnw,
C              xsw, xse, and xne are assumed zero. IFD.ne.1 means a nine
C              point operator.
C              (Note that nine point operators may be generated on the
C              coarser grids even in the case of a five point operator.)
C
C     IU       Number of relaxation sweeps to be performed
C              on a coarse grid before interpolation to a fine grid.
C
C     ID       Number of relaxation sweeps to be performed on a fine
C              grid before the problem is transferred to a coarse grid.
C
C     IM       Unused
C
C     ISTRT    Indicator for number of multigrid cycles to be performed
C              if ISTRT.gt.0 the algorithm will begin on the coarsest grid
C              with a full multigrid cycle (however,
C              without cubic interpolation). It will continue
C              cycling until ISTRT cycles have been performed or until
C              the error criterion is satisfied. If ISTRT.lt.0 the
C              algorithm will begin on the finest grid and will perform
C              -ISTRT cycles unless the error criterion is satisfied first.
C
C     ISETUP   Indicator for whether or not to skip initial setup,
C              initial guess for Q, computation of pointers, and
C              generation of coefficients on the coarser grids.
C            
C                ISETUP=0  => perform setup, including,
C                             - compute workspace pointers
C                             - compute coarse-grid operators
C                             - compute interpolation transfer operators
C
C              If ISETUP.NE.0 then we have the following cases:
C         
C                ISETUP=1  => skip computation of workspace pointers
C                             - compute coarse-grid operators
C                             - compute interpolation operators      
C                ISETUP=2  => skip all setup 
C
C              It is important to note the interaction that ISETUP 
C              has with ISTRT.  Specifically, with ISETUP.NE.0 we have
C            
C                 ISTRT > 0 => FMG so residual is coarsened appropriately
C                 ISTRT < 0 => n-cycle so nothing extra to do here.
C
C              An example of usage would be a constant coefficient 
C              time-dependent problem where ISETUP=0 on the first time step 
C              and ISETUP=2 on subsequent time steps.
C
C     BMG_IOFLAG   Logical array of I/O and debugging flags. Admittedly
C                  this is not well documented yet, but look in the include
C                  file "BMG_SER_parameters.h" for hints (the names are verbose.)
C
C     IRELAX   Indicator for the relaxation algorithm:
C 
C              IRELAX=1  => Colored Gauss-Seidel relaxation:
C                        -  red-black on the finest grid for 5-point stencils
C                        -  four color otherwise
C
C              IRELAX=2  => red-black x-line Gauss-Seidel relaxation
C
C              IRELAX=3  => red-black y-line line Gauss-Seidel relaxation
C
C              IRELAX=4  => red-black x-line Gauss-Seidel relaxation
C                           red-black y-line Gauss-Seidel relaxation
C
C     IRELAX_SYM  Symmetric cycles require that the relaxation ordering
C                 be reversed between the restriction and interpolation
C                 stage of the cycle.   This is important when BOXMG 
C                 is being used as a preconditioner.
C              
C                    IRELAX_SYM = BMG_SER_RELAX_NONSYM
C                    => use the stated relaxation ordering (nonsymmetric)
C
C                    IRELAX_SYM = BMG_SER_RELAX_SYM
C                    => make the cycle symmetric.
C
C                 as defined in BMG_SER_parameters.h
C
C     IVW      Indicator for type of cycle. IVW=1 means v-cycles will
C              be performed. IVW=2 means w-cycles will be performed, etc.
C
C     MCYCL    Currently disabled!
C
C     NFm      Maximum storage for a vector on all grids. This should be
C              computed by BMG2_SER_SymStd_SETUP_space and is checked in boxmg.
C
C     NCm      Maximum storage for a vector on all coarse grids.  This should
C              be computed in BMG2_SER_SymStd_SETUP_space and is checked in boxmg.
C
C     NSOm     See INPUT/OUTPUT
C
C     NSORm    See WORK ARRAYS.
C
C     NCIm     See WORK ARRAYS.
C
C     NOGm     Maximum number of grids that can be supported by the 
C              pointer array IGRD. It should be set in calling program.
C
C     NXYc     Limit for number of points on coarsest grid.  Coarsening
C              occurs until the number of x or y unknowns minus one is 
C              less than or equal to NXYc.  NXYc should be 3 or greater 
C              or the code will abort.  This parameter is included because 
C              efficiency is a function of machine.   On vector machines, 
C              it frequently pays to take NXYc larger than on scalar 
C              machines.
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     SO       User defined real array that contains the coefficient
C              matrix. See above description for format.  If SETUP is
C              specified then coarse-grid operators are created and returned. 
C
C     NSOm     Dimension of SO, the calling program should have computed
C              this with a call to BMG2_SER_SymStd_SETUP_space.  Clearly,
C            
C              IFD .EQ. 1  =>    NSOm = 3*NFm + 2*NCm 
C              IFD .NE. 1  =>    NSOm = 5*NFm 
C  
C              BMG2_SER_SymStd_SOLVE_boxmg verifies that NSOm is large enough.
C
C     CI       CI is a real array which should be dimensioned
C              to (NCIm) in calling program. It is used to contain
C              the interpolation coefficients.
C
C     NCIm     Dimension of CI, set in calling program.  NCIm = 8*NCm is 
C              enough.  BOXMG checks to see if NCIm is large enough.
C
C     QF       The user defined array that contains the right hand side.
C              It is set to dimension NFm.
C
C     Q        The user defined array that contains The solution vector.
C              It is set to dimension NFm.
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   WORK ARRAYS:
C  --------------------
C
C     ABD      User declared two dimensional real array, which
C              is used to store the coefficient matrix for the coarsest
C              grid. It is then used by the LINPACK routine DPBSL. It
C              should be dimensioned to (NCBWm,NCUm).
C
C     BBD      User declared real array of dimension NCUm for use in the
C              the LINPACK routine DPBSL.
C
C     NCBWm    Maximum first subscript of ABD, which needs to be > or =
C              number of grid points plus 2 in x direction on coarsest
C              grid. BOXMG checks to see if it is large enough.
C
C     NCUm     Maximum second subscript of ABD, which needs to be > or =
C              (x-dimension)*(y-dimension) on coarsest grid. BOXMG
C              checks to see if it is large enough.
C
C     IGRD     A work array. IGRD(.,1),...,IGRD(.,6)
C              are used to store starting locations for the arrays on each
C              grid. IGRD(.,7),...,IGRD(.,9) are used as scratch integer
C              arrays needed by BOXMG. IGRD should be dimensioned to
C              IGRD(NOGm,9) in the calling program.
C
C     SOR      SOR is a real array, which is set to dimension NSORm.
C              It is used as workspace for the tridiagonal factorization 
C              and solves used in the case of line relaxation.  
C
C     NSORm    Dimension of SOR, set in calling program.  Actual requirements
C              need to be cleaned up!
C
C              IRELAX .EQ. 1   =>   NSORm = 0
C              IRELAX .EQ. 2   =>   NSORm = 2*NFm 
C              IRELAX .EQ. 3   =>   NSORm = 4*NFm  (this should be fixed)
C              IRELAX .EQ. 4   =>   NSORm = 4*NFm
C
C              BOXMG checks to see if NSORm is large enough.
C
C =======================================================================
C  --------------------
C   REFERENCES:
C  --------------------
C
C     Dendy, J. E. Jr., "Black Box Multigrid", Journal of Computational
C     Physics, Vol. 48, pp. 366-386, 1982
C     
C     Dendy, J. E. Jr., "Black Box Multigrid for Nonsymmetric Problems",
C     Applied Mathematics and Computation, Vol. 13, pp. 261-283, 1983
C
C     Dendy, J. E. Jr., "Black Box Multigrid for Systems", Applied
C     Mathematics and Computation, Vol. 19, pp. 57-74, 1986
C     
C     Dendy, J. E. Jr., "Two Multigrid Methods for Three Dimensional
C     Problems with Discountinuos and Anisotropic Coefficients", SIAM
C     Journal of Scientific and Satatistical Computing, Vol. 8, No. 2,
C     September 1987
C     
C     Dendy, J. E. Jr., "Black Box Multigrid for Periodic and Singular
C     Problems", Applied Mathematics and Computation, Vol. 25, pp. 1-10,
C     1988
C
C ==========================================================================
C  --------------------
C   Error Codes:
C  --------------------
C
C     -------------------------------------------------------------------
C     CODE    CONDITION                                 ORIGIN
C     -------------------------------------------------------------------
C      1     NOG=0                                     SETUP_parts
C      2     NOG<0                                     SETUP_parts
C      3     initial residual=0                        SOLVE_boxmg
C      4     min coarse grid dimension < 3             SETUP_space
C      5     computed number of grids too small        SETUP_space
C      6     computed number of grids < 1              SETUP_space
C      7     fatal error in SETUP_PtrGrid              SETUP_PtrGrid
C      8     IRELAX out of range                       SETUP_PtrGrid
C      9     number of grids < 2                       ncycle
C     10     updown parameter out of range             updown
C     11     NOG=1                                     SETUP_relax
C     12     NOG=0                                     SETUP_relax
C     13     NOG<0                                     SETUP_relax
C     14     fatal error in SETUP_PtrWork              SETUP_PtrWork
C     15     memory allocation mode unspecified        SETUP_PtrWork
C     16     KF-1 .ne. KG in SETUP_cg_ITLI             SETUP_cg_ITLI
C     17     Cholesky decomposition failed, look for 
C            the return code of the LAPACK routine in
C            BMG_iPARMS(id_BMG2_SER_Ext_Err_Code)          SETUP_cg_LU
C     18     NStncl .ne. 3,5                           SETUP_cg_LU
C     19     KF-1 .ne. KG in SETUP_interp_OI           SETUP_interp_OI
C     20     Coarse grid solve failed,  look for the
C            return code of the LAPACK routine in
C            BMG_iPARMS(id_BMG2_SER_Ext_Err_Code)          SOLVE_cg
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_stencils.h'
      INCLUDE 'BMG_SER_workspace.h'
      INCLUDE 'BMG_SER_parameters.h'

C ----------------------------
C     Argument Declarations
C 

      INTEGER  NCIm, NSOm, NSORm, NCBWm, NCUm, NCm, NFm, NOGm, Nxm, Nym

      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS), IGRD(NOGm,9)
      REAL*RKIND   ABD(NCBWm*NCUm), BBD(NCUm), 
     &         BMG_rPARMS(NBMG_SER_rPARMS),
     &         CI(NCIm), SO(NSOm), SOR(NSORm), Q(NFm), QF(NFm), RES(NFm)   
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  I, iBC, ID, IFD, IM,
     &         IRELAX, IRELAX_SYM, ISETUP, ISTOP, ISTRT, IU, IVW,
     &         K, KC, KF, KFMG, MCYC, MCYCL,
     &         NC, NCBW, NCI, NCYC, NCU, NF, NOG, NSO, NSOR, NSORv,
     &         NStncl, NXYc, Nx, Nxc, Nxf, Ny, Nyc, Nyf,
     &         p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC

      REAL*RKIND   RES_L2, RES_L2_0, REL_RES_L2, TOL
      REAL*8 T, T1, T2
      LOGICAL  NCYCLE_FLAG

C =========================================================================


      IF (BMG_iPARMS(id_BMG2_SER_Err_Code) .ne. iZERO) THEN
         RETURN
      ENDIF

      T  = dZERO
      T1 = dZERO
      T2 = dZERO

      IF ( BMG_IOFLAG(iBMG2_SER_BUG_PARAMETERS) ) THEN
         CALL BMG2_SER_SymStd_DUMP_parms( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &             )
      ENDIF

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      IFD    = BMG_iPARMS(id_BMG2_SER_STENCIL)
      IBC    = BMG_iPARMS(id_BMG2_SER_BC) 
      ISETUP = BMG_iPARMS(id_BMG2_SER_SETUP)
      NXYc   = BMG_iPARMS(id_BMG2_SER_CG_MIN_DIM)

      IRELAX = BMG_iPARMS(id_BMG2_SER_RELAX)
      IRELAX_SYM = BMG_iPARMS(id_BMG2_SER_RELAX_SYM )

      ID = BMG_iPARMS(id_BMG2_SER_NRELAX_DOWN)
      IU = BMG_iPARMS(id_BMG2_SER_NRELAX_UP)
      IM = BMG_iPARMS(id_BMG2_SER_NRELAX_FG)

      IF ( BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS).EQ.BMG_SER_N_CYCLE ) THEN
         ISTRT = - BMG_iPARMS(id_BMG2_SER_MAX_ITERS)
      ELSEIF ( BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS)
     &        .EQ.BMG_SER_FMG_CYCLE ) THEN 
         ISTRT = BMG_iPARMS(id_BMG2_SER_MAX_ITERS)
      ENDIF

      IVW   = BMG_iPARMS(id_BMG2_SER_NCYCLE_TYPE)
      MCYCL = BMG_iPARMS(id_BMG2_SER_FMG_NNCYCLE)

      ISTOP = BMG_iPARMS(id_BMG2_SER_STOP_TEST)
      TOL   = BMG_rPARMS(id_BMG2_SER_STOP_TOL)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> BEGIN:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      IF ( ISETUP.EQ.BMG_SER_SETUP_ptrs_opers .OR. 
     &     ISETUP.EQ.BMG_SER_SETUP_only ) THEN

         CALL BMG2_SER_SymStd_SETUP_PtrGrid( 
     &             Nxm, Nym, BMG_iPARMS,
     &             NOGm, NFm, NCm, NCIm, NSOm, NSORm,
     &             NCBWm, NCUm,
     &             NOG, NF, NC, NCI, NSO, NSOR,
     &             NCBW, NCU, IGRD,
     &             BMG_IOFLAG
     &             )
         
         IF (BMG_iPARMS(id_BMG2_SER_Err_Code) .ne. iZERO) THEN
            RETURN
         ENDIF

         ! Output      
         IF( BMG_IOFLAG(iBMG2_SER_OUT_WSPACE_SIZE) ) THEN
            WRITE (*,260) 
     &            'Storage for a vector on all grids = ', NF
         ENDIF

         ! Output
         IF( BMG_IOFLAG(iBMG2_SER_OUT_WSPACE_SIZE) ) THEN
            WRITE (*,270) 
     &            'Storage for a vector on all coarse grids = ', NC
         ENDIF

      ELSE
         
         NF   = NFm
         NC   = NCm
         NSO  = NSOm
         NSOR = NSORm
         NCI  = NCIm
         
         NCBW = NCBWm
         NCU  = NCUm

      ENDIF

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>> END:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>> BEGIN:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      IF ( ISETUP.EQ.BMG_SER_SETUP_ptrs_opers .OR. 
     &     ISETUP.EQ.BMG_SER_SETUP_opers .OR. 
     &     ISETUP.EQ.BMG_SER_SETUP_only ) THEN 

         ! Start timer
         CALL BMG_timer(T1) 

         !
         ! Construct coarse-grid and interpolation operators.  In 
         ! addition, setup any necessary components for relaxation.
         !
         CALL BMG2_SER_SymStd_SETUP_parts( 
     &             IFD, IBC, IRELAX, SO, NSO, SOR, NSOR, CI, NCI,
     &             QF, RES, NF,  !! Part of workspace to be zeroed !!
     &             ISETUP, ISTRT,  !! Temporary for periodicity
     &             ABD, BBD, NCBW, NCU, IGRD, NOGm, NOG, BMG_IOFLAG,
     &             BMG_iPARMS)

         IF (BMG_iPARMS(id_BMG2_SER_Err_Code) .ne. iZERO) THEN
            RETURN
         ENDIF

         ! Compute the setup time
         CALL BMG_timer(T2)
         T=T2-T1

         ! Output the setup time
         IF( BMG_IOFLAG(iBMG2_SER_OUT_TIME_SETUP) )  THEN
            WRITE (*,220) '(2D) SETUP TIME =', T
         ENDIF
         
      ENDIF

      !
      !  No solve just the setup, thanks!  
      !
      IF (ISETUP.EQ.BMG_SER_SETUP_only) RETURN

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN: MULTIGRID CYCLING <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      ! Start the timer
      CALL BMG_timer(T1)

C ----------------------------------------------------------------
C     Initial residual:
C ----------------------------------------------------------------

      !
      ! Get pointers for the finest grid (GRID LEVEL = NOG)
      !
      CALL BMG2_SER_SymStd_GET_pointers(
     &          NOG, IGRD, NOGm, Nx, Ny, 
     &          p_U, p_SO, p_SOR, p_CI
     &          )

      !
      ! Set NStncl for the finest grid
      !
      IF ( IFD.NE.1 ) THEN
         NStncl=5
      ELSE
         NStncl=3
      ENDIF
      
      !
      ! Compute the initial residual
      !
      CALL BMG2_SER_SymStd_residual( 
     &          NOG, SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &          Nx, Ny, NOG, IFD, NStncl, IBC, 
     &          IRELAX, IRELAX_SYM, BMG_SER_UP
     &          )
      !
      CALL BMG2_SER_SymStd_UTILS_norm_l2(
     &          RES(p_U), Nx, Ny, RES_L2_0 
     &          )

      !
      !  Output
      !
      IF ( BMG_IOFLAG(iBMG2_SER_OUT_ITERATIONS) ) THEN
         WRITE (*,405) '*** THE INITIAL RESIDUAL (l2-NORM) = ', RES_L2_0
      ENDIF

      IF ( RES_L2_0 .EQ. rZERO ) THEN
         IF ( BMG_IOFLAG(iBMG2_SER_WARN_ZERO_RESIDUAL) ) THEN
            WRITE(*,505) '*** WARNING: boxmg.f '
            WRITE(*,505) '    Zero initial residual !!!! '
            WRITE(*,505) '    D''ooh! Initial guess is the solution!'
         ENDIF
         RETURN
      ENDIF

C ----------------------------------------------------------------
C     Direct solve in the degenerate case of one grid!
C ----------------------------------------------------------------

      IF ( NOG.EQ.1 ) THEN

         !
         ! Solve on the coarsest grid
         !
         IF ( IBC.EQ.BMG_SER_BCs_definite ) THEN
            CALL BMG2_SER_SymStd_SOLVE_cg(
     &                Q(p_U), QF(p_U), Nx, Ny,
     &                ABD, BBD, NCBW, NCU,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )
            IF (BMG_iPARMS(id_BMG2_SER_Err_Code) .ne. iZERO) THEN
               RETURN
            ENDIF


         ELSE
            CALL BMG2_SER_PerSymStd_SOLVE_cg(
     &                Q(p_U), QF(p_U), Nx, Ny, 
     &                ABD, BBD, NCBW, NCU, IBC 
     &                )
         ENDIF
         !
         ! Compute the final residual
         !
         CALL BMG2_SER_SymStd_residual( 
     &             NOG, SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &             Nx, Ny, NOG, IFD, NStncl, IBC, 
     &             IRELAX, IRELAX_SYM, BMG_SER_UP
     &             )
         !
         CALL BMG2_SER_SymStd_UTILS_norm_l2(
     &             RES(p_U), Nx, Ny, RES_L2 
     &             )
         !
         IF( BMG_IOFLAG(iBMG2_SER_BUG_RES_CG_SOLVE) ) THEN 
            WRITE (*,404) '*** THE FINAL RESIDUAL (l2-NORM) = ', RES_L2
         ENDIF
         !
         TOL = RES_L2
         ! Jump to final I/O
         GOTO 150
         !
      ENDIF

C --------------------------------
C     Cycle Preliminaries:
C -------------------------------- 

      KF   = NOG                ! finest grid index for cycling
      KC   = 1                  ! coarsest grid index
      KFMG = 1                  ! FMG current finest grid
      MCYC = ABS(ISTRT)         ! Maximum number of multigrid cycles
      NCYC = 1                  ! multigrid cycle counter
      NCYCLE_FLAG = .TRUE.      ! default to do n-cycles
      
      !
      !  Number of temporary vectors in SOR
      !
      IF ( iRELAX.EQ.BMG_SER_GS_RB_point 
     &    .OR. iRELAX.EQ.BMG_SER_GS_RB_x_lines ) THEN
         NSORv = 2
      ELSE 
         NSORv = 4
      ENDIF

C ----------------------------------------------------------------
C     Perform multigrid F-cycle:
C     (if necessary, then continue with n-cycles)
C ----------------------------------------------------------------

      IF ( ISTRT.GT.0 ) THEN
      
         ! Coarsen the right hand side
         DO K=KF-1, 1, -1
            CALL BMG2_SER_SymStd_GET_pointers(
     &                K, IGRD, NOGm, Nxc, Nyc,
     &                p_UC, p_SOC, p_SORC, p_CIC
     &                )
            CALL BMG2_SER_SymStd_GET_pointers( 
     &                K+1, IGRD, NOGm, Nxf, Nyf,
     &                p_U, p_SO, p_SOR, p_CI
     &                )
            IF ( IBC.EQ.BMG_SER_BCs_definite ) THEN
               CALL BMG2_SER_SymStd_restrict(
     &                   K+1, K, QF(p_U), QF(p_UC), CI(p_CIC),
     &                   Nxf, Nyf, Nxc, Nyc, IBC
     &                   )
               ! Copy right hand side into RES(i,j)
               DO i=0, Nxf*Nyf-1
                  RES(p_U+i)=QF(p_U+i)
               END DO
               !
            ELSE
C
C           There should be a periodic restriction here but
C           I can't quite figure out how to make it work.
C
C               CALL BMG2_SER_PerSymStd_residual(
C     &                   K+1, K, iONE,
C     &                   SO(p_SO), SO(p_SOC), QF(p_U), QF(p_UC), 
C     &                   Q(p_U), Q(p_UC), CI(p_CIC), SOR(p_SOR), 
C     &                   Nxc, Nyc, Nx, Ny, TAUERR, NStncl, NSORv, IBC
C     &                   )

            ENDIF
         END DO

         ! FMG current finest grid index is set to the coarsest grid
         KFMG=KC     
         ! FMG: solve exactly on the coarsest grid
         CALL BMG2_SER_SymStd_GET_pointers(
     &             KFMG, IGRD, NOGm, Nx, Ny,
     &             p_U, p_SO, p_SOR, p_CI
     &             )
         IF ( IBC.EQ.BMG_SER_BCs_definite ) THEN
            CALL BMG2_SER_SymStd_SOLVE_cg(
     &                Q(p_U), QF(p_U), Nx, Ny,
     &                ABD, BBD, NCBW, NCU,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )
            IF (BMG_iPARMS(id_BMG2_SER_Err_Code) .ne. iZERO) THEN
               RETURN
            ENDIF
         ELSE
            CALL BMG2_SER_PerSymStd_SOLVE_cg( 
     &                          Q(p_U), QF(p_U), Nx, Ny,
     &                          ABD, BBD, NCBW, NCU, IBC
     &                          )
         ENDIF

         RES_L2=rZERO
         IF ( BMG_IOFLAG(iBMG2_SER_OUT_ITERATIONS) ) THEN
            WRITE (*,230) KFMG, RES_L2
         ENDIF
            

 100     CONTINUE    ! >>>>>>>> LOOP BOUNDARY: multigrid F-cycle

            KFMG=KFMG+1
            !
            ! FMG: How many points in the stencil
            !
            IF ( KFMG.EQ.KF .AND. IFD.EQ.1 ) THEN
               NStncl=3
            ELSE
               NStncl=5
            ENDIF
            ! 
            ! FMG: Interpolate
            !                       
            CALL BMG2_SER_SymStd_GET_pointers(
     &                KFMG, IGRD, NOGm, Nxf, Nyf, 
     &                p_U, p_SO, p_SOR, p_CI
     &                )
            CALL BMG2_SER_SymStd_GET_pointers(
     &                KFMG-1, IGRD, NOGm, Nxc, Nyc,
     &                p_UC, p_SOC, p_SORC, p_CIC
     &                )
            
            IF ( IBC.EQ.BMG_SER_BCs_definite ) THEN
               CALL BMG2_SER_SymStd_interp_add(
     &                   KFMG-1, KFMG, Q(p_U), Q(p_UC),
     &                   RES(p_U), SO(p_SO), CI(p_CIC),
     &                   Nxc, Nyc, Nxf, Nyf, NStncl, IBC
     &                   )
            ELSE
               CALL BMG2_SER_PerSymStd_interp_add(
     &                   KFMG-1, KFMG, Q(p_U), Q(p_UC), SOR(p_SOR), 
     &                   CI(p_CIC), Nxc, Nyc, Nxf, Nyf, NSORv, IBC  
     &                   )
            ENDIF

            !
            ! FMG; Perform an n-cycle 
            !
            IF ( BMG_IOFLAG(iBMG2_SER_OUT_ITERATIONS) ) THEN 
               WRITE(*,410) '*** Performing n-cycle:  Coarsest =', KC,
     &                      'Finest =', KFMG
            ENDIF
            CALL BMG2_SER_SymStd_ncycle( 
     &                KC, KFMG, KF, IFD, IBC, IU, ID, IVW,
     &                IRELAX, IRELAX_SYM, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                ABD, BBD, NCBW, NCU, 
     &                IGRD, NOGm, RES_L2,
     &                BMG_iPARMS              
     &                )
            IF (BMG_iPARMS(id_BMG2_SER_Err_Code) .ne. iZERO) THEN
               RETURN
            ENDIF

            IF ( KFMG.EQ.KF ) THEN
               !
               ! Compute the final residual
               !
               CALL BMG2_SER_SymStd_residual( 
     &                   KF, SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &                   Nxf, Nyf, KF, IFD, NStncl, IBC,
     &                   IRELAX, IRELAX_SYM, BMG_SER_UP
     &                   )
               CALL BMG2_SER_SymStd_UTILS_norm_l2( 
     &                   RES(p_U), Nxf, Nyf, RES_L2 
     &                   )
               
               ! Note this is still vulnerable to very small RES_L2_0
               IF ( RES_L2_0.EQ.rZERO ) THEN
                  WRITE(*,505) '*** WARNING : boxmg.f '
                  WRITE(*,505) '    Initial residual has become zero!!!'
                  
                  CALL BMG2_SER_SymStd_ErrTrap(BMG_iPARMS,3)
                  RETURN
                  
               ELSE
                  REL_RES_L2 = RES_L2/RES_L2_0
               ENDIF

               !
               ! Output
               !
               IF ( ( ISTOP.EQ.BMG_SER_STOP_REL_RES_L2 )
     &            .AND. BMG_IOFLAG(iBMG2_SER_OUT_ITERATIONS) ) THEN
                  WRITE (*,400) '*** ITERATION =', NCYC, 
     &                          '*** RELATIVE RESIDUAL =', REL_RES_L2
               ELSE IF ( ( ISTOP.EQ.BMG_SER_STOP_ABS_RES_L2 ) 
     &                 .AND. BMG_IOFLAG(iBMG2_SER_OUT_ITERATIONS) ) THEN 
                  WRITE (*,400) '*** ITERATION =', NCYC, 
     &                          '*** ABSOLUTE RESIDUAL =', RES_L2
               ENDIF
               !
               IF ( ( ISTOP.EQ.BMG_SER_STOP_REL_RES_L2 ) 
     &            .AND. ( REL_RES_L2.LT.TOL )     ) THEN
                  ! FMG: converged
                  !      - set TOL and skip n-cycles
                  TOL = REL_RES_L2
                  NCYCLE_FLAG=.FALSE.
               ELSE IF ( ( ISTOP.EQ.BMG_SER_STOP_ABS_RES_L2 ) 
     &                 .AND. ( RES_L2.LT.TOL )         ) THEN
                  ! FMG: converged
                  ! - set TOL and skip n-cycles
                  TOL = RES_L2
                  NCYCLE_FLAG=.FALSE.
               ELSE IF ( NCYC.GE.MCYC ) THEN
                  ! FMG: cycle limit reached 
                  ! - set TOL and skip n-cycles
                  IF ( ISTOP.EQ.BMG_SER_STOP_REL_RES_L2 ) THEN
                     TOL = -REL_RES_L2
                  ELSE
                     TOL = -RES_L2
                  ENDIF
                  NCYCLE_FLAG=.FALSE.
               ELSE 
                  ! FMG: failed to converge
                  !      - increase counter and continue with n-cycles
                  NCYC=NCYC+1
               ENDIF
               !
            ELSE
               ! FMG: Continue the F-cycle
               GOTO 100
            ENDIF
            
 120     CONTINUE    ! <<<<<<<< LOOP BOUNDARY: multigrid F-cycle

      ENDIF

C ----------------------------------------------------------------
C     Perform multigrid n-cycle(s):
C ----------------------------------------------------------------

      IF ( NCYCLE_FLAG ) THEN

         CALL BMG2_SER_SymStd_GET_pointers( 
     &             KF, IGRD, NOGm, Nxf, Nyf,
     &             p_U, p_SO, p_SOR, p_CI
     &             )

 130     CONTINUE    ! >>>>>>>> LOOP BOUNDARY: multigrid n-cycles

           !
           ! Perform a multigrid n-cycle
           !
           CALL BMG2_SER_SymStd_ncycle( 
     &               KC, KF, KF, IFD, IBC, IU, ID, IVW,
     &               IRELAX, IRELAX_SYM, BMG_IOFLAG,
     &               Q, QF, RES, NF, NC,
     &               SO, NSO, SOR, NSOR, CI, NCI,
     &               ABD, BBD, NCBW, NCU, 
     &               IGRD, NOGm, RES_L2,
     &               BMG_iPARMS                 
     &               )

            IF (BMG_iPARMS(id_BMG2_SER_Err_Code) .ne. iZERO) THEN
               RETURN
            ENDIF

           !
           ! Compute the final residual
           !
           CALL BMG2_SER_SymStd_residual( 
     &               KF, SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &               Nxf, Nyf, KF, IFD, NStncl, IBC, 
     &               IRELAX, IRELAX_SYM, BMG_SER_UP
     &               )
           CALL BMG2_SER_SymStd_UTILS_norm_l2( 
     &               RES(p_U), Nxf, Nyf, RES_L2 
     &               )

           ! Note this is still vulnerable to very small RES_L2_0
           IF ( RES_L2_0.EQ.rZERO ) THEN

              WRITE(*,505) '*** FATAL ERROR: boxmg.f '
              WRITE(*,505) '    Initial residual has become zero!!!'

              CALL BMG2_SER_SymStd_ErrTrap(BMG_iPARMS,3)
              RETURN

           ELSE
              REL_RES_L2 = RES_L2/RES_L2_0
           ENDIF

           !
           ! Output
           !
           IF ( ( ISTOP.EQ.BMG_SER_STOP_REL_RES_L2 )
     &        .AND. BMG_IOFLAG(iBMG2_SER_OUT_ITERATIONS) ) THEN
              WRITE (*,400) '*** ITERATION =', NCYC, 
     &                      '*** RELATIVE RESIDUAL =', REL_RES_L2
           ELSE IF ( ( ISTOP.EQ.BMG_SER_STOP_ABS_RES_L2 )
     &             .AND. BMG_IOFLAG(iBMG2_SER_OUT_ITERATIONS) ) THEN 
              WRITE (*,400) '*** ITERATION =', NCYC, 
     &                      '*** ABSOLUTE RESIDUAL =', RES_L2
           ENDIF

           !
           ! Check Convergence 
           !
           IF ( ( ISTOP.EQ.BMG_SER_STOP_REL_RES_L2 )
     &        .AND. ( REL_RES_L2.LT.TOL )     ) THEN
              ! 
              ! n-cycles have converged in the RELATIVE RESIDUAL
              ! - set TOL and return
              !
              TOL = REL_RES_L2
           ELSE IF ( ( ISTOP.EQ.BMG_SER_STOP_ABS_RES_L2 ) 
     &             .AND. ( RES_L2.LT.TOL )         ) THEN
              ! 
              ! n-cycles have converged in the ABSOLUTE RESIDUAL
              ! - set TOL and return
              !
              TOL = REL_RES_L2
           ELSE IF ( NCYC.GE.MCYC ) THEN
              !
              ! n-cycles failed to converge, cycle limit reached
              ! - set TOL and return
              !
              IF ( ISTOP.EQ.BMG_SER_STOP_REL_RES_L2 ) THEN
                 TOL = -REL_RES_L2
              ELSE
                 TOL = -RES_L2
              ENDIF
           ELSE 
              !
              ! increase cycle count and continue n-cycling
              !
              NCYC=NCYC+1
              GOTO 130
           ENDIF

 140     CONTINUE    ! <<<<<<<< LOOP BOUNDARY: multigrid n-cycles

      ENDIF

 150  CONTINUE       ! <<<<<<<< FINAL I/O
      !

      BMG_rPARMS(id_BMG2_SER_STOP_TOL)  = TOL

      !
      !  Compute the solve time and update the total time
      !
      CALL BMG_timer(T2)
      T=T+T2-T1

      ! Output the multigrid cycling time
      IF ( BMG_IOFLAG(iBMG2_SER_OUT_TIME_CYCLING) ) THEN
         WRITE(*,240) '(2D) MULTIGRID CYCLING TIME =', T2-T1
      ENDIF

      IF ( BMG_IOFLAG(iBMG2_SER_OUT_TIME_TOTAL) ) THEN
         WRITE(*,240) '(2D) TOTAL TIME = ', T
      ENDIF


C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END: MULTIGRID CYCLING <<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      
C -------------------------------------
C     Output
C -------------------------------------         
      
      IF ( BMG_IOFLAG(iBMG_SER_OUT_SOLUTION) ) THEN
         ! print the solution on the finest grid
         ! in the current working directory.
         CALL BMG2_SER_SymStd_GET_pointers(
     &             KF, IGRD, NOGm, Nx, Ny,
     &             p_U, p_SO, p_SOR, p_CI
     &             )
         CALL BMG2_SER_SymStd_DUMP_vector( 
     &             BMG_IOFLAG, Q, Nx, Ny, KF, NOG,
     &             '', 'Q-solution', .FALSE.
     &             )
      ENDIF

      IF ( BMG_IOFLAG(iBMG_SER_OUT_RHS) ) THEN  
         ! print QF on the first coarse grid
         ! in the current working directory.
         CALL BMG2_SER_SymStd_GET_pointers(
     &             KF-1, IGRD, NOGm, Nx, Ny,
     &             p_UC, p_SOC, p_SORC, p_CIC
     &             )
         CALL BMG2_SER_SymStd_DUMP_vector( 
     &             BMG_IOFLAG, QF(p_UC), Nx, Ny, KF-1, NOG,
     &             '', 'QF-source', .FALSE.
     &             )
      ENDIF

C =======================================================================

 220  FORMAT (/,2X,A,1X,F9.3,/)
 230  FORMAT (' LEVEL',I2,' RESIDUAL NORM= ',1P,E10.3)
 240  FORMAT (/,2X,A,1X,F12.3,/)
 260  FORMAT (/,/,2X,A,1X,I5)
 270  FORMAT (2X,A,1X,I5,/)

 400  FORMAT (1X,A,1X,I2,4X,A,1X,1P,E16.9)
 404  FORMAT (1X,A,1X,1P,E16.9,/)
 405  FORMAT (/,1X,A,1X,1P,E16.9,/)
 410  FORMAT (1X,A,1X,I2,4X,A,1X,I2)

 500  FORMAT (/,2X,A)
 505  FORMAT (5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)

C ============================

      RETURN
      END
