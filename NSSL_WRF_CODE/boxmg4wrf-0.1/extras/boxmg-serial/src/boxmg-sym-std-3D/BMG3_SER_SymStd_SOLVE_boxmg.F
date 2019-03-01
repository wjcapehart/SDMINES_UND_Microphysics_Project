      SUBROUTINE BMG3_SER_SymStd_SOLVE_boxmg(
     &                Nxm, Nym, Nzm,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm, 
     &                SO, NSOm, SOR, NSORm, CI, NCIm, 
     &                ABD, BBD, NCBWm, NCUm, 
     &                IGRD, NOGm, NOG,
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                BMG_rWORK_PL, NBMG_rWORK_PLm
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_SOLVE_boxmg is a black box matrix solver. It takes a
C     matrix defined by the user from a given (fine) grid and constructs
C     coarser grids and their associated coefficient matrices. BMG3D
C     then performs the multigrid algorithm and returns the solution
C     vec- tor. The matrix is set by a user written subroutine. The
C     difference stencil at (i,j,k) spans three planes. The k-plane has
C     difference stencil
C
C                        xnw xn xne
C                        xw  xo xe
C                        xsw xs xse.
C
C     The (k+1)-plane has difference stencil
C
C                        xtnw xtn xtne
C                        xtw  xt  xte
C                        xtsw xts xtse.
C
C     The(k-1)-plane has difference stencil
C
C                        xbnw xbn xbne
C                        xbw  xb  xbe
C                        xbsw xbs xbse.
C
C     Here 
C           xo   = SO(i,j,k,ko)
C           xw   = - SO(i,j,k,kw)
C           xs   = - SO(i,j,k,ks)
C           xb   = - SO(i,j,k,kb)
C           xsw  = - SO(i,j,k,ksw)
C           xnw  = - SO(i,j+1,k,knw)
C           xbw  = - SO(i,j,k,kbw)
C           xbnw = - SO(i,j+1,k,kbnw)
C           xbn  = - SO(i,j+1,k,kbn)
C           xbne = - SO(i+1,j+1,k,kbne)
C           xbe  = - SO(i+1,j,k,kbe)
C           xbse = - SO(i+1,j,k,kbse)
C           xbs  = - SO(i,j,k,kbs)
C           xbsw = - SO(i,j,k,kbsw)
C
C     The remaining coefficents:
C
C     xn, xne,xe,xse, xt, xtw, xtnw, xtn, xtne, xte, xtse, xts, and xtsw 
C
C     are specified by symmetry. In the case of a seven point operator,
C     only xo, xw, xs, and xb need to be set. The difference scheme is
C     assumed to be positive definite.
C     
C     The fictitious points 
C
C         (1,j,k), j=1,jj, k=1,kk, (ii,j,k), j=1,jj, k=1,kk, 
C         (i,1,k), i=1,ii, k=1,kk, (i,jj,k), i=1,ii, k=1,kk, 
C         (i,j,1), i=1,ii, j=1,jj, (i,j,kk), i=1,ii, j=1,jj,
C
C     are assumed for ease of programming. However, it is also assumed
C     that the user has set any difference coefficient referring
C     to these points to zero. For example, if the above difference
C     template is centered at (2,2), then xb, xbw, xbnw, xbn, xbne, xbe,
C     xbse, xbs xnw, xw, xsw, xs, xse, xtnw, xtw, xtsw, xts, xtse
C     would be assumed zero.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     NXM      x-dimension of the grid, excluding fictitious points
C     NYM      y-dimension of the grid, excluding fictitious points
C     NZM      z-dimension of the grid, excluding fictitious points
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
C     IFD      Indicator for difference scheme. IFD=1 means a seven
C              point operator. I.e., in the template given above, xnw
C              xne, xse, xsw, xbw, xbnw, xbn, xbne, xbe, xbse, xbs, xbsw,
C              xtw, xtnw, xtn, xtne, xte, xtse, xts, xtsw.
C              are assumed zero. IFD.ne.1. means a 27 point
C              operator. (Note that 27 point operators will be
C              generated on the coarser grids even in the case of a seven
C              point operator.
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
C              -ISTRT cycles unless the error criterion is satisfied
C              first.
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
C              IRELAX=5  => alternating red-black plane Gauss-Seidel relaxation
C                        -  ordering is xy planes, yz planes, xz planes
C                        -  for symmetric relaxation this is reversed
C                           as necessary when moving "up" or "down" levels
C
C     IRELAX_SYM  Symmetric cycles require that the relaxation ordering
C                 be reversed between the restriction and interpolation
C                 stage of the cycle.   This is important when BMG3D
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
C              computed by BMG3_SER_SymStd_SETUP_space and is checked in boxmg.
C
C     NCm      Maximum storage for a vector on all coarse grids.  This should
C              be computed in BMG3_SER_SymStd_SETUP_space and is checked in boxmg.
C
C     NSOm     See INPUT/OUTPUT
C
C     NSORm    See INPUT/OUTPUT
C
C     NCIm     See INPUT/OUTPUT
C     
C     NOGm     Maximum number of grids that can be supported by the 
C              pointer array IGRD. It should be set in calling program.
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
C              IFD .EQ. 1  =>    NSOm = 4*NFm + 10*NCm 
C              IFD .NE. 1  =>    NSOm = 14*NFm 
C
C              BMG3_SER_SymStd_SOLVE_boxmg verifies that NSOm is large enough.
C
C     CI       CI is a real array which should be dimensioned
C              to (NCIm) in calling program. It is used to contain
C              the interpolation coefficients.
C
C     NCIm     Dimension of CI, set in calling program.  NCIm = 26*NCm is 
C              enough.  BOXMG checks to see if NCIm is large enough.
C
C     QF       The user defined array that contains the right hand side.
C              It must be dimensioned to at least NFMAX.
C
C     Q        the user defined array that contains the solution vector.
C              It must be dimensioned to at least NFMAX.
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C
C =======================================================================
C  --------------------
C   WORK ARRAYS:
C  --------------------
C
C     ABD      User declared two dimensional real array, which
C              is used to store the coefficient matrix for the coarsest
C              grid. It is then used by the LAPACK routines. It
C              should be dimensioned to (NCBWm,NCUm).
C
C     BBD      User declared real array of dimension NCUm for use in the
C              the LAPACK routines
C
C     NCBWm    Maximum first subscript of ABD, which needs to be > or =
C              (x-dimension+1)*(y-dimension+2) on coarsest grid. 
C
C     NCUm     Maximum second subscript of ABD, which needs to be > or =
C              (x-dimension+1)*(y-dimension+1)*(z-dimension+1) on
C              coarsest grid.
C
C     IGRD     A work array. IGRD should be dimensioned to
C              IGRD(NOGm,NBMG_SER_pIGRD) in the calling program.
C
C     SOR      SOR is a real array, which should be dimensioned to
C              (NSORm) in the calling program. It is used to store
C              residuals, reciprocals of SO(.,1) if IRELAX=1, and lu
C              decompositions if IRELAX.gt.1.
C
C     NSORm    Dimension of SOR, set in calling program. NSORm=2*NFMAX
C              is enough for IRELAX=1 or 2. For IRELAX=3 or 4, 4*NFMAX
C              is enough. 
C
C
C =======================================================================
C  --------------------
C   REFERENCES:
C  --------------------
C
C     Dendy, J. E. Jr., "Two Multigrid Methods for Three Dimensional
C     Problems with Discountinuos and Anisotropic Coefficients", SIAM
C     Journal of Scientific and Satatistical Computing, Vol. 8, No. 2,
C     September 1987
C     
C ==========================================================================
C  --------------------
C   Error Codes:
C  --------------------
C
C     -------------------------------------------------------------------
C     CODE    DESCRIPTION                           ORIGIN
C     -------------------------------------------------------------------
C      1     initial residual has become zero     SOLVE_boxmg
C      2     NOG = 0                              SETUP_parts
C      3     NOG < 0                              SETUP_parts
C      4     NStncl out of range                  SETUP_cg_LU
C      5     Cholesky decomposition failed, look
C            in BMG_iPARMS(id_BMG3_SER_Ext_Err_Code)
C            for the return code of the LAPACK
C            routine                              SETUP_cg_LU
C      6    NOG = 1                               SETUP_relax
C      7    NOG = 0                               SETUP_relax
C      8    NOG < 0                               SETUP_relax
C      9    fatal setup error                     SETUP_PtrGrid
C     10    IRELAX out of range                   SETUP_PtrGrid
C     11    min coarse grid dim < 3 ( NXYZc )     SETUP_space
C     12    min coarse grid dim < 3 ( NXYc )      SETUP_space
C     13    IRELAX out of range                   SETUP_space
C     14    computed number of grids too small    SETUP_space
C     15    computed number of grids < 1          SETUP_space
C     16    fatal error (check output)            SETUP_PtrWork
C     17    memory allocation mode unspecified    SETUP_PtrWork
C     18    KCF = KC                              ncycle
C     19    updown out of range                   updown
C     20    coarse grid solve failed, look in
C           BMG_iPARMS(id_BMG3_SER_Ext_Err_Code)
C           for the return code of the LAPACK
C           routine                               SOLVE_cg
C     21    inconsistent # of stencil entries     COPY_SO_xy
C     22    inconsistent # of stencil entries     COPY_SO_yz
C     23    inconsistent # of stencil entries     COPY_SO_xz
C     24    invalid cg_operator type              SETUP_parts
C     25    invalid cg_operator construction      SETUP_parts
C
C ==========================================================================

      IMPLICIT NONE

C ------------------------------------------------
C     Includes
C
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_stencils.h'
      INCLUDE 'BMG_SER_workspace.h'
      INCLUDE 'BMG_SER_parameters.h'

C ------------------------------------------------
C     Argument Declarations
C
      INTEGER  NBMG_iWORK_PLm, NBMG_rWORK_PLm, NOG, NOGm

      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS),
     &         BMG_iWORK_PL(NBMG_iWORK_PLm), NCIm, NSOm, NSORm,
     &         IGRD(NOGm,NBMG_SER_pIGRD), NCBWm, NCUm, NCm, NFm,
     &         nxm, nym, nzm
      REAL*RKIND   abd(NCBWm*NCUm), bbd(NCUm), 
     &         BMG_rPARMS(NBMG_SER_rPARMS), 
     &         BMG_rWORK_PL(NBMG_rWORK_PLm), ci(NCIm), q(NFm), 
     &         qf(NFm), RES(NFm), SO(NSOm), sor(NSORm), tol
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)

C ------------------------------------------------
C     Local Declarations
C
      INTEGER i, iBC, id, ifd, ii, iic, iif, im, irelax,
     &        IRELAX_SYM, ISETUP, ISTOP, istrt, istrt2, iu, ivw,
     &        jj, jjc, jjf, k, kc, kfmg, kf, kk, kkc, kkf,
     &        mcyc, mcycl, ncyc, NStncl,
     &        p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC,
     &        NC, NCBW, NCI, NCU, NF, NSO, NSOR, NXYZc,
     &        NBMG_iWORK_PL, NBMG_rWORK_PL

      REAL*8  REL_RES_L2, RES_L2, RES_L2_0
      REAL*8 t, t1, t2
      LOGICAL NCYCLE_FLAG

C ==========================================================================

      IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
         RETURN
      END IF


C -------------------
C     Zero Times:
C -------------------

      t=dZERO
      t1=dZERO
      t2=dZERO

      IF ( BMG_IOFLAG(iBMG3_SER_BUG_PARAMETERS) ) THEN
         CALL BMG3_SER_SymStd_DUMP_parms( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &             )
      ENDIF

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      BMG_iPARMS(id_BMG3_SER_BC) = BMG_SER_BCs_definite
      BMG_iPARMS(id_BMG2_SER_BC) = BMG_SER_BCs_definite

      IFD    = BMG_iPARMS(id_BMG3_SER_STENCIL)
      IBC    = BMG_iPARMS(id_BMG3_SER_BC)
      ISETUP = BMG_iPARMS(id_BMG3_SER_SETUP)
      NXYZc  = BMG_iPARMS(id_BMG3_SER_CG_MIN_DIM)

      IRELAX = BMG_iPARMS(id_BMG3_SER_RELAX)
      IRELAX_SYM = BMG_iPARMS(id_BMG3_SER_RELAX_SYM )

      ID = BMG_iPARMS(id_BMG3_SER_NRELAX_DOWN)
      IU = BMG_iPARMS(id_BMG3_SER_NRELAX_UP)
      IM = BMG_iPARMS(id_BMG3_SER_NRELAX_FG)

      IF (BMG_iPARMS(id_BMG3_SER_CYCLE_CLASS).EQ.BMG_SER_N_CYCLE) THEN
         ISTRT = - BMG_iPARMS(id_BMG3_SER_MAX_ITERS)
      ELSEIF ( BMG_iPARMS(id_BMG3_SER_CYCLE_CLASS)
     &        .EQ.BMG_SER_FMG_CYCLE ) THEN 
         ISTRT = BMG_iPARMS(id_BMG3_SER_MAX_ITERS)
      ENDIF

      IVW   = BMG_iPARMS(id_BMG3_SER_NCYCLE_TYPE)
      MCYCL = BMG_iPARMS(id_BMG3_SER_FMG_NNCYCLE)

      ISTOP = BMG_iPARMS(id_BMG3_SER_STOP_TEST)
      TOL   = BMG_rPARMS(id_BMG3_SER_STOP_TOL)

C ================================================

      IF ( BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS).EQ.BMG_SER_N_CYCLE ) THEN
         ISTRT2 = - BMG_iPARMS(id_BMG2_SER_MAX_ITERS)
      ELSEIF ( BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS)
     &        .EQ.BMG_SER_FMG_CYCLE) THEN 
         ISTRT2 = BMG_iPARMS(id_BMG2_SER_MAX_ITERS)
      ENDIF

C ================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> BEGIN:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      IF ( ISETUP.EQ.BMG_SER_SETUP_ptrs_opers .OR. 
     &     ISETUP.EQ.BMG_SER_SETUP_only ) THEN

         CALL BMG3_SER_SymStd_SETUP_PtrGrid( 
     &             Nxm, Nym, Nzm, BMG_iPARMS, BMG_iWORK_PL, 
     &             NOGm, NFm, NCm, NSOm, NSORm, NCIm, 
     &             NCBWm, NCUm, NBMG_iWORK_PLm, NBMG_rWORK_PLm,
     &             NOG, NF, NC, NSO, NSOR, NCI, NCBW, NCU, IGRD,
     &             NBMG_iWORK_PL, NBMG_rWORK_PL, BMG_IOFLAG
     &             )
         IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF
         

         ! Output      
         IF( BMG_IOFLAG(iBMG3_SER_OUT_WSPACE_SIZE) ) THEN
            WRITE (*,260) 
     &            'Storage for a vector on all grids = ', NF
         ENDIF

         ! Output
         IF( BMG_IOFLAG(iBMG3_SER_OUT_WSPACE_SIZE) ) THEN
            WRITE (*,270) 
     &            'Storage for a vector on all coarse grids = ', NC
         ENDIF

      ELSE

         !
         ! No allocation setup so we use dimensions provided by user.
         !

         NF   = NFm
         NC   = NCm
         NSO  = NSOm
         NSOR = NSORm
         NCI  = NCIm
         
         NCBW = NCBWm
         NCU  = NCUm

         NBMG_iWORK_PL = NBMG_iWORK_PLm
         NBMG_rWORK_PL = NBMG_rWORK_PLm

      ENDIF

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>> END:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>> BEGIN:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

 89   CONTINUE  ! >>>>>>>>> BEGIN: Component Setup

      IF ( ISETUP.EQ.BMG_SER_SETUP_ptrs_opers .OR. 
     &     ISETUP.EQ.BMG_SER_SETUP_opers .OR. 
     &     ISETUP.EQ.BMG_SER_SETUP_only ) THEN 

         !  Start the timer
         CALL BMG_timer(t1)

         !
         ! Construct coarse-grid and interpolation operators.  In 
         ! addition, setup any necessary components for relaxation.
         !
         CALL BMG3_SER_SymStd_SETUP_parts( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             SO, NSOm, SOR, NSORm, CI, NCIm,
     &             ABD, BBD, NCBWm, NCUm, IGRD, NOGm, NOG,
     &             BMG_iWORK_PL, NBMG_iWORK_PL,
     &             BMG_rWORK_PL, NBMG_rWORK_PL
     &             )
         IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF


         ! Compute the setup time
         CALL BMG_timer(t2)
         t=t+t2-t1

         ! Output the setup time
         IF( BMG_IOFLAG(iBMG3_SER_OUT_TIME_SETUP) )  THEN
            WRITE (*,220) '(3D) SETUP TIME =', T
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

 95   CONTINUE

      ! Start the timer
      CALL BMG_timer(t1)

C ----------------------------------------------------------------
C     Initial residual:
C ----------------------------------------------------------------

      !
      ! Get pointers for the finest grid (GRID LEVEL = NOG)
      !
      CALL BMG3_SER_SymStd_GET_pointers( 
     &          NOG, IGRD, NOGm, p_U, p_SO, p_SOR, 
     &          p_CI, II, JJ, KK 
     &          )

      !
      ! Set NStncl for the finest grid
      !
      IF ( IFD.NE.1 ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF
      
      !
      ! Compute the initial residual
      !
      CALL BMG3_SER_SymStd_residual( 
     &          KF, KF, IFD, Q(p_U), QF(p_U), SO(p_SO), RES(p_U),
     &          II, JJ, KK, NStncl
     &          )
      !
      CALL BMG3_SER_SymStd_UTILS_norm_l2(
     &          RES(p_U), II, JJ, KK, RES_L2_0 
     &          )

      !
      !  Output
      !
      IF ( BMG_IOFLAG(iBMG3_SER_OUT_ITERATIONS) ) THEN
         WRITE (*,405) '*** THE INITIAL RESIDUAL (l2-NORM) = ', RES_L2_0
      ENDIF

      IF ( RES_L2_0 .EQ. dZERO ) THEN
         IF ( BMG_IOFLAG(iBMG3_SER_WARN_ZERO_RESIDUAL) ) THEN
            WRITE(*,505) '*** WARNING: bmg3d.f '
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
         IF ( BMG_IOFLAG(iBMG3_SER_OUT_ITERATIONS) ) THEN
            WRITE(*,*) '*** PERFORMING EXACT SOLVE ON FINE GRID '
         ENDIF

         CALL BMG3_SER_SymStd_SOLVE_cg(
     &             Q(p_U), QF(p_U), II, JJ, KK,
     &             ABD, BBD, NCBW, NCU,
     &             BMG_IOFLAG, BMG_iPARMS
     &             )
         IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         !
         ! Compute the final residual
         !
         CALL BMG3_SER_SymStd_residual( 
     &             NOG, NOG, IFD, Q(p_U), QF(p_U), SO(p_SO),
     &             RES(p_U), II, JJ, KK, NStncl
     &             )
         !
         CALL BMG3_SER_SymStd_UTILS_norm_l2(
     &             RES(p_U), II, JJ, KK, RES_L2 
     &             )
         !
         IF( BMG_IOFLAG(iBMG3_SER_BUG_RES_CG_SOLVE) ) THEN 
            WRITE (*,405) '*** THE FINAL RESIDUAL (l2-NORM) = ', RES_L2
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

      KF = NOG                  ! finest grid index
      KC = 1                    ! coarsest grid index
      KFMG = 1                  ! FMG current finest grid
      MCYC = ABS(ISTRT)         ! Maximum number of multigrid cycles
      NCYC = 1                  ! multigrid cycle counter
      NCYCLE_FLAG = .TRUE.      ! default to do n-cycles

C ----------------------------------------------------------------
C     Perform multigrid F-cycle:
C     (if necessary, then continue with n-cycles)
C ----------------------------------------------------------------

      IF ( ISTRT.GT.0 ) THEN
      
         ! Coarsen the right hand side
         DO k=KF-1, 1, -1
            CALL BMG3_SER_SymStd_GET_pointers( 
     &                k+1, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, IIF, JJF, KKF 
     &                )
            CALL BMG3_SER_SymStd_GET_pointers( 
     &                k, IGRD, NOGm,
     &                p_UC, p_SOC, p_SORC, p_CIC, IIC, JJC, KKC 
     &                )
            ! restrict the right hand side 
            CALL BMG3_SER_SymStd_restrict(
     &                k+1, k, QF(p_U), QF(p_UC), CI(p_CIC),
     &                IIF, JJF, KKF, IIC, JJC, KKC
     &                )
            ! Copy right hand side into RES(i,j,k) for interpolation
            DO i=0, IIF*JJF*KKF-1
               RES(p_U+i)=QF(p_U+i)
            END DO
         END DO

         ! FMG current finest grid index is set to the coarsest grid
         KFMG=KC     
         ! FMG: solve exactly on the coarsest grid
         CALL BMG3_SER_SymStd_GET_pointers(
     &             KFMG, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI, II, JJ, KK 
     &             )
         CALL BMG3_SER_SymStd_SOLVE_cg( 
     &             q(p_U), qf(p_U), II, JJ, KK, 
     &             abd, bbd, NCBWm, NCUm,
     &             BMG_IOFLAG, BMG_iPARMS 
     &             )
         IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         RES_L2=dZERO
         IF ( BMG_IOFLAG(iBMG3_SER_OUT_ITERATIONS) ) THEN
            WRITE (*,230) KFMG, RES_L2
         ENDIF

 100     CONTINUE    ! >>>>>>>> LOOP BOUNDARY: multigrid F-cycle

            KFMG=KFMG+1
            !
            ! FMG: How many points in the stencil
            !
            IF ( KFMG.EQ.KF .AND. IFD.EQ.1 ) THEN
               NStncl=4
            ELSE
               NStncl=14
            ENDIF
            ! 
            ! FMG: Interpolate
            !
            CALL BMG3_SER_SymStd_GET_pointers( 
     &                KFMG, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, ii, jj, kk 
     &                )
            CALL BMG3_SER_SymStd_GET_pointers( 
     &                KFMG-1, IGRD, NOGm,
     &                p_UC, p_SOC, p_SORC, p_CIC, iic, jjc, kkc 
     &                )
            CALL BMG3_SER_SymStd_interp_add( 
     &                k-1, k, Q(p_U) ,Q(p_UC),
     &                SO(p_SO), RES(p_U), CI(p_CIC),
     &                IIC, JJC, KKC, II, JJ, KK, NStncl
     &                )

            !
            ! FMG; Perform an n-cycle 
            !
            IF ( BMG_IOFLAG(iBMG3_SER_OUT_ITERATIONS) ) THEN 
               WRITE(*,410) '*** Performing n-cycle:  Coarsest =', KC,
     &                      'Finest =', KFMG
            ENDIF
            CALL BMG3_SER_SymStd_ncycle(
     &                KC, KFMG, KF, IFD, IU, ID, IVW,
     &                IRELAX, IRELAX_SYM,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm, 
     &                SO, NSOm, SOR, NSORm, CI, NCIm,
     &                ABD, BBD, NCBWm, NCUm, 
     &                IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )

            IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF
            


            IF ( KFMG.EQ.KF ) THEN
               !
               !  Compute the final residual
               !
               CALL BMG3_SER_SymStd_residual( 
     &                   KF, KF, IFD, Q(p_U), QF(p_U), SO(p_SO),
     &                   RES(p_U), II, JJ, KK, NStncl
     &                   )
               CALL BMG3_SER_SymStd_UTILS_norm_l2( 
     &                   RES(p_U), II, JJ, KK, RES_L2 
     &                   )

               ! Note this is still vulnerable to very small RES_L2_0
               IF ( RES_L2_0.EQ.dZERO ) THEN
                  IF (BMG_IOFLAG(iBMG3_SER_OUT_STOP_ERROR)) THEN
                    WRITE(*,505)'*** FATAL ERROR: bmg3d.f '
                    WRITE(*,505)'    Initial residual has become zero!'
                  END IF
                  
                  CALL BMG3_SER_SymStd_ErrTrap(BMG_iPARMS,1)
                  RETURN
                  
               ELSE
                  REL_RES_L2 = RES_L2/RES_L2_0
               ENDIF

               !
               !  Output 
               !
               IF ( ( ISTOP.EQ.BMG_SER_STOP_REL_RES_L2 )
     &            .AND. BMG_IOFLAG(iBMG3_SER_OUT_ITERATIONS) ) THEN
                  WRITE (*,400) '*** ITERATION =', NCYC, 
     &                          '*** RELATIVE RESIDUAL = ', REL_RES_L2
               ELSE IF ( ( ISTOP.EQ.BMG_SER_STOP_ABS_RES_L2 ) 
     &                 .AND. BMG_IOFLAG(iBMG3_SER_OUT_ITERATIONS) ) THEN
                  WRITE (*,400) '*** ITERATION =', NCYC, 
     &                          '*** ABSOLUTE RESIDUAL = ', RES_L2
               ENDIF
               !
               !  Check Convergence
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
                  !      - set TOL and skip n-cycles
                  TOL = RES_L2
                  NCYCLE_FLAG=.FALSE.
               ELSE IF ( NCYC.GE.MCYC ) THEN
                  ! FMG: cycle limit reached 
                  !      - set TOL and skip n-cycles
                  IF ( ISTOP.EQ.BMG_SER_STOP_REL_RES_L2 ) THEN
                     TOL = -REL_RES_L2
                  ELSE
                     TOL = -RES_L2
                  ENDIF
                  NCYCLE_FLAG=.FALSE.
               ELSE 
                  ! FMG: didn't converge
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

         
         CALL BMG3_SER_SymStd_GET_pointers(
     &             KF, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI, II, JJ, KK 
     &             )
         
 130     CONTINUE   ! >>>>>>>> LOOP BOUNDARY: multigrid n-cycles


            !
            ! Perform a multigrid n-cycle
            !
            CALL BMG3_SER_SymStd_ncycle(
     &                KC, KF, KF, IFD, IU, ID, IVW,
     &                IRELAX, IRELAX_SYM,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm, 
     &                SO, NSOm, SOR, NSORm, CI, NCIm,
     &                ABD, BBD, NCBWm, NCUm, 
     &                IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )

            IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
            !  Compute the final residual
            !
            CALL BMG3_SER_SymStd_residual( 
     &                KF, KF, IFD, Q(p_U), QF(p_U), SO(p_SO),
     &                RES(p_U), II, JJ, KK, NStncl
     &                )
            CALL BMG3_SER_SymStd_UTILS_norm_l2( 
     &                RES(p_U), II, JJ, KK, RES_L2 
     &                )

            ! Note this is still vulnerable to very small RES_L2_0
            IF ( RES_L2_0.EQ.dZERO ) THEN
               IF (BMG_IOFLAG(iBMG3_SER_OUT_STOP_ERROR)) THEN
                  WRITE(*,505) '*** FATAL ERROR: bmg3d.f '
                  WRITE(*,505) '    Initial residual has become zero!'
               END IF

               CALL BMG3_SER_SymStd_ErrTrap(BMG_iPARMS,1)
               RETURN

            ELSE
               REL_RES_L2 = RES_L2/RES_L2_0
            ENDIF

            !
            !  Output 
            !
            IF ( ( ISTOP.EQ.BMG_SER_STOP_REL_RES_L2 )
     &           .AND. BMG_IOFLAG(iBMG3_SER_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** ITERATION =', NCYC, 
     &                       '*** RELATIVE RESIDUAL = ', REL_RES_L2
            ELSE IF ( ( ISTOP.EQ.BMG_SER_STOP_ABS_RES_L2 )
     &              .AND. BMG_IOFLAG(iBMG3_SER_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** ITERATION =', NCYC, 
     &                      '*** ABSOLUTE RESIDUAL = ', RES_L2
            ENDIF

            !
            ! Check Convergence 
            !
            IF ( ( ISTOP.EQ.BMG_SER_STOP_REL_RES_L2 )
     &         .AND. ( REL_RES_L2.LT.TOL )     ) THEN
               ! 
               ! n-cycles have converged in the RELATIVE RESIDUAL
               ! - set TOL and return
               !
               TOL = REL_RES_L2
            ELSE IF ( ( ISTOP.EQ.BMG_SER_STOP_ABS_RES_L2 ) 
     &              .AND. ( RES_L2.LT.TOL )         ) THEN
               ! 
               ! n-cycles have converged in the ABSOLUTE RESIDUAL
               ! - set TOL and return
               !
               TOL = RES_L2
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
      !
      BMG_rPARMS(id_BMG3_SER_STOP_TOL) = TOL
      !
      !  Compute the solve time and update the total time
      !
      CALL BMG_timer(T2)
      T=T+T2-T1

      ! Output the multigrid cycling time
      IF ( BMG_IOFLAG(iBMG3_SER_OUT_TIME_CYCLING) ) THEN
         WRITE(*,240) '(3D) MULTIGRID CYCLING TIME =', T2-T1
      ENDIF

      IF ( BMG_IOFLAG(iBMG3_SER_OUT_TIME_TOTAL) ) THEN
         WRITE(*,240) '(3D) TOTAL TIME = ', T
      ENDIF

C ==========================================================================

 220  FORMAT (/,1X,A,1X,F9.3,/)
 230  FORMAT ('(3D)LEVEL',I2,' RESIDUAL NORM= ',1P,E10.3)
 240  FORMAT (/,1X,A,1X,F12.3,/)
 260  FORMAT (/,/,2X,A,1X,I5)
 270  FORMAT (2X,A,1X,I5,/)

 400  FORMAT (1X,A,1X,I2,4X,A,1X,1P,E16.9)
 405  FORMAT (/,1X,A,1X,1P,E16.9,/)
 410  FORMAT (1X,A,1X,I2,4X,A,1X,I2)

 500  FORMAT (/,2X,A)
 505  FORMAT (5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
      
C ==========================

      RETURN
      END
