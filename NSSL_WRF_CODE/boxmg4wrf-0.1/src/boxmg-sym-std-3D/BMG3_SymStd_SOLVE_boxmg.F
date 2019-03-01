      SUBROUTINE BMG3_SymStd_SOLVE_boxmg(
     &                NLx, NLy, NLz, NGx, NGy, NGz, 
     &                iGs, jGs, kGs, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm, 
     &                SO, NSOm, SOR, NSORm, CI, NCIm,
     &                IGRD, NOGm, NOG,
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                BMG_rWORK_PL, NBMG_rWORK_PLm,
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                iWorkMSG, NMSGim, pMSG, pMSGSO, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                NProc, MSG_Buffer, NMSGrm, MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_SOLVE_boxmg is a black box matrix solver. It takes a
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
C                ISTOP = BMG_STOP_REL_RES_L2 
C
C                => iteration has converged if the l2-norm of the
C                   current residual, divided by the l2-norm of the 
C                   initial residual is less than TOL.
C
C              and an absolute residual test
C  
C                ISTOP = BMG_STOP_REL_ABS_L2
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
C                  file "BMG_parameters.h" for hints (the names are verbose.)
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
C                    IRELAX_SYM = BMG_RELAX_NONSYM
C                    => use the stated relaxation ordering (nonsymmetric)
C
C                    IRELAX_SYM = BMG_RELAX_SYM
C                    => make the cycle symmetric.
C
C                 as defined in BMG_parameters.h
C
C     ISEP     Indicator for separable equation such that the SO array can be
C              reduced to 3 1-D arrays in the interior, saving memory bandwidth
C
C               ISEP = BMG_NONSEPARABLE (0): Always use 3d SO arrays
C               ISEP = BMG_SEPARABLE (1): only use 3d SO arrays on boundaries.
C
C
C     IVW      Indicator for type of cycle. IVW=1 means v-cycles will
C              be performed. IVW=2 means w-cycles will be performed, etc.
C
C     MCYCL    Currently disabled!
C
C     NFm      Maximum storage for a vector on all grids. This should be
C              computed by BMG3_SymStd_SETUP_space and is checked in boxmg.
C
C     NCm      Maximum storage for a vector on all coarse grids.  This should
C              be computed in BMG3_SymStd_SETUP_space and is checked in boxmg.
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
C              this with a call to BMG2_SymStd_SETUP_space.  Clearly,
C            
C              IFD .EQ. 1  =>    NSOm = 4*NFm + 10*NCm 
C              IFD .NE. 1  =>    NSOm = 14*NFm 
C
C              BMG3_SymStd_SOLVE_boxmg verifies that NSOm is large enough.
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
C              IGRD(NOGm,NBMG_pIGRD) in the calling program.
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
C            in BMG_iPARMS(id_BMG3_Ext_Err_Code)
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
C           BMG_iPARMS(id_BMG3_Ext_Err_Code)
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
      INCLUDE 'geom_param_fort.h'
      INCLUDE 'mpi_param_fort.h'
      INCLUDE 'MSG.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ------------------------------------------------
C     Argument Declarations
C
      INTEGER  NLx, NLy, NLz, NGx, NGy, NGz, iGs, jGs, kGs,
     &         NFm, NCm, NSOm, NSORm, NCIm, NCBWm, NCUm, NOGm, NOG,
     &         NMSGim, NMSGrm, NProc, MPICOMM, NBMG_MSG_iGRID,
     &         NBMG_iWORK_CSm, NBMG_iWORK_PLm,
     &         NBMG_rWORK_CSm, NBMG_rWORK_PLm

      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  iWorkMSG(NMSGim),
     &         pMSG(NBMG_pMSG,NOGm),
     &         pMSGSO(NBMG_pMSG,NOGm)
      REAL*RKIND   MSG_Buffer(NMSGrm)

      INTEGER  BMG_iWORK_PL(NBMG_iWORK_PLm)
      REAL*RKIND   BMG_rWORK_PL(NBMG_rWORK_PLm)

      INTEGER  BMG_iWORK_CS(NBMG_iWORK_CSm)
      REAL*RKIND   BMG_rWORK_CS(NBMG_rWORK_CSm)

      INTEGER  IGRD(NOGm,NBMG_pIGRD)
      REAL*RKIND CI(NCIm),Q(NFm),QF(NFm),RES(NFm),SO(NSOm),SOR(NSORm)

C ------------------------------------------------
C     Local Declarations
C
      INTEGER i, iBC, id, ifd, iil, iilc, iilf, im, irelax, isep,
     &        IRELAX_SYM, ISETUP, ISTOP, istrt, istrt2, iu, ivw,
     &        jjl, jjlc, jjlf, k, kc, kfmg, kf, kkl, kklc, kklf,
     &        mcyc, mcycl, ncyc, NStncl, NMSGi, NMSGr, p_ABD, p_BBD,
     &        p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC,
     &        IIG, IIGC, IIGF, JJG, JJGC, JJGF, KKG, KKGC, KKGF,
     &        NC, NCBW, NCI, NCU, NF, NSO, NSOR, NXYZc,
     &        NBMG_iWORK_CS, NBMG_iWORK_PL, NBMG_rWORK_CS, 
     &        NBMG_rWORK_PL


      INTEGER  iGs_c, iGs_f, jGs_c, jGs_f, kGs_c, kGs_f

      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK,
     &         NProcI, NProcJ, NProcK,
     &         p_NLx, p_NLy, p_NLz,
     &         p_NLx_kg, p_NLy_kg, p_NLz_kg,
     &         p_ProcCoord, p_ProcGrid, pSI_MSG

      INTEGER  BMG_MSG_MyProc, MPI_IERR

      REAL*8  REL_RES_L2, RES_L2, RES_L2_0
      REAL*8 T, T1, T2, TOL
      REAL*8  TT1, TT2
      LOGICAL NCYCLE_FLAG, INIT_CHECK

C ==========================================================================

      IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
         RETURN
      END IF

C -------------------
C     Zero Times:
C -------------------

      t=dZERO
      t1=dZERO
      t2=dZERO

      IF ( BMG_IOFLAG(iBMG3_BUG_PARAMETERS) ) THEN
         CALL BMG3_SymStd_DUMP_parms( 
     &                    BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &                    )
      ENDIF

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------
C --------------------------------------------
C     MPI Parallel Data Decomposition:
C --------------------------------------------
      
      !
      ! Processor Grid:
      !
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_NLx       = BMG_MSG_pGRID(ip_BMG_MSG_NLx_Grid)
      p_NLy       = BMG_MSG_pGRID(ip_BMG_MSG_NLy_Grid)
      p_NLz       = BMG_MSG_pGRID(ip_BMG_MSG_NLz_Grid)
      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)
      
C ----------------------------------------
C     MSG workspace pointers:
C ----------------------------------------

      p_NLx_kg = 1
      p_NLy_kg = p_NLx_kg + NProcI*NOGm
      p_NLz_kg = p_NLy_kg + NProcJ*NOGm

      pSI_MSG  = p_NLz_kg + NProcK*NOGm  !!NB: initalization of pSI_MSG

C ----------------------------------------
C     Multigrid Cycle parameters:
C ----------------------------------------

      IFD    = BMG_iPARMS(id_BMG3_STENCIL)
      IBC    = BMG_iPARMS(id_BMG3_BC)
      ISETUP = BMG_iPARMS(id_BMG3_SETUP)

      IRELAX = BMG_iPARMS(id_BMG3_RELAX)
      IRELAX_SYM = BMG_iPARMS(id_BMG3_RELAX_SYM )

      ISEP = BMG_iPARMS(id_BMG3_SEPARABLE)

      ID = BMG_iPARMS(id_BMG3_NRELAX_DOWN)
      IU = BMG_iPARMS(id_BMG3_NRELAX_UP)
      IM = BMG_iPARMS(id_BMG3_NRELAX_FG)

      IF (BMG_iPARMS(id_BMG3_CYCLE_CLASS).EQ.BMG_N_CYCLE) THEN
         ISTRT = - BMG_iPARMS(id_BMG3_MAX_ITERS)
      ELSEIF (BMG_iPARMS(id_BMG3_CYCLE_CLASS).EQ.BMG_FMG_CYCLE) THEN 
         ISTRT = BMG_iPARMS(id_BMG3_MAX_ITERS)
      ENDIF

      IVW   = BMG_iPARMS(id_BMG3_NCYCLE_TYPE)
      MCYCL = BMG_iPARMS(id_BMG3_FMG_NNCYCLE)

      ISTOP = BMG_iPARMS(id_BMG3_STOP_TEST)
      TOL   = BMG_rPARMS(id_BMG3_STOP_TOL)

C ================================================

      IF (BMG_iPARMS(id_BMG2_CYCLE_CLASS).EQ.BMG_N_CYCLE) THEN
         ISTRT2 = - BMG_iPARMS(id_BMG2_MAX_ITERS)
      ELSEIF (BMG_iPARMS(id_BMG2_CYCLE_CLASS).EQ.BMG_FMG_CYCLE) THEN 
         ISTRT2 = BMG_iPARMS(id_BMG2_MAX_ITERS)
      ENDIF

C ================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C ------------------------------------------------------------------- 
C     Enable MSG communications checking to make sure all initialized
C -------------------------------------------------------------------

      CALL MPI_INITIALIZED(INIT_CHECK, MPI_IERR)
      IF( .NOT. INIT_CHECK ) THEN
         PRINT *,' MPI NOT INITIALIZED!'
         RETURN
      ENDIF
      
      IF( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR. 
     &    ISETUP.EQ.BMG_SETUP_only ) THEN

         CALL MSG_set_comm_parent(MPICOMM)
         CALL MSG_enable(BMG_MSG_MyProc,NProc)
         BMG_MSG_MyProc = MSG_myproc()
         NProc  = MSG_nproc()
            
      ELSE
         
         IF( MPICOMM .EQ. MSG_COMM_PARENT) THEN
            BMG_MSG_MyProc = MSG_myproc()
            NProc  = MSG_nproc()
         ELSE
            IF( MyProc.EQ.1 ) PRINT *,' MSG NOT INITIALIZED!'
         ENDIF

      ENDIF

C ---------------------------------
C     Initialize timers
C ---------------------------------

      T  = dZERO
      T1 = dZERO
      T2 = dZERO

C ----------------------------------
C     Sanity checks!
C ----------------------------------

      !!
      !!   We need to add way more some sanity checks here  !!!
      !!   Maybe a subroutine to check all parameter values??
      !!   

      !
      ! Processor indexing
      !
      IF ( BMG_MSG_MyProc.NE.BMG_MSG_iGRID(id_BMG_MSG_MyProc) ) THEN 
         WRITE(*,*) '**** ERROR: BMG3_SymStd_SOLVE_boxmg.f ... '
         WRITE(*,*) ' MyProc =         ', MyProc
         WRITE(*,*) ' BMG_MSG_MyProc = ', BMG_MSG_MyProc
         STOP
      ENDIF


C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> BEGIN:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      IF ( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR. 
     &     ISETUP.EQ.BMG_SETUP_only .OR. 
     &     ISETUP.EQ.BMG_SETUP_ptrs ) THEN

         TT1 = MPI_Wtime()

         CALL BMG3_SymStd_SETUP_PtrGrid( 
     &             NLx, NLy, NLz, NGx, NGy, NGz,
     &             iGs, jGs, kGs, 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,    
     &             MyProc, NProc, NProcI, NProcJ, NProcK, 
     &             MPICOMM, BMG_iWORK_PL, BMG_iWORK_CS, 
     &             NOGm, NFm, NCm, NSOm, NSORm, NCIm,
     &             NBMG_iWORK_PLm, NBMG_rWORK_PLm,
     &             NBMG_iWORK_CSm, NBMG_rWORK_CSm, NMSGim, NMSGrm,
     &             NOG, NF, NC, NSO, NSOR, NCI, IGRD,
     &             NBMG_iWORK_PL, NBMG_rWORK_PL,
     &             NBMG_iWORK_CS, NBMG_rWORK_CS, NMSGi, NMSGr
     &             )
         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF


         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_SETUP_PTR_GRID) =
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_PTR_GRID) + TT2 - TT1

         ! Output      
         IF( BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE) ) THEN
            WRITE (*,260) 
     &            'Storage for a vector on all grids = ', NF
         ENDIF

         ! Output
         IF( BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE) ) THEN
            WRITE (*,270) 
     &            'Storage for a vector on all coarse grids = ', NC
         ENDIF

         !
         ! Setup of the MSG arrays (pSI_MSG initialized above).
         !
         TT1 = MPI_Wtime()


         CALL BMG3_SymStd_SETUP_MSG( 
     &             pMSG, pMSGSO, iWorkMSG, NMSGi, pSI_MSG, 
     &             IGRD, NOG, NOGm, NProc, MyProc, 
     &             iWorkMSG(p_NLx_kg), iWorkMSG(p_NLy_kg),
     &             iWorkMSG(p_NLz_kg), 
     &             BMG_MSG_iGRID(p_NLx), BMG_MSG_iGRID(p_NLy),
     &             BMG_MSG_iGRID(p_NLz),
     &             BMG_MSG_iGRID(p_ProcGrid),
     &             NProcI, NProcJ, NProcK,
     &             MPICOMM 
     &             )

         TT2 = MPI_Wtime()

         BMG_rPARMS(id_BMG3_TIME_SETUP_MSG) = 
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_MSG) + TT2 - TT1
     
      ELSE

         !
         !  Maximum dimensions are the actual dimensions.
         !

         NF   = NFm
         NC   = NCm
         NSO  = NSOm
         NSOR = NSORm
         NCI  = NCIm

         NBMG_iWORK_CS = NBMG_iWORK_CSm
         NBMG_rWORK_CS = NBMG_rWORK_CSm
 
         NMSGi = NMSGim
         NMSGr = NMSGrm

      ENDIF

      !
      ! Create local pointers for LU solve
      !
      IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN

         p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
         p_BBD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CU)
         
         NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
         NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)

      ENDIF

      !
      !  No components just setup pointers, thanks!
      !
      IF ( ISETUP.EQ.BMG_SETUP_ptrs ) RETURN

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>> END:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>> BEGIN:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------


 89   CONTINUE  ! >>>>>>>>> BEGIN: Component Setup

      IF ( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR. 
     &     ISETUP.EQ.BMG_SETUP_opers .OR. 
     &     ISETUP.EQ.BMG_SETUP_only ) THEN

         !  Start the timer
         t1 = MPI_Wtime()

         !
         ! Construct coarse-grid and interpolation operators.  In 
         ! addition, setup any necessary components for relaxation.
         !
         TT1 = MPI_Wtime()
         
         CALL BMG3_SymStd_SETUP_parts( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm, NOG,
     &             BMG_iWORK_PL, NBMG_iWORK_PL,
     &             BMG_rWORK_PL, NBMG_rWORK_PL,
     &             BMG_iWORK_CS, NBMG_iWORK_CS,
     &             BMG_rWORK_CS, NBMG_rWORK_CS,
     &             iWorkMSG, NMSGi, pSI_MSG, pMSG, pMSGSO,
     &             MSG_Buffer, NMSGr,
     &             MyProcI, MyProcJ, MyProcK,
     &             NProcI, NProcJ, NProcK, NProc,
     &             BMG_MSG_iGRID(p_ProcGrid),
     &             BMG_MSG_iGRID(p_ProcCoord),
     &             iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &             iWorkMSG(p_NLx_kg), iWorkMSG(p_NLy_kg),
     &             iWorkMSG(p_NLz_kg),
     &             MPICOMM
     &             )
         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_SETUP_PARTS) =
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_PARTS) + TT2 - TT1

         ! Compute the setup time
         t2 = MPI_Wtime()
         t=t+t2-t1
         BMG_rPARMS(id_BMG3_TIME_SETUP_TOTAL) =
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_TOTAL) + t2 - t1

         
         ! Output the setup time
         IF( BMG_IOFLAG(iBMG3_OUT_TIME_SETUP) )  THEN
            WRITE (*,220) '(3D) SETUP TIME =', T
         ENDIF

      ENDIF

      !
      !  No solve just the setup, thanks!  
      !
      IF (ISETUP.EQ.BMG_SETUP_only) RETURN

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN: MULTIGRID CYCLING <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

 95   CONTINUE

      ! Start the timer
      t1 = MPI_Wtime()

C ----------------------------------------------------------------
C     Initial residual:
C ----------------------------------------------------------------

      !
      ! Get pointers for the finest grid (GRID LEVEL = NOG)
      !
      CALL BMG3_SymStd_GET_pointers( 
     &          NOG, IGRD, NOGm, 
     &          p_U, p_SO, p_SOR, p_CI, 
     &          IIL, JJL, KKL, IIG, JJG, KKG
     &          )

C      FIXME: The arguments passed in should be correct, but this 
C             is a reasonable place to double check
C
C      iGs = IGRD(NOG,idL_BMG_Icoord)
C      jGs = IGRD(NOG,idL_BMG_Jcoord)
C      kGs = IGRD(NOG,idL_BMG_Kcoord)

      !
      ! Set NStncl for the finest grid
      !
      IF ( IFD.NE.1 ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF
     
      !
      ! if requested, sync the initial guess
      !
      IF ( BMG_iPARMS(id_BMG3_SYNC_INITIAL_GUESS) 
     &     .eq. BMG_SYNC_INITIAL_GUESS ) THEN
         
         CALL BMG3_SymStd_UTILS_update_ghosts(
     &             NOG, Q(p_U), IIL, JJL, KKL, iWorkMSG, pMSG,
     &             MSG_Buffer, MPICOMM
     &             )

      END IF

      CALL BMG3_SymStd_residual(
     &          NOG, NOGm, IFD,
     &          Q(p_U), QF(p_U), SO(p_SO), RES(p_U), 
     &          IIL, JJL, KKL, NStncl,
     &          iWorkMSG, NMSGi, pMSG, 
     &          MSG_Buffer, NMSGr, MPICOMM
     &          )

      CALL BMG3_SymStd_UTILS_norm_l2( 
     &           RES(p_U), IIL, JJL, KKL, RES_L2_0, MPICOMM 
     &           )

      !
      !  Output
      !
      IF ( BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
         WRITE (*,405) '*** THE INITIAL RESIDUAL (l2-NORM) = ', RES_L2_0
      ENDIF

      IF ( RES_L2_0 .EQ. rZERO ) THEN
         IF ( BMG_IOFLAG(iBMG3_WARN_ZERO_RESIDUAL) ) THEN
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
         IF ( MyProc.EQ.1 ) THEN 
            WRITE(*,*) '*** Performing direct solve on the fine grid!'
         ENDIF
         !

         IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG3_SymStd_SOLVE_cg_LU( 
     &                q(p_U), qf(p_U), IIL, JJL, KKL,
     &                BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD), 
     &                NCBW, NCU, NOGm,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ELSE
            !
            TT1 = MPI_Wtime()
            !
            
            CALL BMG3_SymStd_SOLVE_cg_boxmg( 
     &                Q(p_U), QF(p_U), NOG,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                IIL, JJL, KKL, IIG, JJG, KKG,
     &                iGs, jGs, kGs,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ENDIF

         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)  
     &       = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1

         !
         ! Compute the final residual
         !
         CALL BMG3_SymStd_residual( 
     &             NOG, NOG, IFD,
     &             Q(p_U), QF(p_U), SO(p_SO), RES(p_U),
     &             IIL, JJL, KKL, NStncl, 
     &             iWorkMSG, NMSGi, pMSG, 
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )
         !
         CALL BMG3_SymStd_UTILS_norm_l2( 
     &             RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM 
     &             )
         !
         IF( BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) .AND. MyProc.EQ.1 ) THEN 
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

            CALL BMG3_SymStd_GET_pointers( k+1, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, 
     &                IILF, JJLF, KKLF, IIGF, JJGF, KKGF 
     &                )
            iGs_f = IGRD(k+1,idL_BMG_Icoord)
            jGs_f = IGRD(k+1,idL_BMG_Jcoord)
            kGs_f = IGRD(k+1,idL_BMG_Kcoord)

            CALL BMG3_SymStd_GET_pointers( k, IGRD, NOGm,
     &                p_UC, p_SOC, p_SORC, p_CIC, 
     &                IILC, JJLC, KKLC, IIGC, JJGC, KKGC  
     &                )
            iGs_c = IGRD(k,idL_BMG_Icoord)
            jGs_c = IGRD(k,idL_BMG_Jcoord)
            kGs_c = IGRD(k,idL_BMG_Kcoord)

            ! restrict the right hand side 

            TT2 = MPI_Wtime()

            CALL BMG3_SymStd_restrict( 
     &                K+1, K,  
     &                QF(p_U), QF(p_UC), CI(p_CIC),
     &                IILF, JJLF, KKLF, IILC, JJLC, KKLC,
     &                iGs_f, jGs_f, kGs_f
     &                )

            TT1 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_restrict) =
     &           BMG_rPARMS(id_BMG3_TIME_restrict) + TT2 - TT1

            ! Copy right hand side into RES(i,j,k) for interpolation
            DO i=0, IILF*JJLF*KKLF-1
               RES(p_U+i)=QF(p_U+i)
            END DO
         END DO


         ! FMG current finest grid index is set to the coarsest grid
         KFMG=KC     
         ! FMG: solve exactly on the coarsest grid
         CALL BMG3_SymStd_GET_pointers(
     &             KFMG, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI, 
     &             IIL, JJL, KKL, IIG, JJG, KKG 
     &             )
         iGs_c = IGRD(KFMG,idL_BMG_Icoord)
         jGs_c = IGRD(KFMG,idL_BMG_Jcoord)
         kGs_c = IGRD(KFMG,idL_BMG_Kcoord)
         
         IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU ) THEN
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG3_SymStd_SOLVE_cg_LU( 
     &                q(p_U), qf(p_U), IIL, JJL, KKL, 
     &                BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD), 
     &                NCBW, NCU, NOGm,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS 
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ELSE
            !
            TT1 = MPI_Wtime()
            !

            CALL BMG3_SymStd_SOLVE_cg_boxmg( 
     &                Q(p_U), QF(p_U), NOG,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                IIL, JJL, KKL, IIG, JJG, KKG,
     &                iGs_c, jGs_c, kGs_c,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            !
            TT2 = MPI_Wtime()
            !
         ENDIF

         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)  
     &       = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1
         
         RES_L2=dZERO
         IF ( BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
            WRITE (*,230) KFMG, RES_L2
         ENDIF

 100     CONTINUE    ! >>>>>>>> LOOP BOUNDARY: multigrid F-cycle


            KFMG=KFMG+1
            ! 
            ! FMG: Interpolate
            !
            CALL BMG3_SymStd_GET_pointers(
     &                KFMG, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, IIL, JJL, KKL,
     &                IIG, JJG, KKG 
     &                )
            iGs_f = IGRD(KFMG,idL_BMG_Icoord)
            jGs_f = IGRD(KFMG,idL_BMG_Jcoord)
            kGs_f = IGRD(KFMG,idL_BMG_Kcoord)

            CALL BMG3_SymStd_GET_pointers(
     &                KFMG-1, IGRD, NOGm,
     &                p_UC, p_SOC, p_SORC, p_CIC, IILC, JJLC, KKLC,
     &                IIGC, JJGC, KKGC  
     &                )
            iGs_c = IGRD(KFMG-1,idL_BMG_Icoord)
            jGs_c = IGRD(KFMG-1,idL_BMG_Jcoord)
            kGs_c = IGRD(KFMG-1,idL_BMG_Kcoord)

            !
            ! Set NStncl for grid KFMG
            !
            IF ( KFMG.LT.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
               NStncl=14
            ELSE
               NStncl=4
            ENDIF
            !


            TT1 = MPI_Wtime()

            CALL BMG3_SymStd_interp_add( 
     &                KFMG-1, KFMG, NOG,
     &                Q(p_U) ,Q(p_UC), RES(p_U),
     &                SO(p_SO), NStncl, CI(p_CIC),
     &                IILC, JJLC, KKLC,
     &                IIL, JJL, KKL,
     &                iGs_f, jGs_f, kGs_f,
     &                iWorkMSG, NMSGi, pMSG,
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )


            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_interp_add) =
     &           BMG_rPARMS(id_BMG3_TIME_interp_add) + TT2 - TT1
            !
            ! FMG; Perform an n-cycle 
            !
            IF ( BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN 
               WRITE(*,410) '*** Performing n-cycle:  Coarsest =', KC,
     &                      'Finest =', KFMG
            ENDIF


            CALL BMG3_SymStd_ncycle( 
     &                KC, KFMG, KF, 
     &                IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm,
     &                SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &                )


            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            IF ( KFMG.EQ.KF ) THEN
               !
               !  Compute the final residual
               !
               CALL BMG3_SymStd_residual(
     &                   KF, KF, IFD,
     &                   Q(p_U), QF(p_U), SO(p_SO), RES(p_U),
     &                   IIL, JJL, KKL, NStncl,
     &                   iWorkMSG, NMSGi, pMSG, 
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )
               CALL BMG3_SymStd_UTILS_norm_l2( 
     &                   RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM
     &                   )

               ! Note this is still vulnerable to very small RES_L2_0
               IF ( RES_L2_0.EQ.dZERO ) THEN
                  IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR).AND.MyProc.EQ.1)
     &                 THEN
                    WRITE(*,505)
     &                    '*** FATAL ERROR: BMG3_SymStd_SOLVE_boxmg.f '
                    WRITE(*,505)'    Initial residual has become zero!'
                  END IF
                  
                  CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,1)
                  RETURN
                  
               ELSE
                  REL_RES_L2 = RES_L2/RES_L2_0
               ENDIF             

               !
               !  Output 
               !
               IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &            .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
                  WRITE (*,400) '*** ITERATION =', NCYC, 
     &                          '*** RELATIVE RESIDUAL = ', REL_RES_L2
               ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 ) 
     &                 .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
                  WRITE (*,400) '*** ITERATION =', NCYC, 
     &                          '*** ABSOLUTE RESIDUAL = ', RES_L2
               ENDIF
               !
               !  Check Convergence
               !
               IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 ) 
     &            .AND. ( REL_RES_L2.LT.TOL )     ) THEN
                  ! FMG: converged
                  !      - set TOL and skip n-cycles
                  TOL = REL_RES_L2
                  NCYCLE_FLAG=.FALSE.
               ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 ) 
     &                 .AND. ( RES_L2.LT.TOL )         ) THEN
                  ! FMG: converged
                  !      - set TOL and skip n-cycles
                  TOL = RES_L2
                  NCYCLE_FLAG=.FALSE.
               ELSE IF ( NCYC.GE.MCYC ) THEN
                  ! FMG: cycle limit reached 
                  !      - set TOL and skip n-cycles
                  IF ( ISTOP.EQ.BMG_STOP_REL_RES_L2 ) THEN
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

      ENDIF  ! ISTRT > 0



C ----------------------------------------------------------------
C     Perform multigrid n-cycle(s):
C ----------------------------------------------------------------

      IF ( NCYCLE_FLAG ) THEN
         
         CALL BMG3_SymStd_GET_pointers(
     &             KF, IGRD, NOGm, p_U, p_SO, 
     &             p_SOR, p_CI, IIL, JJL, KKL,
     &             IIG, JJG, KKG
     &             )
         

 130     CONTINUE   ! >>>>>>>> LOOP BOUNDARY: multigrid n-cycles



            !
            ! Perform a multigrid n-cycle
            !
            CALL BMG3_SymStd_ncycle( 
     &                KC, KF, KF, 
     &                IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm,
     &                SO, NSOm, SOR, NSORm, CI, NCIm, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &                )

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
            !  Compute the final residual
            !
            CALL BMG3_SymStd_residual(
     &                KF, KF, IFD,
     &                Q(p_U), QF(p_U), SO(p_SO), RES(p_U), 
     &                IIL, JJL, KKL, NStncl, 
     &                iWorkMSG, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            CALL BMG3_SymStd_UTILS_norm_l2( 
     &                RES(p_U), IIL, JJL, KKL, RES_L2, MPICOMM
     &                )

            ! Note this is still vulnerable to very small RES_L2_0
            IF ( RES_L2_0.EQ.rZERO ) THEN
               WRITE(*,505) '*** FATAL ERROR: bmg3d.f '
               WRITE(*,505) '    Initial residual has become zero!!!'
               STOP
            ELSE
               REL_RES_L2 = RES_L2/RES_L2_0
            ENDIF

            !
            !  Output 
            !
            IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &           .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** ITERATION =', NCYC, 
     &                       '*** RELATIVE RESIDUAL = ', REL_RES_L2
            ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 )
     &              .AND. BMG_IOFLAG(iBMG3_OUT_ITERATIONS) ) THEN
               WRITE (*,400) '*** ITERATION =', NCYC, 
     &                      '*** ABSOLUTE RESIDUAL = ', RES_L2
            ENDIF

            !
            ! Check Convergence 
            !
            IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &         .AND. ( REL_RES_L2.LT.TOL )     ) THEN
               ! 
               ! n-cycles have converged in the RELATIVE RESIDUAL
               ! - set TOL and return
               !
               TOL = REL_RES_L2
            ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 ) 
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
               IF ( ISTOP.EQ.BMG_STOP_REL_RES_L2 ) THEN
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
      BMG_rPARMS(id_BMG3_STOP_TOL) = TOL
      !
      !  Compute the solve time and update the total time
      !
      t2 = MPI_Wtime()
      T=T+T2-T1
      BMG_rPARMS(id_BMG3_TIME_SOLVE_total) = 
     &     BMG_rPARMS(id_BMG3_TIME_SOLVE_total) + t2 - t1


      ! Output the multigrid cycling time
      IF ( BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING) ) THEN
         WRITE(*,240) '(3D) MULTIGRID CYCLING TIME =', T2-T1

      ENDIF

      IF ( BMG_IOFLAG(iBMG3_OUT_TIME_TOTAL) ) THEN
         WRITE(*,240) '(3D) TOTAL TIME = ', T
      ENDIF

C -------------------------------------
C     Output:
C -------------------------------------         

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
      
C ===========================================

      RETURN
      END
