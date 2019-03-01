      SUBROUTINE BMG2_SymStd_SOLVE_boxmg( 
     &                NLx, NLy, NGx, NGy, iGs, jGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NFm, NCm, 
     &                SO, NSOm, SOR, NSORm, CI, NCIm, 
     &                IGRD, NOGm, NOG, 
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &                MSG_Buffer, NMSGrm, MPICOMM
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
C     BMG2_SymStd_SOLVE_boxmg takes a matrix defined by the user from a
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
C     indeces k* are defined in the include file 'BMG_stencils.h'.  In
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
C                  file "BMG_parameters.h" for hints (the names are verbose.)
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
C                    IRELAX_SYM = BMG_RELAX_NONSYM
C                    => use the stated relaxation ordering (nonsymmetric)
C
C                    IRELAX_SYM = BMG_RELAX_SYM
C                    => make the cycle symmetric.
C
C                 as defined in BMG_parameters.h
C
C     IVW      Indicator for type of cycle. IVW=1 means v-cycles will
C              be performed. IVW=2 means w-cycles will be performed, etc.
C
C     MCYCL    Currently disabled!
C
C     NFm      Maximum storage for a vector on all grids. This should be
C              computed by BMG2_SymStd_SETUP_space and is checked in boxmg.
C
C     NCm      Maximum storage for a vector on all coarse grids.  This should
C              be computed in BMG2_SymStd_SETUP_space and is checked in boxmg.
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
C              this with a call to BMG2_SymStd_SETUP_space.  Clearly,
C            
C              IFD .EQ. 1  =>    NSOm = 3*NFm + 2*NCm 
C              IFD .NE. 1  =>    NSOm = 5*NFm 
C  
C              BMG2_SymStd_SOLVE_boxmg verifies that NSOm is large enough.
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
C            BMG_iPARMS(id_BMG2_Ext_Err_Code)          SETUP_cg_LU
C     18     NStncl .ne. 3,5                           SETUP_cg_LU
C     19     KF-1 .ne. KG in SETUP_interp_OI           SETUP_interp_OI
C     20     Coarse grid solve failed,  look for the
C            return code of the LAPACK routine in
C            BMG_iPARMS(id_BMG2_Ext_Err_Code)          SOLVE_cg
C
C ==========================================================================

      IMPLICIT NONE
     
C -----------------------------
C     Includes
C
      INCLUDE 'geom_param_fort.h'
      INCLUDE 'mpi_param_fort.h'
      INCLUDE 'MSG.h'

#include     "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  NLx, NLy, NGx, NGy, iGs, jGs,
     &         NFm, NCm, NSOm, NSORm, NCIm,
     &         NOGm, NOG,
     &         NMSGim, NMSGrm, MPICOMM,
     &         NBMG_iWORK_CSm, NBMG_rWORK_CSm

      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)


      INTEGER  iWorkMSG(NMSGim),
     &         pMSG(NBMG_pMSG,NOGm),
     &         pMSGSO(NBMG_pMSG,NOGm),
     &         pLS(NBMG_pLS,NOGm)

      INTEGER  IGRD(NOGm,NBMG_pIGRD)
      REAL*RKIND   CI(NCIm), SO(NSOm), SOR(NSORm), Q(NFm), QF(NFm),
     &         RES(NFm), MSG_Buffer(NMSGrm)

      INTEGER  BMG_iWORK_CS(NBMG_iWORK_CSm)
      REAL*RKIND   BMG_rWORK_CS(NBMG_rWORK_CSm)

C ----------------------------
C     Local Declarations
C
      INTEGER  I, IBC, ID, IFD, IIL, IILC, IILF, IM,
     &         IRELAX, IRELAX_SYM, ISETUP, ISTOP, ISTRT, IU, IVW,
     &         JJL, JJLC, JJLF, K, KC, KF, KFMG, MCYC, MCYCL,
     &         NC, NCBW, NCI, NCYC, NCU, NF, NSO, NSOR,
     &         NStncl, NXYc, p_ABD, p_BBD,
     &         p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC,
     &         IIG, IIGC, IIGF, JJG, JJGC, JJGF,
     &         NBMG_iWORK_CS, NBMG_rWORK_CS 

      INTEGER  MyProc, MyProcI, MyProcJ,
     &         NProc, NProcI, NProcJ,
     &         NMSGi, NMSGr, 
     &         p_NLx, p_NLy, p_ProcCoord, p_ProcGrid, pSI_MSG,
     &         p_NLx_kg, p_NLy_kg

      INTEGER  XLINECOMM, YLINECOMM

      INTEGER  BMG_MSG_MyProc, MPI_IERR

      REAL*RKIND   RES_L2, RES_L2_0, REL_RES_L2
      REAL*8   T, T1, T2, TOL
      REAL*8   TT1, TT2
      LOGICAL  NCYCLE_FLAG, INIT_CHECK

C =========================================================================


      IF ( BMG_iPARMS(id_BMG2_Err_Code).NE.iZERO ) THEN
         RETURN
      END IF


C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------
C --------------------------------------------
C     MPI Parallel Data Decomposition:
C --------------------------------------------
      
      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)

      p_NLx       = BMG_MSG_pGRID(ip_BMG_MSG_NLx_Grid)
      p_NLy       = BMG_MSG_pGRID(ip_BMG_MSG_NLy_Grid)
      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)
      
      !
      !  Line Communicators:
      !
      XLINECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_x)
      YLINECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_y)
      
C ----------------------------------------
C     MSG workspace pointers:
C ----------------------------------------

      p_NLx_kg = 1
      p_NLy_kg = p_NLx_kg + NProcI*NOGm

      pSI_MSG  = p_NLy_kg + NProcJ*NOGm  !!NB: initalization of pSI_MSG

C ----------------------------------------
C     Multigrid Cycle parameters:
C ----------------------------------------

      IFD    = BMG_iPARMS(id_BMG2_STENCIL)
      IBC    = BMG_iPARMS(id_BMG2_BC) 
      ISETUP = BMG_iPARMS(id_BMG2_SETUP)
      NXYc   = BMG_iPARMS(id_BMG2_CG_MIN_DIM)

      IRELAX = BMG_iPARMS(id_BMG2_RELAX)
      IRELAX_SYM = BMG_iPARMS(id_BMG2_RELAX_SYM )

      ID = BMG_iPARMS(id_BMG2_NRELAX_DOWN)
      IU = BMG_iPARMS(id_BMG2_NRELAX_UP)
      IM = BMG_iPARMS(id_BMG2_NRELAX_FG)

      IF (BMG_iPARMS(id_BMG2_CYCLE_CLASS).EQ.BMG_N_CYCLE) THEN
         ISTRT = - BMG_iPARMS(id_BMG2_MAX_ITERS)
      ELSEIF (BMG_iPARMS(id_BMG2_CYCLE_CLASS).EQ.BMG_FMG_CYCLE) THEN 
         ISTRT = BMG_iPARMS(id_BMG2_MAX_ITERS)
      ENDIF

      IVW   = BMG_iPARMS(id_BMG2_NCYCLE_TYPE)
      MCYCL = BMG_iPARMS(id_BMG2_FMG_NNCYCLE)

      ISTOP = BMG_iPARMS(id_BMG2_STOP_TEST)
      TOL   = BMG_rPARMS(id_BMG2_STOP_TOL)

C --------------------------------------------
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

      T  = rZERO
      T1 = rZERO
      T2 = rZERO

C ----------------------------------
C     Sanity checks!
C ----------------------------------

      !!
      !!   FIXME: 
      !!
      !!   We need to add way more some sanity checks here  !!!
      !!   Maybe a subroutine to check all parameter values??
      !!   

      !
      ! Processor indexing
      !
      IF ( BMG_MSG_MyProc.NE.BMG_MSG_iGRID(id_BMG_MSG_MyProc) ) THEN 
         WRITE(*,*) '**** ERROR: BMG2_SymStd_SOLVE_boxmg.f ... '
         WRITE(*,*) ' MyProc =         ', MyProc
         WRITE(*,*) ' BMG_MSG_MyProc = ', BMG_MSG_MyProc
         !
         ! FIXME: probably hopeless at this point but we should 
         !        collect the outcome and return an error code.
         !
         STOP
      ENDIF

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> BEGIN:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      IF ( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR. 
     &     ISETUP.EQ.BMG_SETUP_ptrs .OR. 
     &     ISETUP.EQ.BMG_SETUP_only ) THEN

         TT1 = MPI_Wtime()

         CALL BMG2_SymStd_SETUP_PtrGrid( 
     &             NLx, NLy, NGx, NGy, iGs, jGs,
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             MyProc, Nproc, NprocI, NProcJ, 
     &             MPICOMM, BMG_iWORK_CS,
     &             NOGm, NFm, NCm, NSOm, NSORm, NCIm,
     &             NBMG_iWORK_CSm, NBMG_rWORK_CSm, 
     &             NMSGim, NMSGrm,
     &             NOG, NF, NC, NSO, NSOR, NCI, IGRD,
     &             NBMG_iWORK_CS, NBMG_rWORK_CS,
     &             NMSGi, NMSGr
     &             )

         
         IF ( BMG_iPARMS(id_BMG2_Err_Code).NE.iZERO) THEN
            RETURN
         ENDIF

          TT2 = MPI_Wtime()
          BMG_rPARMS(id_BMG2_TIME_SETUP_PTR_GRID) =
     &        BMG_rPARMS(id_BMG2_TIME_SETUP_PTR_GRID) + TT2 - TT1

        ! Output
         IF (MyProc.eq.1) THEN
            IF( BMG_IOFLAG(iBMG2_OUT_WSPACE_SIZE) ) THEN
               WRITE (*,260) 
     &              'Storage for a vector on all grids = ', NF
            ENDIF
            
            ! Output
            IF( BMG_IOFLAG(iBMG2_OUT_WSPACE_SIZE) ) THEN
               WRITE (*,270) 
     &              'Storage for a vector on all coarse grids = ', NC
            ENDIF
         ENDIF

         !
         ! Setup of the MSG arrays
         !
         TT1 = MPI_Wtime()

         CALL BMG2_SymStd_SETUP_MSG(
     &             pMSG, pMSGSO, iWorkMSG, NMSGi, pSI_MSG,
     &             IGRD, NOG, NOGm, NProc, MyProc,
     &             iWorkMSG(p_NLx_kg), iWorkMSG(p_NLy_kg),
     &             BMG_MSG_iGRID(p_NLx), BMG_MSG_iGRID(p_NLy),
     &             BMG_MSG_iGRID(p_ProcGrid),
     &             NProcI, NProcJ, 
     &             MPICOMM 
     &             )
         
         TT2 = MPI_Wtime()

         BMG_rPARMS(id_BMG2_TIME_SETUP_MSG) = 
     &        BMG_rPARMS(id_BMG2_TIME_SETUP_MSG) + TT2 - TT1

         !
         ! Setup the data distribution for the line solves
         !
         TT1 = MPI_Wtime()

         CALL BMG2_SymStd_SETUP_LS( 
     &                    iWorkMSG, NMSGi, pMSG, pLS, pSI_MSG,
     &                    BMG_MSG_iGRID(p_ProcGrid),
     &                    NProcI, NProcJ, NOG
     &                    )         
         
         TT2 = MPI_Wtime()
        
         BMG_rPARMS(id_BMG2_TIME_SETUP_LS) = 
     &        BMG_rPARMS(id_BMG2_TIME_SETUP_LS) + TT2 - TT1
        
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
      IF ( BMG_iPARMS(id_BMG2_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN

         p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
         p_BBD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CU)
         
         NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
         NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)

      ENDIF

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>> END:  POINTER SETUP <<<<<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------


C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>> BEGIN:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      IF ( ISETUP.EQ.BMG_SETUP_ptrs_opers .OR. 
     &     ISETUP.EQ.BMG_SETUP_opers      .OR. 
     &     ISETUP.EQ.BMG_SETUP_only ) THEN 

         ! Start timer
         T1 = MPI_Wtime()

         !
         ! Construct coarse-grid and interpolation operators.  In 
         ! addition, setup any necessary components for relaxation.
         !
         TT1 = MPI_Wtime()
         
         CALL BMG2_SymStd_SETUP_parts( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             SO, NSO, SOR, NSOR, CI, NCI, IGRD, NOGm, NOG,
     &             BMG_iWORK_CS, NBMG_iWORK_CS,
     &             BMG_rWORK_CS, NBMG_rWORK_CS,
     &             iWorkMSG, NMSGi, pMSG, pMSGSO, 
     &             MSG_Buffer, NMSGr,
     &             MyProcI, MyProcJ, NProcI, NProcJ, NProc, 
     &             BMG_MSG_iGRID(p_ProcGrid), 
     &             BMG_MSG_iGRID(p_ProcCoord),
     &             iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &             MPICOMM, XLINECOMM, YLINECOMM
     &             )

         !
         ! Save communicators for line relaxation
         !
         BMG_MSG_iGRID(id_BMG_MSG_COMM_x) = XLINECOMM
         BMG_MSG_iGRID(id_BMG_MSG_COMM_y) = YLINECOMM

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG2_TIME_SETUP_PARTS) =
     &        BMG_rPARMS(id_BMG2_TIME_SETUP_PARTS) + TT2 - TT1

         ! Compute the setup time
         t2 = MPI_Wtime()
         t=t+t2-t1
         BMG_rPARMS(id_BMG2_TIME_SETUP_TOTAL) =
     &        BMG_rPARMS(id_BMG2_TIME_SETUP_TOTAL) + t2 - t1


         IF (MyProc.eq.1) THEN
            ! Output the setup time
            IF( BMG_IOFLAG(iBMG2_OUT_TIME_SETUP) )  THEN
               WRITE (*,220) '(2D) SETUP TIME =', T
            ENDIF
         ENDIF

      ENDIF

      !
      !  No solve just the setup, thanks!  
      !
      IF ( ISETUP.EQ.BMG_SETUP_only ) RETURN


C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  COMPONENT SETUP <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------


C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN: MULTIGRID CYCLING <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      ! Start the timer
      T1 = MPI_Wtime()

C --------------------------------
C     Initial residual:
C -------------------------------- 

      !
      ! Get pointers for the finest grid (GRID LEVEL = NOG)
      !
      CALL BMG2_SymStd_GET_pointers( 
     &          NOG, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CI,
     &          IIL, JJL, IIG, JJG 
     &          )

      iGs = IGRD(NOG,idL_BMG_Icoord)
      jGs = IGRD(NOG,idL_BMG_Jcoord)

      !
      ! Set NStncl for the finest grid
      !
      IF ( IFD.NE.1 ) THEN
         NStncl=5
      ELSE
         NStncl=3
      ENDIF

      !
      ! if requested, sync the initial guess
      !
      IF ( BMG_iPARMS(id_BMG2_SYNC_INITIAL_GUESS) 
     &     .eq. BMG_SYNC_INITIAL_GUESS ) THEN
         !
         CALL BMG2_SymStd_UTILS_update_ghosts(
     &             NOG, Q(p_U), IIL, JJL, iWorkMSG, pMSG,
     &             MSG_Buffer, MPICOMM
     &             )
         !
      ENDIF

      !
      ! Compute the initial residual
      !
      CALL BMG2_SymStd_residual( 
     &          NOG, NOG, NOG, RES_L2, 
     &          SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &          IIL, JJL, IFD, NStncl, IRELAX, IRELAX_SYM,
     &          iWorkMSG, NMSGi, pMSG, 
     &          MSG_Buffer, NMSGr, MPICOMM
     &          )


      CALL BMG2_SymStd_UTILS_norm_l2( 
     &          RES(p_U), IIL, JJL, RES_L2_0, MPICOMM )

      !
      !  Output
      !
      IF (MyProc.eq.1) THEN
         IF ( BMG_IOFLAG(iBMG2_OUT_ITERATIONS) ) THEN
            WRITE (*,405) '*** (MPI) THE INITIAL RESIDUAL (l2-NORM) = ', 
     &           RES_L2_0
         ENDIF
      ENDIF

      !
      ! return if the initial residual is zero
      !
      IF (RES_L2_0 .EQ. rZERO) RETURN


C ----------------------------------------------------------------
C     Direct solve in the degenerate case of one grid!
C ----------------------------------------------------------------

      IF ( NOG.EQ.1 ) THEN
         !
         ! Solve on the coarsest grid
         !
         IF ( BMG_iPARMS(id_BMG2_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN
            !
            IF ( MyProc.EQ.1
     &          .AND. BMG_IOFLAG(iBMG2_BUG_RES_CG_SOLVE) ) THEN
               WRITE(*,*) '*** Direct solve on the fine grid!'
               WRITE(*,*) 
            ENDIF
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG2_SymStd_SOLVE_cg_LU(
     &                Q(p_U), QF(p_U), IIL, JJL,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD),
     &                NCBW, NCU, NOGm,
     &                NProcI, NProcJ, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                iWorkMSG, NMSGi, pMSG,
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            !
            TT2 = MPI_Wtime()
         ELSE
            !
            IF ( MyProc.EQ.1
     &          .AND. BMG_IOFLAG(iBMG2_BUG_RES_CG_SOLVE) ) THEN
               WRITE(*,*) '*** Serial BoxMG on the fine grid!'
               WRITE(*,*) 
            ENDIF
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG2_SymStd_SOLVE_cg_boxmg( 
     &                Q(p_U), QF(p_U), NOG,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                IIL, JJL, IIG, JJG, iGs, jGs,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                NProcI, NProcJ, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            !
            TT2 = MPI_Wtime()
         ENDIF 

         BMG_rPARMS(id_BMG2_TIME_SOLVE_CG)  
     &           = BMG_rPARMS(id_BMG2_TIME_SOLVE_CG) + TT2 - TT1
            
         IF ( BMG_iPARMS(id_BMG2_Err_Code).NE.iZERO ) THEN
            RETURN
         ENDIF

         !
         ! Compute the final residual 
         !
c$$$
c$$$         CALL BMG2_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, Q(p_U), NOG, NOG,
c$$$     &             IIL, JJL, IIG, JJG, iGs, jGs,
c$$$     &             'final-solution-MPI', .TRUE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProc,
c$$$     &             MyProcI, MyProcJ, MPICOMM
c$$$     &             )
c$$$
         CALL BMG2_SymStd_residual(  
     &             NOG, NOG, NOG, RES_L2, 
     &             SO(p_SO), QF(p_U), Q(p_U), RES(p_U), 
     &             IIL, JJL, IFD, NStncl, IRELAX, IRELAX_SYM,
     &             iWorkMSG, NMSGi, pMSG, 
     &             MSG_Buffer, NMSGr, MPICOMM 
     &             )
c$$$         !
c$$$         CALL BMG2_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, RES(p_U), NOG, NOG,
c$$$     &             IIL, JJL, IIG, JJG, iGs, jGs,
c$$$     &             'final-residual-MPI', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProc,
c$$$     &             MyProcI, MyProcJ, MPICOMM
c$$$     &             )
         !
         CALL BMG2_SymStd_UTILS_norm_l2( 
     &             RES(p_U), IIL, JJL, RES_L2, MPICOMM 
     &             )
         !
         IF (RES_L2_0 .EQ. rZERO ) THEN
            !
            IF(MyProc.EQ.1) THEN
               WRITE(*,505) '*** FATAL ERROR: boxmg.f '
               WRITE(*,505) '    Initial residual has become zero!!!'
            END IF
            !
            CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,3)
            RETURN
            !
         ELSE
            REL_RES_L2 = RES_L2/RES_L2_0
         ENDIF

         IF ( ISTOP.EQ.BMG_STOP_REL_RES_L2 ) THEN
            !
            TOL=REL_RES_L2
            !
            IF ( MyProc.EQ.1 
     &          .AND. ( BMG_IOFLAG(iBMG2_BUG_RES_CG_SOLVE) 
     &                 .OR. BMG_IOFLAG(iBMG2_OUT_ITERATIONS) )
     &           ) THEN
               WRITE (*,405) 
     &              '*** (MPI) FINAL RELATIVE RESIDUAL (l2-NORM) = ',
     &              REL_RES_L2
               !
            ENDIF
         ELSE
            TOL=RES_L2
            IF ( MyProc.EQ.1 
     &           .AND. ( BMG_IOFLAG(iBMG2_BUG_RES_CG_SOLVE)
     &                  .OR. BMG_IOFLAG(iBMG2_OUT_ITERATIONS) )
     &           ) THEN 
               !
               WRITE (*,405) 
     &              '*** (MPI) FINAL ABSOLUTE RESIDUAL (l2-NORM) = ', 
     &              RES_L2
               !
            ENDIF
            !
         ENDIF
         !
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

C ----------------------------------------------------------------
C     Perform multigrid F-cycle:
C     (if necessary, then continue with n-cycles)
C ----------------------------------------------------------------

      IF ( ISTRT.GT.0 ) THEN
      
         ! Coarsen the right hand side
         DO K=KF-1, 1, -1

            CALL BMG2_SymStd_GET_pointers( 
     &                K, IGRD, NOGm,
     &                p_UC, p_SOC, p_SORC, p_CIC, 
     &                IILC, JJLC, IIGC, JJGC 
     &                )
            CALL BMG2_SymStd_GET_pointers( 
     &                K+1, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, 
     &                IILF, JJLF, IIGF, JJGF 
     &                )

            iGs = IGRD(K+1,idL_BMG_Icoord)
            jGs = IGRD(K+1,idL_BMG_Jcoord)

            ! restrict the right hand side 

            TT1 = MPI_Wtime()

            CALL BMG2_SymStd_restrict(
     &                K+1, K, NOGm, 
     &                QF(p_U), QF(p_UC), CI(p_CIC),
     &                IILF, JJLF, IILC, JJLC, iGs, jGs,
     &                iWorkMSG, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )

            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG2_TIME_restrict) =
     &           BMG_rPARMS(id_BMG2_TIME_restrict) + TT2 - TT1

            ! Copy right hand side into RES(i,j)
            DO i=0, IILF*JJLF-1
               RES(p_U+i)=QF(p_U+i)
            END DO

         END DO

         ! FMG current finest grid index is set to the coarsest grid
         KFMG=KC     
         ! FMG: solve exactly on the coarsest grid
         CALL BMG2_SymStd_GET_pointers( 
     &             KFMG, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI, 
     &             IIL, JJL, IIG, JJG 
     &             )

         IF ( BMG_iPARMS(id_BMG2_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG2_SymStd_SOLVE_cg_LU(
     &                Q(p_U), QF(p_U), IIL, JJL,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,            
     &                BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD),
     &                NCBW, NCU, NOGm,
     &                NProcI, NProcJ, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                iWorkMSG, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            !
            TT2 = MPI_Wtime()
         ELSE
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG2_SymStd_SOLVE_cg_boxmg( 
     &                Q(p_U), QF(p_U), NOG,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                IIL, JJL, IIG, JJG, iGs, jGs,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                NProcI, NProcJ, NProc, MyProc,
     &                BMG_MSG_iGRID(p_ProcGrid),
     &                BMG_MSG_iGRID(p_ProcCoord),
     &                iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            !
            TT2 = MPI_Wtime()
         ENDIF 

         BMG_rPARMS(id_BMG2_TIME_SOLVE_CG)  
     &       = BMG_rPARMS(id_BMG2_TIME_SOLVE_CG) + TT2 - TT1
 
         IF (BMG_iPARMS(id_BMG2_Err_Code) .ne. iZERO) THEN
            RETURN
         ENDIF

         RES_L2=rZERO
         IF (MyProc.eq.1) THEN
            IF ( BMG_IOFLAG(iBMG2_OUT_ITERATIONS) ) THEN
               WRITE (*,230) KFMG, RES_L2
            ENDIF
         ENDIF

 100        CONTINUE            ! >>>>>>>> LOOP BOUNDARY: multigrid F-cycle

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
            CALL BMG2_SymStd_GET_pointers(
     &                KFMG, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI,
     &                IIL, JJl, IIG, JJG 
     &                )
            CALL BMG2_SymStd_GET_pointers( 
     &                KFMG-1, IGRD, NOGm, 
     &                p_UC, p_SOC, p_SORC, p_CIC,
     &                IILC, JJLC, IIGC, JJGC 
     &                )
            iGs = IGRD(KFMG-1,idL_BMG_Icoord)
            jGs = IGRD(KFMG-1,idL_BMG_Jcoord)
            !
            TT1 = MPI_Wtime()
            !
            CALL BMG2_SymStd_interp_add(
     &                KFMG-1, KFMG, NOG, 
     &                Q(p_U), Q(p_UC), RES(p_U),
     &                SO(p_SO), NStncl, CI(p_CIC),
     &                IILC, JJLC, IIL, JJL, iGs, jGs,
     &                iWorkMSG, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            !
            TT2 = MPI_Wtime()
            !
             BMG_rPARMS(id_BMG2_TIME_interp_add) =
     &           BMG_rPARMS(id_BMG2_TIME_interp_add) + TT2 - TT1
            !
            ! FMG; Perform an n-cycle 
            !
            IF (MyProc.eq.1) THEN
               IF ( BMG_IOFLAG(iBMG2_OUT_ITERATIONS) ) THEN 
                  WRITE(*,410) '*** Performing n-cycle:  Coarsest =', 
     &                 KC, 'Finest =', KFMG
               ENDIF
            ENDIF

           CALL BMG2_SymStd_ncycle(  
     &               KC, KFMG, KF,
     &               IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, 
     &               BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &               Q, QF, RES, NFm, NCm,
     &               SO, NSOm, SOR, NSORm, CI, NCIm,
     &               IGRD, NOGm, NOG,
     &               BMG_iWORK_CS, NBMG_iWORK_CS,
     &               BMG_rWORK_CS, NBMG_rWORK_CS,
     &               iWorkMSG, NMSGi, pMSG, pLS, 
     &               MSG_Buffer, NMSGr,
     &               BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &               )

            IF (BMG_iPARMS(id_BMG2_Err_Code) .ne. iZERO) THEN
               RETURN
            ENDIF

            IF ( KFMG.EQ.KF ) THEN
               !
               ! Compute the final residual
               !
               CALL BMG2_SymStd_residual( 
     &                   KF, KF, NOG, RES_L2, 
     &                   SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &                   IIL, JJL, IFD, NStncl, IRELAX, IRELAX_SYM,
     &                   iWorkMSG, NMSGi, pMSG, 
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )
               CALL BMG2_SymStd_UTILS_norm_l2( RES(p_U), IIL, JJL, 
     &              RES_L2, MPICOMM )

               IF ( RES_L2_0.EQ.rZERO ) THEN
                  !
                  IF ( BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)
     &                .AND. MyProc.EQ.1
     &               ) THEN
                     !
                     WRITE(*,505) 
     &                    '*** WARNING: BMG2_SymStd_SOLVE_boxmg.f '
                     WRITE(*,505) '   Initial residual has become zero!'
                     !
                  ENDIF
                  !
                  CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,3)
                  RETURN
                  !
               ELSE
                  REL_RES_L2 = RES_L2/RES_L2_0
               ENDIF               

               !
               ! Output
               !
               IF (MyProc.eq.1) THEN
                  IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &                 .AND. BMG_IOFLAG(iBMG2_OUT_ITERATIONS) ) THEN
                     WRITE (*,400) '*** ITERATION =', NCYC, 
     &                    '*** RELATIVE RESIDUAL =', REL_RES_L2
                  ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 ) 
     &                    .AND. BMG_IOFLAG(iBMG2_OUT_ITERATIONS) ) THEN 
                     WRITE (*,400) '*** ITERATION =', NCYC, 
     &                    '*** ABSOLUTE RESIDUAL =', RES_L2
                  ENDIF
               ENDIF
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
                  ! - set TOL and skip n-cycles
                  TOL = RES_L2
                  NCYCLE_FLAG=.FALSE.
               ELSE IF ( NCYC.GE.MCYC ) THEN
                  ! FMG: cycle limit reached 
                  ! - set TOL and skip n-cycles
                  TOL = -RES_L2
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

         CALL BMG2_SymStd_GET_pointers( 
     &             KF, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI,
     &             IIL, JJL, IIG, JJG 
     &             )

 130     CONTINUE  ! >>>>>>>> LOOP BOUNDARY: multigrid n-cycles
           !
           ! Perform a multigrid n-cycle
           !
           CALL BMG2_SymStd_ncycle(  
     &               KC, KF, KF,
     &               IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, 
     &               BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &               Q, QF, RES, NFm, NCm,
     &               SO, NSOm, SOR, NSORm, CI, NCIm,
     &               IGRD, NOGm, NOG,
     &               BMG_iWORK_CS, NBMG_iWORK_CS,
     &               BMG_rWORK_CS, NBMG_rWORK_CS,
     &               iWorkMSG, NMSGi, pMSG, pLS, 
     &               MSG_Buffer, NMSGr,
     &               BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &               )
           !
           IF (BMG_iPARMS(id_BMG2_Err_Code) .ne. iZERO) THEN
              RETURN
           ENDIF
           !
           ! Compute the final residual
           !
           CALL BMG2_SymStd_residual( 
     &               KF, KF, NOG, RES_L2, 
     &               SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &               IIL, JJL, IFD, NStncl, IRELAX, IRELAX_SYM,
     &               iWorkMSG, NMSGi, pMSG,
     &               MSG_Buffer, NMSGr, MPICOMM
     &               )
           CALL BMG2_SymStd_UTILS_norm_l2( 
     &               RES(p_U), IIL, JJL, RES_L2, MPICOMM
     &               )
           !
           IF (RES_L2_0 .EQ. rZERO ) THEN
              !
              IF(MyProc.EQ.1) THEN
                 WRITE(*,505) '*** FATAL ERROR: boxmg.f '
                 WRITE(*,505) '    Initial residual has become zero!!!'
              END IF
              !
              CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,3)
              RETURN
              !
           ELSE
              REL_RES_L2 = RES_L2/RES_L2_0
           ENDIF

           !
           ! Output
           !
           IF (MyProc.eq.1) THEN
              IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &             .AND. BMG_IOFLAG(iBMG2_OUT_ITERATIONS) ) THEN
                 WRITE (*,400) '*** ITERATION =', NCYC, 
     &                '*** RELATIVE RESIDUAL =', REL_RES_L2
              ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 )
     &                .AND. BMG_IOFLAG(iBMG2_OUT_ITERATIONS) ) THEN 
                 WRITE (*,400) '*** ITERATION =', NCYC, 
     &                '*** ABSOLUTE RESIDUAL =', RES_L2
              ENDIF
           ENDIF
           !
           ! Check Convergence 
           !
           IF ( ( ISTOP.EQ.BMG_STOP_REL_RES_L2 )
     &        .AND. ( REL_RES_L2.LT.TOL )     ) THEN
              ! 
              ! n-cycles have converged in the RELATIVE RESIDUAL
              ! - set TOL and return
              !
              TOL = REL_RES_L2
           ELSE IF ( ( ISTOP.EQ.BMG_STOP_ABS_RES_L2 ) 
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
              TOL = -RES_L2
           ELSE 
              !
              ! increase cycle count and continue n-cycling
              !
              NCYC=NCYC+1
              GOTO 130
           ENDIF

 140     CONTINUE    ! <<<<<<<< LOOP BOUNDARY: multigrid n-cycles

      ENDIF

 150  CONTINUE

      !
      !
      !
      BMG_rPARMS(id_BMG2_STOP_TOL)  = TOL

      !
      !  Compute the solve time and update the total time
      !
      T2 = MPI_Wtime()
      T=T+T2-T1

      BMG_rPARMS(id_BMG2_TIME_SOLVE_total) = 
     &     BMG_rPARMS(id_BMG2_TIME_SOLVE_total) + t2 - t1

      IF (MyProc.eq.1) THEN
         ! Output the multigrid cycling time
         IF ( BMG_IOFLAG(iBMG2_OUT_TIME_CYCLING) ) THEN
            WRITE(*,240) '(2D) MULTIGRID CYCLING TIME =', T2-T1
         ENDIF
         
         IF ( BMG_IOFLAG(iBMG2_OUT_TIME_TOTAL) ) THEN
            WRITE(*,240) '(2D) TOTAL TIME = ', T
         ENDIF
      ENDIF

      IF ( ISETUP.NE.BMG_SETUP_only ) THEN

         !CALL MSG_disable(ierror)

      ENDIF

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END: MULTIGRID CYCLING <<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      
C -------------------------------------
C     Output
C -------------------------------------         
      
      IF ( BMG_IOFLAG(iBMG_OUT_SOLUTION) ) THEN
         !
         ! print the solution on the finest grid
         !
         CALL BMG2_SymStd_GET_pointers( 
     &             KF, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI,
     &             IILF, JJLF, IIGF, JJGF 
     &             )
         iGs = IGRD(KF,idL_BMG_Icoord)
         jGs = IGRD(KF,idL_BMG_Jcoord)

         CALL BMG2_SymStd_DUMP_vector( 
     &             BMG_IOFLAG, QF(p_U), KF, NOG,
     &             IILF, JJLF, IIGF, JJGF, iGs, jGs,
     &             'Q-solution', .FALSE.,
     &             BMG_MSG_iGRID(p_ProcGrid),
     &             NProcI, NProcJ, NProc,
     &             MyProcI, MyProcJ, MPICOMM
     &             )
         !
      ENDIF

      IF ( BMG_IOFLAG(iBMG_OUT_RHS) ) THEN  
         !
         ! print QF on the first coarse grid
         !
         CALL BMG2_SymStd_GET_pointers(
     &             KF-1, IGRD, NOGm,
     &             p_UC, p_SOC, p_SORC, p_CIC,
     &             IILC, JJLC, IIGC, JJGC 
     &             )
         iGs = IGRD(KF-1,idL_BMG_Icoord)
         jGs = IGRD(KF-1,idL_BMG_Jcoord)

         CALL BMG2_SymStd_DUMP_vector( 
     &             BMG_IOFLAG, QF(p_UC), KF-1, NOG,
     &             IILC, JJLC, IIGC, JJGC, iGs, jGs,
     &             'QF-source', .FALSE.,
     &             BMG_MSG_iGRID(p_ProcGrid),
     &             NProcI, NProcJ, NProc,
     &             MyProcI, MyProcJ, MPICOMM
     &             )
         !
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
