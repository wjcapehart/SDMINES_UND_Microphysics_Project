      SUBROUTINE BMG2_SymStd_SOLVE_pcg ( 
     &                NLx, NLy, NGx, NGy, iGs, jGs,
     &                BMG_PCG_iPARMS, BMG_PCG_rPARMS,
     &                BMG_PCG_RES, MAX_PCG_ITERS, PCG_ITERS,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, NFm, SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
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
C     This subroutine is a preconditioned conjugate gradient code which
C     may use multigrid preconditioning provided by the black box
C     multigrid code, BMG2_SymStd_SOLVE_boxmg.  It may be used for
C     discretizations of second order elliptic partial differential
C     equations that generate, at most, a 9-point stencil on a logically
C     rectangular grid.  It may be applied to other, similarly
C     structured, problems.  See the BMG2_SymStd_SOLVE_boxmg code for
C     additional information.
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
C     BMG_PCG_iPARMS      PCG integer parameters
C
C     * BMG_PCG_iPARMS(2) --  Maximum number of symmetric V-cycles 
C                             to perform as preconditioner.
C
C     * BMG_PCG_iPARMS(3) --  Type of stop test to use for the pcg 
C                             algorithm.  (See PCG_parameters.h)
C
C     BMG_PCG_rPARMS      PCG real parameters
C
C     * BMG_PCG_rPARMS(1) --  Tolerance used to stop PCG iterations.
C
C     BMG_PCG_NORMits     Array of size PCG_iPARMS(1) containing l2
C                         norm of residual after each pcg iteration.
C
C     BMG_iPARMS      See BMG2_SymStd_SOLVE_boxmg code
C     
C     BMG_rPARMS      See BMG2_SymStd_SOLVE_boxmg code.
C     
C     BMG_IOFLAG      See BMG2_SymStd_SOLVE_boxmg code.
C
C     NSOm     See INPUT/OUTPUT
C     
C     NSORm    See WORK ARRAYS.
C     
C     NCIm     See WORK ARRAYS.
C     
C     NOGm     Maximum number of grids that can be supported by the pointer
C              array IGRD. It should be set in calling program.
C
C     NXYc     Limit for number of points on coarsest grid.  Coarsening
C              occurs until the number of x or y unknowns minus one is 
C              less than or equal to NXYc.  NXYc should be 3 or greater or 
C              the code will abort.  This parameter is included because 
C              efficiency is a function of machine.   On vector machines, 
C              it frequently pays to take NXYc larger than on scalar machines.
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
C              BOXMG verifies that NSOm is large enough.
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

      IMPLICIT NONE

C ---------------------------------------------
C     Includes
C ---------------------------------------------

      INCLUDE 'mpif.h'

      INCLUDE 'BMG_PCG_parameters.h'
#include      "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ---------------------------------------------
C     Generic Argument Declarations
C --------------------------------------------- 

      INTEGER  NLx, NLy, NGx, NGy, iGs, jGs

C ---------------------------------------------
C     BOXMG Specific Argument Declarations
C ---------------------------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER   BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND    BMG_rPARMS(NBMG_rPARMS)
      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG)

      !
      ! Workspace: Generic
      !
      INTEGER   NBMG_iWORKm, NBMG_rWORKm
      INTEGER   BMG_iWORK(NBMG_iWORKm)
      REAL*RKIND    BMG_rWORK(NBMG_rWORKm)

      !
      ! Workspace: Coarse-grid Solve
      !
      INTEGER   NBMG_iWORK_CSm, NBMG_rWORK_CSm
      INTEGER   BMG_iWORK_CS(NBMG_iWORK_CSm)
      REAL*RKIND    BMG_rWORK_CS(NBMG_rWORK_CSm)

      !
      ! Workspace: pointers
      !
      INTEGER   BMG_pWORK(NBMG_pWORK)

      !
      ! Solution, RHS, and Stencil
      !
      INTEGER   NFm, NOGm, NSOm
      REAL*RKIND    Q(NFm), QF(NFm), SO(NSOm)

C ---------------------------------------------
C     MPI-Related Argument Declarations
C ---------------------------------------------

      INTEGER  NBMG_MSG_iGRID

      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  MPICOMM, NMSGim, NMSGrm
      INTEGER  iWorkMSG(NMSGim), pMSG(NBMG_pMSG,NOGm),
     &         pMSGSO(NBMG_pMSG,NOGm), pLS(NBMG_pLS*NOGm)
      REAL*RKIND   MSG_Buffer(NMSGrm) 

C ---------------------------------------------
C     PCG Specific Argument Declarations
C ---------------------------------------------

      !
      !  Iterations and convergence history
      !
      INTEGER MAX_PCG_ITERS, PCG_ITERS
      REAL*RKIND  BMG_PCG_RES(0:MAX_PCG_ITERS)

      !
      !  PCG cycle parameters
      !
      INTEGER BMG_PCG_iPARMS(NBMG_PCG_iPARMS)
      REAL*RKIND  BMG_PCG_rPARMS(NBMG_PCG_rPARMS)

C ---------------------------------------------
C     Local Declarations
C ---------------------------------------------

      INTEGER  J, iBC, IFD, IIL, IRELAX, IRELAX_SYM, JJL, NStncl,
     &         NOG, NF, NC, NSO, NCI, NSOR, NCBW, NCU,
     &         p_CI, p_CSO, p_CU, p_iG, p_RES, p_SOR, p_P, p_R, p_Z

      INTEGER  BMG_MSG_MyProc, iter,
     &         MyProc, MyProcI, MyProcJ, 
     &         NProc, NProcI, NProcJ,
     &         NMG_CYCLES, STOP_TEST, OUTPUT, PRECON

      REAL*RKIND   RES_L2, RES_L2_0

      REAL*RKIND tol, derr, malpha, alpha, delta0, delta1, beta, D1MACH
      REAL*8      T1, T2

      EXTERNAL D1MACH

C =========================================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------
      
      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)

      !
      ! PCG parameters
      !
      NMG_CYCLES   = BMG_PCG_iPARMS(id_BMG_PCG_NMG_CYCLES)
      STOP_TEST    = BMG_PCG_iPARMS(id_BMG_PCG_STOP_TEST)
      PRECON       = BMG_PCG_iPARMS(id_BMG_PCG_PRECON)
      OUTPUT       = BMG_PCG_iPARMS(id_BMG_PCG_OUTPUT)

      tol          = BMG_PCG_rPARMS(id_BMG_PCG_STOP_TOL)

      !
      !  BoxMG dimensional parameters
      !
      NOG  = BMG_iPARMS(id_BMG2_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG2_DIM_NF)
      NC   = BMG_iPARMS(id_BMG2_DIM_NC)
      NSO  = BMG_iPARMS(id_BMG2_DIM_NSO)
      NCI  = BMG_iPARMS(id_BMG2_DIM_NCI)
      NSOR = BMG_iPARMS(id_BMG2_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG2_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG2_DIM_NCU)

      !
      !  BoxMG cycle parameters
      !
      IFD          = BMG_iPARMS(id_BMG2_STENCIL)
      IBC          = BMG_iPARMS(id_BMG2_BC) 
      IRELAX       = BMG_iPARMS(id_BMG2_RELAX)
      IRELAX_SYM   = BMG_iPARMS(id_BMG2_RELAX_SYM )

      !
      !  Local pointers: CG
      !
      p_P   = BMG_pWORK(ip_BMG_PCG_P)
      p_R   = BMG_pWORK(ip_BMG_PCG_R)
      p_Z   = BMG_pWORK(ip_BMG_PCG_Z)

      !
      !  Local pointers: BoxMG
      !
      p_RES = BMG_pWORK(ip_RES)
      p_SOR = BMG_pWORK(ip_SOR)
      p_CI  = BMG_pWORK(ip_CI)
      p_CSO = BMG_pWORK(ip_CSO)
      p_CU  = BMG_pWORK(ip_CU)
      p_iG  = BMG_pWORK(ip_iG)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C ---------------------------------    
C     Enable MSG communications
C ---------------------------------

      CALL MSG_set_comm_parent(MPICOMM)
      CALL MSG_enable(BMG_MSG_MyProc,NProc)

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
         WRITE(*,*) '**** ERROR: BMG2_SymStd_SOLVE_boxmg.f ... '
         WRITE(*,*) ' MyProc =         ', MyProc
         WRITE(*,*) ' BMG_MSG_MyProc = ', BMG_MSG_MyProc
         STOP
      ENDIF

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> BEGIN:  PCG ALGORITHM <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      IIL = NLx+2
      JJL = NLy+2

C -------------- peform zeroth iteration as special case -------------

      IF (IFD.NE.1) THEN
         NStncl=5
      ELSE
         NStncl=3
      ENDIF

C -------------------------
C     Perform Setup:
C -------------------------

      IF ( PRECON.EQ.BMG_PCG_PRECON_BMG ) THEN

         IF ( BMG_PCG_iPARMS(id_BMG_PCG_BMG_SETUP)
     &        .EQ. BMG_PCG_BMG_SETUP_all ) THEN

            BMG_iPARMS(id_BMG2_SETUP) = BMG_SETUP_only

            CALL BMG2_SymStd_PRECON_boxmg( 
     &                NLx, NLy, NGx, NGy, iGs, jGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_rWORK(p_Z), BMG_rWORK(p_R), NFm,
     &                SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &                MSG_Buffer, NMSGrm, MPICOMM
     &                )

            IF ( BMG_iPARMS(id_BMG2_Err_Code).NE.iZERO ) THEN
               RETURN
            END IF

            BMG_iPARMS(id_BMG2_SETUP) = BMG_SETUP_none

         ENDIF

      ELSE

         IF ( BMG_PCG_iPARMS(id_BMG_PCG_BMG_SETUP)
     &        .NE. BMG_PCG_BMG_SETUP_none ) THEN

            BMG_iPARMS(id_BMG2_SETUP) = BMG_SETUP_ptrs

            CALL BMG2_SymStd_PRECON_boxmg( 
     &                NLx, NLy, NGx, NGy, iGs, jGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_rWORK(p_Z), BMG_rWORK(p_R), NFm,
     &                SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &                MSG_Buffer, NMSGrm, MPICOMM
     &                )

            IF (BMG_iPARMS(id_BMG2_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            T1 = MPI_Wtime()
            
            CALL BMG2_SymStd_SETUP_fine_stencil(
     &                NOG, SO, NLx+2, NLy+2, NStncl,
     &                iWorkMSG, NMSGim, pMSGSO,
     &                MSG_Buffer, NMSGrm, MPICOMM
     &                )      

            T2 = MPI_Wtime()

            BMG_rPARMS(id_BMG2_TIME_SETUP_FINE_STENCIL)
     &          = BMG_rPARMS(id_BMG2_TIME_SETUP_FINE_STENCIL) + T2 - T1    

            BMG_iPARMS(ID_BMG2_SETUP) = BMG_SETUP_none

         ENDIF

      ENDIF
      
C --------------------------------------- 
C     Calculate the residual R = F - A*x
C ---------------------------------------

      CALL BMG2_SymStd_residual(
     &          NOG, NOG, NOG, RES_L2_0, 
     &          SO, QF, Q, BMG_rWORK(p_R),
     &          IIL, JJL, IFD, NStncl, IRELAX, IRELAX_SYM,
     &          iWorkMSG, NMSGim, pMSG, 
     &          MSG_Buffer, NMSGrm, MPICOMM
     &          )

      CALL BMG2_SymStd_UTILS_norm_l2( 
     &          BMG_rWORK(p_R), IIL, JJL, RES_L2_0, MPICOMM
     &          )

      IF ( MyProc.EQ.1 ) THEN

         IF (  OUTPUT.EQ.BMG_PCG_OUT_ALL .OR. 
     &         OUTPUT.EQ.BMG_PCG_OUT_INIT_FIN_RES .OR.
     &         OUTPUT.EQ.BMG_PCG_OUT_INIT_RES  ) THEN

            PRINT *
            PRINT *,'  ================================================'
            PRINT *
            PRINT *,'  **** Initial Residual (L2 norm) = ',RES_L2_0
            PRINT *
           
         ENDIF

         IF ( OUTPUT.EQ.BMG_PCG_OUT_ALL .OR. 
     &        OUTPUT.EQ.BMG_PCG_OUT_ITS  ) THEN

            PRINT *, '  ==============================================='
            PRINT *, '    Iteration      Stopping Test   Current Error'
            PRINT *, '  ===============================================' 

         ENDIF

      ENDIF

      PCG_ITERS = iZERO      
      BMG_PCG_Res(iZERO) = RES_L2_0

      IF (RES_L2_0 .lt. D1MACH(1) ) RETURN

C -----------------------------------------------------------
C     Compute Z = M^{-1} R where 
C     M^{-1} is BOXMG symmetric V-cycle
C      BMG_PCG_iPARMS(id_BMG_PCG_PRECON) 
C                = BMG_PCG_PRECON_NONE  ->  M = I (IDENTITY)
C      BMG_PCG_iPARMS(id_BMG_PCG_PRECON) 
C                = BMG_PCG_PRECON_DIAG  ->  M = diag(A)
C      BMG_PCG_iPARMS(id_BMG_PCG_PRECON) 
C                = BMG_PCG_PRECON_BMG   ->  M = N BMG NCYCLES
C -----------------------------------------------------------

      IF ( PRECON.eq.BMG_PCG_PRECON_NONE ) then
         !
         CALL BMG2_SymStd_UTILS_dcopy( 
     &             NOG, BMG_rWORK(p_R), BMG_rWORK(p_Z), IIL, JJL,
     &             iONE, iWorkMSG, 
     &             pMSG, MSG_Buffer, MPICOMM
     &             )
         !
      ELSEIF ( PRECON.EQ.BMG_PCG_PRECON_DIAG ) THEN
	 !
         CALL BMG2_SymStd_PRECON_diag(
     &             NOG, SO, BMG_rWORK(p_R), BMG_rWORK(p_Z),
     &             IIL, JJL, NStncl,
     &             iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &             )
	 !
      ELSEIF ( PRECON.EQ.BMG_PCG_PRECON_BMG ) THEN 
         !
         CALL BMG2_SymStd_PRECON_boxmg( 
     &             NLx, NLy, NGx, NGy, iGs, jGs,
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             BMG_rWORK(p_Z), BMG_rWORK(p_R), NFm,
     &             SO, NSOm, NOGm,
     &             BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &             BMG_rWORK, NBMG_rWORKm,
     &             BMG_iWORK_CS, NBMG_iWORK_CSm,
     &             BMG_rWORK_CS, NBMG_rWORK_CSm,
     &             BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &             iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &             MSG_Buffer, NMSGrm, MPICOMM
     &             )
         !
      ENDIF

C ----------------------------------
C     Copy Z into P.
C ----------------------------------
         
      CALL BMG2_SymStd_UTILS_dcopy(
     &          NOG, BMG_rWORK(p_Z), BMG_rWORK(p_P), IIL, JJL,
     &          iONE, iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &          )

      CALL BMG2_SymStd_UTILS_dot_l2(
     &          BMG_rWORK(p_Z), BMG_rWORK(p_R), IIL, JJL,
     &          delta0, MPICOMM
     &          )

C ---------------------------------------
C     Calculate residual norm in either
C     the M-norm or the L2-norm
C ---------------------------------------

      IF ( STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_L2
     &     .OR.STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_L2) THEN
         !
         CALL BMG2_SymStd_UTILS_norm_l2(
     &             BMG_rWORK(p_R), IIL, JJL, RES_L2, MPICOMM
     &             )
         BMG_PCG_Res(1) = RES_L2
         !
      ELSEIF ( STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_M2
     &     .OR.STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_M2) THEN
         !
         BMG_PCG_Res(1) = SQRT(delta0)
         !
      ENDIF

C ------------------------------------------------------
C     Check convergence.  If small, jump to the end!
C ------------------------------------------------------
       
      IF ( STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_L2
     &     .OR. STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_M2 ) THEN
         !
         derr = BMG_PCG_Res(1)
         !
      ELSEIF ( STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_L2
     &     .OR. STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_M2 ) THEN
                                !
         IF (BMG_PCG_Res(1) .gt. 0d0) THEN
            derr = BMG_PCG_Res(1)/BMG_PCG_Res(1)
         ELSE
            derr = 1d0
         ENDIF
         !
      ENDIF
       
      IF ( MyProc.EQ.1 ) THEN

         IF ( OUTPUT.EQ.BMG_PCG_OUT_ALL .OR.
     &          OUTPUT.EQ.BMG_PCG_OUT_ITS ) THEN

            WRITE(*,400) 1, derr, BMG_PCG_Res(1)

         ENDIF
         
      ENDIF
      
      PCG_ITERS = 1

      IF ( derr.LT.tol ) GO TO 200  ! Jump to post processing and return


C ================================
C     Start main loop.
C ================================
                                                                                                        
      DO iter=2, MAX_PCG_ITERS
         
         PCG_ITERS = PCG_ITERS+1

C ---------------------------------------
C
C        Calculate delta0/(Apj,pj) by
C           (2) Zj     = Apj
C           (3) delta1 = (Zj,Pj)
C           (4) alpha  = delta0/delta1
C
C ---------------------------------------

         CALL BMG2_SymStd_UTILS_matvec(
     &             NOG, SO, BMG_rWORK(p_Z), BMG_rWORK(p_P),
     &             IIL, JJL, NOG, IFD, NStncl,
     &             IRELAX, IRELAX_SYM, BMG_UP,
     &             iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &             )
         CALL BMG2_SymStd_UTILS_dot_l2(
     &             BMG_rWORK(p_Z), BMG_rWORK(p_P), IIL, JJL,
     &             delta1, MPICOMM
     &             )

         alpha = delta0/delta1

C ----------------------------------
C        Calculate q <- q + alpha*p
C ----------------------------------
                                                               
         CALL BMG2_SymStd_UTILS_daxpy(
     &             NOG, alpha, BMG_rWORK(p_P), Q, IIL, JJL,
     &             iONE, iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &             )

C ----------------------------------
C        Calculate r <- r - alpha*z
C ----------------------------------
       
         malpha = -1.0*alpha
       
         CALL BMG2_SymStd_UTILS_daxpy(
     &             NOG, malpha, BMG_rWORK(p_Z), BMG_rWORK(p_R),
     &             IIL, JJL, iONE,
     &             iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &             )
                                     
C -----------------------------------------------------------
C        Compute Z = M^{-1} R where
C        M^{-1} is BOXMG symmetric V-cycle
C
C        BMG_PCG_iPARMS(id_BMG_PCG_PRECON)
C               = BMG_PCG_PRECON_NONE  -> M = I (IDENTITY)
C        BMG_PCG_iPARMS(id_BMG_PCG_PRECON)
C               = BMG_PCG_PRECON_DIAG  -> M = diag(A)
C        BMG_PCG_iPARMS(id_BMG_PCG_PRECON)
C               = BMG_PCG_PRECON_BMG   -> M = N BMG NCYCLES
C
C -----------------------------------------------------------
      
         IF ( PRECON.eq.BMG_PCG_PRECON_NONE ) then
            !
            CALL BMG2_SymStd_UTILS_dcopy(
     &                NOG, BMG_rWORK(p_R), BMG_rWORK(p_Z),
     &                IIL, JJL, iONE,
     &                iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &                )
            !
         ELSEIF ( PRECON.eq.BMG_PCG_PRECON_DIAG ) THEN
            !
            CALL BMG2_SymStd_PRECON_diag(
     &                NOG, SO, BMG_rWORK(p_R), BMG_rWORK(p_Z),
     &                IIL, JJL, NStncl,
     &                iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &                )
            !
         ELSEIF ( PRECON.eq.BMG_PCG_PRECON_BMG ) THEN
            !
            CALL BMG2_SymStd_UTILS_rV_zero( BMG_rWORK(p_Z), IIL, JJL )
            !
            CALL BMG2_SymStd_PRECON_boxmg( 
     &                NLx, NLy, NGx, NGy, iGs, jGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_rWORK(p_Z), BMG_rWORK(p_R), NFm,
     &                SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_iWORK_CS, NBMG_rWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &                MSG_Buffer, NMSGrm, MPICOMM
     &                )
            !
         ENDIF

         CALL BMG2_SymStd_UTILS_update_ghosts(
     &             NOG, BMG_rWORK(p_Z), IIL, JJL,
     &             iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &             )

C ------------------------------------------------------
C        Calculate delta1 = <R,Z> = <R,inv(M)R>
C ------------------------------------------------------
       
         CALL BMG2_SymStd_UTILS_dot_l2(
     &             BMG_rWORK(p_R), BMG_rWORK(p_Z), IIL, JJL,
     &             delta1, MPICOMM
     &             )
      
C ---------------------------------------
C        Calculate residual norm in either
C        the M-norm or the L2-norm
C ---------------------------------------
       
         IF ( STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_L2
     &        .OR.STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_L2) THEN
            !
            CALL BMG2_SymStd_UTILS_norm_l2(
     &                BMG_rWORK(p_R), IIL, JJL, RES_L2, MPICOMM
     &                )
            !
            BMG_PCG_Res(PCG_ITERS) = RES_L2
            !
         ELSEIF ( STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_M2
     &           .OR.STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_M2) THEN
            !
            BMG_PCG_Res(PCG_ITERS) = sqrt(delta1)
            !
         ENDIF

C ------------------------------------------------------
C        Check convergence.  If small, jump to the end!
C ------------------------------------------------------
       
         IF ( STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_L2
     &        .OR. STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_M2 ) THEN
            !
            derr = BMG_PCG_Res(PCG_ITERS)
            !
         ELSEIF ( STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_L2
     &        .OR. STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_M2 ) THEN
            !
            derr = BMG_PCG_Res(PCG_ITERS)/BMG_PCG_Res(1)
            !
         ENDIF
       
         IF ( MyProc.EQ.1 ) THEN

            IF ( OUTPUT.EQ.BMG_PCG_OUT_ALL .OR.
     &           OUTPUT.EQ.BMG_PCG_OUT_ITS ) THEN

               WRITE(*,400) PCG_ITERS, derr, BMG_PCG_Res(PCG_ITERS)

            ENDIF

         ENDIF
       
         IF ( derr.lt.tol ) GO TO 200 ! Jump to post processing and return

C ----------------------------------------------
C        Calculate beta.
C ---------------------------------------
       
         beta   = delta1/delta0
         delta0 = delta1
       
C ----------------------------------
C        Calculate p <- z + beta*p
C ----------------------------------
       
         CALL BMG2_SymStd_UTILS_dxpby(
     &             NOG, beta, BMG_rWORK(p_Z), BMG_rWORK(p_P), IIL, JJL,
     &             iONE, iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &             )
                  
         !
      ENDDO
       
 200  CONTINUE      ! <<<<<<<<<<<< OUTSIDE PCG LOOP

      !
      !  Post Processing could be added here
      !
      
      IF ( MyProc.EQ.1 ) THEN
         IF ( OUTPUT.EQ.BMG_PCG_OUT_FIN_RES .OR.
     &        OUTPUT.EQ.BMG_PCG_OUT_INIT_FIN_RES .OR.
     &        OUTPUT.EQ.BMG_PCG_OUT_ALL ) THEN 
            PRINT *,'  ================================================'
            PRINT *
            PRINT *,'  **** Final Residual (L2 norm) = ',
     &           BMG_PCG_RES(PCG_ITERS)
         ENDIF
      ENDIF
      

      IF ( MyProc.EQ.1 .AND. OUTPUT.EQ.BMG_PCG_OUT_ALL ) THEN 
         PRINT *
         PRINT *,'  ================================================'
      ENDIF
       
C =======================================================================

 400  FORMAT ( 6X,I4,9X,1P,E14.7,2x,E14.7 )

C =======================================================================

      RETURN
      END
