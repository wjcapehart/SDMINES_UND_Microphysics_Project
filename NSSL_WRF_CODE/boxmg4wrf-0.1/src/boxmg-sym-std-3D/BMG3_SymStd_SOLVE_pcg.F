      SUBROUTINE BMG3_SymStd_SOLVE_pcg ( 
     &                NLx, NLy, NLz, NGx, NGy, NGz,
     &                iGs, jGs, kGs,
     &                BMG_PCG_iPARMS, BMG_PCG_rPARMS,
     &                BMG_PCG_RES, MAX_PCG_ITERS,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, NFm, SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                BMG_rWORK_PL, NBMG_rWORK_PLm,
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                NProc, MSG_Buffer, NMSGrm, MPICOMM 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This subroutine is a preconditioned conjugate gradient code which
C     may use multigrid preconditioning provided by the black box
C     multigrid code, BMG3_SymStd_SOLVE_boxmg.  It may be used for
C     discretizations of second order elliptic partial differential
C     equations that generate, at most, a 27-point stencil on a logically
C     rectangular grid.  It may be applied to other, similarly
C     structured, problems.  See the BMG3_SymStd_SOLVE_boxmg code for
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
C              IFD .EQ. 1  =>    NSOm = 4*NFm + 10*NCm 
C              IFD .NE. 1  =>    NSOm = 14*NFm 
C
C              BMG3_SymStd_SOLVE_boxmg verifies that NSOm is large enough.
C
C     QF       The user defined array that contains the right hand side.
C              It must be dimensioned to at least NFm.
C
C     Q        the user defined array that contains the solution vector.
C              It must be dimensioned to at least NFm.
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
C              grid. It is then used by the LAPACK routine, which
C              should be dimensioned to (NCBWm,NCUm).
C
C     BBD      User declared real array of dimension NCUm for use in an
C              the LAPACK routine.
C
C     NCBWm    Maximum first subscript of ABD, which needs to be > or =
C              number of grid points plus 2 in x direction on coarsest
C              grid.
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
C     Dendy, J. E. Jr., "Two Multigrid Methods for Three Dimensional
C     Problems with Discountinuos and Anisotropic Coefficients", SIAM
C     Journal of Scientific and Satatistical Computing, Vol. 8, No. 2,
C     September 1987
C     
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'mpif.h'

      INCLUDE 'BMG_PCG_parameters.h'
#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ------------------------------------
C     Argument Declarations
C 
      !
      ! Dimensions, solution, and RHS
      !
      INTEGER NFm, NLx, NLy, NLz, NGx, NGy, NGz, 
     &        iGs, jGs, kGs, NOGm
      REAL*RKIND  Q(NFm), QF(NFm)

      ! 
      ! MPI related arguments
      !
      INTEGER  NBMG_MSG_iGRID,
     &         NBMG_iWORK_CSm, NBMG_iWORK_PLm,
     &         NBMG_rWORK_CSm, NBMG_rWORK_PLm
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID),
     &         BMG_iWORK_CS(NBMG_iWORK_CSm),
     &         BMG_iWORK_PL(NBMG_iWORK_PLm)

      INTEGER MPICOMM, NProc, NMSGim, NMSGrm
      INTEGER iWorkMSG(NMSGim), pMSG(NBMG_pMSG,NOGm), 
     &        pMSGSO(NBMG_pMSG,NOGm), pLS(NBMG_pLS*NOGm)
      REAL*RKIND  MSG_Buffer(NMSGrm)

      !
      ! BoxMG parameters
      !
      INTEGER BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND  BMG_rPARMS(NBMG_rPARMS)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

      !
      ! BoxMG workspace
      !
      INTEGER NBMG_iWORKm, NBMG_rWORKm
      INTEGER BMG_iWORK(NBMG_iWORKm),
     &        BMG_pWORK(NBMG_pWORK)
      REAL*RKIND  BMG_rWORK(NBMG_rWORKm),
     &        BMG_rWORK_CS(NBMG_rWORK_CSm),
     &        BMG_rWORK_PL(NBMG_rWORK_PLm)

      !
      ! BoxMG stencil
      !
      INTEGER NSOm
      REAL*RKIND  SO(NSOm)

      !
      ! PCG 
      !
      INTEGER MAX_PCG_ITERS, PRECON
      INTEGER BMG_PCG_iPARMS(NBMG_PCG_iPARMS)
      REAL*RKIND  BMG_PCG_rPARMS(NBMG_PCG_rPARMS)
      REAL*RKIND  BMG_PCG_RES(MAX_PCG_ITERS)

C ----------------------------
C     Local Declarations
C
      INTEGER  J, IIL, JJL, KKL,
     &         iBC, IFD, IRELAX, IRELAX_SYM, NStncl,
     &         NOG, NF, NC, NSO, NCI, NSOR, NCBW, NCU,
     &         p_CI, p_CSO, p_CU, p_iG, p_RES, p_SOR, 
     &         p_P, pR, pZ, MPI_IERROR

      INTEGER  BMG_MSG_MyProc, MyProc, MyProcI, MyProcJ, MyProcK, 
     &         NProcI, NProcJ, NProcK, NMG_CYCLES, STOP_TEST,
     &         OUTPUT, OUTFREQ

      REAL*8   RES_L2, RES_L2_0

      REAL*RKIND   tol, derr
      REAL*RKIND   malpha, alpha, beta
      REAL*8   delta0, delta1
      REAL*8   TT1, TT2, T1, T2

C ==========================================================================

      !
      ! START TIMER
      !
      t1 = MPI_Wtime()


C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------
      
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

      !
      ! PCG parameters
      !
      NMG_CYCLES   = BMG_PCG_iPARMS(id_BMG_PCG_NMG_CYCLES)
      STOP_TEST    = BMG_PCG_iPARMS(id_BMG_PCG_STOP_TEST)
      PRECON       = BMG_PCG_iPARMS(id_BMG_PCG_PRECON)
      OUTPUT       = BMG_PCG_iPARMS(id_BMG_PCG_OUTPUT)
      OUTFREQ      = BMG_PCG_iPARMS(id_BMG_PCG_OUT_FREQ)

      tol          = BMG_PCG_rPARMS(id_BMG_PCG_STOP_TOL)

      !
      !  BoxMG dimensional parameters
      !
      NOG  = BMG_iPARMS(id_BMG3_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG3_DIM_NF)
      NC   = BMG_iPARMS(id_BMG3_DIM_NC)
      NSO  = BMG_iPARMS(id_BMG3_DIM_NSO)
      NCI  = BMG_iPARMS(id_BMG3_DIM_NCI)
      NSOR = BMG_iPARMS(id_BMG3_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG3_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG3_DIM_NCU)

      !
      !  BoxMG cycle parameters
      !
      IFD          = BMG_iPARMS(id_BMG3_STENCIL)
      IBC          = BMG_iPARMS(id_BMG3_BC) 
      IRELAX       = BMG_iPARMS(id_BMG3_RELAX)
      IRELAX_SYM   = BMG_iPARMS(id_BMG3_RELAX_SYM )

      !
      !  Local pointers: CG
      !
      p_P  = BMG_pWORK(ip_BMG_PCG_P)
      pR   = BMG_pWORK(ip_BMG_PCG_R)
      pZ   = BMG_pWORK(ip_BMG_PCG_Z)

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
         WRITE(*,*) '**** ERROR: BMG3_SymStd_SOLVE_boxmg.f ... '
         WRITE(*,*) ' MyProc =         ', MyProc
         WRITE(*,*) ' BMG_MSG_MyProc = ', BMG_MSG_MyProc
         STOP
      ENDIF

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> BEGIN:  PCG ALGORITHM <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      IIL = NLx+2
      JJL = NLy+2
      KKL = NLz+2

C -------------- peform zeroth iteration as special case -------------

      IF ( IFD.NE.BMG_STENCIL_7pt ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF

C -------------------------
C     Perform Setup:
C -------------------------
    
      IF ( PRECON.eq.BMG_PCG_PRECON_BMG ) THEN

         IF ( BMG_PCG_iPARMS(id_BMG_PCG_BMG_SETUP)
     &        .EQ. BMG_PCG_BMG_SETUP_all ) THEN

            BMG_iPARMS(ID_BMG3_SETUP) = BMG_SETUP_only

            CALL BMG3_SymStd_PRECON_boxmg( 
     &                NLx, NLy, NLz, NGx, NGy, NGz, 
     &                iGs, jGs, kGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_rWORK(pZ), BMG_rWORK(pR), NFm,
     &                SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                BMG_rWORK_PL, NBMG_rWORK_PLm,
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                NProc, MSG_Buffer, NMSGrm, MPICOMM
     &                )

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            BMG_iPARMS(ID_BMG3_SETUP) = BMG_SETUP_none

         ENDIF

      ELSE

         IF ( BMG_PCG_iPARMS(id_BMG_PCG_BMG_SETUP)
     &        .NE. BMG_PCG_BMG_SETUP_none ) THEN

            BMG_iPARMS(ID_BMG3_SETUP) = BMG_SETUP_ptrs

            CALL BMG3_SymStd_PRECON_boxmg( 
     &                NLx, NLy, NLz, NGx, NGy, NGz, 
     &                iGs, jGs, kGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_rWORK(pZ), BMG_rWORK(pR), NFm,
     &                SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                BMG_rWORK_PL, NBMG_rWORK_PLm,
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                NProc, MSG_Buffer, NMSGrm, MPICOMM
     &                )

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            T1 = MPI_Wtime()
            
            CALL BMG3_SymStd_SETUP_fine_stencil(
     &                NOG, SO, 
     &                NLx+2, NLy+2, NLz+2, NStncl,
     &                iWorkMSG, NMSGim, pMSGSO, MSG_Buffer, NMSGrm, 
     &                MPICOMM
     &                )      

            T2 = MPI_Wtime()

            BMG_rPARMS(id_BMG3_TIME_SETUP_FINE_STENCIL)
     &          = BMG_rPARMS(id_BMG3_TIME_SETUP_FINE_STENCIL) + T2 - T1    

            BMG_iPARMS(ID_BMG3_SETUP) = BMG_SETUP_none

         ENDIF

      ENDIF

C --------------------------------------- 
C     Calculate the residual R = F - A*x
C ---------------------------------------

      CALL BMG3_SymStd_residual(
     &                 NOG, NOG, IFD,
     &                 Q, QF, SO, BMG_rWORK(pR),
     &                 IIL, JJL, KKL, NStncl, 
     &                 iWorkMSG, NMSGim, pMSG,
     &                 MSG_Buffer, NMSGrm, MPICOMM
     &                 )

      CALL BMG3_SymStd_UTILS_norm_l2( 
     &                 BMG_rWORK(pR), IIL, JJL, KKL, 
     &                 RES_L2_0, MPICOMM
     &                 )

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


      TT1 = MPI_Wtime()

      IF ( PRECON.eq.BMG_PCG_PRECON_NONE ) then
         !
         CALL BMG3_SymStd_UTILS_dcopy( 
     &        NOG, BMG_rWORK(pR), BMG_rWORK(pZ), 
     &        IIL, JJL, KKL,
     &        iONE, iWorkMSG, 
     &        pMSG, MSG_Buffer, MPICOMM
     &        )
         !
      ELSEIF ( PRECON.EQ.BMG_PCG_PRECON_DIAG ) THEN
	 !
         CALL BMG3_SymStd_PRECON_diag(
     &        NOG, SO, BMG_rWORK(pR), BMG_rWORK(pZ),
     &        IIL, JJL, KKL, NStncl,
     &        iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &        )
	 !
      ELSEIF ( PRECON.EQ.BMG_PCG_PRECON_BMG ) THEN 
         !
         CALL BMG3_SymStd_PRECON_boxmg( 
     &             NLx, NLy, NLz, NGx, NGy, NGz, 
     &             iGs, jGs, kGs, 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             BMG_rWORK(pZ), BMG_rWORK(pR), NFm,
     &             SO, NSOm, NOGm,
     &             BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &             BMG_rWORK, NBMG_rWORKm,
     &             BMG_iWORK_PL, NBMG_iWORK_PLm,
     &             BMG_rWORK_PL, NBMG_rWORK_PLm,
     &             BMG_iWORK_CS, NBMG_iWORK_CSm,
     &             BMG_rWORK_CS, NBMG_rWORK_CSm,
     &             iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &             BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &             NProc, MSG_Buffer, NMSGrm, MPICOMM
     &             )

         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF
         !
      ENDIF

      TT2 = MPI_Wtime()
      BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) = 
     &     BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) + TT2 - TT1

C ----------------------------------
C     Copy Z into P.
C ----------------------------------
         
      CALL BMG3_SymStd_UTILS_dcopy(
     &                 NOG, BMG_rWORK(pZ), BMG_rWORK(p_P), 
     &                 IIL, JJL, KKL, iONE, iWorkMSG, pMSG, 
     &                 MSG_Buffer, MPICOMM
     &                 )

      CALL BMG3_SymStd_UTILS_dot_l2(
     &                 BMG_rWORK(pZ), BMG_rWORK(pR), 
     &                 IIL, JJL, KKL, delta0, MPICOMM
     &                 )

C ---------------------------------------
C     Calculate residual norm in either
C     the M-norm or the L2-norm
C ---------------------------------------

      IF ( STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_L2
     &     .OR.STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_L2) THEN
         !
         CALL BMG3_SymStd_UTILS_norm_l2(
     &        BMG_rWORK(pR), IIL, JJL, KKL, RES_L2, MPICOMM
     &        )
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
         derr = BMG_PCG_Res(1)/BMG_PCG_Res(1)
         !
      ENDIF
       
      IF ( MyProc.EQ.1 ) THEN

         IF ( OUTPUT.EQ.BMG_PCG_OUT_ALL .OR.
     &          OUTPUT.EQ.BMG_PCG_OUT_ITS ) THEN

            WRITE(*,400) 1, derr, BMG_PCG_Res(1)

         ENDIF
         
      ENDIF
  
      IF ( derr.LT.tol ) GO TO 200  ! Jump to post processing and return


C ================================
C     Start main loop.
C ================================
                                                                                   
      DO J=2, MAX_PCG_ITERS

C ---------------------------------------
C
C        Calculate delta0/(Apj,pj) by
C           (2) Zj     = Apj
C           (3) delta1 = (Zj,Pj)
C           (4) alpha  = delta0/delta1
C
C ---------------------------------------

         CALL BMG3_SymStd_UTILS_matvec(
     &                    NOG, SO, BMG_rWORK(pZ), BMG_rWORK(p_P),
     &                    IIL, JJL, KKL, NOG, IFD, NStncl,
     &                    IRELAX, IRELAX_SYM, BMG_UP,
     &                    iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &                    )
         CALL BMG3_SymStd_UTILS_dot_l2(
     &                    BMG_rWORK(pZ), BMG_rWORK(p_P), 
     &                    IIL, JJL, KKL, delta1, MPICOMM
     &                    )

         alpha = delta0/delta1

C ----------------------------------
C        Calculate q <- q + alpha*p
C ----------------------------------
                                                               
         CALL BMG3_SymStd_UTILS_daxpy(
     &                    NOG, alpha, BMG_rWORK(p_P), Q, IIL, JJL, KKL,
     &                    iONE, iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &                    )

C ----------------------------------
C        Calculate r <- r - alpha*z
C ----------------------------------
       
         malpha = -1.0*alpha
       
         CALL BMG3_SymStd_UTILS_daxpy(
     &                    NOG, malpha, BMG_rWORK(pZ), BMG_rWORK(pR),
     &                    IIL, JJL, KKL, iONE, iWorkMSG, pMSG, 
     &                    MSG_Buffer, MPICOMM
     &                    )
                                     
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
      
         TT1 = MPI_Wtime()

         IF ( PRECON.eq.BMG_PCG_PRECON_NONE ) then
            !
            CALL BMG3_SymStd_UTILS_dcopy(
     &                       NOG, BMG_rWORK(pR), BMG_rWORK(pZ),
     &                       IIL, JJL, KKL, iONE,
     &                       iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &                       )
            !
         ELSEIF ( PRECON.eq.BMG_PCG_PRECON_DIAG ) THEN
            !
            CALL BMG3_SymStd_PRECON_diag(
     &                       NOG, SO, BMG_rWORK(pR), BMG_rWORK(pZ),
     &                       IIL, JJL, KKL, NStncl,
     &                       iWorkMSG, pMSG, MSG_Buffer, MPICOMM
     &                       )
            !
         ELSEIF ( PRECON.eq.BMG_PCG_PRECON_BMG ) THEN
            !
            CALL BMG3_SymStd_UTILS_rV_zero( 
     &                NOG, BMG_rWORK(pZ), IIL, JJL, KKL 
     &                )
            !
            CALL BMG3_SymStd_PRECON_boxmg( 
     &                NLx, NLy, NLy, NGx, NGy, NGy, 
     &                iGs, jGs, kGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_rWORK(pZ), BMG_rWORK(pR), NFm,
     &                SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                BMG_rWORK_PL, NBMG_rWORK_PLm,
     &                BMG_iWORK_CS, NBMG_iWORK_CSm,
     &                BMG_rWORK_CS, NBMG_rWORK_CSm,
     &                iWorkMSG, NMSGim, pMSG, pMSGSO, pLS,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                NProc, MSG_Buffer, NMSGrm, MPICOMM
     &                )
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF
            !
         ENDIF

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) =
     &        BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) + TT2 - TT1



C ------------------------------------------------------
C        Calculate delta1 = <R,Z> = <R,inv(M)R>
C ------------------------------------------------------
       
         CALL BMG3_SymStd_UTILS_dot_l2(
     &                    BMG_rWORK(pR), BMG_rWORK(pZ), 
     &                    IIL, JJL, KKL, delta1, MPICOMM
     &                    )
      
C ---------------------------------------
C        Calculate residual norm in either
C        the M-norm or the L2-norm
C ---------------------------------------

         IF ( STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_L2
     &        .OR.STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_L2) THEN
            !
            CALL BMG3_SymStd_UTILS_norm_l2(
     &           BMG_rWORK(pR), IIL, JJL, KKL, RES_L2, MPICOMM
     &           )
            !
            BMG_PCG_Res(J) = RES_L2
            !
         ELSEIF ( STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_M2
     &           .OR.STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_M2) THEN
            !
            BMG_PCG_Res(J) = sqrt(delta1)
            !
         ENDIF

C ------------------------------------------------------
C        Check convergence.  If small, jump to the end!
C ------------------------------------------------------
       
         IF ( STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_L2
     &        .OR. STOP_TEST.EQ.BMG_PCG_STOP_ABS_RES_M2 ) THEN
            !
            derr = BMG_PCG_Res(J)
            !
         ELSEIF ( STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_L2
     &        .OR. STOP_TEST.EQ.BMG_PCG_STOP_REL_RES_M2 ) THEN
            !
            derr = BMG_PCG_Res(J)/BMG_PCG_Res(1)
            !
         ENDIF
       
       
         IF ( MyProc.EQ.1 ) THEN
            !
            IF ( OUTPUT.EQ.BMG_PCG_OUT_ALL .OR.
     &           OUTPUT.EQ.BMG_PCG_OUT_ITS ) THEN
               !
               IF ( MOD(J,OUTFREQ).EQ.0 
     &             .OR. derr.lt.tol .OR. J.EQ.MAX_PCG_ITERS ) THEN
                  WRITE(*,400) J, derr, BMG_PCG_Res(J)
               ENDIF
               !
            ENDIF
            !
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
       
         CALL BMG3_SymStd_UTILS_dxpby(
     &             NOG, beta, BMG_rWORK(pZ), BMG_rWORK(p_P), 
     &             IIL, JJL, KKL, iONE, iWorkMSG, pMSG, MSG_Buffer, 
     &             MPICOMM
     &             )
                  
         !
      ENDDO
       
 200  CONTINUE      ! <<<<<<<<<<<< OUTSIDE PCG LOOP

      T2 = MPI_Wtime()
      BMG_rPARMS(id_BMG3_TIME_PCG_TOTAL) =
     &     BMG_rPARMS(id_BMG3_TIME_PCG_TOTAL) + T2 - T1



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
     &           BMG_PCG_RES(J)
         ENDIF
      ENDIF
      

      IF ( MyProc.EQ.1 .AND. OUTPUT.NE.BMG_PCG_OUT_NONE ) THEN 
         PRINT *
         PRINT *,'  ================================================'
      ENDIF
       
C ==========================================================================

 400  FORMAT ( 6X,I4,9X,1P,E14.7,2x,E14.7 )

C ===========================================

      RETURN
      END
