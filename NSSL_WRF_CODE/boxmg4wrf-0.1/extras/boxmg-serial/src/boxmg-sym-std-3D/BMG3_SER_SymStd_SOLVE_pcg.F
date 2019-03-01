      SUBROUTINE BMG3_SER_SymStd_SOLVE_pcg(
     &                Nx, Ny, Nz,
     &                BMG_SER_PCG_iPARMS, BMG_SER_PCG_rPARMS,
     &                BMG_SER_PCG_RES, MAX_PCG_ITERS, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, NFm, SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                BMG_rWORK_PL, NBMG_rWORK_PLm
     &                )


C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This subroutine is a preconditioned conjugate gradient code which
C     may use multigrid preconditioning provided by the black box
C     multigrid code, BMG3_SER_SymStd_SOLVE_boxmg.  It may be used for
C     discretizations of second order elliptic partial differential
C     equations that generate, at most, a 27-point stencil on a logically
C     rectangular grid.  It may be applied to other, similarly
C     structured, problems.  See the BMG3_SER_SymStd_SOLVE_boxmg code for
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
C     BMG_SER_PCG_iPARMS      PCG integer parameters
C
C     * BMG_SER_PCG_iPARMS(2) --  Maximum number of symmetric V-cycles 
C                             to perform as preconditioner.
C
C     * BMG_SER_PCG_iPARMS(3) --  Type of stop test to use for the pcg 
C                             algorithm.  (See PCG_parameters.h)
C
C     BMG_SER_PCG_rPARMS      PCG real parameters
C
C     * BMG_SER_PCG_rPARMS(1) --  Tolerance used to stop PCG iterations.
C
C     BMG_SER_PCG_NORMits     Array of size PCG_iPARMS(1) containing l2
C                         norm of residual after each pcg iteration.
C
C     BMG_iPARMS      See BMG2_SER_SymStd_SOLVE_boxmg code
C     
C     BMG_rPARMS      See BMG2_SER_SymStd_SOLVE_boxmg code.
C     
C     BMG_IOFLAG      See BMG2_SER_SymStd_SOLVE_boxmg code.
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
C              this with a call to BMG2_SER_SymStd_SETUP_space.  Clearly,
C            
C              IFD .EQ. 1  =>    NSOm = 4*NFm + 10*NCm 
C              IFD .NE. 1  =>    NSOm = 14*NFm 
C
C              BMG3_SER_SymStd_SOLVE_boxmg verifies that NSOm is large enough.
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
      INCLUDE 'BMG_SER_PCG_parameters.h'
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_stencils.h'
      INCLUDE 'BMG_SER_workspace.h'
      INCLUDE 'BMG_SER_parameters.h'

C -------------------------------
C     Argument Declarations
C     
      !
      ! Dimensions, solution, and RHS
      !
      INTEGER  MAX_PCG_ITERS, Nx, Ny, Nz, NFm
      REAL*RKIND   Q(NFm), QF(NFm)

      !
      ! PCG parameters
      !
      INTEGER  BMG_SER_PCG_iPARMS(NBMG_SER_SER_PCG_iPARMS)
      REAL*RKIND   BMG_SER_PCG_rPARMS(NBMG_SER_SER_PCG_rPARMS),
     &         BMG_SER_PCG_Res(MAX_PCG_ITERS) 

      !
      ! BoxMG Specific
      !
      INTEGER  NBMG_iWORKm, NBMG_rWORKm, NSOm, NOG, NOGm,
     &         NBMG_iWORK_PLm, NBMG_rWORK_PLm
      
      INTEGER BMG_iPARMS(NBMG_SER_iPARMS),
     &        BMG_iWORK(NBMG_iWORKm),
     &        BMG_pWORK(NBMG_SER_pWORK),
     &        BMG_iWORK_PL(NBMG_iWORK_PLm)
      REAL*RKIND  BMG_rPARMS(NBMG_SER_rPARMS),
     &        BMG_rWORK(NBMG_rWORKm), SO(NSOm),
     &        BMG_rWORK_PL(NBMG_rWORK_PLm)
      LOGICAL BMG_IOFLAG(NBMG_SER_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  J, II, JJ, KK, T, T1, T2
      INTEGER  NMGCYCLES, STOP_TEST, IFD, IBC
      INTEGER  IRELAX, IRELAX_SYM, NStncl, UPDOWN
      INTEGER  PRECON
      INTEGER  NF, NC, NCI, NSO, NSOR, NCBW, NCU, p_P, p_R, p_Z
      
      REAL*RKIND   tol, derr
      REAL*8       RES_L2_0, RES_L2
      REAL*RKIND   malpha, alpha, delta0, delta1, beta

C ==========================================================================

      T  = rZERO
      T1 = rZERO
      T2 = rZERO

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      !
      !  PCG parameters
      !
      NMGCYCLES    = BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_NMG_CYCLES)
      STOP_TEST    = BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_STOP_TEST)
      PRECON       = BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_PRECON)

      tol          = BMG_SER_PCG_rPARMS(id_BMG_SER_PCG_STOP_TOL)

      !
      !  BoxMG dimensional parameters
      !
      NOG  = BMG_iPARMS(id_BMG3_SER_DIM_NOG)
      NF   = BMG_iPARMS(id_BMG3_SER_DIM_NF)
      NC   = BMG_iPARMS(id_BMG3_SER_DIM_NC)
      NSO  = BMG_iPARMS(id_BMG3_SER_DIM_NSO)
      NCI  = BMG_iPARMS(id_BMG3_SER_DIM_NCI)
      NSOR = BMG_iPARMS(id_BMG3_SER_DIM_NSOR)
      NCBW = BMG_iPARMS(id_BMG3_SER_DIM_NCBW)
      NCU  = BMG_iPARMS(id_BMG3_SER_DIM_NCU)

      !
      !  BoxMG cycle parameters
      !
      IFD          = BMG_iPARMS(id_BMG3_SER_STENCIL)
      IBC          = BMG_iPARMS(id_BMG3_SER_BC) 
      IRELAX       = BMG_iPARMS(id_BMG3_SER_RELAX)
      IRELAX_SYM   = BMG_iPARMS(id_BMG3_SER_RELAX_SYM )

      !
      !  Local pointers
      !
      p_P   = BMG_pWORK(ip_BMG_SER_PCG_P)
      p_R   = BMG_pWORK(ip_BMG_SER_PCG_R)
      p_Z   = BMG_pWORK(ip_BMG_SER_PCG_Z)


C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> BEGIN:  PCG ALGORITHM <<<<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      II = Nx+2
      JJ = Ny+2
      KK = Nz+2

C -------------- peform zeroth iteration as special case -------------

      IF (IFD.NE.1) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF

C -------------------------
C     Perform Setup:
C -------------------------

      IF ( PRECON.eq.BMG_SER_PCG_PRECON_BMG ) THEN

         IF ( BMG_SER_PCG_iPARMS(id_BMG_SER_BMG_iPARMS0_SETUP)
     &        .EQ. BMG_SER_BMG_iPARMS0_SETUP_all ) THEN
 
            BMG_iPARMS(id_BMG3_SER_SETUP) = BMG_SER_SETUP_only

            CALL BMG3_SER_SymStd_PRECON_boxmg( 
     &                Nx, Ny, Nz, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_rWORK(p_Z), BMG_rWORK(p_R), NFm,
     &                SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                BMG_rWORK_PL, NBMG_rWORK_PLm
     &                )

            IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF
         

            BMG_iPARMS(id_BMG3_SER_SETUP) = BMG_SER_SETUP_none

         ENDIF
         
      ENDIF

C --------------------------------------- 
C     Calculate the residual R = F - A*x
C ---------------------------------------

      CALL BMG3_SER_SymStd_residual(
     &          NOG, NOG, IFD, 
     &          Q, QF, SO, BMG_rWORK(p_R), II, JJ, KK, NStncl
     &          )

      CALL BMG3_SER_SymStd_UTILS_norm_l2( 
     &          BMG_rWORK(p_R), II, JJ, KK, RES_L2_0 
     &          )

      WRITE(*,*) 
      WRITE(*,*) 
      WRITE(*,*) ' **** Initial Residual (L2 norm)  = ', RES_L2_0
      WRITE(*,*) 

      WRITE(*,*) '   ================================= '
      WRITE(*,*) '     Iteration       Stopping Test   '
      WRITE(*,*) '   ================================= '


C -----------------------------------------------------------
C     Compute Z = M^{-1} R where 
C     M^{-1} is BOXMG symmetric V-cycle
C      BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_PRECON) 
C                 = BMG_SER_PCG_PRECON_NONE  -> M = I (IDENTITY)
C      BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_PRECON) 
C                 = BMG_SER_PCG_PRECON_DIAG  -> M = diag(A)
C      BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_PRECON) 
C                 = BMG_SER_PCG_PRECON_BMG   -> M = N BMG NCYCLES 
C -----------------------------------------------------------

      IF ( PRECON.eq.BMG_SER_PCG_PRECON_NONE ) THEN
         !
         CALL BMG3_SER_SymStd_UTILS_dcopy(
     &             BMG_rWORK(p_R), BMG_rWORK(p_Z), II, JJ, KK
     &             )
         !
      ELSEIF ( PRECON.eq.BMG_SER_PCG_PRECON_DIAG ) THEN
         !
         CALL BMG3_SER_SymStd_PRECON_diag(
     &             SO, BMG_rWORK(p_R), BMG_rWORK(p_Z),
     &             II, JJ, KK, NStncl
     &              )
         !
      ELSEIF ( PRECON.eq.BMG_SER_PCG_PRECON_BMG ) THEN
         !
         CALL BMG3_SER_SymStd_PRECON_boxmg( 
     &             Nx, Ny, Nz,
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             BMG_rWORK(p_Z), BMG_rWORK(p_R), NFm,
     &             SO, NSOm, NOGm,
     &             BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &             BMG_rWORK, NBMG_rWORKm,
     &             BMG_iWORK_PL, NBMG_iWORK_PLm,
     &             BMG_rWORK_PL, NBMG_rWORK_PLm
     &             )
         IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         !
      ENDIF

C ----------------------------------
C     Copy Z into P.
C ----------------------------------

      CALL BMG3_SER_SymStd_UTILS_dcopy(
     &          BMG_rWORK(p_Z), BMG_rWORK(p_P), II, JJ, KK
     &          )

      CALL BMG3_SER_SymStd_UTILS_dot_l2(
     &          BMG_rWORK(p_Z), BMG_rWORK(p_R), II, JJ, KK, delta0
     &          )

C ---------------------------------------
C     Calculate residual norm in either
C     the M-norm or the L2-norm
C
      IF ( STOP_TEST.EQ.BMG_SER_PCG_STOP_ABS_RES_L2
     &     .OR.STOP_TEST.EQ.BMG_SER_PCG_STOP_REL_RES_L2) THEN
         !
         CALL BMG3_SER_SymStd_UTILS_norm_l2(
     &             BMG_rWORK(p_R), II, JJ, KK, RES_L2
     &             )
         BMG_SER_PCG_Res(1) = RES_L2
         !
      ELSEIF ( STOP_TEST.EQ.BMG_SER_PCG_STOP_ABS_RES_M2
     &     .OR.STOP_TEST.EQ.BMG_SER_PCG_STOP_REL_RES_M2) THEN
         !
         BMG_SER_PCG_Res(1) = SQRT(delta0)
         !
      ENDIF

C ------------------------------------------------------
C     Check convergence.  If small, jump to the end!
C ------------------------------------------------------

      IF ( STOP_TEST.EQ.BMG_SER_PCG_STOP_ABS_RES_L2 
     &     .OR. STOP_TEST.EQ.BMG_SER_PCG_STOP_ABS_RES_M2 ) THEN
         !
         derr = BMG_SER_PCG_Res(1)
         !
      ELSEIF ( STOP_TEST.EQ.BMG_SER_PCG_STOP_REL_RES_L2 
     &     .OR. STOP_TEST.EQ.BMG_SER_PCG_STOP_REL_RES_M2 ) THEN
                                !
         derr = BMG_SER_PCG_Res(1)/BMG_SER_PCG_Res(1)
         !
      ENDIF

      WRITE(*,400) 1, derr

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

         CALL BMG3_SER_SymStd_UTILS_matvec(
     &             BMG_rWORK(p_P), BMG_rWORK(p_Z), SO, 
     &             II, JJ, KK, NStncl 
     &             )
         CALL BMG3_SER_SymStd_UTILS_dot_l2(
     &             BMG_rWORK(p_Z), BMG_rWORK(p_P), II, JJ, KK, delta1
     &             )

         alpha = delta0/delta1

C ----------------------------------
C        Calculate q <- q + alpha*p
C ----------------------------------

         CALL BMG3_SER_SymStd_UTILS_daxpy( 
     &             alpha, BMG_rWORK(p_P), Q, II, JJ, KK
     &             )
         
C ----------------------------------
C        Calculate r <- r - alpha*z
C ----------------------------------
         
         malpha = -1.0*alpha

         CALL BMG3_SER_SymStd_UTILS_daxpy(
     &             malpha, BMG_rWORK(p_Z), BMG_rWORK(p_R), II, JJ, KK
     &             )

C -----------------------------------------------------------
C        Compute Z = M^{-1} R where 
C        M^{-1} is BOXMG symmetric V-cycle
C
C        BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_PRECON) 
C               = BMG_SER_PCG_PRECON_NONE  -> M = I (IDENTITY)
C        BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_PRECON) 
C               = BMG_SER_PCG_PRECON_DIAG  -> M = diag(A)
C        BMG_SER_PCG_iPARMS(id_BMG_SER_PCG_PRECON) 
C               = BMG_SER_PCG_PRECON_BMG   -> M = N BMG NCYCLES
C 
C -----------------------------------------------------------

         IF ( PRECON.eq.BMG_SER_PCG_PRECON_NONE ) then
            !
            CALL BMG3_SER_SymStd_UTILS_dcopy(
     &                BMG_rWORK(p_R), BMG_rWORK(p_Z), II, JJ, KK
     &                )
            !
         ELSEIF ( PRECON.eq.BMG_SER_PCG_PRECON_DIAG ) THEN
            !
            CALL BMG3_SER_SymStd_PRECON_diag(
     &                SO, BMG_rWORK(p_R), BMG_rWORK(p_Z),
     &                II, JJ, KK, NStncl
     &                )
            !
         ELSEIF ( PRECON.eq.BMG_SER_PCG_PRECON_BMG ) THEN
            !
            CALL BMG3_SER_SymStd_UTILS_rV_zero(
     &                BMG_rWORK(p_Z), II, JJ, KK 
     &                )
            !
            CALL BMG3_SER_SymStd_PRECON_boxmg( 
     &                Nx, Ny, Nz,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_rWORK(p_Z), BMG_rWORK(p_R), NFm,
     &                SO, NSOm, NOGm,
     &                BMG_pWORK, BMG_iWORK, NBMG_iWORKm,
     &                BMG_rWORK, NBMG_rWORKm,
     &                BMG_iWORK_PL, NBMG_iWORK_PLm,
     &                BMG_rWORK_PL, NBMG_rWORK_PLm
     &                )
            IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
         ENDIF

C ------------------------------------------------------
C        Calculate delta1 = <R,Z> = <R,inv(M)R>
C ------------------------------------------------------

         CALL BMG3_SER_SymStd_UTILS_dot_l2(
     &             BMG_rWORK(p_R), BMG_rWORK(p_Z), II, JJ, KK, delta1
     &             )
      
C ---------------------------------------
C        Calculate residual norm in either
C        the M-norm or the L2-norm
C ---------------------------------------
    
         IF ( STOP_TEST.EQ.BMG_SER_PCG_STOP_ABS_RES_L2
     &        .OR.STOP_TEST.EQ.BMG_SER_PCG_STOP_REL_RES_L2) THEN
            !
            CALL BMG3_SER_SymStd_UTILS_norm_l2(
     &                BMG_rWORK(p_R), II, JJ, KK, RES_L2 
     &                )
            !
            BMG_SER_PCG_Res(J) = RES_L2
            !
         ELSEIF ( STOP_TEST.EQ.BMG_SER_PCG_STOP_ABS_RES_M2
     &           .OR.STOP_TEST.EQ.BMG_SER_PCG_STOP_REL_RES_M2) THEN
            !
            BMG_SER_PCG_Res(J) = sqrt(delta1)
            !
         ENDIF

C ------------------------------------------------------
C        Check convergence.  If small, jump to the end!
C ------------------------------------------------------

         IF ( STOP_TEST.EQ.BMG_SER_PCG_STOP_ABS_RES_L2 
     &        .OR. STOP_TEST.EQ.BMG_SER_PCG_STOP_ABS_RES_M2 ) THEN
            !
            derr = BMG_SER_PCG_Res(J)
            !
         ELSEIF ( STOP_TEST.EQ.BMG_SER_PCG_STOP_REL_RES_L2 
     &        .OR. STOP_TEST.EQ.BMG_SER_PCG_STOP_REL_RES_M2 ) THEN
            !
            derr = BMG_SER_PCG_Res(J)/BMG_SER_PCG_Res(1)
            !
         ENDIF

         WRITE(*,400) J, derr

         IF ( derr.lt.tol ) GO TO 200 ! Jump to post processing and return

C ----------------------------------------------
C        Calculate beta. 
C ---------------------------------------
     
         beta   = delta1/delta0
         delta0 = delta1

C ----------------------------------
C        Calculate p <- z + beta*p
C ----------------------------------

         CALL BMG3_SER_SymStd_UTILS_dxpby(
     &             beta, BMG_rWORK(p_Z), BMG_rWORK(p_P), II, JJ, KK
     &             )

         !
      ENDDO

 200  CONTINUE      ! <<<<<<<<<<<< OUTSIDE PCG LOOP

      !
      !  Post Processing could be added here
      !

      WRITE(*,*) '   ================================= '

C ==========================================================================

 400  FORMAT ( 8X,I4,9X,1P,E14.7 )

C ============================

      RETURN
      END
