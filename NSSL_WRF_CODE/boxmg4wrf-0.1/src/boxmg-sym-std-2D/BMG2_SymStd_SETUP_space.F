      SUBROUTINE BMG2_SymStd_SETUP_space( 
     &                NLx, NLy, NGx, NGy, iGs, jGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                MyProc, Nproc, NprocI, NProcJ,
     &                MPICOMM,
     &                NOG, NF, NC, NCI, NSO, NSOR, NCBW, NCU,
     &                NBMG_iWORK_CS, NBMG_rWORK_CS,
     &                NBMG_SER_iWORK, NBMG_SER_rWORK,
     &                NOGm_SER, NMSGi, NMSGr
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SETUP_space is used to compute the space requirements
C     for the various arrays that are used by BOXMG.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C   -------------------------------
C    Local Fine-Grid Dimensions:
C   -------------------------------
C
C     NLx       Number of points in x-direction (excluding ghost points)
C     NLy       Number of points in y-direction (excluding ghost points)
C
C   ------------------------
C    BMG_iPARMS:
C   ------------------------
C
C     * Coarsening Criteria:
C
C     NLXYc     Minimum number of coarse grid points in either direction
C
C     * Discretization:
C
C     iSTNCL   Discrete operator index:
C              iSTNCL .EQ. 1 => a 5  point discretization
C              iSTNCL .NE. 1 => a 9  point discretization
C
C     * Smoothing
C
C     iRELAX   Relaxation index (refer to BOXMG).
C
C     * Boundary Conditions:
C
C     iBC      Boundary condition index (refer to BOXMG,BOXMGP)
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C     NOG       Number of grids needed for the given (Nx,Ny) based on NXYc
C
C     NF        Storage for arrays Q and QF on all grids
C
C     NC        Storage for Q and QF on coarser grids
C               - Q the solution vector on all 2D grids
C               - QF the source vector on all 2D grids
C
C     NSO       Storage for the array SO
C               - SO holds the stencil on all 2D grids.
C
C     NSOR      Storage for the array SOR: 
C               - SOR holds the current residual, the reciprocal of 
C                 the central stencil weight, and the LU decomposition
C                 if iRELAX > 1.
C     NCI       Storage for the array CI
C               - CI holds the interpolation weights on all 2D grids
C
C     NCBW      First dimension of ABD, which is set to the bandwidth
C               of the coarsest grid stencil (including the diagonal).
C     NCU       Second dimension of ABD, the dimension of the coarse grid.
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C     NLx_c     Number of points in the x-direction on a coarser grid
C     NLy_c     Number of points in the y-direction on a coarser grid
C
C     kg        Grid index
C
C ==========================================================================

      IMPLICIT   NONE

C ---------------------------
C     Includes:
C ---------------------------

      INCLUDE 'mpif.h'

#include      "BMG_constants.h"
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

      INCLUDE 'BMG_SER_parameters.h'
      INCLUDE 'BMG_SER_workspace-partial.h'

C ---------------------------
C     Argument Declarations:
C ---------------------------

      !
      ! BoxMG Parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  NBMG_iWORK_CS, NBMG_rWORK_CS

      !
      ! BoxMG Dimensions
      !
      INTEGER  iGs, jGs, MPICOMM, MyProc, NGx, NGy, NLx, NLy,
     &         NProc, NProcI, NProcJ
      INTEGER  NC, NCBW, NCU, NCI, NOG, NF, NSO, NSOR, NMSGi, NMSGr

C ---------------------------
C     Local Declarations:
C
      INTEGER  BMG_SER_iPARMS(NBMG_SER_iPARMS),
     &         BMG_SER_pWORK(NBMG_SER_pWORK),
     &         NBMG_SER_iWORK, NBMG_SER_iWORK_PL, 
     &         NBMG_SER_rWORK, NBMG_SER_rWORK_PL,
     &         NFm_SER, NOGm_SER, NSOm_SER,
     &         pSI_SER, pSR_SER
      REAL*RKIND   BMG_SER_rPARMS(NBMG_SER_rPARMS)
      LOGICAL  BMG_SER_InWORK(NBMG_SER_InWork),
     &         BMG_SER_IOFLAG(NBMG_SER_IOFLAG)
      
      INTEGER  i, iBC, ierr, iRELAX, iSTNCL, kg, NFSO, NLXYc,
     &         NGx_c, NGy_c, NLx_c, NLy_c

C ==========================================================================

C ==========================================================================
C --------------------------------------------
C     Unpack Parameters
C --------------------------------------------

      NLXYc    = BMG_iPARMS(id_BMG2_CG_MIN_DIM)
      IStncl   = BMG_iPARMS(id_BMG2_STENCIL)
      IRELAX   = BMG_iPARMS(id_BMG2_RELAX)
      IBC      = BMG_iPARMS(id_BMG2_BC)

C --------------------------------------------
C     Check coarsening criteria:
C --------------------------------------------

      IF(NLXYc.LT.3)THEN
         WRITE(*,500) 'Minimum coarse grid dimension is less than 3'
         WRITE(*,520) 'HAVE: NXYc = ', NLXYc

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,4)
         RETURN

      ENDIF

C ==========================================================================

C ==========================================================================
C -----------------------------------
C     Compute the number of grids:
C -----------------------------------

      CALL BMG2_SymStd_SETUP_nog( 
     &          NLx, NLy, NLXYc, NOG, iGs, jGs, MPICOMM 
     &          )

C --------------------------------------------
C     Dimensions of the coarsest grid:
C --------------------------------------------

C     !!! NEED proper error handling here !!!
C
C     !!! ANY processor needs to be able to halt execution 
C         and issue the error message !!!
      
      IF ( NOG.LT.BMG_iPARMS(id_BMG2_MIN_NOG)) THEN
C         WRITE(*,500) 'The computed number of grids is too small!'
C         WRITE(*,510) 'HAVE: BMG_iPARMS(id_BMG2_MIN_NOG)    = ', 
C     &                 BMG_iPARMS(id_BMG2_MIN_NOG)
C         WRITE(*,510) 'HAVE: BMG_iPARMS(id_BMG2_CG_MIN_DIM) = ',
C     &                 BMG_iPARMS(id_BMG2_CG_MIN_DIM)
C         WRITE(*,520) 'HAVE: NOG = ', NOG
C         STOP
      ELSEIF ( NOG.GT.1 ) THEN
         NGx_c = (NGx-1)/2**(NOG-1)+1
         NGy_c = (NGy-1)/2**(NOG-1)+1
      ELSEIF ( NOG.EQ.1 ) THEN
         NGx_c = NGx
         NGy_c = NGy
      ELSE
C         WRITE(*,500) 'The computed number of grids is less than 1!'
C         WRITE(*,510) 'HAVE: Fine grid dimension NGx = ', NGx
C         WRITE(*,510) 'HAVE: Fine grid dimension NGy = ', NGy
C         WRITE(*,520) 'HAVE: NOG = ', NOG
C         STOP
      ENDIF

C --------------------------------------------
C     Storage on all grids: NF
C --------------------------------------------

      NF=(NLx+2)*(NLy+2)
      DO 30 kg=1, NOG-1
         NF=NF+(3+(NLx-1)/2**kg)*(3+(NLy-1)/2**kg)
 30   CONTINUE

C ------------------------------------------
C     Storage on the coarse grids: NC
C ------------------------------------------

      NC = NF-(NLx+2)*(NLy+2)

c -----------------------------------------------
c     Storage for all stencils (one coefficient)
c -----------------------------------------------

      NFSO=0
      DO 35 kg=0, NOG-1
         NFSO=NFSO+(5+(NLx-1)/2**kg)*(5+(NLy-1)/2**kg)
 35   continue
  
C ----------------------------------------
C     Interpolation operator:
C ----------------------------------------

      IF (NC.EQ.0) THEN
         NCI=2
      ELSE
         NCI=8*NC
      ENDIF
 
C -----------------------------------------
C     Stencil:
C -----------------------------------------

      IF ( ISTNCL.EQ.1 ) THEN
         NSO = 3*NFSO + 2*NC 
      ELSE
         NSO = 5*NFSO
      ENDIF


C ----------------------------------------
C     Storage for the workspace SOR:
C ----------------------------------------

      IF ( iRELAX.EQ.BMG_GS_RB_point 
     &    .OR. iRELAX.EQ.BMG_GS_RB_x_lines ) THEN
         NSOR = 2*NF
      ELSE 
         NSOR = 4*NF
      ENDIF

C ----------------------------------------
C     Solving on the coarsest grid:
C ----------------------------------------

      NBMG_iWORK_CS = 0
      NBMG_rWORK_CS = 0

      !
      ! Allocate space for pointers
      !
      NBMG_iWORK_CS = NBMG_iWORK_CS + NBMG_iWORK_CS_ptrs

      !
      !  Direct solve dimensions
      !    
      IF ( BMG_iPARMS(id_BMG2_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN

         NCU  = NGx_c*NGy_c
         NCBW = NGx_c+2       ! for non-periodic 
         
         NBMG_rWORK_CS = NBMG_rWORK_CS + NCU + NCBW*NCU

      ELSE ! BMG_CG_SOLVE_BOXMG

         !
         ! Allocate space for BMG_pWORK_CS
         !
         NBMG_iWORK_CS = NBMG_iWORK_CS + NBMG_SER_pWORK

         !
         ! Allocate space for BMG_SER_iPARMS, BMG_SER_rPARMS
         !
         NBMG_iWORK_CS = NBMG_iWORK_CS + NBMG_SER_iPARMS
         NBMG_rWORK_CS = NBMG_rWORK_CS + NBMG_SER_rPARMS

         !
         ! Synchronize the default parameters
         !
         CALL BMG2_SymStd_SETUP_cg_parms(
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             BMG_SER_iPARMS, BMG_SER_rPARMS, BMG_SER_IOFLAG,
     &             NGx_c+2, NGy_c+2, NOG
     &             )

         !
         ! We just need dimensional information at this stage
         ! (assume we want consistent parameters in the serial)
         !
         BMG_SER_iPARMS(id_BMG2_SER_POINTERS) = BMG_SER_USE_pointers

         !
         ! Allocate multigrid related arrays in the workspace
         !
         DO i=1, NBMG_SER_InWORK
            BMG_SER_InWORK(i) = .FALSE.
         ENDDO
         
         BMG_SER_InWork(i_InWORK_SO)  = .TRUE.
         BMG_SER_InWork(i_InWORK_U)   = .TRUE.
         BMG_SER_InWork(i_InWORK_Q)   = .TRUE.
         BMG_SER_InWork(i_InWORK_RES) = .TRUE.

         !
         ! Initialize irrelevant data
         !
         NOGm_SER = 1 
         NFm_SER = 1
         NSOm_SER = 1

         pSR_SER = 1
         pSI_SER = 1

         !
         ! Get the workspace estimate ...
         !
         CALL BMG2_SER_SymStd_SETUP_PtrWork( 
     &             NGx_c, NGy_c, BMG_SER_iPARMS, 
     &             NOGm_SER, NFm_SER, NSOm_SER, 
     &             NBMG_SER_iWORK, NBMG_SER_rWORK,
     &             BMG_SER_pWORK, BMG_SER_InWork,
     &             pSR_SER, pSI_SER, BMG_SER_IOFLAG
     &             )

         NOGm_SER = BMG_SER_iPARMS(id_BMG2_SER_DIM_NOG)

         !
         ! Accumulate the space requirements
         !
         NBMG_iWORK_CS = NBMG_iWORK_CS + NBMG_SER_iWORK

         NBMG_rWORK_CS = NBMG_rWORK_CS + NBMG_SER_rWORK

      ENDIF

C -------------------------------------------
C     Storage for MSG communications
C -------------------------------------------

      !
      !  Storage for NLx, NLy, NLz arrays for all levels
      ! 
      NMSGi = (NProcI+NProcJ)*NOG

      !
      !  Add storage for MSG workspace (Not a sharp bound)
      !
      NMSGi = NMSGi + NProc+2 + NOG*(16*(NLx+NLy+6)+18*NProc+26) 

      !
      ! Add storage for MSGSO workspace
      !
      NMSGi = NMSGi + NProc+2 + NOG*(24*(NLx+NLy+4)+18*NProc+26)

      !
      ! Add storage for Line Solves
      !
      NMSGi = NMSGi + 4*NOG*NProc

      !
      !  Need to fix the this bound!!
      !
      !  - MSG real buffer space 
      !  - Workspace for coarse-grid solve communication.
      !
      NMSGr = max(  2*max(NLx+6, NLy+6), 
     &              10*(5*NGx_c*NGy_c+5),
     &              4*NLy*NProcI + 5*NLx, 4*NLx*NProcJ + 5*NLy )

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_space',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END

