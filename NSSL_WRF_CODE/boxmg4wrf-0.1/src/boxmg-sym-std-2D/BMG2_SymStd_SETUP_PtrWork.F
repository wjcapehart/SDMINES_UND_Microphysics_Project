      SUBROUTINE BMG2_SymStd_SETUP_PtrWork( 
     &                NLx, NLy, NGx, NGy, iGs, jGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &                MPI_MyProc, NProc, NprocI, NprocJ, 
     &                MPICOMM,
     &                NOGm, NFm, NSOm, 
     &                NBMG_iWORKm, NBMG_rWORKm, 
     &                NBMG_iWORK_CSm, NBMG_rWORK_CSm,
     &                BMG_pWORK, BMG_InWORK, pSR, pSI 
     &                )


C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SETUP_PtrWork is used to compute the pointers into the
C     user's work arrays for a direct call to BMG2_SymStd_SOLVE_boxmg.
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
C    Coarsening Criteria:
C   ------------------------
C
C     NXYc      Minimum number of coarse grid points in either direction
C
C   ------------------------
C    Discretization:
C   ------------------------
C
C     iSTNCL   Discrete operator index:
C              iSTNCL .EQ. 1 => a 5  point discretization
C              iSTNCL .NE. 1 => a 9  point discretization
C
C   ------------------------
C    Smoothing:
C   ------------------------
C     
C     iRELAX   Relaxation index (refer to BMG2_SymStd_SOLVE_boxmg).
C
C   ------------------------
C    Boundary Conditions:
C   ------------------------
C     
C     iBC      Boundary condition index (refer to BOXMG,BOXMGP)
C
C   ------------------------
C    Dimensions:
C   ------------------------
C   
C     NOG       Number of grids needed for the given (Nx,Ny)
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
C
C     NCI       Storage for the array CI
C               - CI holds the interpolation weights on all 2D grids
C
C     NCBW      First dimension of ABD, which is set to the bandwidth
C               of the coarsest grid stencil (including the diagonal).
C
C     NCU       Second dimension of ABD, the dimension of the coarse grid.
C
C   ------------------------
C     Workspace:
C   ------------------------
C
C     InWork    Logical array that controls which arrays are contained in rWORK
C
C     NInWork   Dimension of InWork
C
C     NBMG_pWORK    Dimension of pointer array BMG_pWORK
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C     BMG_pWORK     Integer array of (fake) pointers to arrays in rWORK
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C     Nxc       Number of points in the x-direction on a coarser grid
C     Nyc       Number of points in the y-direction on a coarser grid
C
C     kg        Grid index
C
C     NIGRD     INTEGER workspace for grid indexing
C
C ==========================================================================

      IMPLICIT NONE

C ---------------------------
C     Includes
C
      INCLUDE  'BMG_parameters.h'
      INCLUDE  'BMG_workspace.h'
#include      "BMG_constants.h"

C ---------------------------
C     Argument Declarations:
C
      INTEGER  NFm, NBMG_iWORKm, NBMG_iWORK_CSm, NOGm, 
     &         NBMG_rWORKm, NBMG_rWORK_CSm, NSOm, MPICOMM

      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  iGs, jGs

      INTEGER  NC, NCBW, NCU, NCI, NF, NGx, NGy, NLx, NLy, MPI_MyProc,
     &         NOG, NProc, NSO, NSOR, pSI, pSR, BMG_pWORK(NBMG_pWORK),
     &         NMSGi, NMSGr, NProcI, NProcJ

      LOGICAL  BMG_InWORK(NBMG_InWORK)

C ---------------------------
C     Local Declarations:
C
      INTEGER  iBC, iRELAX, iSTNCL, NIGRD,
     &         NBMG_iWORK, NBMG_iWORK_CS,
     &         NBMG_rWORK, NBMG_rWORK_CS,
     &         NBMG_SER_iWORK, NBMG_SER_rWORK,
     &         NOGm_SER,
     &         NXYc, pSI_S, pSR_S, MPI_IERROR
      INTEGER  MyProc

      LOGICAL  FAIL_NOG, FAIL_Q, FAIL_RES, FAIL_SO, FAIL_U,
     &         FAIL_WORK, FAIL_WORK_CS

C ==========================================================================

C ==========================================================================
C --------------------------------------------
C     Unpack Parameters
C --------------------------------------------
      
      NXYc   = BMG_iPARMS(id_BMG2_CG_MIN_DIM)
      IStncl = BMG_iPARMS(id_BMG2_STENCIL)
      IRELAX = BMG_iPARMS(id_BMG2_RELAX)
      IBC    = BMG_iPARMS(id_BMG2_BC)

      MyProc = MPI_MyProc + 1

C ==========================================================================

C ==========================================================================
C --------------------------------------------
C    Compute the space requirements:
C --------------------------------------------

      CALL BMG2_SymStd_SETUP_space( 
     &          NLx, NLy, NGx, NGy, iGs, jGs,
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          MyProc, NProc, NprocI, NProcJ, 
     &          MPICOMM,
     &          NOG, NF, NC, NCI, NSO, NSOR, NCBW, NCU,
     &          NBMG_iWORK_CS, NBMG_rWORK_CS,
     &          NBMG_SER_iWORK, NBMG_SER_rWORK,
     &          NOGm_SER, NMSGi, NMSGr
     &          )

C --------------------------------------------
C     Store the space requirements:
C --------------------------------------------

      BMG_iPARMS(id_BMG2_DIM_NOG)  = NOG
      BMG_iPARMS(id_BMG2_DIM_NF)   = NF
      BMG_iPARMS(id_BMG2_DIM_NC)   = NC
      BMG_iPARMS(id_BMG2_DIM_NCI)  = NCI
      BMG_iPARMS(id_BMG2_DIM_NSO)  = NSO
      BMG_iPARMS(id_BMG2_DIM_NSOR) = NSOR
      BMG_iPARMS(id_BMG2_DIM_NCBW) = NCBW
      BMG_iPARMS(id_BMG2_DIM_NCU)  = NCU

      BMG_iPARMS(id_BMG2_DIM_NMSGi) = NMSGi
      BMG_iPARMS(id_BMG2_DIM_NMSGr) = NMSGr

C --------------------------------------------
C    Grid Indexing and pointers:
C --------------------------------------------

      NIGRD = NBMG_pIGRD*NOG

C --------------------------------------------
C ==========================================================================


C ==========================================================================
C ------------------------------------------------------
C     Compute the pointers into the real work space:
C ------------------------------------------------------

      pSI_S = pSI
      pSR_S = pSR

      !
      ! Include the the PCG workspace vector P in RWORK
      !
      IF (BMG_InWORK(i_InWORK_PCG_P)) THEN
         BMG_pWORK(ip_BMG_PCG_P) = pSR
         pSR = BMG_pWORK(ip_BMG_PCG_P) + NF
      ELSE 
         BMG_pWORK(ip_BMG_PCG_P) = iZERO
      ENDIF


      !
      ! Include the the PCG workspace vector R in RWORK
      !
      IF (BMG_InWORK(i_InWORK_PCG_R)) THEN
         BMG_pWORK(ip_BMG_PCG_R) = pSR
         pSR = BMG_pWORK(ip_BMG_PCG_R) +  NF
      ELSE 
         BMG_pWORK(ip_BMG_PCG_R) = iZERO
      ENDIF

      !
      ! Include the the PCG workspace vector Z in RWORK
      !
      IF (BMG_InWORK(i_InWORK_PCG_Z)) THEN
         BMG_pWORK(ip_BMG_PCG_Z) = pSR
         pSR = BMG_pWORK(ip_BMG_PCG_Z) + NF
      ELSE 
         BMG_pWORK(ip_BMG_PCG_Z) = iZERO
      ENDIF

      !
      ! Include the stencil, SO, in RWORK
      !
      IF (BMG_InWORK(i_InWORK_SO)) THEN 
         BMG_pWORK(ip_SO) = pSR
         pSR = BMG_pWORK(ip_SO) + NSO
      ELSE
         BMG_pWORK(ip_SO) = iZERO
      ENDIF

      !
      ! Include the solution, U, in RWORK
      !      
      IF (BMG_InWORK(i_InWORK_U)) THEN
         BMG_pWORK(ip_U) = pSR
         pSR = BMG_pWORK(ip_U) + NF
      ELSE
         BMG_pWORK(ip_U) = iZERO
      ENDIF

      !
      ! Include the right hand side, QF, in RWORK
      !
      IF (BMG_InWORK(i_InWORK_Q)) THEN
         BMG_pWORK(ip_Q) = pSR
         pSR = BMG_pWORK(ip_Q) + NF
      ELSE 
         BMG_pWORK(ip_Q) = iZERO
      ENDIF

      !
      ! Include the residual, RES, in RWORK
      !
      IF (BMG_InWORK(i_InWORK_RES)) THEN
         BMG_pWORK(ip_RES) = pSR
         pSR = BMG_pWORK(ip_RES) + NF
      ELSE 
         BMG_pWORK(ip_RES) = iZERO
      ENDIF

      !
      ! Storage of SOR and Interpolation 
      !
      BMG_pWORK(ip_SOR) = pSR
      BMG_pWORK(ip_CI)  = BMG_pWORK(ip_SOR) + NSOR

      pSR = BMG_pWORK(ip_CI) + NCI

      !
      ! Not used in the MPI version
      !
      BMG_pWORK(ip_CSO) = -1
      BMG_pWORK(ip_CU)  = -1

      !      BMG_pWORK(ip_CSO) + NCBW*NCU
      !
      !      pSR = BMG_pWORK(ip_CU) + NCU

      !
      ! the buffer for MSG/MPI communication
      !      
      BMG_pWORK(ip_MSG_BUF) = pSR
      pSR = BMG_pWORK(ip_MSG_BUF) + NMSGr

      NBMG_rWORK = pSR - pSR_S

C -------------------------------------------------------
C     Compute the pointers into the integer work space:
C -------------------------------------------------------

      BMG_pWORK(ip_iG) = pSI

      BMG_pWORK(ip_MSG) = BMG_pWORK(ip_iG) + NIGRD

      pSI = BMG_pWORK(ip_MSG) + NMSGi

      NBMG_iWORK = pSI - pSI_S

C --------------------------------------------------------
C     Verify the allocation:
C --------------------------------------------------------   

      IF ( BMG_iPARMS(id_BMG2_POINTERS).EQ.BMG_NO_pointers ) THEN
         !
         !  Verify allocation for standard f77-style workspace
         !
         FAIL_NOG     = (NOG.GT.NOGm)
         FAIL_SO      = (.NOT.BMG_InWORK(i_InWORK_SO)).AND.(NSO.GT.NSOm)
         FAIL_U       = (.NOT.BMG_InWORK(i_InWORK_U)).AND.(NF.GT.NFm)
         FAIL_Q       = (.NOT.BMG_InWORK(i_InWORK_Q)).AND.(NF.GT.NFm)
         FAIL_RES     = (.NOT.BMG_InWORK(i_InWORK_RES)).AND.(NF.GT.NFm)
         FAIL_WORK    = (pSR-1.GT.NBMG_rWORKm).OR.(pSI-1.GT.NBMG_iWORKm)
         FAIL_WORK_CS = (NBMG_rWORK_CS.GT.NBMG_rWORK_CSm)
     &                  .OR. (NBMG_iWORK_CS.GT.NBMG_iWORK_CSm)

         !
         ! WRONG: needs to gather error results and fail if any 
         !        processor fails.  Then control must be returned
         !        to the main code.
         !
         IF ( FAIL_NOG .OR. FAIL_SO .OR. FAIL_U .OR. FAIL_Q 
     &        .OR. FAIL_RES .OR. FAIL_WORK .OR. FAIL_WORK_CS ) THEN 

            IF ( MyProc.EQ.1 ) THEN

            WRITE(*,*) 
            WRITE(*,*) '*** ERROR: MEMORY ALLOCATION OUT OF BOUNDS ***'
            WRITE(*,*) 

            IF ( FAIL_NOG ) THEN 
               WRITE(*,*) 'There are not enough grids !!!! '
               WRITE(*,*) '  - required number of grids, NOG = ', NOG
               WRITE(*,*) '  - maximum number of grids, NOGm = ', NOGm
            ENDIF

            IF ( FAIL_SO ) THEN
               WRITE(*,*) 'The memory allocation for SO failed !!'
               WRITE(*,*) '  - required storage, NSO  = ', NSO
               WRITE(*,*) '  - maximum storage,  NSOm = ', NSOm
               WRITE(*,*) 
            ENDIF

            IF ( FAIL_U ) THEN
               WRITE(*,*) 'The memory allocation for U failed !!'
               WRITE(*,*) '  - required storage, NF  = ', NF
               WRITE(*,*) '  - maximum storage,  NFm = ', NFm
               WRITE(*,*) 
            ENDIF

            IF ( FAIL_Q ) THEN
               WRITE(*,*) 'The memory allocation for Q failed !!'
               WRITE(*,*) '  - required storage, NF  = ', NF
               WRITE(*,*) '  - maximum storage,  NFm = ', NFm
               WRITE(*,*) 
            ENDIF

            IF ( FAIL_RES ) THEN
               WRITE(*,*) 'The memory allocation for RES failed !!'
               WRITE(*,*) '  - required storage, NF  = ', NF
               WRITE(*,*) '  - maximum storage,  NFm = ', NFm
               WRITE(*,*) 
            ENDIF

            IF ( FAIL_WORK ) THEN
               WRITE(*,*) 'Workspace allocation failed !!'
               WRITE(*,*) 
               WRITE(*,*) ' - required REAL allocation, pSR-1 = ',
     &                    pSR - 1
               WRITE(*,*) ' - required REAL allocation, NBMG_rWORK = ',
     &                    NBMG_rWORK
               WRITE(*,*) ' - actual allocation, NBMG_rWORKm = ',
     &                    NBMG_rWORKm
     &                    
               WRITE(*,*) 
               WRITE(*,*) ' - required INTEGER allocation, pSI-1  = ',
     &                    pSI - 1
               WRITE(*,*) ' - required INTEGER allocation,',
     &                    ' NBMG_iWORK = ',
     &                    NBMG_iWORK
               WRITE(*,*) ' - actual allocation, NBMG_iWORKm = ',
     &                    NBMG_iWORKm
               WRITE(*,*)
            ENDIF

            IF ( FAIL_WORK_CS ) THEN
               WRITE(*,*) 'Insufficient Workspace allocated: ',
     &                    'coarse solve !!'
               WRITE(*,*) 
               WRITE(*,*) ' - required REAL allocation = ',
     &                    NBMG_rWORK_CS
               WRITE(*,*) ' - actual REAL allocation   = ',
     &                    NBMG_rWORK_CSm
               WRITE(*,*) 
               WRITE(*,*) ' - required INTEGER allocation = ',
     &                    NBMG_iWORK_CS
               WRITE(*,*) ' - actual INTEGER allocation   = ',
     &                    NBMG_iWORK_CSm
               WRITE(*,*)
            ENDIF


            WRITE(*,*) ' *** CHECK WORKSPACE POINTERS *** '
            WRITE(*,*)
            WRITE(*,*) 'Initial pSI =   ', pSI_S
            WRITE(*,*) 'Initial pSR =   ', pSR_S
            WRITE(*,*)
            WRITE(*,*) 'BMG_pWORK(ip_BMG_PCG_P)  ',
     &                 BMG_pWORK(ip_BMG_PCG_P)
            WRITE(*,*) 'BMG_pWORK(ip_BMG_PCG_R)  ', 
     &                 BMG_pWORK(ip_BMG_PCG_R)
            WRITE(*,*) 'BMG_pWORK(ip_BMG_PCG_Z)  ', 
     &                 BMG_pWORK(ip_BMG_PCG_Z)
            WRITE(*,*) 'BMG_pWORK(ip_SO)         ',
     &                 BMG_pWORK(ip_SO)
            WRITE(*,*) 'BMG_pWORK(ip_U)          ',
     &                 BMG_pWORK(ip_U)
            WRITE(*,*) 'BMG_pWORK(ip_Q)          ',
     &                 BMG_pWORK(ip_Q)
            WRITE(*,*) 'BMG_pWORK(ip_RES)        ',
     &                 BMG_pWORK(ip_RES)
            WRITE(*,*) 'BMG_pWORK(ip_SOR)        ',
     &                 BMG_pWORK(ip_SOR)
            WRITE(*,*) 'BMG_pWORK(ip_CI)         ',
     &                 BMG_pWORK(ip_CI)
            WRITE(*,*) 'BMG_pWORK(ip_CSO)        ',
     &                 BMG_pWORK(ip_CSO)
            WRITE(*,*) 'BMG_pWORK(ip_CU)         ',
     &                 BMG_pWORK(ip_CU)
            WRITE(*,*) 'BMG_pWORK(ip_MSG_BUF)    ',
     &                 BMG_pWORK(ip_MSG_BUF)
            WRITE(*,*) 
            WRITE(*,*) 'Required REAL allocation = ', NBMG_rWORK
            WRITE(*,*) 'Actual REAL Allocation   = ', NBMG_rWORKm
            WRITE(*,*) 
            WRITE(*,*) 'BMG_pWORK(ip_iG)         ',
     &                 BMG_pWORK(ip_iG)
            WRITE(*,*) 'BMG_pWORK(ip_MSG)        ',
     &                 BMG_pWORK(ip_MSG)
            WRITE(*,*) 'BMG_pWork(ip_MSGSO)      ',
     &                 BMG_pWORK(ip_MSGSO)
            WRITE(*,*) 'IWORK(BMG_pWORK(ip_LS)   ',
     &                 BMG_pWORK(ip_LS)
            WRITE(*,*)
            WRITE(*,*) 'Required INTEGER Allocation  = ', NBMG_iWORK
            WRITE(*,*) 'Actual INTEGER Allocation    = ', NBMG_iWORKm
            WRITE(*,*)
            WRITE(*,*)

            WRITE(*,*) 'Required REAL allocation for coarse solve = ',
     &                 NBMG_rWORK_CS
            WRITE(*,*) 'Actual REAL allocation           = ',
     &                 NBMG_rWORK_CSm
            WRITE(*,*) 
            WRITE(*,*) 'Required INTEGER allocation for coarse solve = '
     &                 , NBMG_iWORK_CS
            WRITE(*,*) 'Actual INTEGER Allocation        = ',
     &                 NBMG_iWORK_CSm
            WRITE(*,*)
            WRITE(*,*) 

            WRITE(*,101) '****************************************'
            WRITE(*,102) 'CHECK PARAMETER VALUES'
            WRITE(*,103) '****************************************'
            
            WRITE(*,100) 'Vector on all grids, NF', NF 
            WRITE(*,100) 'Vector on coarse grids, NC', NC
            WRITE(*,100) 'Discrete operator on all grids, NSO', NSO
            WRITE(*,100) 'Workspace, at least, NSOR = 2*NF', NSOR
            WRITE(*,100) 'Interpolation operator coarse grids, NCI', NCI
            
            WRITE(*,100) 'Vector on the coarsest grid, NCU', NCU
            WRITE(*,100) 'Coarsest grid operator bandwidth, NCBW', NCBW
 
         ENDIF

         CALL MPI_FINALIZE(MPI_IERROR)
         STOP

         ENDIF
         !
      ELSEIF ( BMG_iPARMS(id_BMG2_POINTERS)
     &        .EQ. BMG_USE_pointers      ) THEN

         !
         !  Return workspace dimensions for CRAY-style pointers/malloc
         !
         NBMG_iWORKm = NBMG_iWORK
         NBMG_rWORKm = NBMG_rWORK

         NBMG_iWORK_CSm = NBMG_iWORK_CS
         NBMG_rWORK_CSm = NBMG_rWORK_CS

      ELSE

         !
         !  WRONG:  We need to verify that BMG_iPARMS memory 
         !          mode is the same on all processors and valid.
         !
         IF ( MyProc.EQ.1 ) THEN

            WRITE(*,500) '*** MEMORY ALLOCATION MODE UNSPECIFIED ***'

            WRITE(*,*) 'BMG_iPARMS(id_BMG2_POINTERS) = ', 
     &           BMG_iPARMS(id_BMG2_POINTERS)
            WRITE(*,*)
            WRITE(*,*) 'Admissible values are: '
            WRITE(*,*) 
            WRITE(*,*) 'BMG_USE_pointers = ', BMG_USE_pointers
            WRITE(*,*) 'BMG_NO_pointers  = ', BMG_NO_pointers
            
         ENDIF

         CALL MPI_FINALIZE(MPI_IERROR)
         STOP

      ENDIF

C -------------------------------------------------------
C ==========================================================================

 100  FORMAT (8X,A,T24,I10)
      
 101  FORMAT (/,2X,A)
 102  FORMAT (10X,A)
 103  FORMAT (2X,A,/)
 104  FORMAT (/)

 500  FORMAT (/,1X,'FATAL ERROR: BMG2_SymStd_SETUP_PtrWork.f',//,1X,A,/)

C ============================

      RETURN
      END


