      SUBROUTINE BMG3_SymStd_SETUP_PtrWork( 
     &                NLx, NLy, NLz, NGx, NGy, NGz, iGs, jGs, kGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &                MyProc, NProc, NProcI, NProcJ, NProcK,
     &                MPICOMM,
     &                NOGm, NFm, NSOm, 
     &                NBMG_iWORKm, NBMG_rWORKm,
     &                NBMG_iWORK_PLm, NBMG_rWORK_PLm,
     &                NBMG_iWORK_CSm, NBMG_rWORK_CSm,
     &                BMG_pWORK, BMG_InWork, pSR, pSI
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
C     NLz       Number of points in z-direction (excluding ghost points)
C
C   ------------------------
C    Discretization:
C   ------------------------
C
C     iSTNCL   Discrete operator index:
C              iSTNCL .EQ. 1 => a 7  point discretization
C              iSTNCL .NE. 1 => a 27  point discretization
C
C   ------------------------
C    Smoothing:
C   ------------------------
C     
C     iRELAX   Relaxation index (refer to BMG3_SymStd_SOLVE_boxmg).
C
C   ------------------------
C    Boundary Conditions:
C   ------------------------
C     
C     iBC      Boundary condition index (refer to BMG3_SymStd_SOLVE_boxmg).
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
      INCLUDE  'BMG_workspace.h'
#include    "BMG_constants.h"
      INCLUDE  'BMG_parameters.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  NFm, NBMG_iWORKm, NBMG_iWORK_CSm, NBMG_iWORK_PLm, NOGm,
     &         NBMG_rWORKm, NBMG_rWORK_CSm, NBMG_rWORK_PLm, NSOm

      INTEGER  BMG_pWORK(NBMG_pWORK), 
     &         NC, NCBW, NCU, NCI, NOG, NF, NSO,
     &         NSOR, NLx, NLy, NLz, NGx, NGy, NGz, pSI, pSR

      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  iGs, jGs, kGs
      INTEGER  MPICOMM, MPI_IERROR, MyProc,
     &         NProc, NProcI, NProcJ, NProcK 

      LOGICAL  BMG_InWORK(NBMG_InWORK)

C ---------------------------
C     Local Declarations:
C
      INTEGER  NCG_TEMP, NIGRD3, NBMG_iWORK, NBMG_iWORK_CS,
     &         NBMG_iWORK_PL, NBMG_RWORK, NBMG_rWORK_CS,
     &         NBMG_rWORK_PL, NMSGi, NMSGr, pSR_S, pSI_S,
     &         NBMG_SER_iWORK, NBMG_SER_rWORK,
     &         NBMG_SER_iWORK_PL, NBMG_SER_rWORK_PL, 
     &         NOGm_SER

      LOGICAL  FAIL_NOG, FAIL_Q, FAIL_RES, FAIL_SO, FAIL_U,
     &         FAIL_WORK, FAIL_WORK_CS, FAIL_WORK_PL

C ==========================================================================

C ==========================================================================
C --------------------------------------------
C    Compute the space requirements:
C --------------------------------------------

      CALL BMG3_SymStd_SETUP_space( 
     &          NLx, NLy, NLz, NGx, NGy, NGz, iGs, jGs, kGs,
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          MyProc, NProc, NProcI, NProcJ, NProcK, 
     &          MPICOMM,
     &          NOG, NF, NC, NSO, NSOR, NCI, NCBW, NCU, 
     &          NCG_TEMP, NBMG_iWORK_PL, NBMG_rWORK_PL,
     &          NBMG_iWORK_CS, NBMG_rWORK_CS,
     &          NBMG_SER_iWORK, NBMG_SER_rWORK,
     &          NBMG_SER_iWORK_PL, NBMG_SER_rWORK_PL,
     &          NOGm_SER, NMSGi, NMSGr
     &          )  

C --------------------------------------------
C     Store the space requirements:
C --------------------------------------------

      BMG_iPARMS(id_BMG3_DIM_NOG)  = NOG
      BMG_iPARMS(id_BMG3_DIM_NF)   = NF
      BMG_iPARMS(id_BMG3_DIM_NC)   = NC
      BMG_iPARMS(id_BMG3_DIM_NCI)  = NCI
      BMG_iPARMS(id_BMG3_DIM_NSO)  = NSO
      BMG_iPARMS(id_BMG3_DIM_NSOR) = NSOR


      BMG_iPARMS(id_BMG3_DIM_NMSGi) = NMSGi
      BMG_iPARMS(id_BMG3_DIM_NMSGr) = NMSGr

C --------------------------------------------
C    Grid Indexing and pointers:
C --------------------------------------------

      NIGRD3 = NBMG_pIGRD*NOG

C --------------------------------------------
C ==========================================================================

C ==========================================================================
C ------------------------------------------------------
C     Compute the pointers into the real work space:
C ------------------------------------------------------

      pSR_S=pSR
      pSI_S=pSI

      !
      ! Include the the PCG workspace vector P in RWORK
      !
      IF ( BMG_inWORK(i_InWORK_PCG_P) ) THEN
         BMG_pWORK(ip_BMG_PCG_P) = pSR
         pSR = BMG_pWORK(ip_BMG_PCG_P) + NF
      ELSE 
         BMG_pWORK(ip_BMG_PCG_P) = iZERO
      ENDIF


      !
      ! Include the the PCG workspace vector R in RWORK
      !
      IF ( BMG_inWORK(i_InWORK_PCG_R) ) THEN
         BMG_pWORK(ip_BMG_PCG_R) = pSR
         pSR = BMG_pWORK(ip_BMG_PCG_R) +  NF
      ELSE 
         BMG_pWORK(ip_BMG_PCG_R) = iZERO
      ENDIF

      !
      ! Include the the PCG workspace vector Z in RWORK
      !
      IF ( BMG_inWORK(i_InWORK_PCG_Z) ) THEN
         BMG_pWORK(ip_BMG_PCG_Z) = pSR
         pSR = BMG_pWORK(ip_BMG_PCG_Z) + NF
      ELSE 
         BMG_pWORK(ip_BMG_PCG_Z) = iZERO
      ENDIF

      !
      ! Include the stencil, SO, in RWORK
      !
      IF ( BMG_InWORK(i_InWORK_SO) ) THEN 
         BMG_pWORK(ip_SO) = pSR
         pSR = BMG_pWORK(ip_SO) + NSO
      ELSE
         BMG_pWORK(ip_SO) = iZERO
      ENDIF

      !
      ! Include the solution, U, in RWORK
      !      
      IF ( BMG_InWORK(i_InWORK_U) ) THEN
         BMG_pWORK(ip_U) = pSR
         pSR = BMG_pWORK(ip_U) + NF
      ELSE
         BMG_pWORK(ip_U) = iZERO
      ENDIF

      !
      ! Include the right hand side, QF, in RWORK
      !
      IF ( BMG_InWORK(i_InWORK_Q) ) THEN
         BMG_pWORK(ip_Q) = pSR
         pSR = BMG_pWORK(ip_Q) + NF
      ELSE 
         BMG_pWORK(ip_Q) = iZERO
      ENDIF

      !
      ! Include the residual, RES, in RWORK
      !
      IF ( BMG_InWORK(i_InWORK_RES) ) THEN
         BMG_pWORK(ip_RES) = pSR
         pSR = BMG_pWORK(ip_RES) + NF
      ELSE 
         BMG_pWORK(ip_RES) = iZERO
      ENDIF

      BMG_pWORK(ip_SOR) = pSR
      BMG_pWORK(ip_CI)  = BMG_pWORK(ip_SOR) + NSOR

      pSR = BMG_pWORK(ip_CI) + NCI

      !
      !  Not used in MPI version
      !
      BMG_pWORK(ip_CSO) = -1
      BMG_pWORK(ip_CU)  = -1

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
      
      BMG_pWORK(ip_MSG) = BMG_pWORK(ip_iG) + NIGRD3

      pSI = BMG_pWORK(ip_MSG) + NMSGi

      !
      ! Not used
      !
      BMG_pWORK(ip_MSGSO) = -1
      BMG_pWORK(ip_LS)    = -1


      NBMG_iWORK = pSI - pSI_S

C ==========================================================================

C ==========================================================================
C --------------------------------------------------------
C     Memory Allocation:
C --------------------------------------------------------   

      IF ( BMG_iPARMS(id_BMG3_POINTERS).EQ.BMG_NO_pointers ) THEN

         !
         !  Verify allocation for standard f77-style workspace
         !

         FAIL_NOG     = (NOG.GT.NOGm)
         FAIL_SO      = (.NOT.BMG_InWORK(i_InWORK_SO)).AND.(NSO.GT.NSOm)
         FAIL_U       = (.NOT.BMG_InWORK(i_InWORK_U)).AND.(NF.GT.NFm)
         FAIL_Q       = (.NOT.BMG_InWORK(i_InWORK_Q)).AND.(NF.GT.NFm)
         FAIL_RES     = (.NOT.BMG_InWORK(i_InWORK_RES)).AND.(NF.GT.NFm)
         FAIL_WORK    = (NBMG_rWORK.GT.NBMG_rWORKm)
     &                  .OR. (NBMG_iWORK.GT.NBMG_iWORKm)
         FAIL_WORK_PL = (NBMG_rWORK_PL.GT.NBMG_rWORK_PLm)
     &                  .OR. (NBMG_iWORK_PL.GT.NBMG_iWORK_PLm)
         FAIL_WORK_CS = (NBMG_rWORK_CS.GT.NBMG_rWORK_CSm)
     &                  .OR. (NBMG_iWORK_CS.GT.NBMG_iWORK_CSm)

         IF ( FAIL_NOG
     &       .OR. FAIL_SO .OR. FAIL_U .OR. FAIL_Q 
     &       .OR. FAIL_RES .OR. FAIL_WORK
     &       .OR. FAIL_WORK_PL .OR. FAIL_WORK_CS ) THEN 

            IF ( MyProc.EQ.0 ) THEN

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
               WRITE(*,*) 'Insufficient Workspace allocated: 3D !!'
               WRITE(*,*) 
               WRITE(*,*) ' - required REAL allocation = ',
     &                    NBMG_rWORK
               WRITE(*,*) ' - actual REAL allocation   = ',
     &                    NBMG_rWORKm
               WRITE(*,*) 
               WRITE(*,*) ' - required INTEGER allocation = ',
     &                    NBMG_iWORK
               WRITE(*,*) ' - actual INTEGER allocation   = ',
     &                    NBMG_iWORKm
               WRITE(*,*)
            ENDIF

            IF ( FAIL_WORK_PL ) THEN
               WRITE(*,*) 'Insufficient workspace allocated: planes !!'
               WRITE(*,*) 
               WRITE(*,*) ' - required REAL allocation = ',
     &                    NBMG_rWORK_PL
               WRITE(*,*) ' - actual REAL allocation   = ', 
     &                    NBMG_rWORK_PLm
               WRITE(*,*) 
               WRITE(*,*) ' - required INTEGER allocation = ',
     &                    NBMG_iWORK_PL
               WRITE(*,*) ' - actual INTEGER allocation   = ',
     &                    NBMG_iWORK_PLm
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
            WRITE(*,*) 'BMG_pWORK(ip_iG)    ',
     &                 BMG_pWORK(ip_iG)
            WRITE(*,*) 'BMG_pWORK(ip_MSG)   ',
     &                 BMG_pWORK(ip_MSG)
            WRITE(*,*) 'BMG_pWORK(ip_MSGSO) ',
     &                 BMG_pWORK(ip_MSGSO)
            WRITE(*,*)
            WRITE(*,*) 'Required INTEGER Allocation  = ', NBMG_iWORK
            WRITE(*,*) 'Actual INTEGER Allocation    = ', NBMG_iWORKm
            WRITE(*,*)
            WRITE(*,*) 

            WRITE(*,*) 'Required REAL allocation for planes = ',
     &                 NBMG_rWORK_PL
            WRITE(*,*) 'Actual REAL allocation for planes = ',
     &                 NBMG_rWORK_PLm
            WRITE(*,*) 
            WRITE(*,*) 'Required INTEGER Allocation for planes  = ',
     &                 NBMG_iWORK_PL
            WRITE(*,*) 'Actual INTEGER Allocation        = ',
     &                 NBMG_iWORK_PLm
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
            
            WRITE(*,104) 

         ENDIF
            
         CALL MPI_FINALIZE(MPI_IERROR)
         STOP

         ENDIF

      ELSEIF ( BMG_iPARMS(id_BMG3_POINTERS)
     &        .EQ. BMG_USE_pointers      ) THEN

         !
         !  Return workspace dimensions for CRAY-style pointers/malloc
         !

         NBMG_iWORKm = NBMG_iWORK
         NBMG_rWORKm = NBMG_rWORK
         
         NBMG_iWORK_PLm = NBMG_iWORK_PL
         NBMG_rWORK_PLm = NBMG_rWORK_PL

         NBMG_iWORK_CSm = NBMG_iWORK_CS
         NBMG_rWORK_CSm = NBMG_rWORK_CS

      ELSE

         !
         !  BMG_iPARMS should be the same on all processors
         !
         IF ( MyProc.EQ.1 ) THEN
          
            WRITE(*,500) '*** MEMORY ALLOCATION MODE UNSPECIFIED ***'

            WRITE(*,*) 'BMG_iPARMS(id_BMG3_POINTERS) = ', 
     &                 BMG_iPARMS(id_BMG3_POINTERS)
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

 100  FORMAT (8X,A,I8)
      
 101  FORMAT (/,2X,A)
 102  FORMAT (10X,A)
 103  FORMAT (2X,A,/)
 104  FORMAT (/)

 500  FORMAT (/,1X,'FATAL ERROR: BMG3_SymStd_SETUP_PtrWork.f',//,1X,A,/)

C ============================


      RETURN
      END


