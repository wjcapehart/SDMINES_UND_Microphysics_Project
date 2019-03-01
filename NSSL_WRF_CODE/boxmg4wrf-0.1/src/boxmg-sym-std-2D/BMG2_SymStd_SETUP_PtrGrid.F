      SUBROUTINE BMG2_SymStd_SETUP_PtrGrid( 
     &                NLx, NLy, NGx, NGy, iGs, jGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &                MyProc, Nproc, NprocI, NProcJ, 
     &                MPICOMM, BMG_iWORK_CS,
     &                NOGm, NFm, NCm, NSOm, NSORm, NCIm,
     &                NBMG_iWORK_CSm, NBMG_rWORK_CSm, 
     &                NMSGim, NMSGrm,
     &                NOG, NF, NC, NSO, NSOR, NCI, IGRD,
     &                NBMG_iWORK_CS, NBMG_rWORK_CS,
     &                NMSGi, NMSGr
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SETUP_PtrGrid is used to compute the pointers into 
C     the various arrays for each grid within BOXMG.
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
C     NLx         Number of points in x-direction (excluding ghost points)
C     NLy         Number of points in y-direction (excluding ghost points)
C
C   ------------------------
C    Coarsening Criteria:
C   ------------------------
C
C     NXYc       Minimum number of coarse grid points in either direction
C
C   ------------------------
C    Discretization:
C   ------------------------
C
C     iSTNCL     Discrete operator index:
C                iSTNCL .EQ. 1 => a 5  point discretization
C                iSTNCL .NE. 1 => a 9  point discretization
C
C   ------------------------
C    Smoothing:
C   ------------------------
C
C     iRELAX     Relaxation index (refer to BOXMG).
C
C   ------------------------
C    Boundary Conditions:
C   ------------------------
C
C     iBC        Boundary condition index (refer to BOXMG,BOXMGP)
C
C   ---------------------------------------------
C    (PRE)computed dimensions (which should be):
C   ---------------------------------------------
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
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C     IGRD      Grid pointer array for the arrays internal to BOXMG
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
C ==========================================================================

      IMPLICIT   NONE

C ---------------------------
C     Includes:
C
      INCLUDE    'mpif.h'

#include      "BMG_constants.h"
      INCLUDE    'BMG_parameters.h'
      INCLUDE    'BMG_workspace.h'

      INCLUDE    'BMG_SER_workspace-partial.h'
      INCLUDE    'BMG_SER_parameters.h'


C ---------------------------
C     Argument Declarations:
C
      INTEGER  NOGm, NFm, NCm, NSOm, NSORm, NCIm,
     &         NOG, NF, NC, NSO, NSOR, NCI, NCBW, NCU 

      INTEGER  NMSGi, NMSGim, NMSGr, NMSGrm

      INTEGER  NBMG_iWORK_CS, NBMG_iWORK_CSm, 
     &         NBMG_rWORK_CS, NBMG_rWORK_CSm

      INTEGER  iGs, jGs, MyProc, MPICOMM,
     &         NGx, NGy, NLx, NLy,
     &         NProc, NProcI, NProcJ

      INTEGER  BMG_iPARMS(NBMG_iPARMS),
     &         BMG_iWORK_CS(NBMG_iWORK_CSm),
     &         IGRD(NOGm,NBMG_pIGRD) 

      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ---------------------------
C     Local Declarations:
C
      INTEGER  CF, IERR, i, iStncl, iRELAX, iBC, j, kg,
     &         NGxg, NLxg, NGyg, NLyg, NLXYc,
     &         NxNy, NxNySO, MPI_IERROR, pSI, pSR

      INTEGER  NBMG_SER_iWORK, NBMG_SER_rWORK, NOGm_SER

C ==========================================================================

C ==========================================================================
C --------------------------------------------
C     Unpack Parameters
C --------------------------------------------

      NLXYc  = BMG_iPARMS(id_BMG2_CG_MIN_DIM)
      IStncl = BMG_iPARMS(id_BMG2_STENCIL)
      IRELAX = BMG_iPARMS(id_BMG2_RELAX)
      IBC    = BMG_iPARMS(id_BMG2_BC)

C --------------------------------------------
C ==========================================================================

C ==========================================================================
C --------------------------------------------
C    Initialize *_CS Free space pointers
C --------------------------------------------

      BMG_iWORK_CS(ip_BMG_iWORK_CS_SF) = NBMG_iWORK_CS_ptrs + 1
      BMG_iWORK_CS(ip_BMG_rWORK_CS_SF) = 1

C --------------------------------------------
C    Compute the space requirements:
C --------------------------------------------

      CALL BMG2_SymStd_SETUP_space( 
     &          NLx, NLy, NGx, NGy, IGs, jGs,
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          MyProc, Nproc, NprocI, NProcJ,
     &          MPICOMM,
     &          NOG, NF, NC, NCI, NSO, NSOR, NCBW, NCU,
     &          NBMG_iWORK_CS, NBMG_rWORK_CS,
     &          NBMG_SER_iWORK, NBMG_SER_rWORK,
     &          NOGm_SER, NMSGi, NMSGr
     &          )


      IF ( BMG_iPARMS(id_BMG2_Err_Code).NE.iZERO ) THEN
         RETURN
      END IF

C --------------------------------------------
C    Check the actual allocation:
C --------------------------------------------

      IERR = 0
      !
      !  The number of grids
      !
      IF ( NOG.GT.NOGm ) THEN
         WRITE(*,500) 'The number of grids is greater than NOGm'
         WRITE(*,510) 'HAVE: NOGm = ', NOGm
         WRITE(*,520) 'NEED: NOGm = ', NOG
         STOP
      ENDIF

      !
      !  Total storage for a vector on all grids
      !
      IF ( NF.GT.NFm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Total grid space allocation is too small!'
         WRITE(*,510) 'HAVE: NF = ', NFm
         WRITE(*,520) 'NEED: NF = ', NF
      ENDIF

      !
      !  Total storage for a vector on all coarse grids
      !
      IF ( NC.GT.NCm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Total grid space allocation is too small!'
         WRITE(*,510) 'HAVE: NC = ', NCm
         WRITE(*,520) 'NEED: NC = ', NC
      ENDIF

      !
      !  The dimension of the stencil array
      !  -> NSO is IDISO 
      !
      IF ( NSO.GT.NSOm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Stencil allocation is too small!'
         WRITE(*,510) 'HAVE: NSO = ', NSOm
         WRITE(*,520) 'NEED: NSO = ', NSO
      ENDIF

      !
      !  Workspace for line relaxation's factorization
      !  -> NSOR is IDISOR
      !
      IF ( NSOR.GT.NSORm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Factorization allocation is too small!'
         WRITE(*,510) 'HAVE: NSOR = ', NSORm
         WRITE(*,520) 'NEED: NSOR = ', NSOR
      ENDIF

      !
      !  The dimension of the interpolation operator
      !  -> NCI is IDICI
      !
      IF ( NCI.GT.NCIm ) THEN
         IERR = IERR+1
         WRITE(*,500) 'Interpolation allocation is too small!'
         WRITE(*,510) 'HAVE: NCI = ', NCIm
         WRITE(*,520) 'NEED: NCI = ', NCI
      ENDIF

      !
      !  The dimension of the coarse solve INTEGER workspace 
      !  -> NBMG_iWORK_CS
      !
      IF ( NBMG_iWORK_CS.GT.NBMG_iWORK_CSm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Coarse solve INTEGER workspace is too small!'
         WRITE(*,510) 'HAVE: NBMG_iWORK_PLm = ', NBMG_iWORK_CSm
         WRITE(*,520) 'NEED: NBMG_iWORK_PL  = ', NBMG_iWORK_CS
      ENDIF

      !
      !  The dimension of the coarse solve REAL workspace 
      !  -> NBMG_rWORK_CS
      !
      IF ( NBMG_rWORK_CS.GT.NBMG_rWORK_CSm ) THEN
         IERR=IERR+1
         WRITE(*,500) 'Coarse solve REAL workspace is too small!'
         WRITE(*,510) 'HAVE: NBMG_rWORK_CSm = ', NBMG_rWORK_CSm
         WRITE(*,520) 'NEED: NBMG_rWORK_CS  = ', NBMG_rWORK_CS
      ENDIF

      !
      !  If there were errors STOP
      !
      IF ( IERR.GT.0 .AND. MyProc.EQ.1 ) THEN

         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR)) THEN
            WRITE(*,500)
            WRITE(*,510) 'Insufficient memory allocation !'
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,7)
         RETURN

      ENDIF  

C --------------------------------------------
C ==========================================================================

C -------------------------------------------------------
C     Compute pointers for the coarse solve:
C -------------------------------------------------------

      !
      ! Initialize free-space pointer shifts
      !
      pSI = BMG_iWORK_CS(ip_BMG_iWORK_CS_SF)
      pSR = BMG_iWORK_CS(ip_BMG_rWORK_CS_SF)

      IF ( BMG_iPARMS(id_BMG2_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN

         !
         ! Real workspace pointers
         !
         BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO) = pSR
         pSR = pSR + NCBW*NCU
         BMG_iWORK_CS(ip_BMG_iWORK_CS_CU) = pSR
         pSR = pSR + NCU

         !
         ! Integer workspace pointers
         !
         BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW) = NCBW
         BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)  = NCU
      
      ELSE ! CG_SOLVE_BOXMG

         !
         !  Store number of grids
         !
         BMG_iWORK_CS(id_BMG_NOGm_CS) = NOGm_SER

         !
         !  Store workspace dimensions
         !
         BMG_iWORK_CS(id_BMG_iWORK_CS) = NBMG_SER_iWORK
         BMG_iWORK_CS(id_BMG_rWORK_CS) = NBMG_SER_rWORK

         !
         ! BMG_SER_pWORK
         !
         BMG_iWORK_CS(ip_BMG_pWORK_CS) = pSI
         pSI = pSI + NBMG_SER_pWORK

         !
         ! BMG_SER_iPARMS
         !
         BMG_iWORK_CS(ip_BMG_iPARMS_CS) = pSI
         pSI = pSI +  NBMG_SER_iPARMS

         !
         ! BMG_SER_iWORK
         !
         BMG_iWORK_CS(ip_BMG_iWORK_CS) = pSI
         pSI = pSI +  NBMG_SER_iWORK

         !
         ! BMG_SER_rPARMS
         !
         BMG_iWORK_CS(ip_BMG_rPARMS_CS) = pSR
         pSR = pSR +  NBMG_SER_rPARMS

         !
         ! BMG_SER_rWORK
         ! 
         BMG_iWORK_CS(ip_BMG_rWORK_CS)  = pSR
         pSR = pSR + NBMG_SER_rWORK

      ENDIF

      !
      ! Update the free space counters
      !
      BMG_iWORK_CS(ip_BMG_iWORK_CS_SF) = pSI
      BMG_iWORK_CS(ip_BMG_rWORK_CS_SF) = pSR

C ==========================================================================
C --------------------------------------------
C     Initialize the indexing array.
C --------------------------------------------

      DO j=1, NBMG_pIGRD
         DO i=1, NOGm
            IGRD(i,j)= 0
         END DO
      END DO

C -------------------------------------------
C     Finest Grid :
C -------------------------------------------      

      !
      ! Global dimension; finest grid
      !
      NGxg  = NGx+2  ! including ghost points
      NGyg  = NGy+2  ! including ghost points

      IGRD(NOG,idL_BMG_NGx) = NGxg
      IGRD(NOG,idL_BMG_NGy) = NGyg      
      
      !
      ! local Dimension; finest grid
      !
      NLxg = NLx+2
      NLyg = NLy+2
      NxNy = NLxg * NLyg
      NxNySO = (NLxg+1)*(NLyg+1)

      IGRD(NOG,idL_BMG_NLx) = NLxg
      IGRD(NOG,idL_BMG_NLy) = NLyg

      !
      ! location of the local grid in global coordinates
      !
      IGRD(NOG,idL_BMG_Icoord) = IGs
      IGRD(NOG,idL_BMG_Jcoord) = jGs

      !
      ! Pointer: vectors
      !
      IGRD(NOG,ipL_BMG_U)   = 1
      !
      ! Pointer: stencil
      !
      IGRD(NOG,ipL_BMG_SO)  = 1
      !
      ! Pointer: interpolation
      !
      IGRD(NOG,ipL_BMG_CI)   = 1
      !
      ! Pointer: relaxation workspace
      !
      IGRD(NOG,ipL_BMG_SOR)  = 1

C -------------------------------------------
C     Is there only one grid? => RETURN
C -------------------------------------------

      IF ( NOG.EQ.1 ) THEN
         RETURN
      ENDIF

C -------------------------------------------
C     First Coarse Grid:
C -------------------------------------------      

      !
      ! Pointer: vectors
      !
      IGRD(NOG-1,ipL_BMG_U) = IGRD(NOG,ipL_BMG_U) + NxNy
      !
      ! Pointer: stencil
      !
      IF (iSTNCL.EQ.1) THEN
         IGRD(NOG-1,ipL_BMG_SO) = IGRD(NOG,ipL_BMG_SO)+3*NxNySO
      ELSE
         IGRD(NOG-1,ipL_BMG_SO) = IGRD(NOG,ipL_BMG_SO)+5*NxNySO
      ENDIF
      !
      ! Pointer: interpolation
      !
      IGRD(NOG-1,ipL_BMG_CI) = 1
      !
      ! Pointer: relaxation workspace
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point
     &    .OR. IRELAX.EQ.BMG_GS_RB_x_lines ) THEN
         IGRD(NOG-1,ipL_BMG_SOR) = IGRD(NOG,ipL_BMG_SOR)+2*NxNy
      ELSE IF ( IRELAX.EQ.BMG_GS_RB_y_lines 
     &         .OR. IRELAX.EQ.BMG_GS_RB_x_y_lines ) THEN
         IGRD(NOG-1,ipL_BMG_SOR) = IGRD(NOG,ipL_BMG_SOR)+4*NxNy
      ELSE
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,500) 'A parameter is out of range!'
            WRITE(*,520) 'IRELAX is out of range [1,4] = ', IRELAX
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,5)
         RETURN

      ENDIF

C ----------------------------------
C     Remaining Coarse Grids:
C ----------------------------------

      CF = 1

      NLxg = NLx
      NLyg = NLy

      DO kg = NOG-1, 1, -1
         !
         ! Coarsening: previous grid
         !
         CF   = 2*CF
         NGxg  = (NGx-1)/CF + 3
         NGyg  = (NGy-1)/CF + 3

         IGRD(kg,idL_BMG_NGx) = NGxg
         IGRD(kg,idL_BMG_NGy) = NGyg

         if (mod(IGRD(kg+1,idL_BMG_Icoord),2).eq.1) then 
            IGRD(kg,idL_BMG_Icoord) = (IGRD(kg+1,idL_BMG_Icoord)+1)/2
            NLxg = (NLxg+1)/2
         else
            IGRD(kg,idL_BMG_Icoord) = IGRD(kg+1,idL_BMG_Icoord)/2+1
            IF (mod(NLxg,2).eq.1) THEN
               NLxg = (NLxg-1)/2
            ELSE
               NLxg = (NLxg+1)/2
            ENDIF
         endif

         if (mod(IGRD(kg+1,idL_BMG_JCoord),2).eq.1) then
            IGRD(kg,idL_BMG_Jcoord) = (IGRD(kg+1,idL_BMG_Jcoord)+1)/2
            NLyg = (NLyg+1)/2
         else
            IGRD(kg,idL_BMG_Jcoord) = IGRD(kg+1,idL_BMG_Jcoord)/2+1
            IF (mod(NLyg,2).eq.1) THEN
               NLyg = (NLyg-1)/2
            ELSE
               NLyg = (NLyg+1)/2
            ENDIF
         endif

         !
         ! we must check for empty levels
         ! to avoid ideling nodes
         !  


         NxNy   = (NLxg+2)*(NLyg+2) 
         NxNySO = (NLxg+3)*(NLyg+3)
         !
         ! Dimension: previous grid
         !
         IGRD(kg,idL_BMG_NLx) = NLxg+2
         IGRD(kg,idL_BMG_NLy) = NLyg+2


         if (kg.gt.1) then
            !
            ! Pointer: vectors
            !
            IGRD(kg-1,ipL_BMG_U) = IGRD(kg,ipL_BMG_U)+NxNy
            !
            ! Pointer: stencil
            !
            IGRD(kg-1,ipL_BMG_SO) = IGRD(kg,ipL_BMG_SO)+5*NxNySO
            !
            ! Pointer: interpolation
            !
            IGRD(kg-1,ipL_BMG_CI) = IGRD(kg,ipL_BMG_CI)+8*NxNy
            !
            ! Pointer: relaxation workspace
            !
            IF ( IRELAX.EQ.BMG_GS_RB_point 
     &          .OR. IRELAX.EQ.BMG_GS_RB_x_lines ) THEN
               IGRD(kg-1,ipL_BMG_SOR) = IGRD(kg,ipL_BMG_SOR)+2*NxNy
            ELSE IF ( IRELAX.EQ.BMG_GS_RB_y_lines 
     &               .OR. IRELAX.EQ.BMG_GS_RB_x_y_lines ) THEN
               IGRD(kg-1,ipL_BMG_SOR) = IGRD(kg,ipL_BMG_SOR)+4*NxNy
            ELSE
               WRITE(*,500) 'A parameter is out of range!'
               WRITE(*,520) 'IRELAX is out of range [1,4] = ', IRELAX
               STOP
            ENDIF
            !
         ENDIF
      ENDDO




C ---------------------------------------------
C     Check Dimensions
C ---------------------------------------------

C      Need to add code to double check the workspace pointers are 
C      in bounds.

C --------------------------------------------
C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_PtrGrid.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
