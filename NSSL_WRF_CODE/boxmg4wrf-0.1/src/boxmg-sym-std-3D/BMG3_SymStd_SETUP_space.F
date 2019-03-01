      SUBROUTINE BMG3_SymStd_SETUP_space( 
     &                NLx, NLy, NLz, NGx, NGy, NGz, 
     &                iGs, jGs, kGs,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                MyProc, NProc, NProcI, NProcJ, NProcK, 
     &                MPICOMM,
     &                NOG, NF, NC, NSO, NSOR, NCI, NCBW, NCU, 
     &                NCG_TEMP, NBMG_iWORK_PL, NBMG_rWORK_PL,
     &                NBMG_iWORK_CS, NBMG_rWORK_CS,
     &                NBMG_SER_iWORK, NBMG_SER_rWORK,
     &                NBMG_SER_iWORK_PL, NBMG_SER_rWORK_PL,
     &                NOGm_SER, NMSGi, NMSGr
     &                )  

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_SETUP_space is used to compute the storage
C     requirements for the various arrays that are used by 
C     the multigrid solver BMG3_SymStd_SOLVE_boxmg.
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
C     NOG       Number of grids needed for the given (Nx,Ny,Nz)
C
C     NF        Storage for arrays Q and QF on all grids
C     NC        Storage for Q and QF on coarser grids
C               - Q the solution vector on all 3D grids
C               - QF the source vector on all 3D grids
C
C     NSO       Storage for the array SO
C               - SO holds the stencil on all 3D grids.
C
C     NSOR      Storage for the array SOR: 
C               - SOR holds the current residual, the reciprocal of 
C                 the central stencil weight, and the LU decomposition
C                 if iRELAX > 1.

C     NCI       Storage for the array CI
C               - CI holds the 3D interpolation weights on all 3D grids
C
C     NCBW      First dimension of ABD, which is set to the bandwidth
C               of the coarsest grid stencil (including the diagonal).
C     NCU       Second dimension of ABD, the dimension of the coarse grid.
C
C   --------------------------
C    Plane Relaxation:
C   --------------------------
C     
C     NF2       Storage for arrays Q2 and QF2 on all grids.
C
C     NSO3      Storage needed for SO2.
C     NSR3      Storage needed for SOR2.
C     NCI3      Storage needed for CI2.
C     NCBW2     First dimension of ABD2, which is the bandwidth
C               of the 2D coarsest grid stencils (including the diagonal).
C     NCU2      Second dimension of  ABD2, number of unknowns
C               on the largest planes coarsest grid.
C     NCPL2     Third dimension of ABD2, which is the number of 
C               coarse grid 2D planes
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C     Nxc       Number of points in the x-direction on a coarser grid
C     Nyc       Number of points in the y-direction on a coarser grid
C     Nzc       Number of points in the z-direction on a coarser grid
C
C     NPLx      Number of (y,z) planes on a coarser grid
C     NPLy      Number of (x,z) planes on a coarser grid
C     NPLz      Number of (x,y) planes on a coarser grid
C
C     Nx2       -no idea
C     Ny2       -no idea
C
C     NC2       Storage for arrays Q2 and QF2 on coarser grids.
C     
C     NSO2      Storage for SO2 on the largest fine-grid 2D plane
C     NSR2      Storage for SOR2 on the largest fine-grid 2D plane
C     NCI2      Storage for CI2 on the largest fine-grid 2D plane
C
C     kg        Current grid index 
C
C     iTEMP     Temp storage
C
C ==========================================================================

      IMPLICIT NONE

C ---------------------------
C     Includes
C
      INCLUDE 'mpif.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

      INCLUDE 'BMG_SER_parameters.h'
      INCLUDE 'BMG_SER_workspace-partial.h'

C ---------------------------
C     Argument Declarations:
C ---------------------------

      INTEGER  NLx, NLy, NLz, NGx, NGy, NGz, 
     &         iGs, jGs, kGs, MyProc, MPICOMM,
     &         NProc, NProcI, NProcJ, NProcK
      !
      ! BoxMG Parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS), NCG_TEMP,
     &         NBMG_iWORK_CS, NBMG_rWORK_CS,
     &         NBMG_iWORK_PL, NBMG_rWORK_PL

      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      !
      ! BoxMG Dimensions
      !
      INTEGER  NOG, NF, NC,  NSO, NSOR, NCI, NCBW, NCU
      !
      ! MSG communication dimensions
      !
      INTEGER  NMSGi, NMSGr

C ---------------------------
C     Local Declarations:
C ---------------------------

      INTEGER  kg, NBMG_iWORK_OP, NBMG_rWORK_OP,
     &         NCG_TEMP_ZO, NCG_TEMP_YO, NFSO,
     &         NGx_c, NGy_c, NGz_c,
     &         NLx_c, NLy_c, NLz_c,
     &         NLXYZc, NLXYc

      INTEGER  BMG_SER_iPARMS(NBMG_SER_iPARMS),
     &         BMG_SER_pWORK(NBMG_SER_pWORK),
     &         NBMG_SER_iWORK, NBMG_SER_iWORK_PL, 
     &         NBMG_SER_rWORK, NBMG_SER_rWORK_PL,
     &         NFm_SER, NOGm_SER, NSOm_SER,
     &         pSI_SER, pSR_SER
      REAL*RKIND   BMG_SER_rPARMS(NBMG_SER_rPARMS)
      LOGICAL  BMG_SER_InWORK(NBMG_SER_InWork),
     &         BMG_SER_IOFLAG(NBMG_SER_IOFLAG)
      
      INTEGER  i, iBC, ierr, iSTNCL, iRELAX, MPI_IERROR, planesum

C ==========================================================================

C --------------------------------------------
C     Unpack Parameters
C --------------------------------------------

      NLXYZc   = BMG_iPARMS(id_BMG3_CG_MIN_DIM)
      iStncl   = BMG_iPARMS(id_BMG3_STENCIL)
      iRELAX   = BMG_iPARMS(id_BMG3_RELAX)
      iBC      = BMG_iPARMS(id_BMG3_BC)

      NLXYc    = BMG_iPARMS(id_BMG2_CG_MIN_DIM)

C -------------------------------------
C    Check coarsening criteria:
C -------------------------------------

      IF(NLXYZc.LT.3)THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,500) 'Minimum coarse grid dimension is less than 3'
            WRITE(*,520) 'HAVE: NXYZc = ', NLXYZc
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,11)
         RETURN
 
      ENDIF


      IF(NLXYc.LT.3)THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,500) 'Minimum coarse grid dimension is less than 3'
            WRITE(*,520) 'HAVE: NXYc = ', NLXYc
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,12)
         RETURN

      ENDIF

C --------------------------------------------
C    Check relaxation index:
C --------------------------------------------

C      IF ( iRELAX.NE.BMG_RELAX_GS_RB_point .AND. 
C     &           iRELAX.NE.BMG_RELAX_GS_RB_planes_xy_yz_xz ) THEN 
C         WRITE(*,*) 'Fatal Error:  iRELAX can be GS_RB_point or   '
C         WRITE(*,*) '                     GS_RB_planes_xy_yz_xz ...'
C         WRITE(*,*) '  currently iRELAX = ', iRELAX
C         STOP
C      ENDIF

      IF ( iRELAX.NE.BMG_GS_RB_point ) THEN 
         WRITE(*,*) 'Fatal Error:  iRELAX can only be GS_RB_point ... '
         WRITE(*,*) '  currently iRELAX = ', iRELAX
         STOP
      ENDIF

C ==========================================================================

C ==========================================================================
C -------------------------------------
C     Compute the number of 3D grids:
C -------------------------------------

      CALL BMG3_SymStd_SETUP_nog( 
     &                  NLx, NLy, NLz, NLXYZc, NOG, 
     &                  iGs, jGs, kGs, MPICOMM 
     &                  )

C --------------------------------------------
C     Dimensions of the coarsest grid:
C --------------------------------------------

C     !!! NEED proper error handling here !!!

      IF ( NOG.LT.BMG_iPARMS(id_BMG3_MIN_NOG)) THEN
         WRITE(*,500) MyProc, 
     &                'The computed number of grids is too small!'
         WRITE(*,510) MyProc, 
     &                'HAVE: BMG_iPARMS(id_BMG3_MIN_NOG)    = ', 
     &                BMG_iPARMS(id_BMG3_MIN_NOG)
         WRITE(*,510) MyProc, 
     &                'HAVE: BMG_iPARMS(id_BMG3_CG_MIN_DIM) = ',
     &                BMG_iPARMS(id_BMG3_CG_MIN_DIM)
         WRITE(*,520) MyProc, 'HAVE: NOG = ', NOG
         CALL MPI_FINALIZE(MPI_IERROR)
         STOP
      ELSEIF ( NOG.GT.1 ) THEN
         NGx_c = (NGx-1)/2**(NOG-1)+1
         NGy_c = (NGy-1)/2**(NOG-1)+1
         NGz_c = (NGz-1)/2**(NOG-1)+1
      ELSEIF ( NOG.EQ.1 ) THEN
         NGx_c = NGx
         NGy_c = NGy
         NGz_c = NGz
      ELSE
         WRITE(*,500) 'The computed number of grids is less than 1!'
         WRITE(*,510) 'HAVE: Fine grid dimension NGx = ', NGx
         WRITE(*,510) 'HAVE: Fine grid dimension NGy = ', NGy
         WRITE(*,510) 'HAVE: Fine grid dimension NGz = ', NGz
         WRITE(*,520) 'HAVE: NOG = ', NOG
         STOP
      ENDIF    

C --------------------------------------------
C     Storage space for all grids: NF
C --------------------------------------------

      NF=(NLx+2)*(NLy+2)*(NLz+2)
      DO kg=1, NOG-1
         NF = NF + (3+(NLx-1)/2**kg)*(3+(NLy-1)/2**kg)*(3+(NLz-1)/2**kg)
      ENDDO

C ------------------------------------------
C     Storage for the coarse grids (3D):
C ------------------------------------------

      NC = NF - (NLx+2)*(NLy+2)*(NLz+2)

c -----------------------------------------------
c     Storage for all stencils (one coefficient)
c -----------------------------------------------

      NFSO=0
      DO kg=0, NOG-1
         NFSO=NFSO+(5+(NLx-1)/2**kg)*(5+(NLy-1)/2**kg)*(5+(NLz-1)/2**kg)
      ENDDO

C ------------------------------------------
C     Storage for the 3D stencil:
C ------------------------------------------

      IF (iSTNCL.EQ.1) THEN
         NSO = 4*NFSO+10*NC
      ELSE
         NSO = 14*NFSO
      ENDIF

      !PRINT *, NSO, NFSO, NC

C -------------------------------------------
C     Storage for the 3D interpolation:
C -------------------------------------------
      
      NCI = 26*NC

C -------------------------------------------
C     Solving the coarsest grid in 3D
C -------------------------------------------
      
      NBMG_iWORK_CS = 0
      NBMG_rWORK_CS = 0

      !
      ! Allocate space for pointers
      !
      NBMG_iWORK_CS = NBMG_iWORK_CS + NBMG_iWORK_CS_ptrs

      !
      !  Direct solve dimensions
      !    
      IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN

         NCU  = NGx_c*NGy_c*NGz_c
         NCBW = NGx_c*(NGy_c+1)+2
         
         NBMG_rWORK_CS = NBMG_rWORK_CS + NCU + NCBW*NCU

      ELSE ! BMG_CG_SOLVE_BOXMG

         !
         ! Allocate space for BMG_pWORK_CS
         !
         NBMG_iWORK_CS = NBMG_iWORK_CS + NBMG_pWORK

         !
         ! Allocate space for BMG_SER_iPARMS, BMG_SER_rPARMS
         !
         NBMG_iWORK_CS = NBMG_iWORK_CS + NBMG_iPARMS
         NBMG_rWORK_CS = NBMG_rWORK_CS + NBMG_rPARMS

         !
         ! Synchronize the default parameters
         !
         CALL BMG3_SymStd_SETUP_cg_parms(
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             BMG_SER_iPARMS, BMG_SER_rPARMS, BMG_SER_IOFLAG,
     &             NGx_c+2, NGy_c+2, NGz_c+2, NOG
     &             )

         !
         ! We just need dimensional information at this stage
         ! (assume we want consistent parameters in the serial)
         !
         BMG_SER_iPARMS(id_BMG3_SER_POINTERS) = BMG_SER_USE_pointers

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
         CALL BMG3_SER_SymStd_SETUP_PtrWork( 
     &             NGx_c, NGy_c, NGz_c, BMG_SER_iPARMS, 
     &             NOGm_SER, NFm_SER, NSOm_SER, 
     &             NBMG_SER_iWORK, NBMG_SER_rWORK,
     &             NBMG_SER_iWORK_PL, NBMG_SER_rWORK_PL,
     &             BMG_SER_pWORK, BMG_SER_InWork,
     &             pSR_SER, pSI_SER, BMG_SER_IOFLAG
     &             )

         NOGm_SER = BMG_SER_iPARMS(id_BMG3_SER_DIM_NOG)

         !
         ! Accumulate the space requirements
         !
         NBMG_iWORK_CS = NBMG_iWORK_CS 
     &                 + NBMG_SER_iWORK + NBMG_SER_iWORK_PL

         NBMG_rWORK_CS = NBMG_rWORK_CS
     &                 + NBMG_SER_rWORK + NBMG_SER_rWORK_PL

      ENDIF

            
C -------------------------------------------
C     Storage for MSG communications
C -------------------------------------------

      !
      !  Storage for NLx, NLy, NLz arrays for all levels
      ! 
      NMSGi = (NProcI+NProcJ+NProcK)*NOG

      !
      !  Add storage for MSG workspace
      !
      NMSGi = NMSGi + 2*( NProc + 3 
     &                   + 2*( 6*(NLx+5)*(NLy+5)
     &                         + 6*(NLx+5)*(NLz+5)
     &                         + 6*(NLy+5)*(NLz+5) )
     &                   + NOG*( 80 + 20* NProc )
     &                  )

      NMSGi = NMSGi + 2*( NProc + 3 
     &                   + ( 6*(NLx+5)*(NLy+5)
     &                         + 6*(NLx+5)*(NLz+5)
     &                         + 6*(NLy+5)*(NLz+5) )
     &                   + ( 80 + 20* NProc )
     &                  )
      !
      !  Need to fix the this bound!!
      !
      !  - MSG part is probably an incorrect extension of 2D
      !  - Coarse-grid solve
      !

      NMSGr = max(
     &     2*max((NLx+6)*(NLy+6),(NLx+6)*(NLz+6),(NLy+6)*(NLz+6)),
     &     10*(14*NGx_c*NGy_c*NGz_c+7) )

C -------------------------------------------
C     Zero workspace estimates:
C -------------------------------------------

      NBMG_iWORK_OP = 0 
      NBMG_rWORK_OP = 0

      NBMG_iWORK_PL = 0
      NBMG_rWORK_PL = 0 

C ------------------------------------------
C      Workspace for the operator setup:
C ------------------------------------------

      !
      ! Borrow space from plane relaxation work array
      !
      NLy_c = (NLy-1)/2+1
      NLz_c = (NLz-1)/2+1
      
      NCG_TEMP_ZO = 14*(NLx+2)*(NLy + 2)*(NLz_c+2) 
      NCG_TEMP_YO = 14*(NLx+2)*(NLy_c+2)*(NLz_c+2)
      
      NCG_TEMP = NCG_TEMP_YO
      
      NBMG_iWORK_OP = NBMG_iWORK_PL_ptrs
      NBMG_rWORK_OP = NCG_TEMP_YO + NCG_TEMP_ZO

C ------------------------------------------
C     Work space for the Relaxations:
C ------------------------------------------

      IF ( iRELAX.EQ.BMG_GS_RB_point ) THEN

         NSOR = 2*NF

      ELSEIF ( iRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz ) THEN

c$$$         !
c$$$         ! $$$$$ DO NOT USE PLANE RELAXATION YET $$$$$
c$$$         !
c$$$
c$$$         NSOR = 2*NF   ! just to be safe for now, probably can
c$$$                       ! can get rid of this. 
c$$$
c$$$         !
c$$$         ! Allocate space for pointers
c$$$         !
c$$$         NBMG_iWORK_PL = NBMG_iWORK_PL + NBMG_iWORK_PL_ptrs
c$$$
c$$$         !
c$$$         ! Allocation for the pointer array BMG_pWORK_PL
c$$$         !  - one for every plane on every level 
c$$$         !    ( this is an overestimate, but should be easier )
c$$$         !
c$$$         NBMG_iWORK_PL = NBMG_iWORK_PL 
c$$$     &                 + NBMG_pWORK*NOG*((Nx+2)+(Ny+2)+(Nz+2))
c$$$
c$$$         !
c$$$         ! Allocation for the parameter arrays BMG_iPARMS
c$$$         !  - one for each orientation
c$$$         !
c$$$         NBMG_iWORK_PL = NBMG_iWORK_PL + NBMG_iPARMS*NOG*3
c$$$
c$$$         !
c$$$         ! We just need dimensional information at this stage
c$$$         !
c$$$         BMG_iPARMS(id_BMG2_POINTERS)=BMG_USE_pointers
c$$$
c$$$         !
c$$$         ! Initialize irrelevant data
c$$$         !
c$$$         NOGm_2D = 1 
c$$$         NFm_2D = 1
c$$$         NSOm_2D = 1
c$$$
c$$$         pSR_2D = 1
c$$$         pSI_2D = 1
c$$$
c$$$
c$$$         !
c$$$         ! NB: 
c$$$         ! 
c$$$         !
c$$$         ! For the first cut lets put it all in workspace and 
c$$$         ! not worry about the duplication of the solution, RHS 
c$$$         ! and residual space.
c$$$         !
c$$$         ! 
c$$$         ! At some point we'll consider having everything in worspace
c$$$         ! for the first plane on each level and in each orientation
c$$$         ! However, then we'll reuse the solution, RHS, and residual 
c$$$         ! vectors for subsequent planes on the same level and with 
c$$$         ! the same orientation.
c$$$         !
c$$$
c$$$         BMG_InWork_2D(i_InWORK_SO)  = .TRUE.
c$$$         BMG_InWork_2D(i_InWORK_U)   = .TRUE.
c$$$         BMG_InWork_2D(i_InWORK_Q)   = .TRUE.
c$$$         BMG_InWork_2D(i_InWORK_RES) = .TRUE.
c$$$
c$$$         !
c$$$         ! (x,y) - planes
c$$$         !
c$$$         DO kg=2, NOG
c$$$            !
c$$$            !  Set stencil size
c$$$            !
c$$$            IF (kg.NE.NOG .OR. IStncl.NE.1 ) THEN
c$$$               BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_9pt
c$$$            ELSE
c$$$               BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_5pt
c$$$            ENDIF
c$$$            !
c$$$            Nxc = (Nx-1)/2**(NOG-kg)+1
c$$$            Nyc = (Ny-1)/2**(NOG-kg)+1
c$$$            Nzc = (Nz-1)/2**(NOG-kg)+1
c$$$            !
c$$$C            CALL BMG2_SymStd_SETUP_PtrWork( 
c$$$C     &                       Nxc, Nyc, BMG_iPARMS, 
c$$$C     &                       NOGm_2D, NFm_2D, NSOm_2D,
c$$$C     &                       NBMG_iWORK_2D, NBMG_rWORK_2D,
c$$$C     &                       BMG_pWORK_2D, BMG_InWORK_2D,
c$$$C     &                       pSR_2D, pSI_2D 
c$$$C     &                       )
c$$$            NBMG_iWORK_PL = NBMG_iWORK_PL + NBMG_iWORK_2D*Nzc
c$$$            NBMG_rWORK_PL = NBMG_rWORK_PL + NBMG_rWORK_2D*Nzc
c$$$            !
c$$$         ENDDO
c$$$
c$$$         !
c$$$         ! (y,z) - planes
c$$$         !
c$$$         DO kg=2, NOG
c$$$            !
c$$$            !  Set stencil size
c$$$            !
c$$$            IF (kg.NE.NOG .OR. IStncl.NE.1 ) THEN
c$$$               BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_9pt
c$$$            ELSE
c$$$               BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_5pt
c$$$            ENDIF
c$$$            !
c$$$            Nxc = (Nx-1)/2**(NOG-kg)+1
c$$$            Nyc = (Ny-1)/2**(NOG-kg)+1
c$$$            Nzc = (Nz-1)/2**(NOG-kg)+1
c$$$            !
c$$$C            CALL BMG2_SymStd_SETUP_PtrWork( 
c$$$C     &                       Nyc, Nzc, BMG_iPARMS, 
c$$$C     &                       NOGm_2D, NFm_2D, NSOm_2D,
c$$$C     &                       NBMG_iWORK_2D, NBMG_rWORK_2D,
c$$$C     &                       BMG_pWORK_2D, BMG_InWORK_2D,
c$$$C     &                       pSR_2D, pSI_2D 
c$$$C     &                       )
c$$$            NBMG_iWORK_PL = NBMG_iWORK_PL + NBMG_iWORK_2D*Nxc
c$$$            NBMG_rWORK_PL = NBMG_rWORK_PL + NBMG_rWORK_2D*Nxc
c$$$            !
c$$$         ENDDO
c$$$
c$$$         !
c$$$         ! (x,z) - planes
c$$$         !
c$$$         DO kg=2, NOG
c$$$            !
c$$$            !  Set stencil size
c$$$            !
c$$$            IF (kg.NE.NOG .OR. IStncl.NE.1 ) THEN
c$$$               BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_9pt
c$$$            ELSE
c$$$               BMG_iPARMS(id_BMG2_STENCIL) = BMG_STENCIL_5pt
c$$$            ENDIF
c$$$            !
c$$$            Nxc = (Nx-1)/2**(NOG-kg)+1
c$$$            Nyc = (Ny-1)/2**(NOG-kg)+1
c$$$            Nzc = (Nz-1)/2**(NOG-kg)+1
c$$$            !
c$$$C            CALL BMG2_SymStd_SETUP_PtrWork( 
c$$$C     &                       Nxc, Nzc, BMG_iPARMS, 
c$$$C     &                       NOGm_2D, NFm_2D, NSOm_2D,
c$$$C     &                       NBMG_iWORK_2D, NBMG_rWORK_2D,
c$$$C     &                       BMG_pWORK_2D, BMG_InWORK_2D,
c$$$C     &                       pSR_2D, pSI_2D 
c$$$C     &                       )
c$$$            NBMG_iWORK_PL = NBMG_iWORK_PL + NBMG_iWORK_2D*Nyc
c$$$            NBMG_rWORK_PL = NBMG_rWORK_PL + NBMG_rWORK_2D*Nyc
c$$$            !
c$$$         ENDDO

      ENDIF

C ------------------------------------------
C     Dimension of shared workspace:
C ------------------------------------------

      NBMG_iWORK_PL = MAX( NBMG_iWORK_OP, NBMG_iWORK_PL )
      NBMG_rWORK_PL = MAX( NBMG_rWORK_OP, NBMG_rWORK_PL )

C ==========================================================================

 500  FORMAT (/,'[',I4,']',
     &       'FATAL ERROR: BMG3_SymStd_SETUP_space.f',/,5X,A)
 510  FORMAT (5X,'[',I4,']',A,I7)
 520  FORMAT (5X,'[',I4,']',A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
