      SUBROUTINE BMG3_SER_SymStd_SETUP_parts( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                ABD, BBD, NABD1, NABD2, IGRD, NOGm, NOG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_SETUP_parts.f loops over all grids calling the
C     routines to construct the operator induced interpolation and the
C     coarse grid operators.  It also factors the operator on the
C     coarsest grid.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C -----------------------------
C     Includes
C
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_stencils.h'
      INCLUDE 'BMG_SER_workspace.h'
      INCLUDE 'BMG_SER_parameters.h'

C ---------------------------
C    Argument Declarations:
C
      INTEGER  NABD1, NABD2, NBMG_iWORK_PL, NBMG_rWORK_PL, 
     &         NCI, NOGm, NSO, NSOR

      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS),
     &         BMG_iWORK_PL(NBMG_iWORK_PL),
     &         IGRD(NOGm*NBMG_SER_pIGRD), NOG
      REAL*RKIND   abd(nabd1*nabd2), bbd(nabd2),
     &         BMG_rPARMS(NBMG_SER_rPARMS), 
     &         BMG_rWORK_PL(NBMG_rWORK_PL), CI(NCI), SO(NSO),
     &         SOR(NSOR)
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)

C ---------------------------
C    Local Declarations:
C
      INTEGER  i, IFD, IBC, IRELAX, K, NStncl, NStncl_CG,
     &         Nx, Nxc, Ny, Nyc, Nz, Nzc, 
     &         p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC,
     &         p_CGTEMP_yo, p_CGTEMP_zo

C ==========================================================================

      !
      ! Unpack
      !
      IFD    = BMG_iPARMS(id_BMG3_SER_STENCIL)
      IBC    = BMG_iPARMS(id_BMG3_SER_BC)
      IRELAX = BMG_iPARMS(id_BMG3_SER_RELAX)

      !
      ! Coarse-grid operators are always 27 points.
      !
      NStncl_CG=14

      !
      ! Sanity check: NOG
      !
      IF ( NOG.EQ.0 ) THEN
         !
         IF (BMG_IOFLAG(iBMG3_SER_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'There are no grids?'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF
         
         CALL BMG3_SER_SymStd_ErrTrap(BMG_iPARMS,2)
         RETURN
         !
      ELSE IF ( NOG.LT.0 ) THEN
         !
         IF (BMG_IOFLAG(iBMG3_SER_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The number of grids is negative!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF
         
         CALL BMG3_SER_SymStd_ErrTrap(BMG_iPARMS,3)
         RETURN
         !
      ENDIF

      !
      ! Sanity check: CG type
      !
      IF ( 
     &   BMG_iPARMS(id_BMG3_SER_CG_TYPE).NE.BMG_SER_CG_ITLI_IzIyIx
     &   .AND.
     &   BMG_iPARMS(id_BMG3_SER_CG_TYPE).NE.BMG_SER_CG_ITLI 
     &   ) THEN
         !
         IF (BMG_IOFLAG(iBMG3_SER_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The Coarse-Grid Operator Type is invalid!'
            WRITE(*,520) 
     &           'HAVE: BMG_iPARMS(id_BMG3_SER_CG_TYPE) =',
     &           BMG_iPARMS(id_BMG3_SER_CG_TYPE)
         END IF
         
         CALL BMG3_SER_SymStd_ErrTrap(BMG_iPARMS,24)
         RETURN
         !
      ENDIF

      !
      ! Sanity check: CG construction method
      !
      IF ( 
     &     BMG_iPARMS(id_BMG3_SER_CG_CONSTRUCT)
     &        .NE.BMG_SER_CG_CONS_explicit 
     & 
     &   .AND.
     & 
     &     BMG_iPARMS(id_BMG3_SER_CG_CONSTRUCT)
     &        .NE.BMG_SER_CG_CONS_block
     &
     &   ) THEN
         !
         IF (BMG_IOFLAG(iBMG3_SER_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The CG-operator construction is invalid!'
            WRITE(*,520)
     &           'HAVE: BMG_iPARMS(id_BMG3_SER_CG_CONSTRUCT) = ',
     &           BMG_iPARMS(id_BMG3_SER_CG_CONSTRUCT)
         END IF
         
         CALL BMG3_SER_SymStd_ErrTrap(BMG_iPARMS,25)
         RETURN
         !
      ENDIF

      !
      ! Fine Grid stencil size
      !
      IF ( IFD.NE.BMG_SER_STENCIL_7pt ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF

      !
      ! Fine Grid dimensions (Nx,Ny,Nz)
      ! 
      CALL BMG3_SER_SymStd_GET_pointers(
     &          NOG, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &          )

      !
      ! Interactively examine the fine grid stencil
      !
      IF ( BMG_IOFLAG(iBMG3_SER_BUG_STENCIL_FG) ) THEN
         CALL BMG3_SER_SymStd_DUMP_stencil( 
     &             BMG_IOFLAG, SO(p_SO), Nx, Ny, Nz, 
     &             IFD, NStncl, NOG, NOG
     &             )
      ENDIF

      !
      ! Zero workspace
      ! 
      DO i=NStncl*Nx*Ny*Nz+1, NSO
         SO(i)=rZERO
      ENDDO

      DO i=1, NSOR
         SOR(i)=rZERO
      ENDDO

      DO i=1, NCI
         CI(i)=rZERO
      ENDDO
      
      DO i=1, NABD1*NABD2
         ABD(i)=rZERO
      ENDDO
      
      DO i=1, NABD2
         BBD(i) = rZERO
      ENDDO

      !
      !   Zero out temporary workspace for operator construction
      !
      DO i=NBMG_iWORK_PL_ptrs+1, NBMG_iWORK_PL
         BMG_iWORK_PL(i) = iZERO
      ENDDO
      
      DO i=1, NBMG_rWORK_PL
         BMG_rWORK_PL(i) = rZERO
      ENDDO

C ------------------------------------------
C     Setup operators on all coarse grids.
C ------------------------------------------

      IF ( NOG.GT.1 ) THEN  ! NOG=1 => there isn't a coarse grid


         IF ( 
     &      BMG_iPARMS(id_BMG3_SER_CG_TYPE).EQ.BMG_SER_CG_ITLI_IzIyIx
     &      ) THEN

            !
            ! This approximation is only formed explicitly, so 
            ! we are ignoring the BMG_iPARMS(id_BMG3_SER_CG_CONSTRUCT).
            ! A warning or flag should be set earlier.
            !

            !
            ! Loop over grids
            !
            DO K = NOG, 2, -1

               !
               ! Determine the number of points in the stencil
               !
               IF ( K.NE.NOG .OR. IFD.NE.BMG_SER_STENCIL_7pt ) THEN
                  NStncl=14
               ELSE
                  NStncl=4
               ENDIF

               CALL BMG3_SER_SymStd_GET_pointers(
     &                   K, IGRD, NOGm,
     &                   p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &                   )
               CALL BMG3_SER_SymStd_GET_pointers( 
     &                   K-1, IGRD, NOGm,
     &                   p_UC, p_SOC, p_SORC, p_CIC, Nxc, Nyc, Nzc 
     &                   )

               !
               ! Temporary interpolation based on tensor product 
               ! approach is used to compute cg.
               !
               p_CGTEMP_zo=BMG_iWORK_PL(ip_BMG_SER_CGTEMP_zo)
               p_CGTEMP_yo=BMG_iWORK_PL(ip_BMG_SER_CGTEMP_yo)
               CALL BMG3_SER_SymStd_SETUP_ITLI_Izyx( 
     &                   K, K-1,
     &                   SO(p_SO), SO(p_SOC), SOR(p_SOR), CI(p_CIC),
     &                   Nx, Ny, Nz, Nxc, Nyc, Nzc, NOG,
     &                   ifd, NStncl, irelax,
     &                   BMG_rWORK_PL(p_CGTEMP_zo),
     &                   BMG_rWORK_PL(p_CGTEMP_yo)
     &                   ) 

               !
               ! Interactively examine the coarse grid stencil
               !
               IF ( BMG_IOFLAG(iBMG3_SER_BUG_STENCIL_CG) ) THEN
                  CALL BMG3_SER_SymStd_DUMP_stencil( 
     &                      BMG_IOFLAG, SO(p_SOC), Nxc, Nyc, Nzc, 
     &                      IFD, NStncl_CG, K-1, NOG
     &                      )
               ENDIF


               !
               ! Operator Induced Interpolation
               !
               p_CGTEMP_yo=BMG_iWORK_PL(ip_BMG_SER_CGTEMP_yo)
               CALL BMG3_SER_SymStd_SETUP_interp_OI(
     &                   K, K-1, so(p_SO), so(p_SOC), ci(p_CIC), 
     &                   Nx, Ny, Nz, Nxc, Nyc, Nzc, NOG, 
     &                   ifd, NStncl, irelax,
     &                   BMG_rWORK_PL(p_CGTEMP_yo)
     &                   )

               !
               ! Interactively examine the restriction operator
               !
               IF ( BMG_IOFLAG(iBMG3_SER_BUG_RESTRICT) ) THEN
                  CALL BMG3_SER_SymStd_DUMP_restrict( 
     &                      BMG_IOFLAG, CI(p_CIC), K, K-1, 
     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
     &                      )
               ENDIF

            ENDDO

         ELSE IF ( 
     &             BMG_iPARMS(id_BMG3_SER_CG_TYPE)
     &                .EQ.BMG_SER_CG_ITLI 
     & 
     &             .AND. 
     &             
     &             BMG_iPARMS(id_BMG3_SER_CG_CONSTRUCT)
     &                .EQ.BMG_SER_CG_CONS_explicit 
     &
     &           ) THEN

            ! ---------------------------------------------------------
            !  Explicit construction of the "honest" I^{T}LI operator.
            ! ---------------------------------------------------------

            !
            ! Loop over grids
            !
            DO K = NOG, 2, -1

               !
               ! Determine the number of points in the stencil
               !
               IF ( K.NE.NOG .OR. IFD.NE.BMG_SER_STENCIL_7pt ) THEN
                  NStncl=14
               ELSE
                  NStncl=4
               ENDIF

               CALL BMG3_SER_SymStd_GET_pointers(
     &                   K, IGRD, NOGm,
     &                   p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &                   )
               CALL BMG3_SER_SymStd_GET_pointers( 
     &                   K-1, IGRD, NOGm,
     &                   p_UC, p_SOC, p_SORC, p_CIC, Nxc, Nyc, Nzc 
     &                   )

               !
               ! Operator Induced Interpolation
               !
               p_CGTEMP_yo=BMG_iWORK_PL(ip_BMG_SER_CGTEMP_yo)
               CALL BMG3_SER_SymStd_SETUP_interp_OI(
     &                   K, K-1, so(p_SO), so(p_SOC), ci(p_CIC), 
     &                   Nx, Ny, Nz, Nxc, Nyc, Nzc, NOG, 
     &                   ifd, NStncl, irelax,
     &                   BMG_rWORK_PL(p_CGTEMP_yo)
     &                   )

               !
               ! Interactively examine the restriction operator
               !
               IF ( BMG_IOFLAG(iBMG3_SER_BUG_RESTRICT) ) THEN
                  CALL BMG3_SER_SymStd_DUMP_restrict( 
     &                      BMG_IOFLAG, CI(p_CIC), K, K-1, 
     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
     &                      )
               ENDIF


               IF ( IFD.NE.1 .OR. k.LT.NOG ) then
                  CALL BMG3_SER_SymStd_SETUP_ITLI27_ex(
     &                      SO(p_SO), SO(p_SOC), CI(p_CIC), 
     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
     &                      )
               ELSE
                  CALL BMG3_SER_SymStd_SETUP_ITLI07_ex(
     &                      SO(p_SO), SO(p_SOC), CI(p_CIC), 
     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
     &                      )
               ENDIF

               !
               ! Interactively examine the coarse grid stencil
               !
               IF ( BMG_IOFLAG(iBMG3_SER_BUG_STENCIL_CG) ) THEN
                  CALL BMG3_SER_SymStd_DUMP_stencil( 
     &                      BMG_IOFLAG, SO(p_SOC), Nxc, Nyc, Nzc, 
     &                      IFD, NStncl_CG, K-1, NOG
     &                      )
               ENDIF

            ENDDO

         ELSE IF ( 
     &             BMG_iPARMS(id_BMG3_SER_CG_TYPE)
     &                .EQ.BMG_SER_CG_ITLI 
     & 
     &             .AND. 
     &             
     &             BMG_iPARMS(id_BMG3_SER_CG_CONSTRUCT)
     &                .EQ.BMG_SER_CG_CONS_block
     &
     &           ) THEN

            ! ---------------------------------------------------------
            !  Block construction of the "honest" I^{T}LI operator.
            ! ---------------------------------------------------------

            !
            ! Loop over grids
            !
            DO K = NOG, 2, -1

               !
               ! Determine the number of points in the stencil
               !
               IF ( K.NE.NOG .OR. IFD.NE.BMG_SER_STENCIL_7pt ) THEN
                  NStncl=14
               ELSE
                  NStncl=4
               ENDIF

               CALL BMG3_SER_SymStd_GET_pointers(
     &                   K, IGRD, NOGm,
     &                   p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &                   )
               CALL BMG3_SER_SymStd_GET_pointers( 
     &                   K-1, IGRD, NOGm,
     &                   p_UC, p_SOC, p_SORC, p_CIC, Nxc, Nyc, Nzc 
     &                   )

               !
               ! Operator Induced Interpolation
               !
               p_CGTEMP_yo=BMG_iWORK_PL(ip_BMG_SER_CGTEMP_yo)
               CALL BMG3_SER_SymStd_SETUP_interp_OI(
     &                   K, K-1, so(p_SO), so(p_SOC), ci(p_CIC), 
     &                   Nx, Ny, Nz, Nxc, Nyc, Nzc, NOG, 
     &                   ifd, NStncl, irelax,
     &                   BMG_rWORK_PL(p_CGTEMP_yo)
     &                   )

               !
               ! Interactively examine the restriction operator
               !
               IF ( BMG_IOFLAG(iBMG3_SER_BUG_RESTRICT) ) THEN
                  CALL BMG3_SER_SymStd_DUMP_restrict( 
     &                      BMG_IOFLAG, CI(p_CIC), K, K-1, 
     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
     &                      )
               ENDIF


               IF ( IFD.NE.1 .OR. k.LT.NOG ) then
                  CALL BMG3_SER_SymStd_SETUP_ITLI27_bl(
     &                      SO(p_SO), SO(p_SOC), CI(p_CIC), 
     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
     &                      )
               ELSE
                  CALL BMG3_SER_SymStd_SETUP_ITLI07_bl(
     &                      SO(p_SO), SO(p_SOC), CI(p_CIC), 
     &                      Nx, Ny, Nz, Nxc, Nyc, Nzc
     &                      )
               ENDIF

               !
               ! Interactively examine the coarse grid stencil
               !
               IF ( BMG_IOFLAG(iBMG3_SER_BUG_STENCIL_CG) ) THEN
                  CALL BMG3_SER_SymStd_DUMP_stencil( 
     &                      BMG_IOFLAG, SO(p_SOC), Nxc, Nyc, Nzc, 
     &                      IFD, NStncl_CG, K-1, NOG
     &                      )
               ENDIF

            ENDDO

         ELSE IF ( 
     &      BMG_iPARMS(id_BMG3_SER_CG_TYPE).EQ.BMG_SER_CG_USER
     &      ) THEN

            
            !
            !  Set error trap: USER option not supported yet.
            !

         ELSE
            
            !
            !  Set error trap: Unsupported CG (TYPE,CONSTRUCTION)
            !

         ENDIF


         !
         !  Re-zero workspace
         !
         DO i=1, NBMG_rWORK_PL
            BMG_rWORK_PL(i) = rZERO
         ENDDO

         !
         !  Relaxation
         ! 
         CALL BMG3_SER_SymStd_SETUP_relax( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             IFD, IBC, IRELAX, NOGm, NOG, IGRD, 
     &             SO, NSO, SOR, NSOR, CI, NCI,
     &             BMG_iWORK_PL, NBMG_iWORK_PL,
     &             BMG_rWORK_PL, NBMG_rWORK_PL
     &             )
         IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF
         

      ENDIF

C ------------------------------------------
C     Setup the coarse grid solve
C ------------------------------------------

      IF ( NOG.EQ. 1 ) THEN 
         !
         ! Determine the number of points in the stencil
         !
         IF ( IFD.NE.BMG_SER_STENCIL_7pt ) THEN
            NStncl=14
         ELSE
            NStncl=4
         ENDIF
      ELSE
         !
         ! There has been at least one coarsening.
         !
         NStncl=14
         !
      ENDIF
      !
      ! Obtain pointers for the stencil on the coarsest grid
      ! 
      K=1          ! Coarsest-grid index is always 1
      CALL BMG3_SER_SymStd_GET_pointers(
     &          K, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &          )

      !
      ! Interactively examine the stencil on the coarsest grid
      !
      IF ( BMG_IOFLAG(iBMG3_SER_BUG_STENCIL_CG1) ) THEN
         CALL BMG3_SER_SymStd_DUMP_stencil( 
     &             BMG_IOFLAG, SO(p_SO), Nx, Ny, Nz, 
     &             IFD, NStncl, K, NOG
     &             )
      ENDIF

      !
      ! Setup and factor the matrix for direct solve on the coarsest grid.
      !
      CALL BMG3_SER_SymStd_SETUP_cg_LU( 
     &          SO(p_SO), Nx, Ny, Nz, NStncl, abd, nabd1, nabd2,
     &          BMG_IOFLAG, BMG_iPARMS 
     &          )
      IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
         RETURN
      END IF
      

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SER_SymStd_SETUP_parts.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
