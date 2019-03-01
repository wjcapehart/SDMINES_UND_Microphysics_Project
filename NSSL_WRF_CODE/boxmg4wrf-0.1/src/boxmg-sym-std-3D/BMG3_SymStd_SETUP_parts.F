      SUBROUTINE BMG3_SymStd_SETUP_parts( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                SO, NSO, SOR, NSOR, CI, NCI, IGRD, NOGm, NOG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWork, NMSGi, pSI_MSG, pMSG, pMSGSO,
     &                BUFFER, NMSGr,
     &                MyProcI, MyProcJ, MyProcK, NProcI, NProcJ, NProcK, 
     &                NProc, ProcGrid, ProcCoord, LocArrSize, 
     &                DimX, DimY, DimZ, MPICOMM 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_SETUP_parts.f loops over all grids calling the
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
      INCLUDE 'mpif.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ---------------------------
C    Argument Declarations:
C
      INTEGER  NBMG_iWORK_CS, NBMG_iWORK_PL,
     &         NBMG_rWORK_CS, NBMG_rWORK_PL, 
     &         NCI, NOGm, NSO, NSOR, NOG,
     &         NMSGi, NMSGr, pSI_MSG

      INTEGER  iWork(NMSGi), pMSG(NBMG_pMSG,NOG), pMSGSO(NBMG_pMSG,NOG) 

      INTEGER  MPICOMM, NProc
      INTEGER  MyProcI, MyProcJ, MyProcK, NProcI, NProcJ, NProcK
      INTEGER  ProcGrid(NProcI,NProcJ,NProcK), ProcCoord(2,NProc)
      INTEGER  LocArrSize(3,*)

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_iWORK_CS(NBMG_iWORK_CS),
     &         BMG_iWORK_PL(NBMG_iWORK_PL), IGRD(NOG,NBMG_pIGRD)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS), BMG_rWORK_CS(NBMG_rWORK_CS),
     &         BMG_rWORK_PL(NBMG_rWORK_PL), CI(NCI), SO(NSO),
     &         SOR(NSOR), BUFFER(NMSGr)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  DimX(NProcI,NOGm), DimY(NProcJ,NOGm), DimZ(NPRocK,NOGm)

C ---------------------------
C    Local Declarations:
C
      INTEGER  iGs, jGs, kGs, MPI_MyProc, MyProc, MPI_IERROR
    

      INTEGER  i, IFD, IBC, IRELAX, K, NCBW, NCU, NStncl, NStncl_CG,
     &         NLx, NLy, NLz, NLx_c, NLy_c, NLz_c, 
     &         NGx, NGy, NGz, NGx_c, NGy_c, NGz_c,
     &         p_ABD, p_BBD, p_SOR1, p_SOR2, p_CI, p_CIC,
     &         p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC

      REAL*8   T1, T2

C ==========================================================================

      !
      ! Unpack
      !
      IFD    = BMG_iPARMS(id_BMG3_STENCIL)
      IBC    = BMG_iPARMS(id_BMG3_BC)
      IRELAX = BMG_iPARMS(id_BMG3_RELAX)

      MyProc = ProcGrid(MyProcI,MyProcJ,MyProcK)
      MPI_MyProc = MyProc - 1

      !
      ! Sanity check
      !

      IF ( NOG.EQ.0 ) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,500) 'There are no grids?'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,1)
         RETURN

      ELSE IF ( NOG.LT.0 ) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,500) 'The number of grids is negative!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,2)
         RETURN

      ENDIF

      !
      ! Sanity check: CG type
      !
      IF ( 
     &   BMG_iPARMS(id_BMG3_CG_TYPE).NE.BMG_CG_ITLI_IzIyIx
     &   .AND.
     &   BMG_iPARMS(id_BMG3_CG_TYPE).NE.BMG_CG_ITLI 
     &   ) THEN
         !
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The Coarse-Grid Operator Type is invalid!'
            WRITE(*,520) 
     &           'HAVE: BMG_iPARMS(id_BMG3_CG_TYPE) =',
     &           BMG_iPARMS(id_BMG3_CG_TYPE)
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,24)
         RETURN
         !
      ENDIF

      !
      ! Sanity check: CG construction method
      !
      IF ( 
     &     BMG_iPARMS(id_BMG3_CG_CONSTRUCT)
     &        .NE.BMG_CG_CONS_explicit 
     & 
     &   .AND.
     & 
     &     BMG_iPARMS(id_BMG3_CG_CONSTRUCT)
     &        .NE.BMG_CG_CONS_block
     &
     &   ) THEN
         !
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'The CG-operator construction is invalid!'
            WRITE(*,520)
     &           'HAVE: BMG_iPARMS(id_BMG3_CG_CONSTRUCT) = ',
     &           BMG_iPARMS(id_BMG3_CG_CONSTRUCT)
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,25)
         RETURN
         !
      ENDIF

      !
      ! Fine Grid stencil size
      !
      IF ( IFD.NE.BMG_STENCIL_7pt ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF

      !
      ! Coarse-grid operators are always 27 points.
      !
      NStncl_CG=14

      !
      ! Fine Grid dimensions (Nx,Ny,Nz)
      ! 
      CALL BMG3_SymStd_GET_pointers( 
     &          NOG, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CI, 
     &          NLx, NLy, NLz, NGx, NGy, NGz
     &          )


      T1 = MPI_Wtime()

      CALL BMG3_SymStd_SETUP_fine_stencil(
     &          NOG, SO(p_SO), 
     &          NLx, NLy, NLz, NStncl,
     &          iWork, NMSGi, pMSGSO, BUFFER, NMSGr, 
     &          MPICOMM
     &          )      

      T2 = MPI_Wtime()
      BMG_rPARMS(id_BMG3_TIME_SETUP_FINE_STENCIL) = 
     &  BMG_rPARMS(id_BMG3_TIME_SETUP_FINE_STENCIL) + T2 - T1    

      !
      ! Interactively examine the fine grid stencil
      !
      IF ( BMG_IOFLAG(iBMG3_BUG_STENCIL_FG) ) THEN
         CALL BMG3_SymStd_DUMP_stencil( 
     &                    BMG_IOFLAG, SO(p_SO), NStncl, NOG, NOG,
     &                    NLx, NLy, NLz, NGx, NGy, NGz,
     &                    IGRD(NOG,idL_BMG_ICoord), 
     &                    IGRD(NOG,idL_BMG_JCoord),
     &                    IGRD(NOG,idL_BMG_KCoord),
     &                    ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &                    MyProcI, MyProcJ, MyProcK, MPICOMM
     &                    )
      ENDIF

      !
      ! Zero workspace
      ! 
      p_SOC=p_SO+NStncl*(NLx+1)*(NLy+1)*(NLz+1)
      DO i=p_SOC, NSO
         SO(i)=rZERO
      ENDDO

      DO i=1, NSOR
         SOR(i)=rZERO
      ENDDO

      DO i=1, NCI
         CI(i)=rZERO
      ENDDO
      
      !
      DO i=NBMG_iWORK_PL_ptrs+1, NBMG_iWORK_PL
         BMG_iWORK_PL(i) = iZERO
      ENDDO
      
      DO i=1, NBMG_rWORK_PL
         BMG_rWORK_PL(i) = rZERO
      ENDDO

      !
      DO i=NBMG_iWORK_CS_ptrs+1, NBMG_iWORK_CS
         BMG_iWORK_CS(i) = iZERO
      ENDDO

      DO i=1, NBMG_rWORK_CS
         BMG_rWORK_CS(i) = rZERO
      ENDDO

C ------------------------------------------
C     Setup operators on all coarse grids.
C ------------------------------------------

      IF ( NOG.GT.1 ) THEN  ! NOG=1 => there isn't a coarse grid

         IF ( 
     &      BMG_iPARMS(id_BMG3_CG_TYPE).EQ.BMG_CG_ITLI_IzIyIx
     &      ) THEN
            ! WRITE(*,*) 'Coarse grid construction uses successive'
            ! WRITE(*,*) 'semi-coarsening.'
            !
            ! This approximation is only formed explicitly, so 
            ! we are ignoring the BMG_iPARMS(id_BMG3_CG_CONSTRUCT).
            ! A warning or flag should be set earlier.
            !
         !
         ! Loop over grids
         !
         DO K = NOG, 2, -1

            !
            ! Determine the number of points in the stencil
            !
            IF ( K.NE.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
               NStncl=14
            ELSE
               NStncl=4
            ENDIF

            CALL BMG3_SymStd_GET_pointers( 
     &                 K, IGRD, NOGm,
     &                 p_U, p_SO, p_SOR, p_CI, 
     &                 NLx, NLy, NLz, NGx, NGy, NGz
     &                 )
            iGs = IGRD(K,idL_BMG_ICoord)
            jGs = IGRD(K,idL_BMG_JCoord)
            kGs = IGRD(K,idL_BMG_KCoord)


            CALL BMG3_SymStd_GET_pointers( 
     &                 K-1, IGRD, NOGm,
     &                 p_UC, p_SOC, p_SORC, p_CIC, 
     &                 NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c
     &                 )


            p_SOR1 = p_SOR
            p_SOR2 = p_SOR + NLx*NLy*NLz
            !
            ! Temporary interpolation based on tensor product 
            ! approach is used to compute cg.
            !

            T1 = MPI_Wtime()
            CALL BMG3_SymStd_SETUP_ITLI_Izyx( 
     &                K, K-1, 
     &                SO(p_SO), SO(p_SOC), SOR(p_SOR1), SOR(p_SOR2), 
     &                CI(p_CIC), NLx, NLy, NLz, NLx_c, NLy_c, NLz_c, 
     &                NOG, IFD, NStncl, IRELAX,
     &                BMG_rWORK_PL(BMG_iWORK_PL(ip_BMG_CGTEMP_zo)),
     &                BMG_rWORK_PL(BMG_iWORK_PL(ip_BMG_CGTEMP_yo)),
     &                NOGm, IGRD, iWork, NMSGi, pSI_MSG, pMSG, pMSGSO,
     &                BUFFER, NMSGr, NProc, ProcGrid, 
     &                MyProcI, MyProcJ, MyProcK, NProcI, NProcJ, NProcK, 
     &                BMG_IOFLAG, DimX, DimY, DimZ, MPICOMM
     &                )
            T2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_SETUP_CG_ITLI) = 
     &           BMG_rPARMS(id_BMG3_TIME_SETUP_CG_ITLI) + T2 - T1

c$$$            CALL MPI_FINALIZE(MPI_IERROR)
c$$$            STOP
c$$$
               !           
               ! Operator Induced Interpolation
               !
            T1 = MPI_Wtime()   
            CALL BMG3_SymStd_SETUP_interp_OI(
     &                  K, K-1, so(p_SO), so(p_SOC), ci(p_CIC), 
     &                  NLx, NLy, NLz, NLx_c, NLy_c, NLz_c, 
     &                  NOG, IFD, NStncl, IRELAX,
     &                  BMG_rWORK_PL(BMG_iWORK_PL(ip_BMG_CGTEMP_yo)),
     &                  NOGm, IGRD, iWork, NMSGi, pMSG,
     &                  BUFFER, NMSGr, MyProc, MPICOMM
     &                  )
            T2 = MPI_Wtime()
            
            BMG_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI) = 
     &           BMG_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI) + T2 - T1


            !
            ! Interactively examine the coarse grid stencil
            !
            IF ( BMG_IOFLAG(iBMG3_BUG_STENCIL_CG) ) THEN
               !
               CALL BMG3_SymStd_DUMP_stencil( 
     &              BMG_IOFLAG, SO(p_SOC), NStncl_CG, K-1, NOG,
     &              NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c,
     &              IGRD(K-1,idL_BMG_ICoord), 
     &              IGRD(K-1,idL_BMG_JCoord),
     &              IGRD(K-1,idL_BMG_KCoord),
     &              ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &              MyProcI, MyProcJ, MyProcK, MPICOMM
     &              )
               !
            ENDIF

            !
            ! Interactively examine the restriction
            !
            IF ( BMG_IOFLAG(iBMG3_BUG_RESTRICT) ) THEN
               !
               CALL BMG3_SymStd_DUMP_restrict( 
     &              BMG_IOFLAG, CI(p_CIC), K, K-1, NOG,
     &              NLx, NLy, NLz, NGx, NGy, NGz,
     &              NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c,
     &              IGRD(K,idL_BMG_ICoord), 
     &              IGRD(K,idL_BMG_JCoord),
     &              IGRD(K,idL_BMG_KCoord),
     &              IGRD(K-1,idL_BMG_ICoord), 
     &              IGRD(K-1,idL_BMG_JCoord),
     &              IGRD(K-1,idL_BMG_KCoord),
     &              ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &              MyProcI, MyProcJ, MyProcK, MPICOMM
     &              )
               !
            ENDIF

         ENDDO

         ELSE IF ( 
     &             BMG_iPARMS(id_BMG3_CG_TYPE)
     &                .EQ.BMG_CG_ITLI 
     & 
     &             .AND. 
     &             
     &             BMG_iPARMS(id_BMG3_CG_CONSTRUCT)
     &                .EQ.BMG_CG_CONS_explicit 
     &
     &           ) THEN
             ! WRITE(*,*) 'Coarse grid construction uses standard'
             ! WRITE(*,*) 'coarsening and explicit construction of'
             ! WRITE(*,*) 'I{^T}LI.'
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
            IF ( K.NE.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
               NStncl=14
            ELSE
               NStncl=4
            ENDIF

            CALL BMG3_SymStd_GET_pointers( 
     &                 K, IGRD, NOGm,
     &                 p_U, p_SO, p_SOR, p_CI, 
     &                 NLx, NLy, NLz, NGx, NGy, NGz
     &                 )
            iGs = IGRD(K,idL_BMG_ICoord)
            jGs = IGRD(K,idL_BMG_JCoord)
            kGs = IGRD(K,idL_BMG_KCoord)


            CALL BMG3_SymStd_GET_pointers( 
     &                 K-1, IGRD, NOGm,
     &                 p_UC, p_SOC, p_SORC, p_CIC, 
     &                 NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c
     &                 )


            p_SOR1 = p_SOR
            p_SOR2 = p_SOR + NLx*NLy*NLz
            !
            ! Temporary interpolation based on tensor product 
            ! approach is used to compute cg.
            !
            T1 = MPI_Wtime()

c$$$            CALL MPI_FINALIZE(MPI_IERROR)
c$$$            STOP
c$$$
               !           
               ! Operator Induced Interpolation
               !
            T1 = MPI_Wtime()   
            CALL BMG3_SymStd_SETUP_interp_OI(
     &                  K, K-1, so(p_SO), so(p_SOC), ci(p_CIC), 
     &                  NLx, NLy, NLz, NLx_c, NLy_c, NLz_c, 
     &                  NOG, IFD, NStncl, IRELAX,
     &                  BMG_rWORK_PL(BMG_iWORK_PL(ip_BMG_CGTEMP_yo)),
     &                  NOGm, IGRD, iWork, NMSGi, pMSG,
     &                  BUFFER, NMSGr, MyProc, MPICOMM
     &                  )
            T2 = MPI_Wtime()
            
            BMG_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI) = 
     &           BMG_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI) + T2 - T1


            T2 = MPI_Wtime()
            IF (IFD.NE.1 .OR. K.LT.NOG) THEN            

               CALL BMG3_SymStd_SETUP_ITLI27_ex(
     &              K, K-1,SO(p_SO), SO(p_SOC), CI(p_CIC),
     &              NLx, NLy, NLz,NLx_c, NLy_C, NLz_c,
     &              iGs, jGs, kGS, NOG, NOGm, IGRD,
     &              iWORK, NMSGi, pSI_MSG, pMSGSO,
     &              BUFFER, NMSGr, NProc, ProcGrid, MyProcI, MyProcJ,
     &              MyProcK, NProcI, NProcJ, NProcK, DimX, DimY, DimZ,
     &              MPICOMM
     &              )
            ELSE

               CALL BMG3_SymStd_SETUP_ITLI07_ex(
     &              K, K-1, SO(p_SO), SO(p_SOC), CI(p_CIC),
     &              NLx, NLy, NLz, NLx_c, NLy_c, NLz_c,
     &              iGs, jGs, kGS, NOG, NOGm, IGRD,
     &              iWORK, NMSGi, pMSGSO,
     &              BUFFER, NMSGr, NProc, ProcGrid, MyProcI, MyProcJ,
     &              MyProcK, NProcI, NProcJ, NProcK, DimX, DimY, DimZ,
     &              MPICOMM
     &              )

            ENDIF 
  
            T2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_SETUP_CG_ITLI) = 
     &           BMG_rPARMS(id_BMG3_TIME_SETUP_CG_ITLI) + T2 - T1  

            !
            ! Interactively examine the coarse grid stencil
            !
            IF ( BMG_IOFLAG(iBMG3_BUG_STENCIL_CG) ) THEN
               !
               CALL BMG3_SymStd_DUMP_stencil( 
     &              BMG_IOFLAG, SO(p_SOC), NStncl_CG, K-1, NOG,
     &              NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c,
     &              IGRD(K-1,idL_BMG_ICoord), 
     &              IGRD(K-1,idL_BMG_JCoord),
     &              IGRD(K-1,idL_BMG_KCoord),
     &              ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &              MyProcI, MyProcJ, MyProcK, MPICOMM
     &              )
               !
            ENDIF

            !
            ! Operator Induced Interpolation
            !
            
            BMG_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI) = 
     &           BMG_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI) + T2 - T1

            !
            ! Interactively examine the restriction
            !
            IF ( BMG_IOFLAG(iBMG3_BUG_RESTRICT) ) THEN
               !
               CALL BMG3_SymStd_DUMP_restrict( 
     &              BMG_IOFLAG, CI(p_CIC), K, K-1, NOG,
     &              NLx, NLy, NLz, NGx, NGy, NGz,
     &              NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c,
     &              IGRD(K,idL_BMG_ICoord), 
     &              IGRD(K,idL_BMG_JCoord),
     &              IGRD(K,idL_BMG_KCoord),
     &              IGRD(K-1,idL_BMG_ICoord), 
     &              IGRD(K-1,idL_BMG_JCoord),
     &              IGRD(K-1,idL_BMG_KCoord),
     &              ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &              MyProcI, MyProcJ, MyProcK, MPICOMM
     &              )
               !
            ENDIF

         ENDDO

         ELSE IF ( 
     &             BMG_iPARMS(id_BMG3_CG_TYPE)
     &                .EQ.BMG_CG_ITLI 
     & 
     &             .AND. 
     &             
     &             BMG_iPARMS(id_BMG3_CG_CONSTRUCT)
     &                .EQ.BMG_CG_CONS_block
     &
     &           ) THEN
            ! WRITE(*,*) 'Coarse grid construction uses standard'
            ! WRITE(*,*) 'coarsening and block construction of'
            ! WRITE(*,*) 'I^{T}LI.'

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
            IF ( K.NE.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
               NStncl=14
            ELSE
               NStncl=4
            ENDIF

            CALL BMG3_SymStd_GET_pointers( 
     &                 K, IGRD, NOGm,
     &                 p_U, p_SO, p_SOR, p_CI, 
     &                 NLx, NLy, NLz, NGx, NGy, NGz
     &                 )
            iGs = IGRD(K,idL_BMG_ICoord)
            jGs = IGRD(K,idL_BMG_JCoord)
            kGs = IGRD(K,idL_BMG_KCoord)


            CALL BMG3_SymStd_GET_pointers( 
     &                 K-1, IGRD, NOGm,
     &                 p_UC, p_SOC, p_SORC, p_CIC, 
     &                 NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c
     &                 )

c$$$            CALL MPI_FINALIZE(MPI_IERROR)
c$$$            STOP
c$$$
               !           
               ! Operator Induced Interpolation
               !
            T1 = MPI_Wtime()   
            CALL BMG3_SymStd_SETUP_interp_OI(
     &                  K, K-1, so(p_SO), so(p_SOC), ci(p_CIC), 
     &                  NLx, NLy, NLz, NLx_c, NLy_c, NLz_c, 
     &                  NOG, IFD, NStncl, IRELAX,
     &                  BMG_rWORK_PL(BMG_iWORK_PL(ip_BMG_CGTEMP_yo)),
     &                  NOGm, IGRD, iWork, NMSGi, pMSG,
     &                  BUFFER, NMSGr, MyProc, MPICOMM
     &                  )
            T2 = MPI_Wtime()
            
            BMG_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI) = 
     &           BMG_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI) + T2 - T1

            T1 = MPI_Wtime()
                CALL BMG3_SymStd_SETUP_ITLI_bl(
     &                K, K-1, SO(p_SO), SO(p_SOC), CI(p_CIC),
     &                NLx, NLy, NLz, NLx_c, NLy_c, NLz_c,
     &                iGs, jGs, kGs, NOG, IFD, NStncl, NOGm, IGRD,
     &                iWORK, NMSGi, pSI_MSG, pMSGSO,
     &                BUFFER, NMSGr, NProc, ProcGrid, MyProcI, MyProcJ,
     &                MyProcK, NProcI, NprocJ, NProcK, DimX, DimY, DimZ,
     &                MPICOMM
     &                )
  
            T2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_SETUP_CG_ITLI) = 
     &           BMG_rPARMS(id_BMG3_TIME_SETUP_CG_ITLI) + T2 - T1  

            !
            ! Interactively examine the coarse grid stencil
            !
            IF ( BMG_IOFLAG(iBMG3_BUG_STENCIL_CG) ) THEN
               !
               CALL BMG3_SymStd_DUMP_stencil( 
     &              BMG_IOFLAG, SO(p_SOC), NStncl_CG, K-1, NOG,
     &              NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c,
     &              IGRD(K-1,idL_BMG_ICoord), 
     &              IGRD(K-1,idL_BMG_JCoord),
     &              IGRD(K-1,idL_BMG_KCoord),
     &              ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &              MyProcI, MyProcJ, MyProcK, MPICOMM
     &              )
               !
            ENDIF
            
            BMG_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI) = 
     &           BMG_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI) + T2 - T1

            !
            ! Interactively examine the restriction
            !
            IF ( BMG_IOFLAG(iBMG3_BUG_RESTRICT) ) THEN
               !
               CALL BMG3_SymStd_DUMP_restrict( 
     &              BMG_IOFLAG, CI(p_CIC), K, K-1, NOG,
     &              NLx, NLy, NLz, NGx, NGy, NGz,
     &              NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c,
     &              IGRD(K,idL_BMG_ICoord), 
     &              IGRD(K,idL_BMG_JCoord),
     &              IGRD(K,idL_BMG_KCoord),
     &              IGRD(K-1,idL_BMG_ICoord), 
     &              IGRD(K-1,idL_BMG_JCoord),
     &              IGRD(K-1,idL_BMG_KCoord),
     &              ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &              MyProcI, MyProcJ, MyProcK, MPICOMM
     &              )
               !
            ENDIF

         ENDDO
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
         T1 = MPI_Wtime()

         CALL BMG3_SymStd_SETUP_relax( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             IFD, IBC, IRELAX, NOGm, NOG, IGRD, 
     &             SO, NSO, SOR, NSOR, CI, NCI,
     &             BMG_iWORK_PL, NBMG_iWORK_PL,
     &             BMG_rWORK_PL, NBMG_rWORK_PL,
     &             iWork, NMSGi, pMSG, BUFFER, NMSGr,
     &             NProc, MyProc, MyProcI, MyProcJ, MyProcK,
     &             MPICOMM
     &             )
         IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
            RETURN
         END IF

         T2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_SETUP_RELAX) =
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_RELAX) + T2 - T1

      ENDIF

C ------------------------------------------
C     Setup the coarse grid solve
C ------------------------------------------

      IF ( NOG.EQ. 1 ) THEN 
         !
         ! Determine the number of points in the stencil
         !
         IF ( IFD.NE.BMG_STENCIL_7pt ) THEN
            NStncl=14
         ELSE
            NStncl=4
         ENDIF
         
      ELSE
         
         NStncl=14

      ENDIF
      !
      ! Obtain pointers for the stencil on the coarsest grid
      ! 
      K=1          ! Coarsest-grid index is always 1
      CALL BMG3_SymStd_GET_pointers( 
     &            K, IGRD, NOGm,
     &            p_U, p_SO, p_SOR, p_CI, 
     &            NLx_c, NLy_c, NLz_c, NGx_c, NGy_c,NGz_c
     &            )

      iGs = IGRD(K,idL_BMG_ICoord)
      jGs = IGRD(K,idL_BMG_JCoord)
      kGs = IGRD(K,idL_BMG_KCoord)

      !
      ! Interactively examine the stencil on the coarsest grid
      !
      IF ( BMG_IOFLAG(iBMG3_BUG_STENCIL_CG1) ) THEN
         CALL BMG3_SymStd_DUMP_stencil( 
     &             BMG_IOFLAG, SO(p_SO), NStncl, K, NOG,
     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c,
     &             iGs, jGs, kGs,
     &             ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &             MyProcI, MyProcJ, MyProcK, MPICOMM
     &             )
      ENDIF


      IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN
         !
         ! Setup and factor the matrix for direct solve on the coarsest grid.
         !
         p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
         NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
         NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)
         !
         T1 = MPI_Wtime()
      
#if RKIND == 1
           write(0,*) 'BMG3_SymStd_SETUP_parts 1, proc = ',myproc,
     :       'call BMG3_SymStd_SETUP_cg_LU'
#endif 
         CALL BMG3_SymStd_SETUP_cg_LU( 
     &             SO(p_SO), NLx_c, NLy_c, NLz_c, NStncl, 
     &             BMG_rWORK_CS(p_ABD), NCBW, NCU, BUFFER, NMSGr, 
     &             NProcI, NProcJ, NProcK, NProc, MyProc,
     &             ProcGrid, ProcCoord, LocArrSize, MPICOMM,
     &             BMG_IOFLAG, BMG_iPARMS
     &             )

#if RKIND == 1
           write(0,*) 'BMG3_SymStd_SETUP_parts 1, proc = ',myproc,
     :       'after call BMG3_SymStd_SETUP_cg_LU'
#endif 
         T2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_SETUP_CG_LU) =
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_CG_LU) + T2 - T1

      ELSE

         T1 = MPI_Wtime()
#if RKIND == 1
           write(0,*) 'BMG3_SymStd_SETUP_parts 1, proc = ',myproc,
     :       'call BMG3_SymStd_SETUP_cg_boxmg'
         write(0,*) 'setup_parts,proc, ProcGrid: ',myproc,
     :    (ProcGrid(i,1,1),i=1,2),ProcGrid(1,2,1),ProcGrid(2,2,1)
#endif 

         CALL BMG3_SymStd_SETUP_cg_boxmg(
     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c,
     &             iGs, jGs, kGs,
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             SO(p_SO), NStncl, NOG,
     &             BMG_iWORK_CS, NBMG_iWORK_CS,
     &             BMG_rWORK_CS, NBMG_rWORK_CS,
     &             BUFFER, NMSGr, 
     &             NProcI, NProcJ, NProcK, NProc, MyProc,
     &             ProcGrid, ProcCoord, LocArrSize, MPICOMM
     &             )
         
#if RKIND == 1
           write(0,*) 'BMG3_SymStd_SETUP_parts 1, proc = ',myproc,
     :       ' after call BMG3_SymStd_SETUP_cg_boxmg'
#endif 
         T2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_SETUP_CG_LU) =
     &        BMG_rPARMS(id_BMG3_TIME_SETUP_CG_LU) + T2 - T1

      ENDIF

      IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
         RETURN
      END IF


C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SETUP_parts.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
