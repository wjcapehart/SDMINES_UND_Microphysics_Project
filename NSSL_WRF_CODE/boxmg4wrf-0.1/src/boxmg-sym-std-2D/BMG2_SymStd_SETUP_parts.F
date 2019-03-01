      SUBROUTINE BMG2_SymStd_SETUP_parts( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                SO, NSO, SOR, NSOR, CI, NCI, IGRD, NOGm, NOG, 
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWork, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                MyProcI, MyProcJ, NProcI, NProcJ, NProc,
     &                ProcGrid, ProcCoord, LocArrSize, 
     &                MPICOMM, XLINECOMM, YLINECOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SETUP_parts loops over all grids calling the routines
C     to construct the operator induced interpolation and the coarse
C     grid operators.  It also sets up the coarse-grid solve (either
C     an LU factorization, or the parameters and workspace to call the
C     serial version of BoxMG).
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

#include      "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C

      !
      !  BoxMG parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      !
      !  BoxMG dimensions
      !
      INTEGER  NF, NC, NSO, NSOR, NCI, NOGm, NOG,
     &         NBMG_iWORK_CS, NBMG_rWORK_CS

      !
      !  BoxMG internal pointers
      !
      INTEGER  IGRD(NOGm,NBMG_pIGRD)

      !
      !  BoxMG data arrays
      !
      REAL*RKIND   CI(NCI), SO(NSO), SOR(NSOR)

      INTEGER  BMG_iWORK_CS(NBMG_iWORK_CS)
      REAL*RKIND   BMG_rWORK_CS(NBMG_rWORK_CS)

      !
      !  MSG workspace dimensions
      !
      INTEGER  NMSGi, NMSGr
      !
      !  MSG workspace and pointers
      !
      INTEGER  iWork(NMSGi), 
     &         pMSG(NBMG_pMSG,NOG),
     &         pMSGSO(NBMG_pMSG,NOG) 
      REAL*RKIND   MSG_Buffer(NMSGr)
      

      !
      !  BMG_MSG grid stuff
      !
      INTEGER  MPICOMM, XLINECOMM, YLINECOMM, NProc
      INTEGER  MyProcI, MyProcJ, NProcI, NProcJ
      INTEGER  ProcGrid(NProcI,NProcJ), ProcCoord(2,NProc)
      INTEGER  LocArrSize(3,*)

C      INTEGER  NBMG_MSG_iGRID
C      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
C     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)
            
C ----------------------------
C     Local Declarations
C
      !
      !  BoxMG Cycle Parameters
      !
      INTEGER  IBC, IFD, IU, ID, IVW, IRELAX, IRELAX_SYM

      INTEGER  iGs, iGs_c, jGs, jGs_c, MPI_MyProc, MyProc, MPI_IERROR
 
      INTEGER  i, IIC, IIF, JJF, JJC, K, NStncl, NStncl_CG,
     &         NCBW, NCU, IIGF, JJGF, IIGC, JJGC, p_ABD,
     &         p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC

      INTEGER  NLx_c, NLy_c, NGx_c, NGy_c

      REAL*8   TT1, TT2

C ==========================================================================

      !
      ! Unpack
      !
      IFD    = BMG_iPARMS(id_BMG2_STENCIL)
      IBC    = BMG_iPARMS(id_BMG2_BC)
      IRELAX = BMG_iPARMS(id_BMG2_RELAX)

      MyProc = ProcGrid(MyProcI,MyProcJ)
      MPI_MyProc = MyProc - 1

      !
      ! Sanity check
      !
      IF ( NOG.EQ.0 ) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,500) 'There are no grids?'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,1)
         RETURN

      ELSE IF ( NOG.LT.0 ) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,500) 'The number of grids is negative!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF
         
         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,2)
         RETURN

      ENDIF

      !
      ! Fine Grid stencil size
      !
      IF ( IFD.EQ.BMG_STENCIL_5pt ) THEN
         NStncl=3
      ELSE
         NStncl=5
      ENDIF

      !
      ! Coarse-grid operators are always 9 points.
      !
      NStncl_CG=5

      !
      ! synchronize the stencil on the finest grid
      ! by update ghost points within processor subdomains
      !
      CALL BMG2_SymStd_GET_pointers(
     &          NOG, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CI, 
     &          IIF, JJF, IIGF, JJGF 
     &          )
      iGs = IGRD(NOG,idL_BMG_ICoord) 
      jGs = IGRD(NOG,idL_BMG_JCoord)
      
      TT1 = MPI_Wtime()

      CALL BMG2_SymStd_SETUP_fine_stencil (
     &          NOG, SO(p_SO), IIF, JJF, NStncl,
     &          iWork, NMSGi, pMSGSO,
     &          MSG_Buffer, NMSGr, MPICOMM
     &          )
 
      TT2 = MPI_Wtime()
      BMG_rPARMS(id_BMG2_TIME_SETUP_fine_stencil) =
     &     BMG_rPARMS(id_BMG2_TIME_SETUP_fine_stencil) + TT2 - TT1

      !
      ! Interactively examine the fine grid stencil
      !
      IF ( BMG_IOFLAG(iBMG2_BUG_STENCIL_FG) ) THEN
         CALL BMG2_SymStd_DUMP_stencil( 
     &                    BMG_IOFLAG, SO(p_SO), NStncl, NOG, NOG,
     &                    IIF, JJF, IIGF, JJGF, iGs, jGs,
     &                    ProcGrid, NProcI, NProcJ, NProc,
     &                    MyProcI, MyProcJ, MPICOMM
     &                    )
      ENDIF


      !
      ! Zero workspace
      ! 
      DO i=1,NSOR
         SOR(i)=rZERO
      END DO

      DO i=1,NCI
         CI(i)=rZERO
      END DO

      !
      ! Zero coarse-grid solve workspace
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
      
      IF ( NOG.GT.1 ) THEN ! NOG=1 => there isn't a coarse grid
         !
         ! Loop over grids
         !
         DO  K = NOG, 2, -1
            !
            ! Determine the number of points in the stencil
            !
            IF (K.NE.NOG .OR. IFD.NE.1) THEN
               NStncl=5
            ELSE
               NStncl=3
            ENDIF

            !
            ! (fake) memory pointers
            !
            CALL BMG2_SymStd_GET_pointers(
     &                K, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, 
     &                IIF, JJF, IIGF, JJGF 
     &                )
            iGs = IGRD(K,idL_BMG_ICoord)
            jGs = IGRD(K,idL_BMG_JCoord)

            CALL BMG2_SymStd_GET_pointers( 
     &                K-1, IGRD, NOGm, 
     &                p_UC, p_SOC, p_SORC, p_CIC,
     &                IIC, JJC, IIGC, JJGC
     &                )
            iGs_c = IGRD(K-1,idL_BMG_ICoord)
            jGs_c = IGRD(K-1,idL_BMG_JCoord)

            IF ( IBC.EQ.BMG_BCs_definite ) THEN

               ! ---------------------------------------------------------
               !  Setup for the definite case
               ! ---------------------------------------------------------

               !
               ! Operator Induced Interpolation
               !
               TT1 = MPI_Wtime()

               CALL BMG2_SymStd_SETUP_interp_OI(
     &                   K, K-1, SO(p_SO), CI(p_CIC),
     &                   IIF, JJF, IIC, JJC, NOG, IFD, NStncl,
     &                   BMG_iPARMS, BMG_IOFLAG,
     &                   NOGm, IGRD, iWork, NMSGi, pMSG, 
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )


               TT2 = MPI_Wtime()
               BMG_rPARMS(id_BMG2_TIME_SETUP_interp_OI) =
     &              BMG_rPARMS(id_BMG2_TIME_SETUP_interp_OI) + TT2 - TT1
               
               IF ( BMG_IOFLAG(iBMG2_BUG_RESTRICT) )THEN
                  !
                  CALL BMG2_SymStd_DUMP_restrict( 
     &                 BMG_IOFLAG, CI(p_CIC), K, K-1, NOG,
     &                 IIF, JJF, IIGF, JJGF,
     &                 IIC, JJC, IIGC, JJGC,
     &                 iGs, jGs, iGs_c, jGs_c,
     &                 ProcGrid, NProcI, NProcJ, NProc,
     &                 MyProcI, MyProcJ, MPICOMM
     &                 )
                  !
               ENDIF
               !
               ! Coarse-grid operator I^{T}*L*I
               !
               TT1 = MPI_Wtime()

               IF ( BMG_iPARMS(id_BMG2_CG_CONSTRUCT)
     &            .EQ.BMG_CG_CONS_explicit
     &            ) THEN
                  !
                  ! Coarse-grid operator I^{T}*L*I: explicit construction
                  !
                  CALL BMG2_SymStd_SETUP_ITLI_ex( 
     &                      K, K-1, SO(p_SO), SO(p_SOC), CI(p_CIC),
     &                      IIF, JJF, IIC, JJC, iGs, jGs, 
     &                      NOG, IFD, NStncl, BMG_iPARMS, BMG_IOFLAG,
     &                      iWork, NMSGi, pMSGSO,
     &                      MSG_Buffer, NMSGr, MPICOMM
     &                      )

               ELSE IF ( BMG_iPARMS(id_BMG2_CG_CONSTRUCT)
     &                 .EQ.BMG_CG_CONS_block
     &                 ) THEN
                  !
                  ! Coarse-grid operator I^{T}*L*I:  block construction
                  !
                  CALL BMG2_SymStd_SETUP_ITLI_bl( 
     &                      K, K-1, SO(p_SO), SO(p_SOC), CI(p_CIC),
     &                      IIF, JJF, IIC, JJC, iGs, jGs,
     &                      NOG, IFD, NStncl, BMG_iPARMS, BMG_IOFLAG,
     &                      iWork, NMSGi, pMSGSO,
     &                      MSG_Buffer, NMSGr, MPICOMM
     &                      )
               
               ENDIF

               TT2 = MPI_Wtime()
               BMG_rPARMS(id_BMG2_TIME_SETUP_cg_ITLI) =
     &              BMG_rPARMS(id_BMG2_TIME_SETUP_cg_ITLI) + TT2 - TT1
               
               !
               ! Interactively examine the coarse grid stencil
               !
               IF ( BMG_IOFLAG(iBMG2_BUG_STENCIL_CG) ) THEN
                  !
                  CALL BMG2_SymStd_DUMP_stencil( 
     &                 BMG_IOFLAG, SO(p_SOC), NStncl_CG, K-1, NOG,
     &                 IIC, JJC, IIGC, JJGC, iGs, jGs, 
     &                 ProcGrid, NProcI, NProcJ, NProc,
     &                 MyProcI, MyProcJ, MPICOMM
     &                 )
                  !
               ENDIF

            ELSE
               
               WRITE(*,*) 'Periodic BCs not currently available'

            ENDIF
            !
         END DO
         ! 
         ! Relaxation 
         !
         IF ( IBC.EQ.BMG_BCs_definite ) THEN
            !

            TT1 = MPI_Wtime()

            CALL BMG2_SymStd_SETUP_relax(
     &                IFD, IRELAX, SO, NSO, SOR, NSOR, 
     &                IGRD, NOGm, NOG, BMG_iPARMS, BMG_IOFLAG,
     &                iWork, pMSG, NMSGi, MSG_Buffer, NMSGr,
     &                MyProc, MyProcI, MyProcJ,
     &                MPICOMM, XLINECOMM, YLINECOMM
     &                )

            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG2_TIME_SETUP_relax) =
     &           BMG_rPARMS(id_BMG2_TIME_SETUP_relax) + TT2 - TT1
            

         ENDIF

      ENDIF

C ------------------------------------------
C     Setup the coarse grid solve
C ------------------------------------------

      IF ( NOG.EQ.1 ) THEN 
         !
         ! Determine the number of points in the stencil
         !
         IF ( IFD.NE.1 ) THEN
            NStncl=5
         ELSE
            NStncl=3
         ENDIF

      ELSE

         NStncl=5  ! coarsest grid 

      ENDIF
      
      !
      !  Coarsest grid index is always 1
      !
      K = 1
      CALL BMG2_SymStd_GET_pointers( 
     &          K, IGRD, NOGm,
     &          p_UC, p_SOC, p_SORC, p_CIC,
     &          NLx_c, NLy_c, NGx_c, NGy_c
     &          )

      iGs = IGRD(K,idL_BMG_ICoord)
      jGs = IGRD(K,idL_BMG_JCoord)
     
      !
      ! Interactively examine the stencil on the coarsest grid
      !

      ! WRONG: call is missing

      IF ( BMG_iPARMS(id_BMG2_CG_SOLVER).EQ. BMG_CG_SOLVE_LU ) THEN
         !
         ! Setup and factor the matrix for direct solve on the coarsest grid.
         !
         p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
         NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
         NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)
         !
         TT1 = MPI_Wtime()
         !
         CALL BMG2_SymStd_SETUP_cg_LU( 
     &             SO(p_SOC), NLx_c, NLy_c, NStncl, 
     &             BMG_rWORK_CS(p_ABD), NCBW, NCU, MSG_Buffer, NMSGr, 
     &             NProcI, NProcJ, NProc, MyProc,
     &             ProcGrid, ProcCoord, LocArrSize, MPICOMM,
     &             BMG_IOFLAG, BMG_iPARMS
     &             )
         !
         TT2 = MPI_Wtime()
         !
      ELSE
         !
         TT1 = MPI_Wtime()
         !
         CALL BMG2_SymStd_SETUP_cg_boxmg(
     &             NLx_c, NLy_c, NGx_c, NGy_c, iGs, jGs,
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &             SO(p_SOC), NStncl, NOG,
     &             BMG_iWORK_CS, NBMG_iWORK_CS,
     &             BMG_rWORK_CS, NBMG_rWORK_CS,
     &             MSG_Buffer, NMSGr, 
     &             NProcI, NProcJ, NProc, MyProc,
     &             ProcGrid, ProcCoord, LocArrSize, MPICOMM
     &             )
         !
         TT2 = MPI_Wtime()
      ENDIF 

      !
      ! FIXME: SHOULD BE called SETUP_cg_solve
      !
      BMG_rPARMS(id_BMG2_TIME_SETUP_cg_LU) =
     &     BMG_rPARMS(id_BMG2_TIME_SETUP_cg_LU) + TT2 - TT1

c$$$      IF ( MyProc.EQ.1 ) THEN
c$$$         WRITE(*,*) 'BMG err (parts-end) = ',
c$$$     &              BMG_iPARMS(id_BMG2_Err_Code)
c$$$         WRITE(*,*) 'BMG ext err (parts-end) = ',
c$$$     &              BMG_iPARMS(id_BMG2_Ext_Err_Code)
c$$$      ENDIF

      IF ( BMG_iPARMS(id_BMG2_Err_Code).NE.iZERO ) THEN
         RETURN
      END IF

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_parts.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
