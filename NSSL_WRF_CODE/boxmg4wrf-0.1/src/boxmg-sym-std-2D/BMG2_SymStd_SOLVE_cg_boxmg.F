      SUBROUTINE BMG2_SymStd_SOLVE_cg_boxmg(
     &                Q, QF, NOG,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                II, JJ, NGx, NGy, iGs, jGs,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                NProcI, NProcJ, NProc, MyProc,
     &                ProcGrid, ProcCoord, LocArrSize,
     &                WS, NMSGr, MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SOLVE_cg_boxmg solves on the coarsest grid using
C     the serial BoxMG code (BMG2_SER_SymStd_SOLVE_boxmg).
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     QF        Refer to BMG2_SymStd_SOLVE_boxmg
C     
C     II        Number of grid points in x direction, including
C               two fictitious points.
C     JJ        Number of grid points in y direction, including
C               two fictitious points.
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
C     Q         Refer to BMG2_SymStd_SOLVE_boxmg
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================
#if RKIND == 4
#define FLOAT_PRECISION MPI_REAL
#else
#define FLOAT_PRECISION MPI_DOUBLE_PRECISION
#endif

      IMPLICIT NONE

C ----------------------------
C     Includes
C 
      INCLUDE 'BMG_stencils.h'
#include      "BMG_constants.h"
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

      INCLUDE 'BMG_SER_workspace-partial.h'
      INCLUDE 'BMG_SER_parameters.h'

      INCLUDE 'mpif.h'

C ----------------------------
C     Argument Declarations
C ----------------------------
 
      !
      !  Global/Local indexing
      !
      INTEGER II, JJ,
     &        NGx, NGy,
     &        iGs, jGs

      INTEGER NOG
      !
      !  Solution and RHS
      !
      REAL*RKIND Q(II,JJ), QF(II,JJ)

      !
      !  Processor grid
      !
      INTEGER NProcI, NProcJ, NProc, MyProc
      INTEGER ProcCoord(2,NProc), ProcGrid(NProcI,NProcJ)

      !
      !  MPI REAL buffer space
      !
      INTEGER NMSGr
      REAL*RKIND  WS(NMSGr)

      !
      !  MPI communicator
      !
      INTEGER  MPICOMM

      !
      !    
      ! 
      INTEGER LocArrSize(3,*) 

      !
      !
      !
      INTEGER BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND  BMG_rPARMS(NBMG_rPARMS)
      LOGICAL BMG_IOFLAG(NBMG_IOFLAG)

      !
      !  Coarse-grid solve workspace
      !
      INTEGER  NBMG_iWORK_CS, NBMG_rWORK_CS
      INTEGER  BMG_iWORK_CS(NBMG_iWORK_CS)
      REAL*RKIND   BMG_rWORK_CS(NBMG_rWORK_CS)

C ----------------------------
C     Local Declarations
C ----------------------------

      INTEGER i, i1, i2, ibw, j, j1, k, k1, kl, KK_t, n
      INTEGER IIG, JJG, proc, KKMAX, P1, P2
      INTEGER P1SUM, P2SUM, PXP1, PYP2, p_WS
      INTEGER INT_TEMP, INT_TEMP1, INT_TEMP2
      INTEGER LARGESTNODES, IERR, MPI_IERROR

      INTEGER MyProcI, MyProcJ

      INTEGER i_WS, iGs_ws, jGs_ws, NLx_ws, NLy_ws,
     &        Proc_ws, ProcI_ws, ProcJ_ws

      INTEGER p_pWORK, p_iPARMS, p_iWORK, p_iWORK_PL, 
     &        p_rPARMS, p_rWORK, p_rWORK_PL

      INTEGER NBMG_SER_iWORK, NBMG_SER_iWORK_PL, 
     &        NBMG_SER_rWORK, NBMG_SER_rWORK_PL

      INTEGER NOGm_SER, NOG_SER, NFm_SER, NSOm_SER
      INTEGER NF, NC, NCI, NSO, NSOR, NCBW, NCU,
     &        p_CI, p_CSO, p_CU, p_iGRD, p_Q,
     &        p_RES, p_SO, p_SOR, p_U
      INTEGER pSI, pSR

      LOGICAL BMG_SER_IOFLAG(NBMG_SER_IOFLAG),
     &        BMG_SER_InWORK(NBMG_SER_InWORK)

C =========================================================================

C ------------------------------------------------
C     Unpack 
C ------------------------------------------------

c$$$      MyProcI = ProcCoord(1,MyProc)
c$$$      MyProcJ = ProcCoord(2,MyProc)
c$$$
c$$$      CALL BMG2_SymStd_DUMP_vector( 
c$$$     &          BMG_IOFLAG, QF, iONE, NOG,
c$$$     &          II, JJ, NGx, NGy,
c$$$     &          iGs, jGs, 'QF-cg', .FALSE.,
c$$$     &          ProcGrid, NProcI, NProcJ, NProc,
c$$$     &          MyProcI, MyProcJ, MPICOMM
c$$$     &          )
c$$$

C  -----------------------------------------------
C     Calculate the global number of points
C ------------------------------------------------

      !
      ! Global number in x
      !
      IIG=2
      do i=1,NProcI
         IIG = IIG + LocArrSize(1,ProcGrid(i,1)) - 2
      enddo
      
      !
      ! Global number in y
      !
      JJG=2
      do j=1,NProcJ
         JJG = JJG + LocArrSize(2,ProcGrid(1,j)) - 2
      enddo
      
C ------------------------------------------------
C     Find the largest local array
C ------------------------------------------------

      !
      ! Loop over local x dimensions
      !
      INT_TEMP = LocArrSize(1,1) - 2
      INT_TEMP1 = INT_TEMP
      DO proc=2, NProcI*NProcJ
         INT_TEMP = LocArrSize(1,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP1) THEN
            INT_TEMP1 = INT_TEMP
         END IF
      END DO
      
      !
      ! Loop over local y dimensions
      !
      INT_TEMP = LocArrSize(2,1) - 2
      INT_TEMP2 = INT_TEMP
      DO proc=2, NProcI*NProcJ
         INT_TEMP = LocArrSize(2,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP2) THEN
            INT_TEMP2 = INT_TEMP
         END IF
      END DO
      
      !
      ! Conservative: take largest from each
      !
      LARGESTNODES = INT_TEMP1 * INT_TEMP2 + 5
          
C ------------------------------------------------
C     Copy all information into the buffer
C ------------------------------------------------

      WS(1) = MyProc
      WS(2) = II
      WS(3) = JJ
      WS(4) = iGs
      WS(5) = jGs
      INT_TEMP = 6

      DO J = 2, JJ - 1
         DO I = 2, II - 1
            WS(INT_TEMP) = QF(I,J)
            INT_TEMP = INT_TEMP + 1
         END DO
      END DO

C ------------------------------------------------
C     Send/Receive information to/from everybody
C ------------------------------------------------

      IF (BMG_iPARMS(id_BMG2_CG_COMM) .eq. BMG_CG_ALLGATHER) THEN
         
         CALL MPI_ALLGATHER( 
     &        WS(1), LARGESTNODES, FLOAT_PRECISION,
     &        WS(LARGESTNODES+1), LARGESTNODES, FLOAT_PRECISION,
     &        MPICOMM, IERR
     &        ) 
      
      ELSE

         CALL MPI_GATHER( 
     &        WS(1), LARGESTNODES, FLOAT_PRECISION,
     &        WS(LARGESTNODES+1), LARGESTNODES, FLOAT_PRECISION,
     &        iZERO, MPICOMM, IERR
     &        ) 
         
      ENDIF

C ------------------------------------------------
C     Assemble the global right hand side
C ------------------------------------------------

      ! 
      ! Setup local pointers
      !
      p_pWORK  = BMG_iWORK_CS(ip_BMG_pWORK_CS)

      p_iPARMS = BMG_iWORK_CS(ip_BMG_iPARMS_CS)
      p_rPARMS = BMG_iWORK_CS(ip_BMG_rPARMS_CS)
         
      p_iWORK  = BMG_iWORK_CS(ip_BMG_iWORK_CS)
      p_rWORK  = BMG_iWORK_CS(ip_BMG_rWORK_CS)
      
      !
      ! Setup workspace dimensions
      !
      NBMG_SER_iWORK = BMG_iWORK_CS(id_BMG_iWORK_CS)
      NBMG_SER_rWORK = BMG_iWORK_CS(id_BMG_rWORK_CS)
      
      NBMG_SER_iWORK_PL = BMG_iWORK_CS(id_BMG_iWORK_PL_CS)
      NBMG_SER_rWORK_PL = BMG_iWORK_CS(id_BMG_rWORK_PL_CS)
      
      !
      ! Set the maximum number of grids
      !
      NOGm_SER = BMG_iWORK_CS(id_BMG_NOGm_CS)

      !
      !  Local pointers for the SERIAL call
      !
      p_SO  = BMG_iWORK_CS(p_pWORK+ip_SO-1)
      p_U   = BMG_iWORK_CS(p_pWORK+ip_U-1)
      p_Q   = BMG_iWORK_CS(p_pWORK+ip_Q-1)
      p_CI  = BMG_iWORK_CS(p_pWORK+ip_CI-1)
      p_RES = BMG_iWORK_CS(p_pWORK+ip_RES-1)
      p_SOR = BMG_iWORK_CS(p_pWORK+ip_SOR-1)
      p_CSO = BMG_iWORK_CS(p_pWORK+ip_CSO-1)
      p_CU  = BMG_iWORK_CS(p_pWORK+ip_CU-1)
         
      p_iGRD  = BMG_iWORK_CS(p_pWORK+ip_iG-1)

      !
      ! Initialize irrelevant data
      !
      NFm_SER = 1
      NSOm_SER = 1

      !
      !  Local dimensional parameters for the serial call
      !
      NOG_SER  = BMG_iWORK_CS(p_iPARMS+id_BMG2_SER_DIM_NOG-1)

      NF       = BMG_iWORK_CS(p_iPARMS+id_BMG2_SER_DIM_NF-1)
      NC       = BMG_iWORK_CS(p_iPARMS+id_BMG2_SER_DIM_NC-1)
      NSO      = BMG_iWORK_CS(p_iPARMS+id_BMG2_SER_DIM_NSO-1)
      NCI      = BMG_iWORK_CS(p_iPARMS+id_BMG2_SER_DIM_NCI-1)
      NSOR     = BMG_iWORK_CS(p_iPARMS+id_BMG2_SER_DIM_NSOR-1)
      NCBW     = BMG_iWORK_CS(p_iPARMS+id_BMG2_SER_DIM_NCBW-1)
      NCU      = BMG_iWORK_CS(p_iPARMS+id_BMG2_SER_DIM_NCU-1)


      !
      ! Synchronize the default parameter values
      !
      CALL BMG2_SymStd_SETUP_cg_parms(
     &          BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &          BMG_iWORK_CS(p_iPARMS), 
     &          BMG_rWORK_CS(p_rPARMS), 
     &          BMG_SER_IOFLAG,
     &          NGx, NGy, NOG
     &          )

      !
      !  Zero the intial guess:
      !
      DO i = 1, NF
         BMG_rWORK_CS(p_U+i-1) = rZERO
      ENDDO


      IF ( BMG_iPARMS(id_BMG2_CG_COMM) .eq. BMG_CG_ALLGATHER .OR. 
     &     MyProc .eq. iONE ) THEN

         !
         !  Copy WS blocks into the SERIAL RHS
         !
         DO proc=1, NProcI*NProcJ
            
            i_WS = proc * LARGESTNODES + 1

            Proc_ws = WS(i_WS)
            i_WS = i_WS + 1
                        
            NLx_ws = WS(i_WS)
            NLy_ws = WS(i_WS+1)
            i_WS = i_WS + 2

            iGs_ws = WS(i_WS)
            jGs_ws = WS(i_WS+1)
            i_WS = i_WS + 2

            IF ( NLx_ws .NE. LocArrSize(1,Proc_ws)
     &         .OR.  NLy_ws .NE. LocArrSize(2,Proc_ws)  ) THEN
               !
               IERR = IERR + 1
               !
               WRITE(*,*) 'Error: LocArrSize is inconsitent ... '
               WRITE(*,*) 
               WRITE(*,*) ' Proc_ws = ', Proc_ws
               WRITE(*,*) 
               WRITE(*,*) ' NLx_ws = ', NLx_ws
               WRITE(*,*) ' NLy_ws = ', NLy_ws
               WRITE(*,*) 
               WRITE(*,*) ' NLx_ws = ', NLx_ws
               WRITE(*,*) ' NLy_ws = ', NLy_ws
            ENDIF

            ProcI_ws = ProcCoord(1,Proc_ws)
            ProcJ_ws = ProcCoord(2,Proc_ws)

            !
            ! Copy the current block of WS into the 
            ! correct block of the SERIAL RHS.
            !
            CALL BMG2_SymStd_COPY_cg_WS_RHS(
     &                WS(i_WS), NLx_ws-2, NLy_ws-2,
     &                BMG_rWORK_CS(p_Q), NGx, NGy,
     &                iGs_ws, jGs_ws
     &                )


            !
         ENDDO

         !
         !  Override defaults (.FALSE.) for debugging:
         !
         BMG_SER_IOFLAG(iBMG2_SER_BUG_STENCIL_FG)   = .FALSE.
         BMG_SER_IOFLAG(iBMG2_SER_BUG_STENCIL_CG)   = .FALSE.
         BMG_SER_IOFLAG(iBMG2_SER_BUG_STENCIL_CG1)  = .FALSE.

         BMG_SER_IOFLAG(iBMG2_SER_OUT_TIME_SETUP)   = .FALSE.
         BMG_SER_IOFLAG(iBMG2_SER_OUT_TIME_CYCLING) = .FALSE.
         BMG_SER_IOFLAG(iBMG2_SER_OUT_TIME_TOTAL)   = .FALSE.

         !
         ! Call SERIAL BoxMG to solve
         !
         BMG_iWORK_CS(p_iPARMS+id_BMG2_SER_SETUP-1) = BMG_SER_SETUP_none
c$$$
c$$$         CALL BMG2_SER_SymStd_DUMP_vector( 
c$$$     &             BMG_SER_IOFLAG, BMG_rWORK_CS(p_Q),
c$$$     &             NGx, NGy, NOG_SER, NOGm_SER,
c$$$     &             'output', 'residual-serial-start', .FALSE.
c$$$     &             )

         CALL BMG2_SER_SymStd_SOLVE_boxmg( 
     &             NGx-2, NGy-2,
     &             BMG_iWORK_CS(p_iPARMS), BMG_rWORK_CS(p_rPARMS),
     &             BMG_SER_IOFLAG,
     &             BMG_rWORK_CS(p_U), BMG_rWORK_CS(p_Q),
     &             BMG_rWORK_CS(p_RES), NF, NC, 
     &             BMG_rWORK_CS(p_SO), NSO,
     &             BMG_rWORK_CS(p_SOR), NSOR,
     &             BMG_rWORK_CS(p_CI), NCI,
     &             BMG_rWORK_CS(p_CSO),
     &             BMG_rWORK_CS(p_CU), NCBW, NCU,
     &             BMG_iWORK_CS(p_iGRD), NOGm_SER, NOG_SER
     &             )

c$$$         CALL BMG2_SER_SymStd_DUMP_vector( 
c$$$     &             BMG_SER_IOFLAG, BMG_rWORK_CS(p_U),
c$$$     &             NGx, NGy, NOG_SER, NOGm_SER,
c$$$     &             'output', 'solution-serial-return', .FALSE.
c$$$     &             )

      ENDIF

C ------------------------------------------------
C     Extract the local data 
C ------------------------------------------------
      
      IF ( BMG_iPARMS(id_BMG2_CG_COMM).EQ.BMG_CG_GATHER_SCATTER) THEN

         !
         ! Broadcast the solution vector
         !
         CALL MPI_Bcast(
     &            BMG_rWORK_CS(p_U), NGx*NGy, FLOAT_PRECISION, 
     &            iZERO, MPICOMM, IERR
     &            )
         
      END IF

      CALL BMG2_SymStd_COPY_cg_rV_G_L(
     &                BMG_rWORK_CS(p_U), NGx, NGy,
     &                Q, II, JJ, iGs, jGs,
     &                ProcCoord, NProcI, NProcJ,
     &                Nproc, MyProc
     &                )


C ========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SOLVE_cg_boxmg.f',/,5X,A)
 510  FORMAT (5X,A,1X,I3)

C ===========================================

      RETURN
      END
