      SUBROUTINE BMG2_SymStd_ncycle( 
     &                KC, KCF, KF,
     &                IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm, NOG,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWorkMSG, NMSGi, pMSG, pLS, 
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &               )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   BMG2_SymStd_ncycle performs a single multigrid n-cycle.  
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
      !  BoxMG Cycle Parameters
      !
      INTEGER  KC, KCF, KF, IFD, IU, ID, IVW, IRELAX, IRELAX_SYM

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
      REAL*RKIND   CI(NCI), Q(NF), QF(NF), RES(NF), SO(NSO), SOR(NSOR)

      INTEGER  BMG_iWORK_CS(NBMG_iWORK_CS)
      REAL*RKIND   BMG_rWORK_CS(NBMG_rWORK_CS)

      !
      !  MSG workspace dimensions
      !
      INTEGER  NMSGi, NMSGr
      !
      !  MSG workspace and pointers
      !
      INTEGER  iWorkMSG(NMSGi), 
     &         pMSG(NBMG_pMSG,NOG),
     &         pLS(NBMG_pLS,NOG)
      REAL*RKIND   MSG_Buffer(NMSGr)
      

      !
      !  BMG_MSG grid stuff
      !
      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)
      
C ----------------------------
C     Local Declarations
C
      INTEGER  MyProc, MyProcI, MyProcJ, NProc, NProcI, NProcJ, 
     &         MPICOMM, XLINECOMM, YLINECOMM, MPI_IERROR, NStncl

      INTEGER  iGs, jGs, NLx, NLy, NGx, NGy

      INTEGER  p_ProcGrid, p_ProcCoord

      INTEGER  i, j, k, p_ABD, p_BBD, p_CI, p_SO, p_SOR, p_U,
     &         NCBW, NCU
      REAL*RKIND   RES_L2
      REAL*8 TT1, TT2

C =========================================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------
C --------------------------------------------
C     MPI Parallel Data Decomposition:
C --------------------------------------------

      !
      ! MPI Communicator
      !
      MPICOMM   = BMG_MSG_iGRID(id_BMG_MSG_COMM)

      !
      ! MPI Communicators for line relaxation
      !
      XLINECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_x)
      YLINECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_y)

      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      k = KCF   ! set the current grid index

      IF (KCF.EQ.KC) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,*) '*** FATAL ERROR: BMG2_SymStd_ncycle.f ****'
            WRITE(*,*) 'The current finest grid = ', KCF
            WRITE(*,*) 'The coarest grid        = ', KC
            WRITE(*,*) 'ERROR; minimum number of grids is 2 !'
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,9)
         RETURN

      ENDIF

      ! set the n-cycle counter to zero on all grids
      DO i = KC, KCF
         IGRD(i,idL_BMG_IVW)=iZERO
      END DO

C -----------------------------------------------------------
C     Begin the n-cycle
C     (this is also the n-cycle recursive loop boundary)
C -----------------------------------------------------------

 140  CONTINUE     

         IF ( k.EQ.KC .AND. IGRD(k,idL_BMG_IVW).EQ.iZERO ) THEN
            !
            ! Solve on the coarsest grid
            !
            CALL BMG2_SymStd_GET_pointers( 
     &                KC, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI,
     &                NLx, NLy, NGx, NGy 
     &                )
            !
            iGs = IGRD(KC,idL_BMG_Icoord)
            jGs = IGRD(KC,idL_BMG_Jcoord)
            !
            !
            IF ( BMG_iPARMS(id_BMG2_CG_SOLVER).EQ.BMG_CG_SOLVE_LU ) THEN
               !
               p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
               p_BBD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CU)
               !
               NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
               NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)
               !
               TT1 = MPI_Wtime()
               !
               CALL BMG2_SymStd_SOLVE_cg_LU(
     &                   Q(p_U), QF(p_U), NLx, NLy,
     &                   BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                   BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD),
     &                   NCBW, NCU, NOGm,
     &                   NProcI, NProcJ, NProc, MyProc,
     &                   BMG_MSG_iGRID(p_ProcGrid),
     &                   BMG_MSG_iGRID(p_ProcCoord),
     &                   iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                   iWorkMSG, pMSG, NMSGi,
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )
               !
               TT2 = MPI_Wtime()
            ELSE
               !
               TT1 = MPI_Wtime()
               !
c$$$               CALL BMG2_SymStd_DUMP_vector( 
c$$$     &              BMG_IOFLAG, QF(p_U), k, NOG,
c$$$     &              NLx, NLy, NGx, NGy, iGs, jGs,
c$$$     &              'QF-solve_cg_boxmg-MPI', .FALSE.,
c$$$     &              BMG_MSG_iGRID(p_ProcGrid),
c$$$     &              NProcI, NProcJ, NProc,
c$$$     &              MyProcI, MyProcJ, MPICOMM
c$$$     &              )

               CALL BMG2_SymStd_SOLVE_cg_boxmg( 
     &                   Q(p_U), QF(p_U), NOG,
     &                   BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                   NLx, NLy, NGx, NGy, iGs, jGs,
     &                   BMG_iWORK_CS, NBMG_iWORK_CS,
     &                   BMG_rWORK_CS, NBMG_rWORK_CS,
     &                   NProcI, NProcJ, NProc, MyProc,
     &                   BMG_MSG_iGRID(p_ProcGrid),
     &                   BMG_MSG_iGRID(p_ProcCoord),
     &                   iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )

c$$$               CALL BMG2_SymStd_DUMP_vector( 
c$$$     &              BMG_IOFLAG, Q(p_U), k, NOG,
c$$$     &              NLx, NLy, NGx, NGy, iGs, jGs,
c$$$     &              'Q-solve_cg_boxmg-MPI', .FALSE.,
c$$$     &              BMG_MSG_iGRID(p_ProcGrid),
c$$$     &              NProcI, NProcJ, NProc,
c$$$     &              MyProcI, MyProcJ, MPICOMM
c$$$     &              )

               !
               TT2 = MPI_Wtime()
               !
            ENDIF 

            IF ( BMG_iPARMS(id_BMG2_Err_Code).NE.iZERO ) THEN
               RETURN
            END IF

            BMG_rPARMS(id_BMG2_TIME_SOLVE_CG) =
     &           BMG_rPARMS(id_BMG2_TIME_SOLVE_CG) + TT2 - TT1

            IGRD(k,idL_BMG_IVW) = IVW
            
            IF( BMG_IOFLAG(iBMG2_BUG_RES_CG_SOLVE)
     &          .AND. BMG_iPARMS(id_BMG2_CG_SOLVER).EQ.BMG_CG_SOLVE_LU
     &        ) THEN
               NStncl = 5  ! always the same on coarse grid if NOG > 1
               !
               CALL BMG2_SymStd_residual(
     &                   K, KF, NOG, RES_L2,
     &                   SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &                   NLx, NLy, IFD, NStncl, IRELAX, IRELAX_SYM,
     &                   iWorkMSG, NMSGi, pMSG,
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )
               !
               CALL BMG2_SymStd_UTILS_norm_l2(
     &                   RES(p_U), NLx, NLy, RES_L2, MPICOMM 
     &                   )
               !
               IF ( MyProc.EQ.1 ) THEN
                  WRITE (*,100) K, RES_L2
               ENDIF
               !
            ENDIF
            !
            ! continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_IVW).LT.IVW ) THEN
            !
            ! Move down
            !
            CALL BMG2_SymStd_updown( 
     &                k, KF, NOG, BMG_DOWN,
     &                IFD, IU, ID, IRELAX, IRELAX_SYM, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                iWorkMSG, NMSGi, pMSG, pLS,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID 
     &                )

            IGRD(k,idL_BMG_IVW) = IGRD(k,idL_BMG_IVW) + 1 
            k=k-1
            !
            ! continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_IVW).EQ.IVW ) THEN
            !
            ! Move up
            !
            IGRD(k,idL_BMG_IVW)=iZERO   ! reset n-cycle counter to zero
            k=k+1
            !
            CALL BMG2_SymStd_updown( 
     &                k, KF, NOG, BMG_UP,
     &                IFD, IU, ID, IRELAX, IRELAX_SYM, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC, 
     &                SO, NSO, SOR, NSOR, CI, NCI, 
     &                IGRD, NOGm,
     &                iWorkMSG, NMSGi, pMSG, pLS,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID 
     &                )

            !
            ! Are we on the finest grid? 
            ! 
            IF (k.EQ.KCF) THEN 
               ! n-cycle is done
               RETURN
            ELSE
               ! continue the n-cycle
               GOTO 140  
            ENDIF
            !
         ENDIF

C =========================================================================

 100     FORMAT (' LEVEL',I2,' RESIDUAL NORM= ',1P,E12.5)

C ==============================

         RETURN
         END
