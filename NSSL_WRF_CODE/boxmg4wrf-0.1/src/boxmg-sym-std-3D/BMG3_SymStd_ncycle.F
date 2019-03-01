      SUBROUTINE BMG3_SymStd_ncycle( 
     &                KC, KCF, NOG, 
     &                IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                BMG_iWORK_CS, NBMG_iWORK_CS,
     &                BMG_rWORK_CS, NBMG_rWORK_CS,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   BMG3_SymStd_ncycle performs a single multigrid n-cycle.  
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
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  KC, KCF, NOG, IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, ISEP,
     &         NF, NC, NSO, NSOR, NCI, NOGm, NMSGi, NMSGr,
     &         NBMG_MSG_iGRID, NBMG_iWORK_CS, NBMG_iWORK_PL,
     &         NBMG_rWORK_CS, NBMG_rWORK_PL

      INTEGER  BMG_iPARMS(NBMG_iPARMS),
     &         BMG_iWORK_CS(NBMG_iWORK_CS),
     &         BMG_iWORK_PL(NBMG_iWORK_PL), 
     &         iWorkMSG(NMSGi), pMSG(NBMG_pMSG,NOGm), 
     &         pMSGSO(NBMG_pMSG,NOGm), IGRD(NOGm,NBMG_pIGRD),
     &         BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS),
     &         BMG_rWORK_CS(NBMG_rWORK_CS),
     &         BMG_rWORK_PL(NBMG_rWORK_PL),
     &         CI(NCI), Q(NF), QF(NF), RES(NF),
     &         SO(NSO), SOR(NSOR), MSG_Buffer(NMSGr)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, NProc, NProcI, NProcJ, 
     &         NProcK, MPICOMM, MPI_IERROR, NStncl

      INTEGER  iGs, jGs, kGs, NLx, NLy, NLz, NGx, NGy, NGz

      INTEGER  p_ProcGrid, p_ProcCoord

      INTEGER  i, j, k, p_ABD, p_BBD, p_CI, p_SO, p_SOR, p_U,
     &         NCBW, NCU
      REAL*8   RES_L2
      REAL*8  TT1, TT2

C ==========================================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------
C --------------------------------------------
C     MPI Parallel Data Decomposition:
C --------------------------------------------
      !
      ! MPI Communicators
      !
      MPICOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM)
      !
      ! Processor Grid:
      !
      NProc  = BMG_MSG_iGRID(id_BMG_MSG_NProc)
      NProcI = BMG_MSG_iGRID(id_BMG_MSG_NProcI)
      NProcJ = BMG_MSG_iGRID(id_BMG_MSG_NProcJ)
      NProcK = BMG_MSG_iGRID(id_BMG_MSG_NProcK)

      MyProc  = BMG_MSG_iGRID(id_BMG_MSG_MyProc)
      MyProcI = BMG_MSG_iGRID(id_BMG_MSG_MyProcI)
      MyProcJ = BMG_MSG_iGRID(id_BMG_MSG_MyProcJ)
      MyProcK = BMG_MSG_iGRID(id_BMG_MSG_MyProcK)

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      k = KCF   ! set the current grid index

      IF (KCF.EQ.KC) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,*) '*** FATAL ERROR: BMG3_SymStd_ncycle.f ****'
            WRITE(*,*) 'The current finest grid = ', KCF
            WRITE(*,*) 'The coarest grid        = ', KC
            WRITE(*,*) 'ERROR; minimum number of grids is 2 !'
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,18)
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
         
         IF ( K.EQ.KC .AND. IGRD(k,idL_BMG_IVW).EQ.iZERO ) THEN
            !
            ! Solve on the coarsest grid
            !
            CALL BMG3_SymStd_GET_pointers( 
     &                KC, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, 
     &                NLx, NLy, NLz, NGx, NGy, NGz 
     &                )
            !
            iGs = IGRD(KC,idL_BMG_Icoord)
            jGs = IGRD(KC,idL_BMG_Jcoord)
            kGs = IGRD(KC,idL_BMG_Kcoord)
            !
            IF ( BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU ) THEN
               !
               p_ABD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CSO)
               p_BBD = BMG_iWORK_CS(ip_BMG_iWORK_CS_CU)
               !
               NCBW  = BMG_iWORK_CS(id_BMG_iWORK_CS_NCBW)
               NCU   = BMG_iWORK_CS(id_BMG_iWORK_CS_NCU)
               !
               TT1 = MPI_Wtime()
               !
               CALL BMG3_SymStd_SOLVE_cg_LU( 
     &                   q(p_U), qf(p_U), NLx, NLy, NLz, 
     &                   BMG_rWORK_CS(p_ABD), BMG_rWORK_CS(p_BBD), 
     &                   NCBW, NCU, NOGm,
     &                   NProcI, NProcJ, NProcK, NProc, MyProc,
     &                   BMG_MSG_iGRID(p_ProcGrid),
     &                   BMG_MSG_iGRID(p_ProcCoord),
     &                   iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                   MSG_Buffer, NMSGr, MPICOMM, BMG_iPARMS
     &                   )
               !
               TT2 = MPI_Wtime()
               !
            ELSE
               !
               TT1 = MPI_Wtime()
               !
               CALL BMG3_SymStd_SOLVE_cg_boxmg( 
     &                   Q(p_U), QF(p_U), NOG,
     &                   BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                   NLx, NLy, NLz, NGx, NGy, NGz,
     &                   iGs, jGs, kGs,
     &                   BMG_iWORK_CS, NBMG_iWORK_CS,
     &                   BMG_rWORK_CS, NBMG_rWORK_CS,
     &                   NProcI, NProcJ, NProcK, NProc, MyProc,
     &                   BMG_MSG_iGRID(p_ProcGrid),
     &                   BMG_MSG_iGRID(p_ProcCoord),
     &                   iWorkMSG(pMSG(ipL_MSG_LocalArraySize,1)),
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )

               !
               TT2 = MPI_Wtime()
               !
            ENDIF

            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            BMG_rPARMS(id_BMG3_TIME_SOLVE_CG)
     &          = BMG_rPARMS(id_BMG3_TIME_SOLVE_CG) + TT2 - TT1
            
            !
            IGRD(k,idL_BMG_IVW) = IVW
            IF( BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE)
     &          .AND. BMG_iPARMS(id_BMG3_CG_SOLVER).EQ.BMG_CG_SOLVE_LU
     &        ) THEN
               NStncl = 14   ! always the same on coarse grid if NOG > 1
               CALL BMG3_SymStd_residual( 
     &                   KC, NOG, IFD, 
     &                   q(p_U), qf(p_U), so(p_SO), RES(p_U),
     &                   NLx, NLy, NLz, NStncl,
     &                   iWorkMSG, NMSGi, pMSG,
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )               
               CALL BMG3_SymStd_UTILS_norm_l2( 
     &                   RES(p_U), NLx, NLy, NLz, RES_L2, MPICOMM
     &                   )
               !
               IF ( MyProc.EQ.1 ) THEN 
                  WRITE (*,100) k, RES_L2
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
            CALL BMG3_SymStd_updown( 
     &                k, KCF, NOG, BMG_DOWN, 
     &                IFD, IU, ID, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID 
     &                )
            !
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            IGRD(k,idL_BMG_IVW) = IGRD(k,idL_BMG_IVW) + 1 
            k=k-1
            !
            ! Continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_IVW).EQ.IVW ) THEN
            !
            ! Move up
            !
            IGRD(k,idL_BMG_IVW) = iZERO   ! reset n-cycle counter to zero
            k=k+1
            !
            CALL BMG3_SymStd_updown(
     &                k, KCF, NOG, BMG_UP,
     &                IFD, IU, ID, IRELAX, IRELAX_SYM, ISEP,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iWorkMSG, NMSGi, pMSG, pMSGSO,
     &                MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID 
     &                )
            !
            IF (BMG_iPARMS(id_BMG3_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            !
            ! Are we on the finest grid? 
            ! 
            IF (K.EQ.KCF) THEN 
               ! n-cycle is done
               RETURN
            ELSE
               ! continue the n-cycle
               GOTO 140  
            ENDIF
            !
         ENDIF

C ==========================================================================

 100     FORMAT (1X,'(3D)LEVEL ',I2,' RESIDUAL NORM = ',1P,E12.5)

C ===========================================

         RETURN
         END
