      SUBROUTINE BMG2_SymStd_updown( 
     &                K, KF, NOG, UPDOWN,
     &                IFD, IU, ID, IRELAX, IRELAX_SYM, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &                Q, QF, RES, NF, NC, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                iWorkMSG, NMSGi, pMSG, pLS,
     &                MSG_Buffer, NMSGr, 
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_updown performs all the necessary tasks associated
C     with moving "down" to a coarser grid or moving "up" to a finer
C     grid.  It is necessary because f77 does not support recursion.
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
      INCLUDE 'MSG.h'

#include      "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

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
      INTEGER  KC, KCF, KF, IFD, IU, ID, IVW, IRELAX, IRELAX_SYM,
     &         UPDOWN

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
     &         MPICOMM, XLINECOMM, YLINECOMM, MPI_IERROR

      INTEGER  iGs, jGs, NLx, NLxc, NLy, NLyc, NGx, NGxc, NGy, NGyc

      INTEGER  p_ProcGrid, p_ProcCoord

      INTEGER  i, k, NStncl, 
     &         p_CI, p_CIC, p_SO, p_SOC,
     &         p_SOR, p_SORC, p_U, p_UC

      REAL*RKIND   RES_L2
      REAL*8     TT1, TT2

C ==========================================================================

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

      !
      ! Processor Coordinates
      !
      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      !
      ! How many points in the stencil?
      !
      IF ( K.EQ.KF .AND. IFD.EQ.1 ) THEN
         NStncl=3
      ELSE
         NStncl=5
      ENDIF

      IF (UPDOWN.EQ.BMG_DOWN) THEN

         ! Get pointers for grids k and k-1
         CALL BMG2_SymStd_GET_pointers(
     &             K, iGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI,
     &             NLx, NLy, NGx, NGy 
     &             )
         iGs = IGRD(k,idL_BMG_Icoord)
         jGs = IGRD(k,idL_BMG_Jcoord)

         CALL BMG2_SymStd_GET_pointers( 
     &             K-1, iGRD, NOGm,
     &             p_UC, p_SOC, p_SORC, p_CIC,
     &             NLxc, NLyc, NGxc, NGyc 
     &             )
         
         ! Relaxation
         DO i=1, ID

            TT1 = MPI_Wtime()

            CALL BMG2_SymStd_relax ( 
     &                K, KF, NOG, BMG_DOWN, 
     &                SO(p_SO), QF(p_U), Q(p_U), RES(p_U), SOR(p_SOR), 
     &                NLx, NLy, iGs, jGs, BMG_iPARMS,
     &                IFD, NStncl, IRELAX, IRELAX_SYM,
     &                iWorkMSG, NMSGi, pMSG,
     &                MSG_Buffer, NMSGr,
     &                iWorkMSG(pLS(ipL_LS_XDataDist,k)),
     &                iWorkMSG(pLS(ipL_LS_YDataDist,k)),
     &                MPICOMM, XLINECOMM, YLINECOMM
     &                )
            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG2_TIME_relax) =
     &           BMG_rPARMS(id_BMG2_TIME_relax) + TT2 - TT1
 
            IF( BMG_IOFLAG(iBMG2_BUG_RES_RELAX) ) THEN
               !
               CALL BMG2_SymStd_residual(
     &                   K, KF, NOG, RES_L2,
     &                   SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &                   NLx, NLy, IFD, NStncl, IRELAX, IRELAX_SYM,
     &                   iWorkMSG, NMSGi, pMSG,
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )
               CALL BMG2_SymStd_UTILS_norm_l2( 
     &                   RES(p_U), NLx, NLy, RES_L2, MPICOMM 
     &                   )
               !
               IF ( MyProc.EQ.1 ) THEN
                  WRITE (*,100) K, RES_L2
               ENDIF
               !
            ENDIF

         END DO
         
         ! Residual calculation (only necessary when we're 
         ! not in debuging mode
         IF ( .not. BMG_IOFLAG(iBMG2_BUG_RES_RELAX) ) THEN
            CALL BMG2_SymStd_residual(
     &                K, KF, NOG, RES_L2,
     &                SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &                NLx, NLy, IFD, NStncl, IRELAX, IRELAX_SYM,
     &                iWorkMSG, NMSGi, pMSG,
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
         ENDIF

         TT1 = MPI_Wtime()

         ! Restriction
         CALL BMG2_SymStd_restrict( 
     &             K, K-1, NOG,
     &             RES(p_U), QF(p_UC), CI(p_CIC),
     &             NLx, NLy, NLxc, NLyc, iGs, jGs,
     &             iWorkMSG, NMSGi, pMSG, 
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )
         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG2_TIME_restrict) =
     &        BMG_rPARMS(id_BMG2_TIME_restrict) + TT2 - TT1


         ! Zero the initial guess
         CALL BMG2_SymStd_UTILS_rV_zero( Q(p_UC), NLxc, NLyc )

      ELSE IF (UPDOWN.EQ.BMG_UP) THEN

         ! Interpolate and Correct 
         CALL BMG2_SymStd_GET_pointers( 
     &             k, iGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI,
     &             NLx, NLy, NGx, NGy 
     &             )
         iGs = IGRD(k,idL_BMG_Icoord)
         jGs = IGRD(k,idL_BMG_Jcoord)

         CALL BMG2_SymStd_GET_pointers( 
     &             k-1, iGRD, NOGm,
     &             p_UC, p_SOC, p_SORC, p_CIC,
     &             NLxc, NLyc, NGxc, NGyc 
     &             )

         TT1 = MPI_Wtime()
c$$$
c$$$         CALL BMG2_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, Q(p_UC), k-1, NOG,
c$$$     &             NLxc, NLyc, NGxc, NGyc,
c$$$     &             IGRD(k-1,idL_BMG_Icoord),
c$$$     &             IGRD(k-1,idL_BMG_Jcoord),
c$$$     &             'UP-Q-pre-interp-MPI', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProc,
c$$$     &             MyProcI, MyProcJ, MPICOMM
c$$$     &             )

         CALL BMG2_SymStd_interp_add( 
     &             K-1, K, NOG,
     &             Q(p_U), Q(p_UC), RES(p_U),
     &             SO(p_SO), NStncl, CI(p_CIC),
     &             NLxc, NLyc, NLx, NLy, iGs, jGs, 
     &             iWorkMSG, NMSGi, pMSG,
     &             MSG_Buffer, NMSGr, MPICOMM                
     &             )

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG2_TIME_interp_add) =
     &        BMG_rPARMS(id_BMG2_TIME_interp_add) + TT2 - TT1

c$$$         
c$$$         CALL BMG2_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, Q(p_U), k, NOG,
c$$$     &             NLx, NLy, NGx, NGy, iGs, jGs,
c$$$     &             'UP-Q-post-interp-MPI', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProc,
c$$$     &             MyProcI, MyProcJ, MPICOMM
c$$$     &             )
c$$$
c$$$         CALL BMG2_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, QF(p_U), k, NOG,
c$$$     &             NLx, NLy, NGx, NGy, iGs, jGs,
c$$$     &             'UP-QF-post-interp-MPI', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProc,
c$$$     &             MyProcI, MyProcJ, MPICOMM
c$$$     &             )

         ! Relaxation
         DO i=1, IU

            TT1 = MPI_Wtime()   

            CALL BMG2_SymStd_relax ( 
     &                K, KF, NOG, BMG_UP,
     &                SO(p_SO), QF(p_U), Q(p_U), RES(p_U), SOR(p_SOR), 
     &                NLx, NLy, iGs, jGs, BMG_iPARMS,
     &                IFD, NStncl, IRELAX, IRELAX_SYM,
     &                iWorkMSG, NMSGi, pMSG,
     &                MSG_Buffer, NMSGr,
     &                iWorkMSG(pLS(ipL_LS_XDataDist,k)),
     &                iWorkMSG(pLS(ipL_LS_YDataDist,k)),
     &                MPICOMM, XLINECOMM, YLINECOMM
     &                )

            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG2_TIME_relax) =
     &           BMG_rPARMS(id_BMG2_TIME_relax) + TT2 - TT1
            
            IF( BMG_IOFLAG(iBMG2_BUG_RES_RELAX) ) THEN
               !
               CALL BMG2_SymStd_residual(
     &                   K, KF, NOG, RES_L2,
     &                   SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &                   NLx, NLy, IFD, NStncl, IRELAX, IRELAX_SYM,
     &                   iWorkMSG, NMSGi, pMSG,
     &                   MSG_Buffer, NMSGr, MPICOMM
     &                   )
               CALL BMG2_SymStd_UTILS_norm_l2( 
     &                   RES(p_U), NLx, NLy, RES_L2, MPICOMM 
     &                   )
               !
               IF ( MyProc.EQ.1 ) THEN
                  WRITE (*,100) K, RES_L2
               ENDIF
               !
            ENDIF

         END DO

c$$$         CALL BMG2_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, RES(p_U), k, NOG,
c$$$     &             NLx, NLy, NGx, NGy, iGs, jGs,
c$$$     &             'UP-RES-post-smooth-up-MPI', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProc,
c$$$     &             MyProcI, MyProcJ, MPICOMM
c$$$     &             )

      ELSE

         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,*) '*** UPDOWN out of range: '
            WRITE(*,*) 'HAVE: UPDOWN = ', UPDOWN
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,10)
         RETURN

      ENDIF

C ==========================================================================

 100  FORMAT (' LEVEL',I2,' RESIDUAL NORM= ',1P,E12.5)

C =============================

      RETURN
      END

