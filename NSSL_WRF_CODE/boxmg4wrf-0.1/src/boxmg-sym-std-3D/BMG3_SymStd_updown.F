      SUBROUTINE BMG3_SymStd_updown ( 
     &                K, KF, NOG, UPDOWN,
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


C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_updown performs all the necessary tasks associated
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
C -----------------------------

      INCLUDE 'mpif.h'

#include    "BMG_constants.h" 
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

C ----------------------------
C     Argument Declarations
C ----------------------------

      !
      !  BoxMG parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      !
      !  BoxMG Cycle Parameters
      !
      INTEGER  k, KC, KCF, KF, IFD, IU, ID, IVW,
     &         IRELAX, IRELAX_SYM, UPDOWN, ISEP

      !
      !  BoxMG dimensions
      !
      INTEGER  NF, NC, NSO, NSOR, NCI, NOGm, NOG,
     &         NBMG_iWORK_CS, NBMG_iWORK_PL,
     &         NBMG_rWORK_CS, NBMG_rWORK_PL

      !
      !  BoxMG internal pointers
      !
      INTEGER  IGRD(NOGm,NBMG_pIGRD)

      !
      !  BoxMG data arrays
      !
      REAL*RKIND   CI(NCI), Q(NF), QF(NF), RES(NF), SO(NSO), SOR(NSOR)

      !
      !  BoxMG plane relaxation workspace
      !
      INTEGER BMG_iWORK_PL(NBMG_iWORK_PL)
      REAL*RKIND  BMG_rWORK_PL(NBMG_rWORK_PL)

      !
      !  MSG workspace dimensions
      !
      INTEGER  NMSGi, NMSGr
      !
      !  MSG workspace and pointers
      !
      INTEGER  iWorkMSG(NMSGi), 
     &         pMSG(NBMG_pMSG,NOG),
     &         pMSGSO(NBMG_pMSG,NOG)
      REAL*RKIND   MSG_Buffer(NMSGr)
 
      !
      !  BMG_MSG grid stuff
      !
      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)
      
C ----------------------------
C     Local Declarations
C ----------------------------

      INTEGER  I, J, NLx, NLy, NLz, NLx_c, NLy_c, NLz_c, 
     &         NGx, NGy, NGz, NGx_c, NGy_c, NGz_c, 
     &         p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC
      INTEGER  MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM,
     &         iGs, jGs, kGs, NSORv, NStncl
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, MPI_IERROR,
     &         NProc, NProcI, NProcJ, NProcK, p_ProcGrid, p_ProcCoord
      REAL*8   RES_L2
      REAL*8   TT1, TT2

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
      ! MPI Communicatorsn for plane relaxation
      !
      XYPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xy)
      XZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_xz)
      YZPLANECOMM = BMG_MSG_iGRID(id_BMG_MSG_COMM_yz)
 
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

      !
      ! Processor Coordinates
      !
      iGs = IGRD(k,idL_BMG_Icoord)
      jGs = IGRD(k,idL_BMG_Jcoord)
      kGs = IGRD(k,idL_BMG_Kcoord)
      
      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
      p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      !
      ! Set NStncl for the finest grid
      !
      IF ( K.LT.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
         NStncl=14
      ELSE
         NStncl=4
      ENDIF

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE 
         NSORv = 2  ! the same for now but this may change
      ENDIF

      IF (UPDOWN.EQ.BMG_DOWN) THEN

         ! Get pointers for grids k and k-1
         CALL BMG3_SymStd_GET_pointers(     
     &               K, IGRD, NOGm,
     &               p_U, p_SO, p_SOR, p_CI, 
     &               NLx, NLy, NLz, NGx, NGy, NGz
     &               )
         CALL BMG3_SymStd_GET_pointers( 
     &               K-1, IGRD, NOGm,
     &               p_UC, p_SOC, p_SORC, p_CIC, 
     &               NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c 
     &               )
        
         ! Relaxation
         DO 10 i=1, ID
            !
            TT1 = MPI_Wtime()

            CALL BMG3_SymStd_relax( 
     &                k, KF, NOG, UPDOWN,
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,
     &                Q, QF, RES, NF, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iGs, jGs, kGs,
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
     &                )

            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_relax) = 
     &           BMG_rPARMS(id_BMG3_TIME_relax) + TT2 - TT1


            !
            IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) .AND. MyProc.EQ.1 ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
            !
 10      CONTINUE

         ! restrict the residual
         TT1 = MPI_Wtime()

         ! calculate the residual
         CALL BMG3_SymStd_residual( 
     &             K, NOG, IFD, 
     &             Q(p_U), QF(p_U), SO(p_SO), RES(p_U),
     &             NLx, NLy, NLz, NStncl,
     &             iWorkMSG, NMSGi, pMSG, 
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) = 
     &        BMG_rPARMS(id_BMG3_TIME_PCG_PRECON) + TT2 - TT1
         
c$$$         CALL BMG3_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, RES(p_U), K, NOG,
c$$$     &             NLx, NLy, NLz, NGx, NGy, NGz,
c$$$     &             IGRD(K,idL_BMG_ICoord),
c$$$     &             IGRD(K,idL_BMG_JCoord),
c$$$     &             IGRD(K,idL_BMG_KCoord),
c$$$     &             'DOWN-RES-post-relax', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProcK, NProc,
c$$$     &             MyProcI, MyProcJ, MyProcK, MPICOMM
c$$$     &             )

c         ! restrict the residual
c         TT1 = MPI_Wtime()

         CALL BMG3_SymStd_restrict( 
     &             K, K-1,  
     &             RES(p_U), QF(p_UC), CI(p_CIC),
     &             NLx, NLy, NLz, NLx_c, NLy_c, NLz_c,
     &             iGs, jGs, kGs
     &             )

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_restrict) = 
     &        BMG_rPARMS(id_BMG3_TIME_restrict) + TT2 - TT1
         

c$$$         CALL BMG3_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, QF(p_UC), K-1, NOG,
c$$$     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy, NGz,
c$$$     &             IGRD(K-1,idL_BMG_ICoord),
c$$$     &             IGRD(K-1,idL_BMG_JCoord),
c$$$     &             IGRD(K-1,idL_BMG_KCoord),
c$$$     &             'DOWN-QF-post-restrict', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProcK, NProc,
c$$$     &             MyProcI, MyProcJ, MyProcK, MPICOMM
c$$$     &             )

         ! Zero the initial guess
         CALL BMG3_SymStd_UTILS_rV_zero( 
     &               k-1, q(p_UC), NLx_c, NLy_c, NLz_c
     &               )


      ELSE IF (UPDOWN.EQ.BMG_UP) THEN

         ! Get pointers for grids k and k-1
         CALL BMG3_SymStd_GET_pointers(      
     &             k, IGRD, NOGm,
     &             p_U, p_SO, p_SOR, p_CI, 
     &             NLx, NLy, NLz, NGx, NGy, NGz
     &             )
         CALL BMG3_SymStd_GET_pointers( 
     &             k-1, IGRD, NOGm,
     &             p_UC, p_SOC, p_SORC, p_CIC, 
     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c 
     &             )
         
c$$$         CALL BMG3_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, Q(p_UC), K-1, NOG,
c$$$     &             NLx_c, NLy_c, NLz_c, NGx_c, NGy, NGz,
c$$$     &             IGRD(K-1,idL_BMG_ICoord),
c$$$     &             IGRD(K-1,idL_BMG_JCoord),
c$$$     &             IGRD(K-1,idL_BMG_KCoord),
c$$$     &             'UP-Q-pre-interp', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProcK, NProc,
c$$$     &             MyProcI, MyProcJ, MyProcK, MPICOMM
c$$$     &             )

         ! Interpolate and Correct 
         TT1 = MPI_Wtime()

         CALL BMG3_SymStd_interp_add( 
     &             k-1, k, NOG,
     &             Q(p_U) ,Q(p_UC), RES(p_U),
     &             SO(p_SO), NStncl, CI(p_CIC),
     &             NLx_c, NLy_c, NLz_c,
     &             NLx, NLy, NLz,
     &             iGs, jGs, kGs,
     &             iWorkMSG, NMSGi, pMSG,
     &             MSG_Buffer, NMSGr, MPICOMM
     &             )

         TT2 = MPI_Wtime()
         BMG_rPARMS(id_BMG3_TIME_interp_add) =
     &        BMG_rPARMS(id_BMG3_TIME_interp_add) + TT2 - TT1

                  
c$$$         CALL BMG3_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, Q(p_U), K, NOG,
c$$$     &             NLx, NLy, NLz, NGx, NGy, NGz,
c$$$     &             IGRD(K,idL_BMG_ICoord),
c$$$     &             IGRD(K,idL_BMG_JCoord),
c$$$     &             IGRD(K,idL_BMG_KCoord),
c$$$     &             'UP-Q-post-interp', .FALSE.,
c$$$     &             BMG_MSG_iGRID(p_ProcGrid),
c$$$     &             NProcI, NProcJ, NProcK, NProc,
c$$$     &             MyProcI, MyProcJ, MyProcK, MPICOMM
c$$$     &             )

         IF ( BMG_IOFLAG(iBMG3_BUG_RES_INTERP) ) THEN
            CALL BMG3_SymStd_residual( 
     &                K, NOG, IFD,
     &                Q(p_U), QF(p_U), SO(p_SO), RES(p_U), 
     &                NLx, NLy, NLz, NStncl,
     &                iWorkMSG, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )
            CALL BMG3_SymStd_UTILS_norm_l2( 
     &                RES(p_U), NLx, NLy, NLz, RES_L2, MPICOMM
     &                )
            WRITE(*,110) K, RES_L2
         ENDIF

         ! Relaxation
         DO 20 i=1, IU
            !
            TT1 = MPI_Wtime()
            
            CALL BMG3_SymStd_relax(  
     &                k, KF, NOG, UPDOWN,
     &                IFD, NStncl, IRELAX, iRELAX_SYM, RES_L2, ISEP,
     &                Q, QF, RES, NF, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iGs, jGs, kGs,
     &                iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr,
     &                BMG_MSG_iGRID, NBMG_MSG_iGRID, BMG_MSG_pGRID,
     &                MPICOMM, XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
     &                )
            TT2 = MPI_Wtime()
            BMG_rPARMS(id_BMG3_TIME_relax) = 
     &           BMG_rPARMS(id_BMG3_TIME_relax) + TT2 - TT1

            !
            IF( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) .AND. MyProc.EQ.1 ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
            !
 20      CONTINUE

      ELSE
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN         
            WRITE(*,*) '*** UPDOWN out of range: '
            WRITE(*,*) 'HAVE: UPDOWN = ', UPDOWN
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,19)
         RETURN

      ENDIF

C ==========================================================================

 100  FORMAT (1X,'(3D)LEVEL ',I2,' RESIDUAL NORM = ',1P,E12.5)
 110  FORMAT (1X,'(3D)LEVEL ',I2,
     &           ' AFTER INTERPOLATION RES_L2 = ',1P,E12.5 )

C ===========================================

      RETURN
      END
