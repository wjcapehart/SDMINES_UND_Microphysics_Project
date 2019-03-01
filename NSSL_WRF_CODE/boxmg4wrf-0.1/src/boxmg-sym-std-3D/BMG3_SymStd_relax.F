      SUBROUTINE BMG3_SymStd_relax( 
     &                KG, KF, NOG, UPDOWN,
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

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_relax performs one relaxation step (either colored
C     point Gauss Seidel, or alternating plane red-black Gauss-Seidel)
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

C ------------------------------------------------
C     Includes
C
#include    "BMG_constants.h"
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ------------------------------------------------
C     Argument Declarations
C 
      INTEGER  NCI, NSO, NSOR, NSORv, NStncl, kg, kf, NBMG_iWORK_PL,
     &         NBMG_rWORK_PL, NOG, NF, NOGm, NMSGi, NMSGr

      INTEGER  NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRID),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER  BMG_iPARMS(NBMG_iPARMS), BMG_iWORK_PL(NBMG_iWORK_PL),
     &         IFD, IGRD(NOGm,24), IRELAX_SYM, UPDOWN, IRELAX, ISEP
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS), BMG_rWORK_PL(NBMG_rWORK_PL),
     &         CI(NCI), q(NF), qf(NF), RES(NF),
     &         SO(NSO), SOR(NSOR), MSG_Buffer(NMSGr)

      REAL*8   RES_L2
      
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  iGs, jGs, kGs,
     &         MPICOMM, iWorkMSG(NMSGi), pMSG(NBMG_pMSG,NOGm)
      INTEGER  XYPLANECOMM, XZPLANECOMM, YZPLANECOMM
      INTEGER  MyProc, MyProcI, MyProcJ, MyProcK, MPI_IERROR,
     &         NProc, NProcI, NProcJ, NProcK


C ----------------------------
C     Local Declarations
C
      INTEGER  NLx, NLy, NLz, NGx, NGy, NGz, 
     &         NLx_f, NLy_f, NLz_f, NGx_f, NGy_f, NGz_f, 
     &         p_CI, p_SO, p_SOR, p_U, p_ProcGrid
      
C ==========================================================================

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>> BEGIN:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

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

      p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)

C --------------------------------------------------------------------
C >>>>>>>>>>>>>>>>>>>> END:  UNPACK PARAMETERS <<<<<<<<<<<<<<<<<<<<<<<
C --------------------------------------------------------------------

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE 
         NSORv = 2  ! the same for now but this may change
      ENDIF

      !
      ! Collect fine-grid dimensions 
      ! ( need to check: KF or NOG )
      !
      CALL BMG3_SymStd_GET_pointers( 
     &            KF, IGRD, NOGm,
     &            p_U, p_SO, p_SOR, p_CI, 
     &            NLx_f, NLy_f, NLz_f, NGx_f, NGy_f, NGz_f
     &            ) 

      CALL BMG3_SymStd_GET_pointers(  
     &            KG, IGRD, NOGm,
     &            p_U, p_SO, p_SOR, p_CI, 
     &            NLx, NLy, NLz, NGx, NGy, NGz
     &            )
      

      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         !
         ! Gauss-Seidel relaxation
         !
         CALL BMG3_SymStd_relax_GS(  
     &             KG, SO(p_SO), QF(p_U), Q(p_U), RES(p_U), SOR(p_SOR),
     &             NLx, NLy, NLz, NGx, NGy, NGz, 
     &             RES_L2, BMG_IOFLAG, NOG, NOGm,
     &             IFD, NStncl, NSORv, IRELAX, IRELAX_SYM, UPDOWN,
     &             iGs, jGs, kGs,
     &             MyProcI, MyProcJ, MyProcK, MyProc,
     &             BMG_MSG_iGRID(p_ProcGrid),
     &             NProcI, NProcJ, NProcK, NProc,
     &             iWorkMSG, NMSGi, pMSG, MSG_Buffer, NMSGr, MPICOMM
     &            ,ISEP)
         !
      ELSEIF ( IRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz ) THEN

         WRITE(*,*) 'Fatal Error: To busy to code planes !!!'
         CALL MPI_FINALIZE(MPI_IERROR)
         STOP

c$$$         !
c$$$         ! Alternating plane relaxation
c$$$         !
c$$$         IF ( ( IRELAX_SYM .EQ.BMG_RELAX_NONSYM )
c$$$     &      .OR. ( UPDOWN.EQ.BMG_DOWN )         )THEN
c$$$            !
c$$$            ! On the way down, or in the non-symmetric case, 
c$$$            ! relax in the order xy, yz, xz
c$$$            !
c$$$            CALL BMG3_SymStd_relax_planes_xy(
c$$$     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
c$$$     &                NLx_f, NLy_f, NLz_f, kg, IGRD, Q, QF, RES, NF,
c$$$     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
c$$$     &                NOGm, NOG, UPDOWN, RES_L2,
c$$$     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xy)),
c$$$     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xy)),
c$$$     &                BMG_iWORK_PL, NBMG_iWORK_PL,
c$$$     &                BMG_rWORK_PL, NBMG_rWORK_PL
c$$$     &                )
c$$$            !
c$$$            CALL BMG3_SymStd_relax_planes_yz(
c$$$     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
c$$$     &                NLx_f, NLy_f, NLz_f, kg, IGRD, Q, QF, RES, NF,
c$$$     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
c$$$     &                NOGm, NOG, UPDOWN, RES_L2,
c$$$     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_yz)),
c$$$     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_yz)),
c$$$     &                BMG_iWORK_PL, NBMG_iWORK_PL,
c$$$     &                BMG_rWORK_PL, NBMG_rWORK_PL
c$$$     &                )
c$$$            !
c$$$            CALL BMG3_SymStd_relax_planes_xz(
c$$$     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
c$$$     &                NLx_f, NLy_f, NLz_f, kg, IGRD, Q, QF, RES, NF,
c$$$     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
c$$$     &                NOGm, NOG, UPDOWN, RES_L2,
c$$$     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xz)),
c$$$     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xz)),
c$$$     &                BMG_iWORK_PL, NBMG_iWORK_PL,
c$$$     &                BMG_rWORK_PL, NBMG_rWORK_PL
c$$$     &                )
c$$$            !
c$$$         ELSE
c$$$            !
c$$$            ! on the way up, in the symmetric case, use
c$$$            ! the opposite order xz, yz, xy
c$$$            !
c$$$            CALL BMG3_SymStd_relax_planes_xz(
c$$$     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
c$$$     &                NLx_f, NLy_f, NLz_f, kg, IGRD, Q, QF, RES, NF,
c$$$     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
c$$$     &                NOGm, NOG, UPDOWN, RES_L2,
c$$$     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xz)),
c$$$     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xz)),
c$$$     &                BMG_iWORK_PL, NBMG_iWORK_PL,
c$$$     &                BMG_rWORK_PL, NBMG_rWORK_PL
c$$$     &                )
c$$$            !
c$$$            CALL BMG3_SymStd_relax_planes_yz(
c$$$     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
c$$$     &                NLx_f, NLy_f, NLz_f, kg, IGRD, Q, QF, RES, NF,
c$$$     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
c$$$     &                NOGm, NOG, UPDOWN, RES_L2,
c$$$     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_yz)),
c$$$     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_yz)),
c$$$     &                BMG_iWORK_PL, NBMG_iWORK_PL,
c$$$     &                BMG_rWORK_PL, NBMG_rWORK_PL
c$$$     &                )
c$$$            !
c$$$            CALL BMG3_SymStd_relax_planes_xy(
c$$$     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
c$$$     &                NLx_f, NLy_f, NLz_f, kg, IGRD, Q, QF, RES, NF,
c$$$     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
c$$$     &                NOGm, NOG, UPDOWN, RES_L2,
c$$$     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xy)),
c$$$     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xy)),
c$$$     &                BMG_iWORK_PL, NBMG_iWORK_PL,
c$$$     &                BMG_rWORK_PL, NBMG_rWORK_PL
c$$$     &                )
c$$$            !
c$$$         ENDIF
         !
      ENDIF

C ==========================================================================

      RETURN
      END









