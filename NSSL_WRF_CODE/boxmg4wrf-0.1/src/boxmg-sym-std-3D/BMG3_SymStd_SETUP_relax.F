      SUBROUTINE BMG3_SymStd_SETUP_relax( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                IFD, IBC, IRELAX, NOGm, NOG, IGRD,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL,
     &                iWork, NMSGi, pMSG, BUFFER, NMSGr, 
     &                NProc, MyProc, MyProcI, MyProcJ, MyProcK,
     &                MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_SETUP_relax.f performs any necessary setup for the
C     chosen relaxation scheme.  In particular, it calls routines that
C     setup the operators for plane relaxation.
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
#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ---------------------------
C    Argument Declarations:
C
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)


      INTEGER  NBMG_iWORK_PL, NBMG_rWORK_PL
      INTEGER  BMG_iWORK_PL(NBMG_iWORK_PL)
      REAL*RKIND   BMG_rWORK_PL(NBMG_rWORK_PL)

      INTEGER  NCI, NMSGi, NMSGr, NOGm, NSO, NSOR
      INTEGER  IBC, IFD, IGRD(24*NOGm), IRELAX, NOG
      REAL*RKIND   CI(NCI), SO(NSO), SOR(NSOR)

      INTEGER  iWork(NMSGi), MPICOMM, MyProc, MyProcI, MyProcJ, MyProcK,
     &         Nproc, pMSG(NBMG_pMSG,NOG)
      REAL*RKIND   BUFFER(NMSGr)

C ---------------------------
C    Local Declarations:
C
      INTEGER  K, MPI_IERROR, NGx, NGy, NGz, NLx, NLy, NLz,
     &         NSORv, NStncl, p_CI, p_SO, p_SOR, p_U

C ==========================================================================
      
      !
      ! Sanity check
      !
      IF ( NOG.EQ.1) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN 
            WRITE(*,500) 'There is only 1 grid!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,6)
         RETURN

      ELSE IF ( NOG.EQ.0 ) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN 
            WRITE(*,500) 'There are no grids?'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF
         
         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,7)
         RETURN
         
      ELSE IF ( NOG.LT.0 ) THEN
         IF (BMG_IOFLAG(iBMG3_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN 
            WRITE(*,500) 'There number of grids is negative!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,8)
         RETURN
         
      ENDIF

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         NSORv = 2
      ELSE 
         NSORv = 2  ! the same for now but this may change
      ENDIF

      IF (IRELAX.EQ.BMG_GS_RB_point) THEN

         !
         ! Loop over grids
         !
         DO K = NOG, 2, -1
            
            IF ( K.NE.NOG .OR. IFD.NE.BMG_STENCIL_7pt ) THEN
               NStncl=14
            ELSE
               NStncl=4
            ENDIF

            CALL BMG3_SymStd_GET_pointers( 
     &                       K, IGRD, NOGm,
     &                       p_U, p_SO, p_SOR, p_CI,
     &                       NLx, NLy, NLz, NGx, NGy, NGz 
     &                       )

            CALL BMG3_SymStd_SETUP_recip(
     &                SO(p_SO), SOR(p_SOR), NLx, NLy, NLz, Nstncl, NSORv
     &                )
            
         ENDDO

      ELSEIF (IRELAX.EQ.BMG_GS_RB_planes_xy_yz_xz) THEN

         WRITE(*,*) 'Fatal Error: To busy to code planes !!!'
         CALL MPI_FINALIZE(MPI_IERROR)
         STOP

c$$$         !
c$$$         ! Collect fine-grid dimensions (Nxf,Nyf,Nzf)
c$$$         !
c$$$         CALL BMG3_SymStd_GET_pointers( NOG, IGRD, NOGm,
c$$$     &               p_U, p_SO, p_SOR, p_CI, Nxf, Nyf, Nzf )
c$$$
c$$$         CALL BMG3_SymStd_SETUP_planes_xy( 
c$$$     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
c$$$     &             Nxf, Nyf, Nzf, IFD, IBC, IGRD, 
c$$$     &             SO, NSO, NOGm, NOG,
c$$$     &             BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xy)),
c$$$     &             BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xy)),
c$$$     &             BMG_iWORK_PL, NBMG_iWORK_PL,
c$$$     &             BMG_rWORK_PL, NBMG_rWORK_PL
c$$$     &             )
c$$$
c$$$         CALL BMG3_SymStd_SETUP_planes_yz( 
c$$$     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
c$$$     &             Nxf, Nyf, Nzf, IFD, IBC, IGRD, 
c$$$     &             SO, NSO, NOGm, NOG,
c$$$     &             BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_yz)),
c$$$     &             BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_yz)),
c$$$     &             BMG_iWORK_PL, NBMG_iWORK_PL,
c$$$     &             BMG_rWORK_PL, NBMG_rWORK_PL
c$$$     &             )
c$$$
c$$$         CALL BMG3_SymStd_SETUP_planes_xz( 
c$$$     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
c$$$     &             Nxf, Nyf, Nzf, IFD, IBC, IGRD, 
c$$$     &             SO, NSO, NOGm, NOG,
c$$$     &             BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xz)),
c$$$     &             BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xz)),
c$$$     &             BMG_iWORK_PL, NBMG_iWORK_PL,
c$$$     &             BMG_rWORK_PL, NBMG_rWORK_PL
c$$$     &             )

      ENDIF
      
C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SETUP_relax.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END


