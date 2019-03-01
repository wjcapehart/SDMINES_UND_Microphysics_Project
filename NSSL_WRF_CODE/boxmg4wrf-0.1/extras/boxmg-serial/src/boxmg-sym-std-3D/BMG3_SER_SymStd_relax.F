      SUBROUTINE BMG3_SER_SymStd_relax( 
     &                kg, kf, IGRD,
     &                Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IFD, NStncl, IRELAX, IRELAX_SYM,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_relax performs one relaxation step (either colored
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
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_workspace.h'
      INCLUDE 'BMG_SER_parameters.h'

C ------------------------------------------------
C     Argument Declarations
C 
      INTEGER  NCI, NSO, NSOR, kg, kf, NBMG_iWORK_PL,
     &         NBMG_rWORK_PL, NOG, NF, NOGm, NStncl

      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS),
     &         BMG_iWORK_PL(NBMG_iWORK_PL),
     &         IFD, IGRD(NOGm,NBMG_SER_pIGRD),
     &         IRELAX_SYM, UPDOWN, IRELAX
      REAL*RKIND   BMG_rPARMS(NBMG_SER_rPARMS),
     &         BMG_rWORK_PL(NBMG_rWORK_PL),
     &         CI(NCI), q(NF), qf(NF), RES(NF),
     &         SO(NSO), SOR(NSOR)
      REAL*8   RES_L2
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  NSORv, Nx, Nxf, Ny, Nyf, Nz, Nzf, p_CI, p_SO, p_SOR, p_U
      
C ==========================================================================

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_SER_GS_RB_point ) THEN
         NSORv = 2
      ELSE 
         NSORv = 2  ! the same for now but this may change
      ENDIF

      !
      ! Collect fine-grid dimensions
      !
      CALL BMG3_SER_SymStd_GET_pointers( 
     &          NOG, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CI, Nxf, Nyf, Nzf 
     &          )

      CALL BMG3_SER_SymStd_GET_pointers( 
     &          kg, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &          )
      

      IF ( IRELAX.EQ.BMG_SER_GS_RB_point ) THEN
         !
         ! Gauss-Seidel relaxation
         !
         CALL BMG3_SER_SymStd_relax_GS( 
     &             kg, SO(p_SO), qf(p_U), q(p_U), RES(p_U), SOR(p_SOR), 
     &             Nx, Ny, Nz, RES_L2, BMG_IOFLAG,
     &             KF, IFD, NStncl, NSORv, IRELAX, IRELAX_SYM, UPDOWN
     &             )
         !
      ELSEIF ( IRELAX.EQ.BMG_SER_GS_RB_planes_xy_yz_xz ) THEN

         !
         ! Alternating plane relaxation
         !
         IF ( ( IRELAX_SYM .EQ.BMG_SER_RELAX_NONSYM )
     &      .OR. ( UPDOWN.EQ.BMG_SER_DOWN )         )THEN
            !
            ! On the way down, or in the non-symmetric case, 
            ! relax in the order xy, yz, xz
            !
            CALL BMG3_SER_SymStd_relax_planes_xy(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xy)),
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xy)),
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )
            !
            CALL BMG3_SER_SymStd_relax_planes_yz(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_yz)),
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_yz)),
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )
            !
            CALL BMG3_SER_SymStd_relax_planes_xz(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xz)),
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xz)),
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )
            !
         ELSE
            !
            ! on the way up, in the symmetric case, use
            ! the opposite order xz, yz, xy
            !
            CALL BMG3_SER_SymStd_relax_planes_xz(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xz)),
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xz)),
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )
            !
            CALL BMG3_SER_SymStd_relax_planes_yz(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_yz)),
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_yz)),
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )
            !
            CALL BMG3_SER_SymStd_relax_planes_xy(
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_iPARMS_PL_xy)),
     &                BMG_iWORK_PL(BMG_iWORK_PL(ip_BMG_pWORK_PL_xy)),
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )
            !
         ENDIF
         !
      ENDIF

C ==========================================================================

      RETURN
      END









