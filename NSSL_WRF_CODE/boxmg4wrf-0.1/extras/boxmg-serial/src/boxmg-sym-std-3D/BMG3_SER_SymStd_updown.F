      SUBROUTINE BMG3_SER_SymStd_updown(
     &                K, KF, UPDOWN, IFD, IU, ID,
     &                IRELAX, IRELAX_SYM,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG, 
     &                Q, QF, RES, NF, NC, SO, NSO,
     &                SOR, NSOR, CI, NCI, IGRD, NOGm, 
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_updown performs all the necessary tasks associated
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
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_stencils.h'
      INCLUDE 'BMG_SER_parameters.h'
      INCLUDE 'BMG_SER_workspace.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  NBMG_iWORK_PL, NBMG_rWORK_PL, NC, NCI,
     &         NF, NOGm, NSO, NSOR

      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS), BMG_iWORK_PL(NBMG_iWORK_PL),
     &         ID, IFD, IGRD(NOGm,NBMG_SER_pIGRD),
     &         IRELAX, IRELAX_SYM, IU, K, KF, UPDOWN
      REAL*RKIND   BMG_rPARMS(NBMG_SER_rPARMS), 
     &         BMG_rWORK_PL(NBMG_rWORK_PL),
     &         CI(NCI), Q(NF), QF(NF), RES(NF),
     &         SO(NSO), SOR(NSOR)
      REAL*8   RES_L2
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, NSORv, NStncl, Nx, Nxc, Ny, Nyc, Nz, Nzc,
     &         p_CI, p_CIC, p_SO, p_SOC, p_SOR, p_SORC, p_U, p_UC

C ==========================================================================

      !
      ! How many points in the stencil?
      !
      IF ( K.EQ.KF .AND. IFD.EQ.1 ) THEN
         NStncl=4
      ELSE
         NStncl=14
      ENDIF

      !
      !  Number of temporary vectors in SOR
      !
      IF ( IRELAX.EQ.BMG_SER_GS_RB_point ) THEN
         NSORv = 2
      ELSE 
         NSORv = 2  ! the same for now but this may change
      ENDIF

      !
      ! Get pointers for grids k and k-1
      !
      CALL BMG3_SER_SymStd_GET_pointers( 
     &          k, IGRD, NOGm,
     &          p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &          )

      CALL BMG3_SER_SymStd_GET_pointers( 
     &          k-1, IGRD, NOGm,
     &          p_UC, p_SOC, p_SORC, p_CIC, Nxc, Nyc, Nzc 
     &          )


      IF (UPDOWN.EQ.BMG_SER_DOWN) THEN

         ! Relaxation
         DO 10 i=1, ID
            CALL BMG3_SER_SymStd_relax( 
     &                k, kf, IGRD, 
     &                q, qf, RES, NF, so, NSO,
     &                sor, NSOR, ci, NCI,
     &                ifd, NStncl, IRELAX, iRELAX_SYM, 
     &                NOGm, KF, UPDOWN, RES_L2,   
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )
            IF( BMG_IOFLAG(iBMG3_SER_BUG_RES_RELAX) ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
 10      CONTINUE

         ! calculate the residual
         CALL BMG3_SER_SymStd_residual(
     &             k, KF, IFD, Q(p_U), QF(p_U), SO(p_SO),
     &             RES(p_U), Nx, Ny, Nz, NStncl
     &             )
         
c$$$         CALL BMG3_SER_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, RES(p_U), Nx, Ny, Nz, K, NOGm,
c$$$     &             'output', 'DOWN-RES-post-relax', .FALSE.
c$$$     &             )


         ! restrict the residual
         CALL BMG3_SER_SymStd_restrict(
     &             k, k-1,  
     &             RES(p_U), QF(p_UC), CI(p_CIC),
     &             Nx, Ny, Nz, Nxc, Nyc, Nzc
     &             )

c$$$         CALL BMG3_SER_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, QF(p_UC), Nxc, Nyc, Nzc, K-1, NOGm,
c$$$     &             'output', 'DOWN-QF-post-restrict', .FALSE.
c$$$     &             )


         ! Zero the initial guess
         CALL BMG3_SER_SymStd_UTILS_rV_zero(
     &             q(p_UC), Nxc, Nyc, Nzc
     &             )

      ELSE IF (UPDOWN.EQ.BMG_SER_UP) THEN

c$$$         CALL BMG3_SER_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, Q(p_UC), Nxc, Nyc, Nzc, K-1, NOGm,
c$$$     &             'output', 'UP-Q-pre-interp', .FALSE.
c$$$     &             )
 
         ! Interpolate and Correct 
         CALL BMG3_SER_SymStd_interp_add(
     &             k-1, k, Q(p_U) ,Q(p_UC),
     &             SO(p_SO), RES(p_U), CI(p_CIC),
     &             Nxc, Nyc, Nzc, Nx, Ny, Nz, NStncl
     &             )
         
c$$$         CALL BMG3_SER_SymStd_DUMP_vector( 
c$$$     &             BMG_IOFLAG, Q(p_U), Nx, Ny, Nz, K, NOGm,
c$$$     &             'output', 'UP-Q-post-interp', .FALSE.
c$$$     &             )

         IF ( BMG_IOFLAG(iBMG3_SER_BUG_RES_INTERP) ) THEN
            CALL BMG3_SER_SymStd_residual(
     &                k, KF, IFD, Q(p_U), QF(p_U), SO(p_SO),
     &                RES(p_U), Nx, Ny, Nz, NStncl
     &                )
            CALL BMG3_SER_SymStd_UTILS_norm_l2( 
     &                RES(p_U), Nx, Ny, Nz, RES_L2
     &                )
            WRITE(*,110) K, RES_L2
         ENDIF

         ! Relaxation
         DO 20 i=1, IU
            CALL BMG3_SER_SymStd_relax(
     &                k, kf, IGRD, 
     &                q, qf, RES, NF, so, NSO,
     &                sor, NSOR, ci, NCI,
     &                ifd, NStncl, IRELAX, iRELAX_SYM, 
     &                NOGm, KF, UPDOWN, RES_L2,   
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )
            IF( BMG_IOFLAG(iBMG3_SER_BUG_RES_RELAX) ) THEN
               WRITE(*,100) K, RES_L2
            ENDIF
 20      CONTINUE

      ELSE
         IF (BMG_IOFLAG(iBMG3_SER_OUT_STOP_ERROR)) THEN         
            WRITE(*,500) '*** UPDOWN out of range: '
            WRITE(*,510) 'HAVE: UPDOWN = ', UPDOWN
         END IF

         CALL BMG3_SER_SymStd_ErrTrap(BMG_iPARMS,19)
         RETURN

      ENDIF

C ==========================================================================

 100  FORMAT (1X,'(3D)LEVEL ',I2,' RESIDUAL NORM = ',1P,E12.5)
 110  FORMAT (1X,'(3D)LEVEL ',I2,
     &           ' AFTER INTERPOLATION RES_L2 = ',1P,E12.5 )

C -----------------------------                                                 
 500    FORMAT (/,'FATAL ERROR: BMG3_SER_SymStd_updown.f',/,5X,A)
 510    FORMAT (5X,A,1X,I2)

C =============================

      RETURN
      END
