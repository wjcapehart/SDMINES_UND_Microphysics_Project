      SUBROUTINE BMG2_SER_SymStd_relax(
     &                K, SO, QF, Q, RES, SOR, NSOR, Nx, Ny, 
     &                KF, IFD, NStncl, IBC, IRELAX, IRELAX_SYM, UPDOWN,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_SymStd_relax performs point or line relaxation, by calling
C     appropriate subroutines based on the value of IRELAX.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     K         index of the current grid
C     KF        index of the finest grid
C
C     Nx        Number of grid points in x direction, including
C               two fictitious points.
C     Ny        Number of grid points in y direction, including
C               two fictitious points.
C
C     SO        Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     QF        Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C     SOR       Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C
C     IFD       Stencil parameter
C     IRELAX    Refer to BMG2_SER_SymStd_SOLVE_boxmg.
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     Q         Refer to BMG2_SER_SymStd_SOLVE_boxmg.
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
      INCLUDE 'BMG_SER_stencils.h'
      INCLUDE 'BMG_SER_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER Nx, Ny, NSOR, NStncl

      INTEGER IBC, IFD, IRELAX, IRELAX_SYM, K, KF, UPDOWN
      REAL*RKIND  Q(Nx*Ny), QF(Nx*Ny), RES(Nx*Ny),
     &        SO(Nx*Ny*NStncl), SOR(NSOR)
      LOGICAL BMG_IOFLAG(NBMG_SER_IOFLAG)
      INTEGER BMG_iPARMS(NBMG_SER_iPARMS)

C ----------------------------
C     Local Declarations
C 
      INTEGER p_SOR, NSORv

C =========================================================================

C     NB: Workspace double time
C
C     RES is used for temporary workspace in the y-lines case.
C     It is included in the call to x-lines for consistency of
C     the calling sequence only.
C

C     -------------------------------------------------------------------
      
      !
      !  Number of temporary vectors in SOR
      !
      IF ( iRELAX.EQ.BMG_SER_GS_RB_point 
     &    .OR. iRELAX.EQ.BMG_SER_GS_RB_x_lines ) THEN
         NSORv = 2
      ELSE 
         NSORv = 4
      ENDIF


      IF ( IRELAX.EQ.BMG_SER_GS_RB_point ) THEN
         !
         CALL BMG2_SER_SymStd_relax_GS( 
     &             K, SO, QF, Q, SOR, Nx, Ny, 
     &             KF, IFD, NStncl, NSORv, IRELAX_SYM, UPDOWN 
     &             )
         !
      ELSE IF ( IRELAX.EQ.BMG_SER_GS_RB_x_lines ) THEN
         !
         p_SOR = (ipL_BMG_SER_LUL1-1)*Nx*Ny + 1
         CALL BMG2_SER_SymStd_relax_lines_x ( 
     &             K, SO, QF, Q, SOR(p_SOR), RES, 
     &             Nx, Ny, KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &             )
         !
      ELSE IF ( IRELAX.EQ.BMG_SER_GS_RB_y_lines ) THEN
         !
         p_SOR = (ipL_BMG_SER_LUL2-1)*Nx*Ny + 1
         CALL BMG2_SER_SymStd_relax_lines_y ( 
     &             K, SO, QF, Q, SOR(p_SOR), RES, 
     &             Nx, Ny, KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &             )
         !
      ELSE IF ( IRELAX.EQ.BMG_SER_GS_RB_x_y_lines ) THEN
         ! 
         IF ( IRELAX_SYM.EQ.BMG_SER_RELAX_NONSYM
     &       .OR.( UPDOWN.EQ.BMG_SER_DOWN
     &            .AND.IRELAX_SYM.EQ.BMG_SER_RELAX_SYM ) ) THEN
            !
            p_SOR = (ipL_BMG_SER_LUL1-1)*Nx*Ny + 1
            CALL BMG2_SER_SymStd_relax_lines_x ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES, 
     &                Nx, Ny, KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &                )
            p_SOR = (ipL_BMG_SER_LUL2-1)*Nx*Ny + 1
            CALL BMG2_SER_SymStd_relax_lines_y ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES, 
     &                Nx, Ny, KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &                )
            !
         ELSE ! on the way up relax in the opposite order
            !
            p_SOR = (ipL_BMG_SER_LUL2-1)*Nx*Ny + 1
            CALL BMG2_SER_SymStd_relax_lines_y ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES, 
     &                Nx, Ny, KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &                )
            p_SOR = (ipL_BMG_SER_LUL1-1)*Nx*Ny + 1
            CALL BMG2_SER_SymStd_relax_lines_x ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES, 
     &                Nx, Ny, KF, IFD, NStncl, IRELAX_SYM, UPDOWN 
     &                )
            !
         ENDIF
         !
      ENDIF

      RETURN
      END
