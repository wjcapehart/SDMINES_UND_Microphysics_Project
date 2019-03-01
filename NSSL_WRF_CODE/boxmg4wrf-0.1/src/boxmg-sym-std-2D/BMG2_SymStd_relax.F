      SUBROUTINE BMG2_SymStd_relax(
     &                K, KF, NOG, UPDOWN, 
     &                SO, QF, Q, RES, SOR, 
     &                II, JJ, iGs, jGs, BMG_iPARMS,
     &                IFD, NStncl, IRELAX, IRELAX_SYM,
     &                iWork, NMSGi, pMSG, BF, NMSGr,
     &                XDATADIST, YDATADIST,
     &                MPICOMM, XLINECOMM, YLINECOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_relax performs point or line relaxation, by calling
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
C     SO        Refer to BMG2_SymStd_SOLVE_boxmg.
C     QF        Refer to BMG2_SymStd_SOLVE_boxmg.
C     SOR       Refer to BMG2_SymStd_SOLVE_boxmg.
C
C     IFD       Stencil parameter
C     IRELAX    Refer to BMG2_SymStd_SOLVE_boxmg.
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     Q         Refer to BMG2_SymStd_SOLVE_boxmg.
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

      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER II, JJ, NOG, NStncl, NMSGi, NMSGr

      INTEGER iWork(NMSGi), pMSG(NBMG_pMSG,NOG)
      INTEGER MPICOMM, XLINECOMM, YLINECOMM
      INTEGER XDATADiST(2,*), YDATADIST(2,*)

      INTEGER KF, IFD, IRELAX, IRELAX_SYM, K, UPDOWN, iGs, jGs
      REAL*RKIND  Q(II*JJ), QF(II*JJ), RES(II*JJ), 
     &        SO((II+1)*(JJ+1)*NStncl), SOR(*),
     &        BF(NMSGr)

      INTEGER  BMG_iPARMS(NBMG_iPARMS)

C ----------------------------
C     Local Declarations
C 
      INTEGER p_SOR

C =========================================================================

C     NB: Workspace double time
C
C     RES is used for temporary workspace in the y-lines case.
C     It is included in the call to x-lines for consistency of
C     the calling sequence only.
C
C     -------------------------------------------------------------------

      IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
         !
         CALL BMG2_SymStd_relax_GS( 
     &             K, SO, QF, Q, SOR, II, JJ, 
     &             KF, IFD, NStncl, IRELAX_SYM, UPDOWN,
     &             iGs, jGs, iWork, NMSGi,
     &             pMSG, BF, NMSGr, MPICOMM
     &             )
         !
      ELSE IF ( IRELAX.EQ.BMG_GS_RB_x_lines ) THEN
         !
         p_SOR = (ipL_BMG_LUL1-1)*II*JJ + 1
         CALL BMG2_SymStd_relax_lines_x ( 
     &             K, SO, QF, Q, SOR(p_SOR), RES,
     &             II, JJ, iGs, jGs, BMG_iPARMS,
     &             NOG, NStncl, IRELAX_SYM, UPDOWN,
     &             XDATADiST, iWork, NMSGi, pMSG, BF, NMSGr, 
     &             MPICOMM, XLINECOMM, YLINECOMM
     &             )
         !
      ELSE IF ( IRELAX.EQ.BMG_GS_RB_y_lines ) THEN
         !
         p_SOR = (ipL_BMG_LUL2-1)*II*JJ + 1
         CALL BMG2_SymStd_relax_lines_y ( 
     &             K, SO, QF, Q, SOR(p_SOR), RES, 
     &             II, JJ, iGs, jGs, BMG_iPARMS,
     &             NOG, NStncl, IRELAX_SYM, UPDOWN,
     &             YDATADIST, iWork, NMSGi, pMSG, BF, NMSGr, 
     &             MPICOMM, XLINECOMM, YLINECOMM
     &             )
         !
      ELSE IF ( IRELAX.EQ.BMG_GS_RB_x_y_lines ) THEN
         ! 
         IF ( IRELAX_SYM.EQ.BMG_RELAX_NONSYM ) THEN
            !
            p_SOR = (ipL_BMG_LUL1-1)*II*JJ + 1
            CALL BMG2_SymStd_relax_lines_x ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES,
     &                II, JJ, iGs, jGs, BMG_iPARMS,
     &                NOG, NStncl, IRELAX_SYM, UPDOWN,
     &                XDATADiST, iWork, NMSGi, pMSG, BF, NMSGr, 
     &                MPICOMM, XLINECOMM, YLINECOMM 
     &                )
            p_SOR = (ipL_BMG_LUL2-1)*II*JJ + 1
            CALL BMG2_SymStd_relax_lines_y ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES, 
     &                II, JJ, iGs, jGs, BMG_iPARMS,
     &                NOG, NStncl, IRELAX_SYM, UPDOWN,
     &                YDATADIST, iWork, NMSGi, pMSG, BF, NMSGr, 
     &                MPICOMM, XLINECOMM, YLINECOMM
     &                )
         ELSEIF (UPDOWN.EQ.BMG_DOWN ) THEN
            !
            p_SOR = (ipL_BMG_LUL1-1)*II*JJ + 1
            CALL BMG2_SymStd_relax_lines_x ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES,
     &                II, JJ, iGs, jGs, BMG_iPARMS,
     &                NOG, NStncl, IRELAX_SYM, UPDOWN,
     &                XDATADiST, iWork, NMSGi, pMSG, BF, NMSGr, 
     &                MPICOMM, XLINECOMM, YLINECOMM
     &                )
            p_SOR = (ipL_BMG_LUL2-1)*II*JJ + 1
            CALL BMG2_SymStd_relax_lines_y ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES, 
     &                II, JJ, iGs, jGs, BMG_iPARMS,
     &                NOG, NStncl, IRELAX_SYM, UPDOWN,
     &                YDATADIST, iWork, NMSGi, pMSG, BF, NMSGr, 
     &                MPICOMM, XLINECOMM, YLINECOMM
     &                )
            !
         ELSE  ! on the way up relax in the opposite order
            !
            p_SOR = (ipL_BMG_LUL2-1)*II*JJ + 1
            CALL BMG2_SymStd_relax_lines_y ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES, 
     &                II, JJ, iGs, jGs, BMG_iPARMS,
     &                NOG, NStncl, IRELAX_SYM, UPDOWN,
     &                YDATADIST, iWork, NMSGi, pMSG, BF, NMSGr, 
     &                MPICOMM, XLINECOMM, YLINECOMM
     &                )
            p_SOR = (ipL_BMG_LUL1-1)*II*JJ + 1
            CALL BMG2_SymStd_relax_lines_x ( 
     &                K, SO, QF, Q, SOR(p_SOR), RES,
     &                II, JJ, iGs, jGs, BMG_iPARMS,
     &                NOG, NStncl, IRELAX_SYM, UPDOWN,
     &                XDATADiST, iWork, NMSGi, pMSG, BF, NMSGr, 
     &                MPICOMM, XLINECOMM, YLINECOMM 
     &                )
            !
         ENDIF
         !
      ENDIF

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_relax.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C =====================

      RETURN
      END
