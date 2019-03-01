      SUBROUTINE BMG2_SER_SymStd_ncycle(
     &                KC, KCF, KF, IFD, IBC, IU, ID, IVW,
     &                IRELAX, IRELAX_SYM, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                ABD, BBD, NABD1, NABD2, 
     &                IGRD, NOGm, RES_L2,
     &                BMG_iPARMS
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   BMG2_SER_SymStd_ncycle performs a single multigrid n-cycle.  
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
      INCLUDE 'BMG_SER_workspace.h'
      INCLUDE 'BMG_SER_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  NABD1, NABD2, NC, NCI, NF, NOGm, NSO, NSOR

      INTEGER  IBC, ID, IFD, IGRD(NOGm,9), IRELAX, IRELAX_SYM,
     &         IU, IVW, KC, KCF, KF
      REAL*RKIND   ABD(NABD1,NABD2), BBD(NABD2), CI(NCI),
     &         Q(NF), QF(NF), RES(NF), RES_L2, SO(NSO), SOR(NSOR)
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, K, Nx, Ny, NStncl, p_CI, p_SO, p_SOR, p_U

C =========================================================================

      k = KCF   ! set the current grid index

      IF (KCF.EQ.KC) THEN
         IF (BMG_IOFLAG(iBMG2_SER_OUT_STOP_ERROR)) THEN
            WRITE(*,*)'*** FATAL ERROR: BMG2_SER_SymStd_ncycle.f ****'
            WRITE(*,*)'The current finest grid = ', KCF
            WRITE(*,*)'The coarest grid        = ', KC
            WRITE(*,*)'ERROR; minimum number of grids for n-cycle is 2!'
         END IF

         CALL BMG2_SER_SymStd_ErrTrap(BMG_iPARMS,9)
         RETURN
         
      ENDIF

      ! set the n-cycle counter to zero on all grids
      DO i = KC, KCF
         IGRD(i,idL_BMG_SER_IVW)=iZERO
      END DO

C -----------------------------------------------------------
C     Begin the n-cycle
C     (this is also the n-cycle recursive loop boundary)
C -----------------------------------------------------------

 140  CONTINUE     
         
         IF ( k.EQ.KC .AND. IGRD(k,idL_BMG_SER_IVW).EQ.iZERO ) THEN
            !
            ! Solve on the coarsest grid
            !
            CALL BMG2_SER_SymStd_GET_pointers(
     &                KC, IGRD, NOGm, Nx, Ny,
     &                p_U, p_SO, p_SOR, p_CI
     &                )

            IF ( IBC.EQ.BMG_SER_BCs_definite ) THEN
               CALL BMG2_SER_SymStd_SOLVE_cg(
     &                   Q(p_U), QF(p_U), Nx, Ny,
     &                   ABD, BBD, NABD1, NABD2,
     &                   BMG_IOFLAG, BMG_iPARMS  
     &                   )
               IF (BMG_iPARMS(id_BMG2_SER_Err_Code) .ne. iZERO) THEN
                  RETURN
               ENDIF
            ELSE
               CALL BMG2_SER_PerSymStd_SOLVE_cg( 
     &                   Q(p_U), QF(p_U), Nx, Ny, 
     &                   ABD, BBD, NABD1,NABD2, IBC 
     &                   )
            ENDIF

            IGRD(k,idL_BMG_SER_IVW) = IVW

            IF( BMG_IOFLAG(iBMG2_SER_BUG_RES_CG_SOLVE) ) THEN 
               !
               ! Residual calculation
               !
               NStncl=5
               CALL BMG2_SER_SymStd_residual(
     &                   KC, SO(p_SO), QF(p_U), Q(p_U), RES(p_U),
     &                   Nx, Ny, KF, IFD, NStncl, IBC, 
     &                   IRELAX, IRELAX_SYM, BMG_SER_DOWN
     &                   )
               CALL BMG2_SER_SymStd_UTILS_norm_l2( 
     &                   RES(p_U), Nx, Ny, RES_L2
     &                   )
               WRITE (*,100) k, RES_L2
               !
            ENDIF
            !
            ! continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_SER_IVW).LT.IVW ) THEN
            !
            ! Move down
            !
            CALL BMG2_SER_SymStd_updown(
     &                k, KF, BMG_SER_DOWN, IFD, IBC, IU, ID,
     &                IRELAX, IRELAX_SYM, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm, RES_L2, BMG_iPARMS
     &                )

            IF (BMG_iPARMS(id_BMG2_SER_Err_Code) .ne. iZERO) THEN
               RETURN
            ENDIF

            IGRD(k,idL_BMG_SER_IVW) = IGRD(k,idL_BMG_SER_IVW) + 1 
            k=k-1
            !
            ! continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_SER_IVW).EQ.IVW ) THEN
            !
            ! Move up
            !
            IGRD(k,idL_BMG_SER_IVW)=iZERO   ! reset n-cycle counter to zero
            k=k+1
            !
            CALL BMG2_SER_SymStd_updown( 
     &                k, KF, BMG_SER_UP, IFD, IBC, IU, ID,
     &                IRELAX, IRELAX_SYM, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC, 
     &                SO, NSO, SOR, NSOR, CI, NCI,
     &                IGRD, NOGm, RES_L2, BMG_iPARMS
     &                )

            IF (BMG_iPARMS(id_BMG2_SER_Err_Code) .ne. iZERO) THEN
               RETURN
            ENDIF

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

C ========================================================================

 100     FORMAT (' LEVEL',I2,' RESIDUAL NORM= ',1P,E10.3)

C ==============================

         RETURN
         END
