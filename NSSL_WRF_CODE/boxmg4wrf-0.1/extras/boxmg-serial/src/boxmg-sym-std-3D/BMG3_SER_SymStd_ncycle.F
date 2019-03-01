      SUBROUTINE BMG3_SER_SymStd_ncycle(
     &                KC, KCF, KF, 
     &                IFD, IU, ID, IVW, IRELAX, IRELAX_SYM, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC, SO, NSO, SOR, NSOR, CI, NCI,
     &                ABD, BBD, NABD1, NABD2, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C   BMG3_SER_SymStd_ncycle performs a single multigrid n-cycle.  
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
      INCLUDE 'BMG_SER_stencils.h'
      INCLUDE 'BMG_SER_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  NABD1, NABD2, NBMG_iWORK_PL, NBMG_rWORK_PL, NC, NCI,
     &         NF, NOGm, NSO, NSOR

      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS), BMG_iWORK_PL(NBMG_iWORK_PL),
     &         ID, IFD, IGRD(NOGm,NBMG_SER_pIGRD),
     &         IRELAX, IRELAX_SYM, IU, IVW, KC, KCF, KF
      REAL*RKIND   ABD(NABD1,NABD2), BBD(NABD2),
     &         BMG_rPARMS(NBMG_SER_rPARMS), BMG_rWORK_PL(NBMG_rWORK_PL),
     &         CI(NCI), Q(NF), QF(NF), RES(NF),
     &         SO(NSO), SOR(NSOR)
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, k, Nx, Ny, Nz, NStncl, p_CI, p_SO, p_SOR, p_U
      REAL*8   RES_L2 

C ==========================================================================

      k = KCF   ! set the current grid index

      IF (KCF.EQ.KC) THEN
         IF (BMG_IOFLAG(iBMG3_SER_OUT_STOP_ERROR)) THEN
            WRITE(*,*) '*** FATAL ERROR: BMG3_SER_SymStd_ncycle.f ****'
            WRITE(*,*) 'The current finest grid = ', KCF
            WRITE(*,*) 'The coarest grid        = ', KC
            WRITE(*,*) 'ERROR; minimum number of grids is 2 !'
         END IF

         CALL BMG3_SER_SymStd_ErrTrap(BMG_iPARMS,18)
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
            CALL BMG3_SER_SymStd_GET_pointers( 
     &                KC, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, Nx, Ny, Nz 
     &                )
            CALL BMG3_SER_SymStd_SOLVE_cg( 
     &                q(p_U), qf(p_U), Nx, Ny, Nz, 
     &                ABD, BBD, NABD1, NABD2,
     &                BMG_IOFLAG, BMG_iPARMS 
     &                )
            IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF
            

            IGRD(k,idL_BMG_SER_IVW) = IVW
            IF( BMG_IOFLAG(iBMG3_SER_BUG_RES_CG_SOLVE) ) THEN
               NStncl = 14      ! NOG > 1 so CG has 27-point stencil
               CALL BMG3_SER_SymStd_residual(
     &                   KC, KF, ifd, 
     &                   q(p_U), qf(p_U), so(p_SO), RES(p_U),
     &                   Nx, Ny, Nz, NStncl
     &                   )
               CALL BMG3_SER_SymStd_UTILS_norm_l2( 
     &                   RES(p_U), Nx, Ny, Nz, RES_L2
     &                   )
               WRITE (*,100) k, RES_L2
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
            CALL BMG3_SER_SymStd_updown(
     &                k, KF, BMG_SER_DOWN, 
     &                IFD, IU, ID, IRELAX, IRELAX_SYM,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC, SO, NSO,
     &                SOR, NSOR, CI, NCI, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )
            IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
               RETURN
            END IF

            IGRD(k,idL_BMG_SER_IVW) = IGRD(k,idL_BMG_SER_IVW) + 1 
            k=k-1
            !
            ! Continue the n-cycle
            !
            GOTO 140
            !
         ELSE IF ( IGRD(k,idL_BMG_SER_IVW).EQ.IVW ) THEN
            !
            ! Move up

            !
            IGRD(k,idL_BMG_SER_IVW) = iZERO   ! reset n-cycle counter to zero
            k=k+1
            !
            CALL BMG3_SER_SymStd_updown(
     &                k, KF, BMG_SER_UP,
     &                IFD, IU, ID, IRELAX, IRELAX_SYM,
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Q, QF, RES, NF, NC, SO, NSO,
     &                SOR, NSOR, CI, NCI, IGRD, NOGm,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )
            IF (BMG_iPARMS(id_BMG3_SER_Err_Code) .ne. iZERO) THEN
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

C ==============================

         RETURN
         END
