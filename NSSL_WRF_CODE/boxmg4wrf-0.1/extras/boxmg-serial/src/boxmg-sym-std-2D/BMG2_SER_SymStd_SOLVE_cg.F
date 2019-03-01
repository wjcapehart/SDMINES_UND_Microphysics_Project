      SUBROUTINE BMG2_SER_SymStd_SOLVE_cg(
     &                       Q, QF, II, JJ, ABD, BBD, NABD1, NABD2,
     &                       BMG_IOFLAG, BMG_iPARMS
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_SymStd_SOLVE_cg does a direct solve on the coarsest grid. it
C     uses the LAPACK routine DPBTRS.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     QF        Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     
C     II        Number of grid points in x direction, including
C               two fictitious points.
C     JJ        Number of grid points in y direction, including
C               two fictitious points.
C
C     ABD       Refer to BMG2_SER_SymStd_SOLVE_boxmg
C
C     NABD1     Refer to BMG2_SER_SymStd_SOLVE_boxmg
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     BBD       Workspace, refer to BMG2_SER_SymStd_SOLVE_boxmg
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C     Q         Refer to BMG2_SER_SymStd_SOLVE_boxmg
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
      INCLUDE  'BMG_SER_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  II, JJ, NABD1, NABD2
      REAL*RKIND   ABD(NABD1,NABD2), BBD(NABD2), Q(II,JJ), QF(II,JJ)
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS)

C ----------------------------
C     Local Declarations
C
      INTEGER  I, I1, I2, J, J1, KK, N, INFO

C =========================================================================

C -------------------------------------------
C     Copy the RHS into LAPACK
C -------------------------------------------

      I1=II-1
      J1=JJ-1

      I2=I1-1
      N=I2*(J1-1)

      KK=0
      DO J=2,J1
         DO I=2,I1
            KK=KK+1
            BBD(KK)=QF(I,J)
         ENDDO
      ENDDO

C -------------------------------------------
C     Solve with LAPACK routine DPBTRS
C -------------------------------------------

#if RKIND == 4
      CALL SPBTRS ('U', KK, I1, 1, ABD, NABD1, BBD, NABD2, INFO) 
#else
      CALL DPBTRS ('U', KK, I1, 1, ABD, NABD1, BBD, NABD2, INFO) 
#endif

      IF (INFO .NE. 0) THEN

         IF (BMG_IOFLAG(iBMG2_SER_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'Coarse grid solve failed!'
            WRITE(*,510) 'INFO = ', INFO
         END IF

         BMG_iPARMS(id_BMG2_SER_Ext_Err_Code) = INFO
         CALL BMG2_SER_SymStd_ErrTrap(BMG_iPARMS,20)
         RETURN

      ENDIF

C -------------------------------------------
C     Copy the solution back from LAPACK
C -------------------------------------------

      KK=0
      DO J=2,J1
         DO I=2,I1
            KK=KK+1
            Q(I,J)=BBD(KK)
         ENDDO
      ENDDO

C =========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SER_SymStd_SOLVE_cg.f',/,5X,A)
 510  FORMAT (5X,A,1X,I3)

C ===========================================

      RETURN
      END
