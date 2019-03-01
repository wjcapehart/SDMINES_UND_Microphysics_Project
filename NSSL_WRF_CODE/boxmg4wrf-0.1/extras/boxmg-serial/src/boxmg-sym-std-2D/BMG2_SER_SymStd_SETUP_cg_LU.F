      SUBROUTINE BMG2_SER_SymStd_SETUP_cg_LU(
     &                SO, II, JJ, NStncl, ABD, NABD1, NABD2,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SER_SymStd_SETUP_cg_LU sets up the matrix on the coarsest grid,
C     and using the linpack routine it forms the l-u decomposition of
C     the matrix.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     SO        Refer to BMG2_SER_SymStd_SOLVE_boxmg
C     II        Number of grid points in x direction, including
C               two fictitious points.
C     JJ        Number of grid points in y direction, including
C               two fictitious points.
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
C     ABD       Refer to BMG2_SER_SymStd_SOLVE_boxmg
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
      INTEGER  II, JJ, NABD1, NABD2, NStncl
      REAL*RKIND   ABD(NABD1,NABD2), SO(II,JJ,NStncl)
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS)

C ----------------------------
C     Local Declarations
C
      INTEGER  I, INFO, I1, I2, J, J1, KK, N

C ==========================================================================

C -------------------------------------------------------
C     Copy the operator on the coarsest grid into ABD 
C -------------------------------------------------------

      I1=II-1
      J1=JJ-1

      I2=I1-1
      N=I2*(J1-1)

      IF ( NStncl.EQ.5 ) THEN
         KK=0
         DO J=2,J1
            DO I=2,I1
               KK=KK+1
               ABD(II,KK)=SO(I,J,KO)
               ABD(I1,KK)=-SO(I,J,KW)
               ABD(3,KK)=-SO(I+1,J,KNW)
               ABD(2,KK)=-SO(I,J,KS)
               ABD(1,KK)=-SO(I,J,KSW)
            ENDDO
         ENDDO
      ELSEIF ( NStncl.EQ.3 ) THEN
         KK=0
         DO J=2,J1
            DO I=2,I1
               KK=KK+1
               ABD(II,KK)=SO(I,J,KO)
               ABD(I1,KK)=-SO(I,J,KW)
               ABD(3,KK)=0
               ABD(2,KK)=-SO(I,J,KS)
               ABD(1,KK)=0
            ENDDO
         ENDDO
      ELSE
         IF (BMG_IOFLAG(iBMG2_SER_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 
            WRITE(*,510) 'NEED: NStncl = 3 or 5 '
            WRITE(*,510) 'HAVE: NStncl = ', NStncl
         END IF
         
         CALL BMG2_SER_SymStd_ErrTrap(BMG_iPARMS,18)
         RETURN

      ENDIF
         
C -------------------------------------------------------
C     Factor using the LAPACK routine DPBTRF
C -------------------------------------------------------

#if RKIND == 4
      CALL SPBTRF('U', KK, I1, ABD, NABD1, INFO) 
#else
      CALL DPBTRF('U', KK, I1, ABD, NABD1, INFO) 
#endif
      
      IF (INFO .NE. 0) THEN

         IF (BMG_IOFLAG(iBMG2_SER_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'Coarse grid Cholesky decomposition failed!'
            WRITE(*,510) 'INFO = ', INFO
         END IF

         BMG_iPARMS(id_BMG2_SER_Ext_Err_Code) = INFO
         CALL BMG2_SER_SymStd_ErrTrap(BMG_iPARMS,17)
         RETURN

      ENDIF

C ==========================================================================

 500    FORMAT (/,'FATAL ERROR: BMG2_SER_SymStd_SETUP_cg_LU.f',/,5X,A)
 510    FORMAT (5X,A,1X,I3)

C ====================

      RETURN
      END
