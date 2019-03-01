      SUBROUTINE BMG3_SER_SymStd_ErrTrap( 
     &                BMG_iPARMS, IERROR
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_SETUP_ErrTrap stores a nonzero error code in
C     BMG_iPARMS for the user to trap.
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
      INCLUDE 'BMG_SER_parameters.h'
#include       "BMG_SER_constants.h"

C ----------------------------
C     Argument Declarations
C 

      INTEGER  IERROR
      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS) 

C ----------------------------
C     Local Declarations
C

C ==========================================================================
      
      !
      ! first test if an error has occured
      !
      IF (IERROR .ne. iZERO) THEN
         !
         ! store the error code in BMG_iPARMS
         !
         BMG_iPARMS(id_BMG3_SER_Err_Code) = IERROR
      ENDIF

C ===========================================

      RETURN
      END
