      SUBROUTINE BMG2_SymStd_SETUP_relax( 
     &                IFD, IRELAX, SO, NSOm, SOR, NSORm, 
     &                IGRD, NOGm, NOG, BMG_iPARMS, BMG_IOFLAG,
     &                iWork, NMSGi, pMSG, BUFFER, NMSGr,
     &                MyProc, MyProcI, MyProcJ,
     &                MPICOMM, XLINECOMM, YLINECOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SETUP_relax performs any necessary setup for the
C     chosen relaxation scheme.  In particular, it performs calls
C     routines to perform the LU factorization of the Tridiagonal
C     systems that arise in line relaxation.
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
#include      "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ---------------------------
C    Argument Declarations:
C
      INTEGER  NMSGi, NMSGr, NSOm, NSORm, NOGm

      INTEGER  IFD, IGRD(NOGm,NBMG_pIGRD), IRELAX, NOG,
     &         iWork(NMSGi), pMSG(NBMG_pMSG,NOG), MPICOMM,
     &         MyProc, MyProcI, MyProcJ, XLINECOMM, YLINECOMM

      REAL*RKIND   SO(NSOm), SOR(NSORm), BUFFER(NMSGr)

      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_IPARMS)

C ---------------------------
C    Local Declarations:
C
      INTEGER  K, NStncl, NLx, NLy, NGx, NGy, p_CI, p_SO, p_SOR, p_U
      INTEGER  ierror

C ==========================================================================

      !
      ! Sanity check
      !
      IF ( NOG.EQ.1) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,500) 'There is only 1 grid!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,11)
         RETURN

      ELSE IF ( NOG.EQ.0 ) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,500) 'There are no grids?'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,12)
         RETURN
          
      ELSE IF ( NOG.LT.0 ) THEN
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,500) 'There number of grids is negative!'
            WRITE(*,520) 'HAVE: NOG = ', NOG
         END IF

         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,13)
         RETURN

      ENDIF

C     
C     create the line-communicators
C
      IF ( IRELAX.NE.BMG_GS_RB_point) THEN
         ! set up both x- and y-line communicators
         CALL MPI_Comm_Split (MPICOMM, MyProcJ, MyProcI, XLINECOMM, 
     &        ierror)
         CALL MPI_Comm_Split (MPICOMM, MyProcI, MyProcJ, YLINECOMM, 
     &        ierror)
      ENDIF



      !
      ! Loop over grids
      !
      DO  K = NOG, 2, -1
         !
         ! Determine the number of points in the stencil
         !
         IF (K.NE.NOG .OR. IFD.NE.1) THEN
            NStncl=5
         ELSE
            NStncl=3
         ENDIF

         !
         ! (fake) memory pointers
         !
         CALL BMG2_SymStd_GET_pointers( K, IGRD, NOGm,
     &                p_U, p_SO, p_SOR, p_CI, NLx, NLy, NGx, NGy )

         IF ( IRELAX.EQ.BMG_GS_RB_point ) THEN
            !
            CALL BMG2_SymStd_SETUP_recip(
     &                SO(p_SO), SOR(p_SOR), NLx, NLy, NStncl
     &                )
         ELSE IF ( IRELAX.EQ.BMG_GS_RB_x_lines ) THEN
            !
            p_SOR = p_SOR + (ipL_BMG_LUL1-1)*NLx*NLy 
            CALL BMG2_SymStd_SETUP_lines_x( 
     &                SO(p_SO), SOR(p_SOR), NLx, NLy, NStncl
     &                )
            !
         ELSE IF ( IRELAX.EQ.BMG_GS_RB_y_lines ) THEN
            !
            p_SOR = p_SOR + (ipL_BMG_LUL2-1)*NLx*NLy
            CALL BMG2_SymStd_SETUP_lines_y(
     &                SO(p_SO), SOR(p_SOR), NLx, NLy, NStncl
     &                )
            !
         ELSE IF ( IRELAX.EQ.BMG_GS_RB_x_y_lines ) THEN
            !
            p_SOR = p_SOR + (ipL_BMG_LUL1-1)*NLx*NLy

            CALL BMG2_SymStd_SETUP_lines_x( 
     &                SO(p_SO), SOR(p_SOR), NLx, NLy, NStncl
     &                )
            p_SOR = p_SOR + (ipL_BMG_LUL2-1)*NLx*NLy 
            CALL BMG2_SymStd_SETUP_lines_y( 
     &                SO(p_SO), SOR(p_SOR), NLx, NLy, NStncl
     &                )
            !
         ENDIF

      ENDDO

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_relax.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C =====================

      RETURN
      END
