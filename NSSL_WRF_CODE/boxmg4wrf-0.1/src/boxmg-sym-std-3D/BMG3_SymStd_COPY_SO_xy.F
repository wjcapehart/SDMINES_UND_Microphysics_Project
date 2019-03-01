      SUBROUTINE BMG3_SymStd_COPY_SO_xy(
     &                       SO, SO_xy, iPL, Nx, Ny, Nz,
     &                       NStncl_3D, NStncl_2D
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_COPY_SO_xy copies the iPL{th}-(x,y) plane of the 3D
C     stencil into SO_xy.
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
#include    "BMG_constants.h"
      INCLUDE    'BMG_stencils.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  Nx, iPL, Ny, Nz, NStncl_2D, NStncl_3D
      REAL*RKIND   SO(Nx,Ny,Nz,NStncl_3D), SO_xy(Nx,Ny,NStncl_2D)

C ---------------------------
C     Local Declarations:
C
      INTEGER  i, j

C ==========================================================================

      IF ( NStncl_3D.EQ.14 .AND. NStncl_2D.EQ.5 ) THEN
         !
         DO j=1, Ny
            DO i=1, Nx
               SO_xy(i,j,ko)  = SO(i,j,iPL,kp)
               SO_xy(i,j,kw)  = SO(i,j,iPL,kpw)
               SO_xy(i,j,ks)  = SO(i,j,iPL,kps)
               SO_xy(i,j,ksw) = SO(i,j,iPL,kpsw)
               SO_xy(i,j,knw) = SO(i,j,iPL,kpnw)
            ENDDO
         ENDDO
         !
      ELSEIF ( NStncl_3D.EQ.4 .AND. NStncl_2D.EQ.3 ) THEN
         !
         DO j=1, Ny
            DO i=1, Nx
               SO_xy(i,j,ko) = SO(i,j,iPL,kp)
               SO_xy(i,j,kw) = SO(i,j,iPL,kpw)
               SO_xy(i,j,ks) = SO(i,j,iPL,kps)
            ENDDO
         ENDDO
         !
      ELSE
         !
         WRITE(*,500) 'Inconsitent number of stencil entries: '
         WRITE(*,510) 'HAVE: NStncl_3D = ', NStncl_3D
         WRITE(*,520) 'HAVE: NStncl_2D = ', NStncl_2D
         STOP
         !
      ENDIF

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_COPY_SO_xy.f',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
