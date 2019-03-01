      SUBROUTINE BMG3_SymStd_COPY_RHS_yz( 
     &                       SO, Q, QF, QF2, iPL, Nx, Ny, Nz, NStncl
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_COPY_RHS_yz creates the right hand side vector, QF,
C     from the iPL{th}-(y,z) plane of the current 3D vector.
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

      IMPLICIT   NONE

C ---------------------------
C     Includes
C
      INCLUDE 'BMG_stencils.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  iPL, Nx, Ny, Nz, NStncl
      REAL*RKIND   Q(Nx,Ny,Nz), QF(Nx,Ny,Nz), QF2(Ny,Nz),
     &         SO(Nx,Ny,Nz,NStncl)

C ---------------------------
C     Local Declarations:
C
      INTEGER  j, k

C ==========================================================================

      IF ( NStncl.EQ.14 ) THEN
         !
         DO k=2, Nz-1
            DO j=2, Ny-1
               !
               QF2(j,k) = QF(iPL,j,k)
     &                  + SO(iPL,j+1,k,kpnw)*Q(iPL-1,j+1,k)
     &                  + SO(iPL,j,k,kpw)*Q(iPL-1,j,k)
     &                  + SO(iPL,j,k,kpsw)*Q(iPL-1,j-1,k)
     &                  + SO(iPL,j+1,k,kbnw)*Q(iPL-1,j+1,k-1)
     &                  + SO(iPL,j,k,kbw)*Q(iPL-1,j,k-1)
     &                  + SO(iPL,j,k,kbsw)*Q(iPL-1,j-1,k-1)
     &                  + SO(iPL,j+1,k+1,kbse)*Q(iPL-1,j+1,k+1)
     &                  + SO(iPL,j,k+1,kbe)*Q(iPL-1,j,k+1)
     &                  + SO(iPL,j,k+1,kbne)*Q(iPL-1,j-1,k+1)
     &                  + SO(iPL+1,j+1,k,kpsw)*Q(iPL+1,j+1,k)
     &                  + SO(iPL+1,j,k,kpw)*Q(iPL+1,j,k)
     &                  + SO(iPL+1,j,k,kpnw)*Q(iPL+1,j-1,k)
     &                  + SO(iPL+1,j+1,k,kbne)*Q(iPL+1,j+1,k-1)
     &                  + SO(iPL+1,j,k,kbe)*Q(iPL+1,j,k-1)
     &                  + SO(iPL+1,j,k,kbse)*Q(iPL+1,j-1,k-1)
     &                  + SO(iPL+1,j+1,k+1,kbsw)*Q(iPL+1,j+1,k+1)
     &                  + SO(iPL+1,j,k+1,kbw)*Q(iPL+1,j,k+1)
     &                  + SO(iPL+1,j,k+1,kbnw)*Q(iPL+1,j-1,k+1)
               !
            ENDDO
         ENDDO
         !
      ELSE
         !
         DO k=2, Nz-1
            DO j=2, Ny-1
               QF2(j,k) = QF(iPL,j,k)
     &                  + SO(iPL,j,k,kpw)*Q(iPL-1,j,k)
     &                  + SO(iPL+1,j,k,kpw)*Q(iPL+1,j,k)
            ENDDO
         ENDDO
         !
      ENDIF

C ==========================================================================

      RETURN
      END

