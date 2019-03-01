      SUBROUTINE BMG3_SymStd_COPY_RHS_xz( 
     &                       SO, Q, QF, QF2, iPL, Nx, Ny, Nz, NStncl
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_COPY_RHS_xz creates the right hand side vector, QF,
C     from the iPL{th}-(x,z) plane of the current 3D vector.
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
      REAL*RKIND   Q(Nx,Ny,Nz), QF(Nx,Ny,Nz), QF2(Nx,Nz),
     &         SO(Nx,Ny,Nz,NStncl)

C ---------------------------
C     Local Declarations:
C
      INTEGER  i, k

C ==========================================================================

      IF ( NStncl.EQ.14 ) THEN
         DO k=2, Nz-1
            DO i=2, Nx-1
               QF2(i,k) = QF(i,iPL,k)
     &                  + SO(i,iPL+1,k,kpnw)*Q(i-1,iPL+1,k)
     &                  + SO(i,iPL+1,k,kps)*Q(i,iPL+1,k)
     &                  + SO(i+1,iPL+1,k,kpsw)*Q(i+1,iPL+1,k)
     &                  + SO(i,iPL+1,k,kbnw)*Q(i-1,iPL+1,k-1)
     &                  + SO(i,iPL+1,k,kbn)*Q(i,iPL+1,k-1)
     &                  + SO(i+1,iPL+1,k,kbne)*Q(i+1,iPL+1,k-1)
     &                  + SO(i,iPL+1,k+1,kbse)*Q(i-1,iPL+1,k+1)
     &                  + SO(i,iPL+1,k+1,kbs)*Q(i,iPL+1,k+1)
     &                  + SO(i+1,iPL+1,k+1,kbsw)*Q(i+1,iPL+1,k+1)
     &                  + SO(i,iPL,k,kpsw)*Q(i-1,iPL-1,k)
     &                  + SO(i,iPL,k,kps)*Q(i,iPL-1,k)
     &                  + SO(i+1,iPL,k,kpnw)*Q(i+1,iPL-1,k)
     &                  + SO(i,iPL,k,kbsw)*Q(i-1,iPL-1,k-1)
     &                  + SO(i,iPL,k,kbs)*Q(i,iPL-1,k-1)
     &                  + SO(i+1,iPL,k,kbse)*Q(i+1,iPL-1,k-1)
     &                  + SO(i,iPL,k+1,kbne)*Q(i-1,iPL-1,k+1)
     &                  + SO(i,iPL,k+1,kbn)*Q(i,iPL-1,k+1)
     &                  + SO(i+1,iPL,k+1,kbnw)*Q(i+1,iPL-1,k+1)
            ENDDO
         ENDDO
      ELSE
         !
         DO k=2, Nz-1
            DO i=2, Nx-1
               QF2(i,k) = QF(i,iPL,k)
     &                  + SO(i,iPL,k,kps)*Q(i,iPL-1,k)
     &                  + SO(i,iPL+1,k,kps)*Q(i,iPL+1,k)
            ENDDO
         ENDDO
         !
      ENDIF

C ==========================================================================

      RETURN
      END

