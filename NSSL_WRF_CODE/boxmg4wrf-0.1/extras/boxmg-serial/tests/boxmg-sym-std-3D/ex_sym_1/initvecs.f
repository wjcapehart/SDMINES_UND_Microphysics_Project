      subroutine initvecs ( u, v, Nx, Ny, Nz )
      
C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Initalize vectors to use in the inner product check of symmetry.
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
      INCLUDE 'BMG_SER_constants.h'

C -----------------------------
C     Argument Declarations
C 
      INTEGER   Nx, Ny, Nz
      REAL*8    u(0:Nx-1,0:Ny-1,0:Nz-1), v(0:Nx-1,0:Ny-1,0:Nz-1)

C -----------------------------
C     Local Declarations
C
      INTEGER   i, j, k

C -----------------------------
C     External Functions
C
      REAL*8   BMG_rand
      EXTERNAL BMG_rand

C ==========================================================================

      DO k=0,Nz-1
         DO j=0,Ny-1
            DO i=0,Nx-1
               IF ( i.eq.0.or.i.eq.Nx-1 .or.
     &              j.eq.0.or.j.eq.Ny-1 .or.
     &              k.eq.0.or.k.eq.Nz-1 ) THEN
                  u(i,j,k) = 0.0
                  v(i,j,k) = 0.0
               ELSE
                  u(i,j,k) = BMG_rand( rZERO )
                  v(i,j,k) = BMG_rand( rZERO )
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      
C ==========================================================================

      return
      end