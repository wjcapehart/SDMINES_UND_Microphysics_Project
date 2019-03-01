      SUBROUTINE BMG3_SymStd_UTILS_dot_l2( 
     &               u, v, Nx, Ny, Nz, l2dot, MPICOMM 
     &               )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_UTILS_dot_l2.f computes the dot product of a grid
C     function u with grid function v.  It assumes that ghost points
C     should be neglected.
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
#if RKIND == 4
#define FLOAT_PRECISION MPI_REAL
#else
#define FLOAT_PRECISION MPI_DOUBLE_PRECISION
#endif

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'mpif.h'
      
#include    "BMG_constants.h"

C ----------------------------
C     Argument Declarations
C 
      INTEGER  Nx, Ny, Nz, MPICOMM
      REAL*RKIND  u(Nx,Ny,Nz), v(Nx,Ny,Nz)
      REAL*8  l2dot
C ----------------------------
C     Local Declarations
C
      INTEGER  i, j, k, ierror
      REAL*8   buffer

C ==========================================================================

      l2dot=dZERO

      DO k = 2, Nz-1
         DO j = 2, Ny-1
            DO i = 2, Nx-1
               l2dot=l2dot+u(i,j,k)*v(i,j,k)
            ENDDO
         ENDDO
      ENDDO

      call MPI_AllReduce(l2dot, buffer,1,    
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPICOMM,ierror)

      l2dot=buffer

C ==========================================================================

      RETURN
      END
