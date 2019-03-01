      SUBROUTINE PUTF(
     &                SO, QF, u, NLx, NLy, iGs, jGs,
     &                NGx, NGy, x, y, deltax, deltay, 
     &                MPI_COMM
     &                )

C ==========================================================================C
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     PUTF creates and stores the discretization in the SO array for use
C     by BoxMG.  This example is Poisson discretized on vertices, but
C     with the source generated from the discrete solution.  Hence, this
C     is not a real example, but handy for debugging as the error in the
C     computed solution will only depend on the tolerance chosen for the
C     iterative solve.
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
C   OUTPUT:
C  --------------------
C
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

C ---------------------------
C    Includes:
C
      INCLUDE  'BMG_stencils.h'
      INCLUDE  'BMG_constants.h'

C ---------------------------
C    Argument Declarations:
C     
      INTEGER   NLx, NLy, NGx, NGy, iGs, jGs
      REAL*8    deltax, deltay, SO(0:NLx+2,0:NLy+2,3),
     &          QF(0:NLx+1,0:NLy+1), u(0:NLx+1,0:NLy+1),
     &          x(0:NLx+1), y(0:NLy+1)

C ---------------------------
C    Local Declarations:
C     
      INTEGER   i, j, MyProc, MPI_COMM, MPI_IERR
      REAL*8    deltascale

C ==========================================================================

      CALL MPI_COMM_RANK( MPI_COMM, MyProc, MPI_IERR )

C -------------------------------------------------
C     Compute 5-point vertex based stencil
C -------------------------------------------------

      !
      ! Scaling of all equations, normally this would be
      ! deltax*deltay, not deltascale**2.  But then again,
      ! this is not a real problem.
      !
      deltascale = min(deltax,deltay)

      DO j = 0, NLy+2
         DO i = 0, NLx+2
            SO(i,j,ks) = 1D0/deltay**2 * deltascale**2
         ENDDO
      ENDDO   

      DO  j = 0, NLy+2
         DO  i = 0, NLx+2
            SO(i,j,kw) = 1D0/deltax**2 * deltascale**2
         ENDDO
      ENDDO
         
      DO  j = 0, NLy+1
         DO  i = 0, NLx+1
            SO(i,j,ko) = SO(i,j,kw)+SO(i,j,ks)+SO(i+1,j,kw)+SO(i,j+1,ks)
         ENDDO
      ENDDO

      DO j = 1, NLy
         DO i=1, NLx
            QF(i,j) =  SO(i,j,ko)  * u(i,j) 
     &                -SO(i,j,kw)  * u(i-1,j)
     &                -SO(i,j,ks)  * u(i,j-1)
     &                -SO(i+1,j,kw)* u(i+1,j)
     &                -SO(i,j+1,ks)* u(i,j+1)
         ENDDO
      ENDDO

C -------------------------------------------------
C     Apply Dirichlet BC by modifying QF (RHS).
C -------------------------------------------------

      IF (iGs.EQ.1) THEN
        DO j = 1, NLy
           QF(1,j) = QF(1,j) + SO(1,j,kw)*u(0,j)
        ENDDO
      ENDIF

      IF (iGs+NLx-1.EQ.NGx) THEN
         DO j = 1, NLy
            QF(NLx,j) = QF(NLx,j) + SO(NLx+1,j,kw)*u(NLx+1,j)
         ENDDO
      ENDIF

      IF (jGs.EQ.1) THEN
         DO i = 1, NLx
            QF(i,1) = QF(i,1) + SO(i,1,ks)*u(i,0)
         ENDDO
      ENDIF

      IF (jGs+NLy-1.EQ.NGy) THEN
         DO i = 1, NLx
            QF(i,NLy) = QF(i,NLy) + SO(i,NLy+1,ks)*u(i,NLy+1)
         ENDDO
      ENDIF

C -------------------------------------------------
C     Zero connections to Global Ghost points
C -------------------------------------------------

      IF (iGs.EQ.1) THEN
         DO j = 1, NLy
            SO(0,j,kw) = 0.0
            SO(0,j,ko) = 0.0
            SO(0,j,ks) = 0.0
            SO(1,j,kw) = 0.0
         ENDDO
         SO(0,0,kw) = 0.0
         SO(0,0,ko) = 0.0
         SO(0,0,ks) = 0.0
         SO(1,0,kw) = 0.0
         SO(0,NLy+1,kw) = 0.0
         SO(0,NLy+1,ko) = 0.0
         SO(0,NLy+1,ks) = 0.0
         SO(1,NLy+1,kw) = 0.0
         SO(0,NLy+2,kw) = 0.0
         SO(0,NLy+2,ko) = 0.0
         SO(0,NLy+2,ks) = 0.0
         SO(1,NLy+2,kw) = 0.0
      ENDIF
      
      IF (iGs+NLx-1.EQ.NGx) THEN
         DO j= 1, NLy+2
            SO(NLx+1,j,kw) = 0.0
            SO(NLx+1,j,ko) = 0.0
            SO(NLx+1,j,ks) = 0.0
            SO(NLx+2,j,kw) = 0.0
            SO(NLx+2,j,ko) = 0.0
            SO(NLx+2,j,ks) = 0.0
         ENDDO

         SO(NLx+1,0,kw) = 0.0
         SO(NLx+1,0,ks) = 0.0
         SO(NLx+1,0,ko) = 0.0

         SO(NLx+2,0,kw) = 0.0
         SO(NLx+2,0,ks) = 0.0
         SO(NLx+2,0,ko) = 0.0

         SO(NLx+1,NLy+1,kw) = 0.0
         SO(NLx+1,NLy+1,ks) = 0.0
         SO(NLx+1,NLy+1,ko) = 0.0

         SO(NLx+1,NLy+2,kw) = 0.0
         SO(NLx+1,NLy+2,ks) = 0.0
         SO(NLx+1,NLy+2,ko) = 0.0

         SO(NLx+2,NLy+1,kw) = 0.0
         SO(NLx+2,NLy+1,ks) = 0.0
         SO(NLx+2,NLy+1,ko) = 0.0

         SO(NLx+2,NLy+2,kw) = 0.0
         SO(NLx+2,NLy+2,ks) = 0.0
         SO(NLx+2,NLy+2,ko) = 0.0
      ENDIF

      IF (jGs.EQ.1) THEN
         DO i = 1, NLx
            SO(i,0,ks) = 0.0
            SO(i,0,ko) = 0.0
            SO(i,0,kw) = 0.0
            SO(i,1,ks) = 0.0
         ENDDO
         SO(0,0,kw) = 0.0
         SO(0,0,ko) = 0.0
         SO(0,0,ks) = 0.0
         SO(0,1,ks) = 0.0
         SO(NLx+1,0,kw) = 0.0
         SO(NLx+1,0,ko) = 0.0
         SO(NLx+1,0,ks) = 0.0
         SO(NLx+1,1,ks) = 0.0
         SO(NLx+2,0,kw) = 0.0
         SO(NLx+2,0,ko) = 0.0
         SO(NLx+2,0,ks) = 0.0
         SO(NLx+2,1,ks) = 0.0
      ENDIF

      IF (jGs+NLy-1.EQ.NGy) THEN
         DO i = 1, NLx+2
            SO(i,NLy+1,ks) = 0.0
            SO(i,NLy+1,kw) = 0.0
            SO(i,NLy+1,ko) = 0.0
            SO(i,NLy+2,ks) = 0.0
            SO(i,NLy+2,kw) = 0.0
            SO(i,NLy+2,ko) = 0.0
         ENDDO
         SO(0,NLy+1,kw) = 0.0
         SO(0,NLy+1,ko) = 0.0
         SO(0,NLy+1,ks) = 0.0

         SO(0,NLy+2,kw) = 0.0
         SO(0,NLy+2,ko) = 0.0
         SO(0,NLy+2,ks) = 0.0

      ENDIF

C ==========================================================================

      RETURN
      END

