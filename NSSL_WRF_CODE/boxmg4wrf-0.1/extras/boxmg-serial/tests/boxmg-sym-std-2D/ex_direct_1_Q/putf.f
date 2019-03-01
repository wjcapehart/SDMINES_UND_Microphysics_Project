      SUBROUTINE PUTF( SO, QF, Q,
     &                 Nx, Ny, hx, hy, kg 
     &               )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Creates the stencil representation of the discretized PDE.
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

C ==========================================================================

C ----------------
C     Includes
C
      INCLUDE  'BMG_SER_constants.h'
      INCLUDE  'BMG_SER_stencils.h'
      
C ----------------
C     Parmeters
C
      REAL*8   rHALF
      PARAMETER ( rHALF = 0.5D0 )

C ----------------
C     Arguments
C
      INTEGER  Nx, Ny, kg
      REAL*8   hx, hy
      REAL*8   QF(Nx+2,Ny+2), SO(Nx+2,Ny+2,3), Q(Nx+2,Ny+2)

C ----------------
C     Local
C
      INTEGER  i, i1, i2, j, j1, j2,
     &         iBEG, iEND, jBEG, jEND,
     &         Nx_g, Ny_g, NSOx_g, NSOy_g
      REAL*8   cm, fs, h2, s, xh, yh

C -----------------
C     External
C
      REAL*8   da, dr
      EXTERNAL da, dr

C -----------------
C     Includes
C
      INCLUDE 'common2.h'

C ==========================================================================

C ----------------------------------------
C     Initalize constants
C ----------------------------------------

      !
      !  Standard halo with depth = [1,1]
      !
      Nx_g = Nx + 2
      Ny_g = Ny + 2

      i1 = Nx + 1
      j1 = Ny + 1

      i2 = Nx
      j2 = Ny

      !
      !  Stencil halo with depth = [1,1]
      !
      NSOx_g = Nx + 2
      NSOy_g = Ny + 2

      !
      !
      !  Grid Spacing
      !
      h2=hx*hy

      xh=hy/hx
      yh=hx/hy

      !
      !  Loop Boundaries: iGs, jGs, kGs
      !
      iBEG = 3
      jBEG = 3


      !
      !  Loop Boundaries: iGf, jGf, kGf
      !
      iEND = i1
      jEND = j1
         
C ----------------------------------------
C     Zero right-hand-side
C ----------------------------------------

      DO j=1, Ny_g
         DO i=1, Nx_g
            !
            Q(i,j)  = rZERO
            QF(i,j) = rZERO
            !
         END DO
      END DO

C ----------------------------------------
C     Zero Stencil
C ----------------------------------------

      DO j=1, NSOy_g
         DO i=1, NSOx_g

            SO(i,j,ko) = rZERO
            SO(i,j,kw) = rZERO
            SO(i,j,ks) = rZERO
            
         END DO
      END DO

C ----------------------------------------
C     Compute the source
C ----------------------------------------

      DO j=2, j1
         DO i=2, i1
            !
            CALL RHS( i,j,hx,hy,fs,s )
            !
            qf(i,j) = fs*h2
            !
         ENDDO
      ENDDO
 
C ----------------------------------------
C     Compute the removal/absorption
C ----------------------------------------

      DO j=2, j1
         DO i=2, i1
            !
            CALL RHS( i,j,hx,hy,fs,s )
            !
            SO(i,j,ko) = s*h2
            !
         ENDDO
      ENDDO

C ----------------------------------------
C     Stencil coefficients
C ----------------------------------------

      !
      !  South
      !
      DO j=jBEG,jEND
         DO i=2, iEND
            SO(i,j,ks)  = da( i,j,hx,hy )*yh
         ENDDO
      ENDDO

      !
      !  West
      !
      DO j=2, jEND
         DO i=iBEG, iEND
            SO(i,j,kw) = dr( i,j,hx,hy )*xh
         ENDDO
      ENDDO

C ----------------------------------------
C     Initial Guess
C ----------------------------------------
  
c$$$      DO j=1, Ny_g
c$$$         DO i=1, Nx_g
c$$$            q(i,j)=0.
c$$$         ENDDO
c$$$      ENDDO

C ---------------------------------------
C     Boundary Conditions
C ---------------------------------------

      !
      !  West boundary
      !
      DO j=2, jEND
         SO(2,j,ks)  = rHALF*SO(2,j,ks)
      ENDDO


      IF ( iBL.EQ.1) THEN
         
         DO j=2, jEND
            SO(2,j,ko) = SO(2,j,ko) + hy
         ENDDO
         
      ENDIF

      !
      !  East boundary
      !
      DO j=2, jEND
         SO(i1,j,ks) = rHALF*SO(i1,j,ks)
      ENDDO

      IF ( iBR.EQ.1) THEN

         DO j=2, jEND
            SO(i1,j,ko) = SO(i1,j,ko) + hy
         ENDDO
         
      ENDIF

      !
      !  Bottom Boundary
      !
      DO i=2, iEND
         SO(i,2,kw)  = rHALF*SO(i,2,kw)
      ENDDO

      IF ( iBB.EQ.1) THEN
         
         DO i=2, iEND
            SO(i,2,ko) = SO(i,2,ko) + hx   
         ENDDO
         
      ENDIF

      !
      !  Top boundary
      !
      DO i=2, iEND
         SO(i,j1,kw) = rHALF*SO(i,j1,kw)
      ENDDO

      IF( iBT.EQ.1) THEN

         DO i=2, iEND
            SO(i,j1,ko) = SO(i,j1,ko) + hx
         ENDDO
         
      ENDIF

      DO j=2, j1
         DO i=2, i1
            cm=rONE
            IF ( i.EQ.2 .OR. i.EQ.i1 ) THEN
               cm = rHALF*cm
            ENDIF
            IF ( j.EQ.2 .OR. j.EQ.j1 ) THEN
               cm = rHALF*cm
            ENDIF
            !
            SO(i,j,ko) = cm*SO(i,j,ko)
     &           + SO(i,j,kw)   + SO(i,j+1,ks)
     &           + SO(i+1,j,kw) + SO(i,j,kps)
            !
            qf(i,j) = cm*qf(i,j)
            !
         ENDDO
      ENDDO

C ==========================================================================

      RETURN
      END


