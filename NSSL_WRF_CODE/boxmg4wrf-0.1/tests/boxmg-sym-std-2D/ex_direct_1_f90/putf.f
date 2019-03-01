      SUBROUTINE PUTF( 
     &                 MyProc, SO, QF, Q,
     &                 NLx, NLy, NGx, NGy, 
     &                 iGs, jGs, hx, hy, kg 
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

C ----------------
C     Includes
C
      INCLUDE  'BMG_constants.h'
      INCLUDE  'BMG_stencils.h'
      
C ----------------
C     Parmeters
C
      REAL*8   rHALF
      PARAMETER ( rHALF = 0.5D0 )

C ----------------
C     Arguments
C
      INTEGER  iGs, jGs, NLx, NLy, NGx, NGy, kg
      REAL*8   hx, hy
      REAL*8   QF(NLx+2,NLy+2), SO(NLx+3,NLy+3,3),
     &         Q(NLx+2,NLy+2)

C ----------------
C     Local
C
      INTEGER  i, i1, i2, iGf, j, j1, j2, jGf, 
     &         iBEG, iEND, jBEG, jEND, MyProc,
     &         NLx_g, NLy_g, NLSOx_g, NLSOy_g,
     &	       is, js
      REAL*8   cm, fs, h2, s, xh, yh

C -----------------
C     External
C
      REAL*8   da, dr, BMG_rand
      EXTERNAL da, dr, BMG_rand

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
      NLx_g = NLx + 2
      NLy_g = NLy + 2

      i1 = NLx + 1
      j1 = NLy + 1

      i2 = NLx
      j2 = NLy

      !
      !  Global indexing
      !
      iGf=iGs+i2-1
      jGf=jGs+j2-1

      !
      !  Stencil halo with depth = [1,2]
      !
      NLSOx_g = NLx + 3
      NLSOy_g = NLy + 3

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
      IF ( iGs.EQ.1 ) THEN
         iBEG = 3
      ELSE
         iBEG = 2
      ENDIF

      IF ( jGs.EQ.1 ) THEN
         jBEG = 3
      ELSE
         jBEG = 2
      ENDIF


      !
      !  Loop Boundaries: iGf, jGf, kGf
      !
      IF ( iGf.EQ.NGx ) THEN
         iEND = i1
      ELSE
         iEND = NLx_g
      ENDIF

      IF ( jGf.EQ.NGy ) THEN
         jEND = j1
      ELSE
         jEND = NLy_g
      ENDIF

         
C ----------------------------------------
C     Zero right-hand-side
C ----------------------------------------

      DO j=1, NLy_g
         DO i=1, NLx_g
            !
            Q(i,j) = rZERO
            QF(i,j)  = rZERO
            !
         END DO
      END DO

C ----------------------------------------
C     Zero Stencil
C ----------------------------------------

      DO j=1, NLSOy_g
         DO i=1, NLSOx_g

            SO(i,j,ko) = rZERO
            SO(i,j,kw) = rZERO
            SO(i,j,ks) = rZERO
            
         END DO
      END DO

C ----------------------------------------
C     Compute the source
C ----------------------------------------

c$$$      DO j=2, j1
c$$$         DO i=2, i1
c$$$            !
!               is = iGs+i-1
!               js = jGs+j-1
!               !
c$$$            CALL RHS( is,js,hx,hy,fs,s )
c$$$            !
c$$$            qf(i,j) = fs*h2
c$$$            !
c$$$         ENDDO
c$$$      ENDDO
 
C ----------------------------------------
C     Compute the removal/absorption
C ----------------------------------------

      DO j=2, j1
         DO i=2, i1
            !
            is = iGs+i-1
            js = jGs+j-1
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
            !
            is = iGs+i-1
            js = jGs+j-1
            !
            SO(i,j,ks)  = da(is,js,hx,hy)*yh
         ENDDO
      ENDDO

      !
      !  West
      !
      DO j=2, jEND
         DO i=iBEG, iEND
            !
            is = iGs+i-1
            js = jGs+j-1
            !
            SO(i,j,kw) = dr(is,js,hx,hy)*xh
         ENDDO
      ENDDO

C ----------------------------------------
C     Initial Guess
C ----------------------------------------

      DO j=1, NLy_g
         DO i=1, NLx_g
            q(i,j)=BMG_rand(rZERO)
         ENDDO
      ENDDO

C ---------------------------------------
C     Boundary Conditions
C ---------------------------------------

      !
      !  West boundary
      !
      IF ( iGs.EQ.1 ) THEN

         DO j=2, jEND
            SO(2,j,ks)  = rHALF*SO(2,j,ks)
         ENDDO


         IF ( iBL.EQ.1) THEN

            DO j=2, jEND
               SO(2,j,ko) = SO(2,j,ko) + hy
            ENDDO

         ENDIF

      ENDIF


      !
      !  East boundary
      !
      IF ( iGf.EQ.NGx ) THEN

         DO j=2, jEND
            SO(i1,j,ks) = rHALF*SO(i1,j,ks)
         ENDDO

         IF ( iBR.EQ.1) THEN

            DO j=2, jEND
               SO(i1,j,ko) = SO(i1,j,ko) + hy
            ENDDO
            
         ENDIF

      ENDIF


      !
      !  Bottom Boundary
      !
      IF ( jGs.EQ.1 ) THEN
         
         DO i=2, iEND
            SO(i,2,kw)  = rHALF*SO(i,2,kw)
         ENDDO

         IF ( iBB.EQ.1) THEN
            
            DO i=2, iEND
               SO(i,2,ko) = SO(i,2,ko) + hx   
            ENDDO
            
         ENDIF
         
      ENDIF


      !
      !  Top boundary
      !
      IF ( jGf.EQ.NGy ) THEN 
         
         DO i=2, iEND
            SO(i,j1,kw) = rHALF*SO(i,j1,kw)
         ENDDO

         IF( iBT.EQ.1) THEN

            DO i=2, iEND
               SO(i,j1,ko) = SO(i,j1,ko) + hx
            ENDDO
            
         ENDIF

      ENDIF


      DO j=2, j1
         DO i=2, i1
            cm=rONE
            IF ( ( i.EQ.2 .AND. iGs.EQ.1 ) 
     &           .OR. ( i.EQ.i1 .AND. iGf.EQ.NGx ) ) THEN
               cm = rHALF*cm
            ENDIF
            IF ( ( j.EQ.2 .AND. jGs.EQ.1 ) 
     &           .OR. ( j.EQ.j1 .AND. jGf.EQ.NGy ) ) THEN
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


