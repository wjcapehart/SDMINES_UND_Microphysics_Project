      SUBROUTINE PUTF( SO, QF, Q,
     &                 NLx, NLy, NLz, NGx, NGy, NGz, 
     &                 iGs, jGs, kGs, hx, hy, hz, kg 
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
      INTEGER  iGs, jGs, kGs, NLx, NLy, NLz, NGx, NGy, NGz, kg
      REAL*8   hx, hy, hz, Q(NLx+2,NLy+2,NLz+2),
     &         QF(NLx+2,NLy+2,NLz+2), SO(NLx+3,NLy+3,NLz+3,4)

C ----------------
C     Local
C
      INTEGER  i, i1, i2, iGf, j, j1, j2, jGf, k, k1, k2, kGf, 
     &         iBEG, iEND, jBEG, jEND, kBEG, kEND, 
     &         NLx_g, NLy_g, NLz_g, NLSOx_g, NLSOy_g, NLSOz_g,
     &	       is, js, kss
      REAL*8   cm, fs, h2, s, xh, yh, zh

C -----------------
C     External
C
      REAL*8   da, dr, dz, BMG_rand
      EXTERNAL da, dr, dz, BMG_rand

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
      NLz_g = NLz + 2

      i1 = NLx + 1
      j1 = NLy + 1
      k1 = NLz + 1

      i2 = NLx
      j2 = NLy
      k2 = NLz

      !
      !  Global indexing
      !
      iGf=iGs+i2-1
      jGf=jGs+j2-1
      kGf=kGs+k2-1

      !
      !  Stencil halo with depth = [1,2]
      !
      NLSOx_g = NLx + 3
      NLSOy_g = NLy + 3
      NLSOz_g = NLz + 3

      !
      !
      !  Grid Spacing
      !
      h2=hx*hy*hz

      xh=hy*hz/hx
      yh=hx*hz/hy
      zh=hx*hy/hz

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

      IF ( kGs.EQ.1 ) THEN
         kBEG = 3
      ELSE
         kBEG = 2
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

      IF ( kGf.EQ.NGz ) THEN
         kEND = k1
      ELSE
         kEND = NLz_g
      ENDIF
         
C ----------------------------------------
C     Zero solution
C ----------------------------------------

      DO k=1, NLz_g
         DO j=1, NLy_g
            DO i=1, NLx_g

               Q(i,j,k)   = rZERO
               QF(i,j,k)  = rZERO

            END DO
         END DO
      END DO

C ----------------------------------------
C     Zero Stencil
C ----------------------------------------

      DO k=1, NLSOz_g
         DO j=1, NLSOy_g
            DO i=1, NLSOx_g

               SO(i,j,k,kp)  = rZERO
               SO(i,j,k,kpw) = rZERO
               SO(i,j,k,kps) = rZERO
               SO(i,j,k,kb)  = rZERO

            END DO
         END DO
      END DO

C ----------------------------------------
C     Compute the source
C ----------------------------------------

!      DO k=2, k1
!         DO j=2, j1
!            DO i=2, i1
!	        !
!               is = iGs+i-1
!               js = jGs+j-1
!               kss = kGs+k-1
!               !
!               CALL RHS(i,j,k,hx,hy,hz,fs,s)
!               !
!               qf(i,j,k)    = fs*h2
!               !
!            ENDDO
!         ENDDO
!      ENDDO

C ----------------------------------------
C     Compute the removal/absorption
C ----------------------------------------

!      DO k=2, k1
!         DO j=2, j1
!            DO i=2, i1
!	        !
!               is = iGs+i-1
!               js = jGs+j-1
!               kss = kGs+k-1
!               !
!               !
!               CALL RHS(i,j,k,hx,hy,hz,fs,s)
!               !
!               SO(i,j,k,kp) = s*h2
!               !
!            ENDDO
!         ENDDO
!      ENDDO

C ----------------------------------------
C     Stencil coefficients
C ----------------------------------------
      !
      !  South
      !
      DO k=2, kEND
         DO j=jBEG,jEND
            DO i=2, iEND
               !
               is = iGs+i-1
               js = jGs+j-1
               kss = kGs+k-1
               !
               SO(i,j,k,kps) = da(i,j,k,hx,hy,hz)*yh
            ENDDO
         ENDDO
      ENDDO

      !
      !  West
      !
      DO k=2, kEND
         DO j=2, jEND
            DO i=iBEG, iEND
               !
               is = iGs+i-1
               js = jGs+j-1
               kss = kGs+k-1
               !
               SO(i,j,k,kpw) = dr(i,j,k,hx,hy,hz)*xh
            ENDDO
         ENDDO
      ENDDO

      !
      !  Bottom
      !
      DO k=kBEG, kEND
         DO j=2, jEND
            DO i=2, iEND
               !
               is = iGs+i-1
               js = jGs+j-1
               kss = kGs+k-1
               !
               SO(i,j,k,kb) = dz(i,j,k,hx,hy,hz)*zh
            ENDDO
         ENDDO
      ENDDO

C ----------------------------------------
C     Initial Guess
C ----------------------------------------

      DO k=2, k1
         DO j=2, j1
            DO i=2, i1
               q(i,j,k)=BMG_rand(rZERO)
            ENDDO
         ENDDO
      ENDDO

C ---------------------------------------
C     Boundary Conditions
C ---------------------------------------

      !
      !  West boundary
      !
      IF ( iGs.EQ.1 ) THEN

         DO k=2, kEND
            DO j=2, jEND
               SO(2,j,k,kps)  = rHALF*SO(2,j,k,kps)
               SO(2,j,k,kb)   = rHALF*SO(2,j,k,kb)
            ENDDO
         ENDDO


         IF ( ibl.EQ.1) THEN

            DO j=2, jEND
               DO k=2, kEND
                  SO(2,j,k,kp) = SO(2,j,k,kp) + hy*hz
               ENDDO
            ENDDO

         ENDIF

      ENDIF


      !
      !  East boundary
      !
      IF ( iGf.EQ.NGx ) THEN

         DO k=2, kEND
            DO j=2, jEND
               SO(i1,j,k,kps) = rHALF*SO(i1,j,k,kps)
               SO(i1,j,k,kb)  = rHALF*SO(i1,j,k,kb)
            ENDDO
         ENDDO

         IF ( ibr.EQ.1) THEN

            DO k=2, kEND
               DO j=2, jEND
                  SO(i1,j,k,kp) = SO(i1,j,k,kp) + hy*hz
               ENDDO
            ENDDO
            
         ENDIF

      ENDIF


      !
      !  South boundary
      !
      IF ( jGs.EQ.1 ) THEN
         
         DO k=2, kEND
            DO i=2, iEND
               SO(i,2,k,kpw)  = rHALF*SO(i,2,k,kpw)
               SO(i,2,k,kb)   = rHALF*SO(i,2,k,kb)
            ENDDO
         ENDDO

         IF ( ibyb.EQ.1 ) THEN

            DO k=2, kEND
               DO i=2, iEND
                  SO(i,2,k,kp) = SO(i,2,k,kp) + hx*hz
               ENDDO
            ENDDO

         ENDIF

      ENDIF


      !
      !  North boundary
      !
      IF ( jGf.EQ.NGy ) THEN

         DO k=2, kEND
            DO i=2, iEND
               SO(i,j1,k,kpw) = rHALF*SO(i,j1,k,kpw)
               SO(i,j1,k,kb)  = rHALF*SO(i,j1,k,kb)
            ENDDO
         ENDDO

         IF (ibyt.EQ.1) THEN

            DO k=2, kEND
               DO i=2, iEND
                  SO(i,j1,k,kp) = SO(i,j1,k,kp) + hx*hz
               ENDDO
            ENDDO
            
         ENDIF


      ENDIF

      !
      !  Bottom Boundary
      !
      IF ( kGs.EQ.1 ) THEN
         
         DO j=2, jEND
            DO i=2, iEND
               SO(i,j,2,kpw)  = rHALF*SO(i,j,2,kpw)
               SO(i,j,2,kps)  = rHALF*SO(i,j,2,kps)
            ENDDO
         ENDDO

         IF ( ibzb.EQ.1) THEN

            DO j=2, jEND
               DO i=2, iEND
                  SO(i,j,2,kp) = SO(i,j,2,kp) + hx*hy
               ENDDO
            ENDDO
            
         ENDIF
         
      ENDIF


      !
      !  Top boundary
      !
      IF ( kGf.EQ.NGz ) THEN 
         
         DO j=2, jEND
            DO i=2, iEND
               SO(i,j,k1,kpw) = rHALF*SO(i,j,k1,kpw)
               SO(i,j,k1,kps) = rHALF*SO(i,j,k1,kps)
            ENDDO
         ENDDO

         IF( ibzt.EQ.1) THEN

            DO j=2, jEND
               DO i=2, iEND
                  SO(i,j,k1,kp) = SO(i,j,k1,kp) + hx*hy
               ENDDO
            ENDDO
            
         ENDIF

      ENDIF




      DO k=2, k1
         DO j=2, j1
            DO i=2, i1
               cm=rONE
               IF ( ( i.EQ.2 .AND. iGs.EQ.1 ) 
     &              .OR. ( i.EQ.i1 .AND. iGf.EQ.NGx ) ) THEN
                  cm = rHALF*cm
               ENDIF
               IF ( ( j.EQ.2 .AND. jGs.EQ.1 ) 
     &              .OR. ( j.EQ.j1 .AND. jGf.EQ.NGy ) ) THEN
                  cm = rHALF*cm
               ENDIF
               IF ( ( k.EQ.2 .AND. kGs.EQ.1 ) 
     &              .OR. ( k.EQ.k1 .AND. kGf.EQ.NGz ) ) THEN
                  cm = rHALF*cm
               ENDIF
               !
               SO(i,j,k,kp) = cm*SO(i,j,k,kp)
     &                      + SO(i,j,k,kpw)+SO(i,j+1,k,kps)
     &                      + SO(i+1,j,k,kpw)+SO(i,j,k,kps)
     &                      + SO(i,j,k,kb)+SO(i,j,k+1,kb)
               !
               qf(i,j,k) = cm*qf(i,j,k)
               !
            ENDDO
         ENDDO
      ENDDO

C ==========================================================================

      RETURN
      END


