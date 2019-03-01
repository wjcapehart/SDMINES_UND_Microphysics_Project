      SUBROUTINE BMG3_SymStd_COPY_cg_WS_RHS(
     &                Q_ws, NLx, NLy, NLz, 
     &                Q, NGxg, NGyg, NGzg,
     &                iGs, jGs, kGs
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C 
C   Copy the RHS for the coarse-grid BoxMG solve.
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

C ----------------------------
C     Includes
C 
      INCLUDE 'BMG_SER_stencils.h'

C ----------------------------
C     Argument Declarations
C
      !
      !  Global/Local indexing
      !
      INTEGER NLx, NLy, NLz,
     &        NGxg, NGyg, NGzg,
     &        iGs, jGs, kGs

      !
      !  Stencil: SER
      !
      REAL*RKIND   Q(NGxg,NGyg,NGzg)
      !
      !  Stencil: WS
      !
      REAL*RKIND   Q_ws(NLx,NLy,NLz)

C ----------------------------
C     Local Declarations
C
      INTEGER  iG, iL, jG, jL, kG, kL

C ==========================================================================

      DO kL=1, NLz
         !
         kG = kGs + kL
         !
         DO jL=1, NLy
            !
            jG = jGs + jL
            !
            DO iL=1, NLx
               !
               iG = iGs + iL
               !
               Q(iG,jG,kG) = Q_ws(iL,jL,kL)
               !
            ENDDO
         ENDDO
      ENDDO

C ==========================================================================

      RETURN
      END
