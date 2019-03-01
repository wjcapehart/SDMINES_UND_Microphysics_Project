      SUBROUTINE BMG3_SymStd_COPY_cg_WS_SO(
     &                SO_ws, NLx, NLy, NLz, 
     &                SO, NGxg, NGyg, NGzg,
     &                i_WS, iGs, jGs, kGs
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C 
C   Copy the stencil for the coarse-grid BoxMG solve after it has
C   been stored in SO_ws during the collection to a single processor.
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
     &        iGs, jGs, kGs, i_WS

      !
      !  Stencil: SER
      !
      REAL*RKIND   SO(NGxg,NGyg,NGzg,14)
      !
      !  Stencil: WS
      !
      REAL*RKIND   SO_ws(14,NLx,NLy,NLz)

C ----------------------------
C     Local Declarations
C
      INTEGER  iG, iL, jG, jL, kG, kL

C ==========================================================================

      DO kL=1, NLz

         kG = kGs + kL

         DO jL=1, NLy

            jG = jGs + jL

            DO iL=1, NLx

               iG = iGs + iL

               SO(iG  ,jG  ,kG  ,kp)   = SO_ws(1,iL,jL,kL)
               SO(iG  ,jG  ,kG  ,kpw)  = SO_ws(2,iL,jL,kL)
               SO(iG+1,jG  ,kG  ,kpnw) = SO_ws(3,iL,jL,kL)
               SO(iG  ,jG  ,kG  ,kps)  = SO_ws(4,iL,jL,kL)
               SO(iG  ,jG  ,kG  ,kpsw) = SO_ws(5,iL,jL,kL)
               SO(iG+1,jG+1,kG  ,kbne) = SO_ws(6,iL,jL,kL)
               SO(iG  ,jG+1,kG  ,kbn)  = SO_ws(7,iL,jL,kL)
               SO(iG  ,jG+1,kG  ,kbnw) = SO_ws(8,iL,jL,kL)
               SO(iG+1,jG  ,kG  ,kbe)  = SO_ws(9,iL,jL,kL)
               SO(iG  ,jG  ,kG  ,kb)   = SO_ws(10,iL,jL,kL)
               SO(iG  ,jG  ,kG  ,kbw)  = SO_ws(11,iL,jL,kL)
               SO(iG+1,jG  ,kG  ,kbse) = SO_ws(12,iL,jL,kL)
               SO(iG  ,jG  ,kG  ,kbs)  = SO_ws(13,iL,jL,kL)
               SO(iG  ,jG  ,kG  ,kbsw) = SO_ws(14,iL,jL,KL)

            ENDDO
         ENDDO
      ENDDO

C ==========================================================================

      RETURN
      END
