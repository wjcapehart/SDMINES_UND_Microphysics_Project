      SUBROUTINE BMG2_SymStd_SETUP_nog( 
     &                NLx, NLy, NLXYc, NOG, iGs, jGs, MPICOMM 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SETUP_nog computes the number of grids (NOG) based on
C     a local processor coarsening limit NLXYc.
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
      
C ---------------------------
C     Includes:
C
      INCLUDE 'mpif.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  iGs, jGs, MPICOMM, NLx, NLy, NLXYc, NOG

C ---------------------------
C     Local Declarations:
C
      INTEGER  iGs_c, jGs_c, kg, NLx_c, NLy_c, ierror

C ==========================================================================

C ==========================================================================
C -----------------------------------
C     Compute the number of grids:
C -----------------------------------

      !
      ! Initilize coarse grid counter (local)
      !
      kg = 1
      !
      ! Initialize coarse grid dimensions (local)
      !
      iGs_c = iGs
      jGs_c = jGs
      NLx_c = NLx
      NLy_c = NLy

 100  CONTINUE     ! >>>>>>>>>>>> LOOP BOUNDARY: computing NOG

         IF ( MOD(iGs_c,2).EQ.1 ) THEN
            iGs_c = (iGs_c+1)/2
            NLx_c = (NLx_c+1)/2
         ELSE
            iGs_c = iGs_c/2+1
            IF ( MOD(NLx_c,2).EQ.1 ) THEN
               NLx_c = (NLx_c-1)/2
            ELSE
               NLx_c = (NLx_c+1)/2
            ENDIF
         ENDIF

         IF ( MOD(jGs_c,2).EQ.1 ) THEN
            jGs_c = (jGs_c+1)/2
            NLy_c = (NLy_c+1)/2
         ELSE
            jGs_c = jGs_c/2+1
            IF ( MOD(NLy_c,2).EQ.1 ) THEN
               NLy_c = (NLy_c-1)/2
            ELSE
               NLy_c = (NLy_c+1)/2
            ENDIF
         ENDIF

         IF ( MIN(NLx_c,NLy_c).GE.NLXYc ) then
            kg=kg+1
            GOTO 100
         endif

 110  CONTINUE     ! >>>>>>>>>>>> LOOP BOUNDARY: computing NOG

      !
      !  Global minimum is stored in NOG
      !
      CALL MPI_Allreduce( kg, NOG, 1,
     &                    MPI_INTEGER, MPI_MIN, MPICOMM, ierror 
     &                   )

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_nog',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
