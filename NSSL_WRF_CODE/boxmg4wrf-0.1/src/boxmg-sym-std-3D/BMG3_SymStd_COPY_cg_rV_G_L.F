      SUBROUTINE BMG3_SymStd_COPY_cg_rV_G_L(
     &                Q_SER, NGx, NGy, NGz,
     &                Q, NLx, NLy, NLz, iGs, jGs, kGs,
     &                ProcCoord, NProcI, NProcJ, NProcK,
     &                Nproc, MyProc
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C 
C   Copy the solution of the coarse-grid solve, back into the
C   MPI-based BoxMG solution vector.
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
C     Argument Declarations
C
      !
      !  Global/Local indexing
      !
      INTEGER NLx, NLy, NLz,
     &        NGx, NGy, NGz,
     &        iGs, jGs, kGs

      !
      !  Global Solution: SER
      !
      REAL*RKIND   Q_SER(NGx,NGy,NGz)
      !
      !  Local Solution: MPI
      !
      REAL*RKIND   Q(NLx,NLy,NLz)

      !
      !  Processor grid
      !
      INTEGER NProcI, NProcJ, NProcK, NProc, MyProc
      INTEGER ProcCoord(3, NProc)

C ----------------------------
C     Local Declarations
C
      INTEGER  iG, iL, jG, jL, kG, kL, 
     &         MyProcI, MyProcJ, MyProcK

      INTEGER  iBEG, iEND, jBEG, jEND, kBEG, kEND


C ==========================================================================

      MyProcI = ProcCoord(1,MyProc)
      MyProcJ = ProcCoord(2,MyProc)
      MyProcK = ProcCoord(3,MyProc)
      
      !
      !  Setup loop boundaries in x
      !      
      IF ( MyProcI.EQ.1 ) THEN
         iBEG = 2
      ELSE
         iBEG = 1
      END IF
      
      IF ( MyProcI.EQ.NProcI ) THEN
         iEND = NLx-1
      ELSE
         iEND = NLx
      END IF

      !
      !  Setup loop boundaries in y
      !
      IF ( MyProcJ.EQ.1 ) THEN
         jBEG = 2
      ELSE
         jBEG = 1
      END IF

      IF ( MyProcJ.EQ.NProcJ ) THEN
         jEND = NLy-1
      ELSE
         jEND = NLy
      ENDIF
         
      !
      !  Setup loop boundaries in z
      !
      IF ( MyProcK.EQ.1 ) THEN
         kBEG = 2
      ELSE
         kBEG = 1
      END IF

      IF ( MyProcK.EQ.NProcK) THEN
         kEND = NLz-1
      ELSE
         kEND = NLz
      ENDIF

      DO kL=kBEG, kEND
         kG = kGs + kL - 1
         DO jL=jBEG, jEND
            jG = jGs + jL - 1 
            DO iL=iBEG, iEND
               iG = iGs + iL - 1
               Q(iL,jL,kL) = Q_SER(iG,jG,kG)
            ENDDO
         ENDDO
      ENDDO
               
C ==========================================================================

      RETURN
      END
