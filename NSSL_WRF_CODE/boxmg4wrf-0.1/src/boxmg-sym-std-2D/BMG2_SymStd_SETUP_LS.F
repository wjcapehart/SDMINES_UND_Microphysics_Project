      SUBROUTINE BMG2_SymStd_SETUP_LS( 
     &                       iWorkMSG, NMSGi, pMSG, pLS, pSI_MSG,
     &                       ProcGrid, NProcI, NProcJ, NOG
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SETUP_LS creates pointers into the integer work space 
C     describing the storage of MSG communication setup
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
      INCLUDE 'mpif.h'
      INCLUDE 'MSG.h'

      INCLUDE 'BMG_workspace.h'

C ---------------------------
C    Argument Declarations:
C
      INTEGER NMSGi, NOG, NProcI, NProcJ
      INTEGER pMSG(NBMG_pMSG,NOG), iWorkMSG(NMSGi)
      INTEGER pLS(NBMG_pLS,NOG), pSI_MSG
      INTEGER ProcGrid(NProcI,NProcJ)

C --------------------------
C     Local Declarations:
C
      INTEGER kg
      
C ==========================================================================


      CALL BMG2_SymStd_SETUP_PtrLS( 
     &                 pLS, pSI_MSG, NProcI, NProcJ, NOG
     &                 )


      DO kg=NOG, 1, -1
         !
         CALL BMG2_SymStd_SETUP_LSGrid(
     &        iWorkMSG(pMSG(ipL_MSG_GlobalCoordLocalData,kg)),
     &        NProcI, NProcJ, ProcGrid,
     &        iWorkMSG(pLS(ipL_LS_XDataDist,kg)),
     &        iWorkMSG(pLS(ipL_LS_YDataDist,kg)))
         !
      ENDDO

      return
      end



