      SUBROUTINE BMG2_SymStd_SETUP_LSGrid(
     &                       GlobalCoordLocalData,
     &                       NProcI, NProcJ, ProcGrid,
     &                       XDataDist, YDataDist
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Setup the line solve grid information.
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
  
      IMPLICIT   NONE

C -----------------------------
C     Includes
C
      INCLUDE    'mpif.h'

      INCLUDE    'BMG_workspace.h'
#include      "BMG_constants.h"

C ---------------------------
C     Argument Declarations:
C
      INTEGER  NProcI, NProcJ
      INTEGER  ProcGrid(NProcI, NProcJ)
      INTEGER  XDataDist(2,NProcI), YDataDist(2,NProcJ)
      INTEGER  GlobalCoordLocalData(2,3,*)

C ---------------------------
C     Local Declarations:
C

      INTEGER  i, j

C ========================================================================

C -----------------------------------------------
C    initialize the Line Solver grid information
C -----------------------------------------------


      DO i=1,NProcI
         XDataDist(1,I) = GlobalCoordLocalData(1,1,ProcGrid(i,1))
         XDataDist(2,I) = GlobalCoordLocalData(2,1,ProcGrid(i,1))
      ENDDO
      XDataDist(1,1) = XDataDist(1,1)+1
      XDataDist(2,NProcI) = XDataDist(2,NProcI)-1
      

      DO j=1,NProcJ
         YDataDist(1,J) = GlobalCoordLocalData(1,2,ProcGrid(1,j))
         YDataDist(2,J) = GlobalCoordLocalData(2,2,ProcGrid(1,j))
      ENDDO
      YDataDist(1,1) = YDataDist(1,1)+1
      YDataDist(2,NProcJ) = YDataDist(2,NprocJ)-1


      RETURN
      END


