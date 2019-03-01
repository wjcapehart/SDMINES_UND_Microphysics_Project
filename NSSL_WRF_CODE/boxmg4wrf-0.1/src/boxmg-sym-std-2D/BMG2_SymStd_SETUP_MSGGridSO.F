      SUBROUTINE BMG2_SymStd_SETUP_MSGGridSO(
     &                NGx, NGy, 
     &                LocalArraySize, GlobalCoordLocalData,
     &                GlobalCoordActData, ActDataStart,
     &                DimX, DimY, ProcGrid, 
     &                NProc, NProcI, NProcJ, NOGm, KG,
     &                MPICOMM 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SETUP_MSGGrid is used to setup a parallel grid for the
C     use with MSG.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C   -------------------------
C    Global Grid Dimensions:
C   -------------------------
C
C     NGX      Number of points in x-direction (excluding ghost points)
C     NGY      Number of points in y-direction (excluding ghost points)
C
C   -------------------------
C    Local Grid Dimensions:
C   -------------------------
C
C     NLX      Number of points in x-direction (excluding ghost points)
C     NLY      Number of points in y-direction (excluding ghost points)
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
C     iGs      Global i coordinate of the local grid's lower left corner.
C     jGs      Global j coordinate of the local grid's lower left corner.
C
C ==========================================================================
  
      IMPLICIT   NONE

C -----------------------------
C     Includes
C
      INCLUDE    'BMG_workspace.h'
#include      "BMG_constants.h"

C ---------------------------
C     Argument Declarations:
C
      INTEGER  NGx, NGy, MPICOMM
      INTEGER  NProc, NProcI, NProcJ, NOGm, KG
 
      INTEGER  LocalArraySize(3,NProc), ActDataStart(3,NProc)
      INTEGER  GlobalCoordLocalData(2,3,NProc)
      INTEGER  GlobalCoordActData(2,3,NProc)
      INTEGER  DimX(NProcI, NOGm), DimY(NProcJ, NOGm)
      INTEGER  ProcGrid(NProcI,NProcJ)

C --------------------------
C     Local Declarations:
C
      INTEGER  I,J, II, JJ, iGs, jGs, IJRank

C ========================================================================

C --------------------------------------------
C    initialize the MSG grid information
C --------------------------------------------


      DO J=1, NProcJ
         DO I=1, NProcI

            iGs = 1
            DO II=1,I-1
               iGs = iGs + DimX(II,KG)
            END DO

            jGs = 1
            DO JJ=1, J-1
               jGs = jGs + DimY(JJ,KG)
            END DO

            IJRank = ProcGrid(I,J)
           
            ActDataStart(1,IJRank) = 1
            ActDataStart(2,IJRank) = 1
            ActDataStart(3,IJRank) = 1


            LocalArraySize(1,IJRank) = DimX(I,KG)+3
            LocalArraySize(2,IJRank) = DimY(J,KG)+3
            LocalArraySize(3,IJRank) = 1


            GlobalCoordLocalData(1,1,IJRank) = iGs+1
            GlobalCoordLocalData(1,2,IJRank) = jGs+1
            GlobalCoordLocalData(1,3,IJRank) = 1

            GlobalCoordLocalData(2,1,IJRank) =
     &           GlobalCoordLocalData(1,1,IJRank)
     &           + LocalArraySize(1,IJRank)-4
            GlobalCoordLocalData(2,2,IJRank) =
     &           GlobalCoordLocalData(1,2,IJRank)
     &           + LocalArraySize(2,IJRank)-4
            GlobalCoordLocalData(2,3,IJRank) = 1


            GlobalCoordActData(1,1,IJRank) = iGs
            GlobalCoordActData(1,2,IJRank) = jGs
            GlobalCoordActData(1,3,IJRank) = 1
            
            GlobalCoordActData(2,1,IJRank) = 
     &           GlobalCoordActData(1,1,IJRank)
     &           + LocalArraySize(1,IJRank)-1
            GlobalCoordActData(2,2,IJRank) = 
     &           GlobalCoordActData(1,2,IJRank)
     &           + LocalArraySize(2,IJRank)-1
            GlobalCoordActData(2,3,IJRank) = 1
            
         END DO
      END DO
      
      RETURN
      END


