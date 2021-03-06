      SUBROUTINE BMG3_SymStd_SETUP_MSGGrid(
     &                NGX, NGY, NGz,
     &                OFFX, OFFY, OFFZ,
     &                LocalArraySize, GlobalCoordLocalData,
     &                GlobalCoordActData, ActDataStart,
     &                DimX, DimY, DimZ, ProcGrid, 
     &                NProc, NProcI, NProcJ, NProcK, NOGm, KG,
     &                MPICOMM 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_SETUP_MSGGrid is used to setup a parallel grid for the
C     use with MSG
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
C     NGZ      Number of points in z-direction (excluding ghost points)
C
C   -------------------------
C    Local Grid Dimensions:
C   -------------------------
C
C     NLX      Number of points in x-direction (excluding ghost points)
C     NLY      Number of points in y-direction (excluding ghost points)
C     NLZ      Number of points in z-direction (excluding ghost points)
C
C   -------------------------------------------------
C    Offsets ( related to number of ghost planes )
C   -------------------------------------------------
C
C     OFFX     Number of offsets in the x-direction
C     OFFY     Number of offsets in the y-direction
C     OFFZ     Number of offsets in the z-direction
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
C     kGs      Global k coordinate of the local grid's lower left corner.
C
C ==========================================================================
  
      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE  'BMG_workspace.h'
#include    "BMG_constants.h"

C ---------------------------
C     Argument Declarations:
C
      INTEGER  NGX, NGY, NGz
      INTEGER  OFFX, OFFY, OFFZ
      INTEGER  MPICOMM
      INTEGER  NProc, NProcI, NProcJ, NProcK, NOGm, KG 

      INTEGER  LocalArraySize(3,NProc), ActDataStart(3,NProc)
      INTEGER  GlobalCoordLocalData(2,3,NProc)
      INTEGER  GlobalCoordActData(2,3,NProc)

      INTEGER  DimX(NprocI,NOGm), DimY(NProcJ,NOGm), DimZ(NProcK,NOGm)
      INTEGER  ProcGrid(NProcI, NProcJ, NProcK)
      
C --------------------------
C     Local Declarations:
C
      INTEGER I, J, K, II, JJ, KK, IJKRank, iGs, jGs, kGs

C ========================================================================

C --------------------------------------------
C    initialize the MSG grid information
C --------------------------------------------

      DO K=1, NProcK
         DO J=1, NProcJ
            DO I=1, NprocI

               iGs = 1
               DO II=1,I-1
                  iGs = iGs + DimX(II,KG)
               END DO
               
               jGs = 1
               DO JJ=1, J-1
                  jGs = jGs + DimY(JJ,KG)
               END DO

               kGs = 1
               DO KK=1,K-1
                  kGs = kGs + DimZ(KK,KG)
               END DO
               
               IJKRank = ProcGrid(I,J,K)
               
               
               LocalArraySize(1,IJKRank) = DimX(I,KG)+2+OFFX
               LocalArraySize(2,IJKRank) = DimY(J,KG)+2+OFFY
               LocalArraySize(3,IJKRank) = DimZ(K,KG)+2+OFFZ
               

               ActDataStart(1,IJKRank) = 1
               ActDataStart(2,IJKRank) = 1
               ActDataStart(3,IJKRank) = 1


               GlobalCoordLocalData(1,1,IJKRank) = iGs+1
               GlobalCoordLocalData(1,2,IJKRank) = jGs+1
               GlobalCoordLocalData(1,3,IJKRank) = kGs+1


               GlobalCoordLocalData(2,1,IJKRank)
     &              = GlobalCoordLocalData(1,1,IJKRank)
     &              + LocalArraySize(1,IJKRank)-3-OFFX
               GlobalCoordLocalData(2,2,IJKRank)
     &              = GlobalCoordLocalData(1,2,IJKRank)
     &              + LocalArraySize(2,IJKRank)-3-OFFY
               GlobalCoordLocalData(2,3,IJKRank)
     &              = GlobalCoordLocalData(1,3,IJKRank)
     &              + LocalArraySize(3,IJKRank)-3-OFFZ

               IF( OFFX.EQ.0 ) THEN
                  IF (iGs .eq. 1) THEN
                     GlobalCoordLocalData(1,1,IJKRank) = 
     &                    GlobalCoordLocalData(1,1,IJKRank) - 1
                  ENDIF

                  IF (NGX .eq. iGs+DimX(I,KG)+1) THEN
                     GlobalCoordLocalData(2,1,IJKRank) = 
     &                    GlobalCoordLocalData(2,1,IJKRank) + 1
                  ENDIF
               ENDIF
         
               IF( OFFY.EQ.0 ) THEN
                  IF (jGs .eq. 1) THEN
                     GlobalCoordLocalData(1,2,IJKRank) = 
     &                    GlobalCoordLocalData(1,2,IJKRank) - 1
                  ENDIF
                  
                  IF (NGY .eq. jGs+DimY(J,KG)+1) THEN
                     GlobalCoordLocalData(2,2,IJKRank) =
     &                    GlobalCoordLocalData(2,2,IJKRank) + 1
                  ENDIF
               ENDIF
               
               IF( OFFZ.EQ.0 ) THEN
                  IF (kGs .eq. 1) THEN
                     GlobalCoordLocalData(1,3,IJKRank) = 
     &                    GlobalCoordLocalData(1,3,IJKRank) - 1
                  ENDIF
                  
                  IF (NGZ .eq. kGs+DimZ(K,KG)+1) THEN
                     GlobalCoordLocalData(2,3,IJKRank) =
     &                    GlobalCoordLocalData(2,3,IJKRank) + 1
                  ENDIF
               ENDIF
               


               GlobalCoordActData(1,1,IJKRank) = iGs
               GlobalCoordActData(1,2,IJKRank) = jGs
               GlobalCoordActData(1,3,IJKRank) = kGs
               
               GlobalCoordActData(2,1,IJKRank) 
     &              = GlobalCoordActData(1,1,IJKRank)
     &              + LocalArraySize(1,IJKRank)-1
               GlobalCoordActData(2,2,IJKRank) 
     &              = GlobalCoordActData(1,2,IJKRank)
     &              + LocalArraySize(2,IJKRank)-1
               GlobalCoordActData(2,3,IJKRank) 
     &              = GlobalCoordActData(1,3,IJKRank)
     &              + LocalArraySize(3,IJKRank)-1
               

            END DO
         END DO
      END DO

C ==========================================================================

      RETURN
      END


