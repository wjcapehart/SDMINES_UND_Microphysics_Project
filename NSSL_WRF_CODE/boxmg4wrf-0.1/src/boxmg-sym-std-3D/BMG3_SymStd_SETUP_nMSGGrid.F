      SUBROUTINE BMG3_SymStd_SETUP_nMSGGrid(
     &                NGX, NGY, NGz, CFx, CFy, CFz, 
     &                LXGP, RXGP, LYGP, RYGP, LZGP, RZGP,
     &                LocalArraySize, GlobalCoordLocalData,
     &                GlobalCoordActData, ActDataStart,
     &                DimX, DimY, DimZ, ProcGrid, 
     &                NProc, NProcI, NProcJ, NProcK, NOGm, KG
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
C     Offsets ( related to number of ghost planes )
C   -------------------------------------------------
C
C     LXGP     Number of ghost points on left in x-direction
C     RXGP     Number of ghost points on right in x-direction
C     LYGP     Number of ghost points on left in y-direction
C     RYGP     Number of ghost points on right in y-direction
C     LZGP     Number of ghost points on left in z-direction
C     RZGP     Number of ghost points on right in z-direction
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

C ---------------------------
C     Includes
C
      INCLUDE  'BMG_workspace.h'
#include    "BMG_constants.h"

C ---------------------------
C     Argument Declarations:
C
      INTEGER  NGX, NGY, NGz
      INTEGER  LXGP, RXGP, LYGP, RYGP, LZGP, RZGP
      INTEGER  CFx, CFy, CFz
      INTEGER  NProc, NProcI, NProcJ, NProcK, NOGm, KG

      INTEGER  LocalArraySize(3,NProc), ActDataStart(3,NProc)
      INTEGER  GlobalCoordLocalData(2,3,NProc)
      INTEGER  GlobalCoordActData(2,3,NProc)
      INTEGER  DimX(NProcI,NOGm), DimY(NProcJ,NOGm), DimZ(NProcK,NOGm)
      INTEGER  ProcGrid(NProcI, NProcJ, NProcK)
      
C --------------------------
C     Local Declarations:
C
      INTEGER  I, J, K, iGs, jGs, kGs, IJKRank, II, JJ, KK

C ==========================================================================

C --------------------------------------------
C    initialize the MSG grid information
C --------------------------------------------


      DO K=1, NProcK
         DO J=1, NProcJ
            DO I=1, NProcI

               iGs = 1
               DO II=1,I-1
                  iGs = iGs + DimX(II,KG-CFx)
               END DO
               
               jGs = 1
               DO JJ=1, J-1
                  jGs = jGs + DimY(JJ,KG-CFy)
               END DO

               kGs = 1
               DO KK=1,K-1
                  kGs = kGs + DimZ(KK,KG-CFz)
               END DO

               IJKRank = ProcGrid(I,J,K)



               LocalArraySize(1,IJKRank) = DimX(I,KG-CFx)+LXGP+RXGP
               LocalArraySize(2,IJKRank) = DimY(J,KG-CFy)+LYGP+RYGP
               LocalArraySize(3,IJKRank) = DimZ(K,KG-CFz)+LZGP+RZGP



               ActDataStart(1,IJKRank) = 1
               ActDataStart(2,IJKRank) = 1
               ActDataStart(3,IJKRank) = 1

     
               GlobalCoordLocalData(1,1,IJKRank) = iGs+LXGP
               GlobalCoordLocalData(1,2,IJKRank) = jGs+LYGP
               GlobalCoordLocalData(1,3,IJKRank) = kGs+LZGP


               GlobalCoordLocalData(2,1,IJKRank)
     &              = GlobalCoordLocalData(1,1,IJKRank)
     &              + LocalArraySize(1,IJKRank)-1-LXGP-RXGP
               GlobalCoordLocalData(2,2,IJKRank)
     &              = GlobalCoordLocalData(1,2,IJKRank)
     &              + LocalArraySize(2,IJKRank)-1-LYGP-RYGP
               GlobalCoordLocalData(2,3,IJKRank)
     &              = GlobalCoordLocalData(1,3,IJKRank)
     &              + LocalArraySize(3,IJKRank)-1-LZGP-RZGP

               
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


