      SUBROUTINE BMG2_SymStd_SETUP_PtrMSG( 
     &                       NLx, NLy, MSG_START, NProc, NOG, pMSG 
     &                       )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SETUP_PtrMSG is used to compute the pointers into the
C     user's integer work array for storage of all necessary information
C     for intra-level communication
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C   --------------------------------------
C    Local Fine-Grid Dimensions:
C   --------------------------------------
C
C     NLx      Number of points in x-direction (excluding ghost points)
C     NLy      Number of points in y-direction (excluding ghost points)
C
C   --------------------------------------
C    Pointer into the integer workspace:
C   --------------------------------------
C
C     MSG_START points to the place in iWORK, where the workspace that is
C               allocated for MSG communication, starts 
C
C   --------------------------------------
C    Dimensions:
C   --------------------------------------
C     
C     NOG       Number of grids needed for the given (NLx,NLy)
C
C   --------------------------------------
C    Workspace:
C   --------------------------------------
C
C     iWork     Integer array that (on exit) contains the data needed to 
C               set up the MSG communication scheme 
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C     IGRD      Grid pointer array for the arrays internal to BOXMG
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C     OFFSET    an offset in the calculation of pointers
C     NSTART    holds the pointer to the beginning of data for a given level
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
      INTEGER  NLx, NLy, MSG_START, NOG, NProc
      INTEGER  pMSG(NBMG_pMSG,NOG), ptwo

C ---------------------------
C     Local Declarations:
C
      INTEGER  OFFSET, NSTART, N

C ========================================================================

      NSTART = MSG_START
      OFFSET = 0
      ptwo = 1

      pMSG(ipL_MSG_ProcGrid,NOG) = NSTART
      OFFSET = OFFSET + NProc
      
      pMSG(ipL_MSG_ProcGridCoord_x,NOG) = NSTART + OFFSET
      OFFSET = OFFSET + 1
      
      pMSG(ipL_MSG_ProcGridCoord_y,NOG) = NSTART + OFFSET
      OFFSET = OFFSET + 1
      

      do N=NOG,1,-1

         pMSG(ipL_MSG_Index,N)                = NSTART + OFFSET
         OFFSET = OFFSET + 8*(NLx/ptwo + NLy/ptwo + 6)

         pMSG(ipL_MSG_LocalArraySize,N)       = NSTART + OFFSET
         OFFSET = OFFSET + 3*NProc

         pMSG(ipL_MSG_Proc,N)                 = NSTART + OFFSET
         OFFSET = OFFSET + 8
 
         pMSG(ipL_MSG_Ipr,N)                  = NSTART + OFFSET
         OFFSET = OFFSET + 17

         pMSG(ipL_MSG_NumAdjProc,N)           = NSTART + OFFSET
         OFFSET = OFFSET + 1
 
         pMSG(ipL_MSG_ActDataStart,N)         = NSTART + OFFSET
         OFFSET = OFFSET + 3*Nproc

         pMSG(ipl_MSG_GlobalCoordLocalData,N) = NSTART + OFFSET
         OFFSET = OFFSET + 2*3*NProc

         pMSG(ipL_MSG_GlobalCoordActData,N)   = NSTART + OFFSET
         OFFSET = OFFSET + 2*3*NProc

         if (N.lt.NOG) then
            pMSG(ipL_MSG_ProcGrid,N) = pMSG(ipL_MSG_ProcGrid,NOG)       
            pMSG(ipL_MSG_ProcGridCoord_x,N) = 
     &           pMSG(ipL_MSG_ProcGridCoord_x,NOG)
            pMSG(ipL_MSG_ProcGridCoord_y,N) = 
     &           pMSG(ipL_MSG_ProcGridCoord_y,NOG)
         endif
            
         NSTART = NSTART + OFFSET
         OFFSET = 0
         ptwo = 2*ptwo

      enddo

      MSG_START = NSTART

C ========================================================================

      RETURN
      END


