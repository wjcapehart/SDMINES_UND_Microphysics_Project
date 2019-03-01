      SUBROUTINE BMG3_SymStd_SETUP_PtrMSG( 
     &                       NLX, NLY, NLZ, OFFX, OFFY, OFFZ, 
     &                       MSG_START, NProc, NOG, NOG_fg, NOG_cg,
     &                       pMSG
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
C     NLz      Number of points in z-direction (excluding ghost points)
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
C   ------------------------------
C    Number Additional Ghost Planes:
C   ------------------------------
C
C     OFFX      Number of additional ghost planes in x-direction
C     OFFY      Number of additional ghost planes in y-direction
C     OFFZ      Number of additional ghost planes in z-direction
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

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE  'BMG_workspace.h'
#include    "BMG_constants.h"

C ---------------------------
C     Argument Declarations:
C
      INTEGER  MSG_START, NLx, NLy, NLz, NOG, NOG_cg, NOG_fg,
     &         NProc, OFFX, OFFY, OFFZ

      INTEGER  pMSG(NBMG_pMSG,NOG)

C ---------------------------
C     Local Declarations:
C
      INTEGER  OFFSET, NSTART, N, ptwo

C ==========================================================================

      NSTART = MSG_START
      OFFSET = 0
      ptwo = 1

      pMSG(ipL_MSG_ProcGrid,NOG) = NSTART
      OFFSET = OFFSET + NProc
      
      pMSG(ipL_MSG_ProcGridCoord_x,NOG) = NSTART + OFFSET
      OFFSET = OFFSET + 1
      
      pMSG(ipL_MSG_ProcGridCoord_y,NOG) = NSTART + OFFSET
      OFFSET = OFFSET + 1
      
      pMSG(ipL_MSG_ProcGridCoord_z,NOG) = NSTART + OFFSET
      OFFSET = OFFSET + 1
 
      DO N=NOG_fg, NOG_cg, -1

         pMSG(ipL_MSG_Index,N)                = NSTART + OFFSET
         OFFSET = OFFSET 
     &          + (4+4*OFFZ)*((NLX-3)/ptwo+3+OFFX)*((NLY-3)/ptwo+3+OFFY)
     &          + (4+4*OFFY)*((NLX-3)/ptwo+3+OFFX)*((NLZ-3)/ptwo+3+OFFZ)
     &          + (4+4*OFFX)*((NLY-3)/ptwo+3+OFFY)*((NLZ-3)/ptwo+3+OFFZ)

         pMSG(ipL_MSG_LocalArraySize,N)       = NSTART + OFFSET
         OFFSET = OFFSET + 3*NProc

         pMSG(ipL_MSG_Proc,N)                 = NSTART + OFFSET
         OFFSET = OFFSET + 26
 
         pMSG(ipL_MSG_Ipr,N)                  = NSTART + OFFSET
         OFFSET = OFFSET + 53

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
            pMSG(ipL_MSG_ProcGridCoord_z,N) = 
     &           pMSG(ipL_MSG_ProcGridCoord_z,NOG)
         endif
            
         NSTART = NSTART + OFFSET
         OFFSET = 0
         ptwo = 2*ptwo

      enddo

      MSG_START = NSTART

C ==========================================================================

      RETURN
      END


