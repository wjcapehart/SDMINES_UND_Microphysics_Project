      SUBROUTINE BMG3_SymStd_update_yo( 
     &                yo, IIF, JJC, KKC, 
     &                pMSGyo, iMSG_Geom, NMSGi, IGRD, BUFFER, NMSGr,
     &                KGF, KGC, NOGm, NProc, MyProc, 
     &                DimX, DimY, DimZ, ProcGrid,
     &                NProcI, NProcJ, NProcK,
     &                MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Update halo (ghost points) in the "yo" workspace array for
C     the series of 1D interpolation approximations in the 
C     ITLI_Izyx approximation of the variational coarse-grid operator.
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

C ------------------------------------------------
C     Includes
C
      INCLUDE 'mpif.h'
      INCLUDE 'MSG.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'

C ------------------------------------------------
C     Argument Declarations
C
      INTEGER MPICOMM, MyProc, NMSGi, Nproc, NOG, NOGm,
     &        IIF, JJC, KKC, NMSGr, KGF, KGC, NProcI,
     &        NProcJ, NProcK

      INTEGER IGRD(NOGm,NBMG_pIGRD), iMSG_Geom(NMSGi)
      
      INTEGER pMSGyo(NBMG_pMSG,2)

      REAL*RKIND  yo(IIF,JJC,KKC,14), BUFFER(NMSGr)

      INTEGER DimX(NProcI,NOGm), DimY(NProcJ,NOGm), DimZ(NProcK,NOGm)
      INTEGER ProcGrid(NProcI, NProcJ, NProcK)


C ------------------------------------------------
C     Local Declarations
C
      INTEGER kpz, ptrn, ierror
      INTEGER LXGP, RXGP, LYGP, RYGP, LZGP, RZGP
      
C ==========================================================================

      LXGP = 1
      RXGP = 1
      LYGP = 1
      RYGP = 1
      LZGP = 1
      RZGP = 1

      CALL  BMG3_SymStd_SETUP_nMSGGrid(
     &           IGRD(KGF,idL_BMG_NGx), IGRD(KGC,idL_BMG_NGy),
     &           IGRD(KGC,idL_BMG_NGz), 
     &           iZERO, iONE, iONE,
     &           LXGP, RXGP, LYGP, RYGP, LZGP, RZGP,
     &           iMSG_Geom(pMSGyo(ipL_MSG_LocalArraySize,1)),
     &           iMSG_Geom(pMSGyo(ipL_MSG_GlobalCoordLocalData,1)),
     &           iMSG_Geom(pMSGyo(ipL_MSG_GlobalCoordActData,1)),
     &           iMSG_Geom(pMSGyo(ipL_MSG_ActDataStart,1)),
     &           DimX, DimY, DimZ, ProcGrid, 
     &           NProc, NProcI, NProcJ, NProcK, NOGm, KGF
     &           )

      iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)) = 0

      CALL MSG_tp_setup (
     &        iMSG_Geom(pMSGyo(ipL_MSG_LocalArraySize,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_ActDataStart,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_GlobalCoordLocalData,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_GlobalCoordActData,1)),
     &        NProc, MyProc,
     &        iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Proc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Ipr,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Index,1)),
     &        1,1,ierror)

      DO kpz=1,9

         ptrn = 1
         CALL MSG_tbdx_send(yo(1,1,1,kpz), BUFFER, 
     &        iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Proc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Ipr,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Index,1)),
     &        ptrn, ierror )
            
         CALL MSG_tbdx_receive(yo(1,1,1,kpz), BUFFER,
     &        iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Proc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Ipr,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Index,1)),
     &        ptrn, ierror )
            
         CALL MSG_tbdx_close(yo(1,1,1,kpz), BUFFER,
     &        iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Proc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Ipr,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Index,1)),
     &        ptrn, ierror )

      ENDDO

      DO kpz=13,14

         ptrn = 1
         CALL MSG_tbdx_send(yo(1,1,1,kpz), BUFFER, 
     &        iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Proc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Ipr,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Index,1)),
     &        ptrn, ierror )
            
         CALL MSG_tbdx_receive(yo(1,1,1,kpz), BUFFER,
     &        iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Proc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Ipr,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Index,1)),
     &        ptrn, ierror )
            
         CALL MSG_tbdx_close(yo(1,1,1,kpz), BUFFER,
     &        iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Proc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Ipr,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Index,1)),
     &        ptrn, ierror )

      ENDDO

C =======================================================================

      LXGP = 2
      RXGP = 0
      LYGP = 1
      RYGP = 1
      LZGP = 1
      RZGP = 1

      CALL  BMG3_SymStd_SETUP_nMSGGrid(
     &         IGRD(KGF,idL_BMG_NGx),   IGRD(KGC,idL_BMG_NGy),
     &         IGRD(KGC,idL_BMG_NGz), 
     &         iZERO, iONE, iONE,
     &         LXGP, RXGP, LYGP, RYGP, LZGP, RZGP,
     &         iMSG_Geom(pMSGyo(ipL_MSG_LocalArraySize,1)),
     &         iMSG_Geom(pMSGyo(ipL_MSG_GlobalCoordLocalData,1)),
     &         iMSG_Geom(pMSGyo(ipL_MSG_GlobalCoordActData,1)),
     &         iMSG_Geom(pMSGyo(ipL_MSG_ActDataStart,1)),
     &         DimX, DimY, DimZ, ProcGrid, 
     &         NProc, NProcI, NProcJ, NProcK, NOGm, KGF
     &         )

      iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)) = 0

      CALL MSG_tp_setup (
     &        iMSG_Geom(pMSGyo(ipL_MSG_LocalArraySize,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_ActDataStart,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_GlobalCoordLocalData,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_GlobalCoordActData,1)),
     &        NProc, MyProc,
     &        iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Proc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Ipr,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Index,1)),
     &        1,1,ierror)
      
       DO kpz=10,12

         ptrn = 1
         CALL MSG_tbdx_send(yo(1,1,1,kpz), BUFFER, 
     &        iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Proc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Ipr,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Index,1)),
     &        ptrn, ierror )
            
         CALL MSG_tbdx_receive(yo(1,1,1,kpz), BUFFER,
     &        iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Proc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Ipr,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Index,1)),
     &        ptrn, ierror )
            
         CALL MSG_tbdx_close(yo(1,1,1,kpz), BUFFER,
     &        iMSG_Geom(pMSGyo(ipL_MSG_NumAdjProc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Proc,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Ipr,1)),
     &        iMSG_Geom(pMSGyo(ipL_MSG_Index,1)),
     &        ptrn, ierror )

      ENDDO

C ==========================================================================

      RETURN
      END
