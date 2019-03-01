      SUBROUTINE BMG2_SymStd_SETUP_MSG(
     &                pMSG, pMSGSO, iMSG_Geom, NMSGi, pSI_MSG,
     &                IGRD, NOG, NOGm, NProc, MyProc, 
     &                DimX, DimY, DimXfine, DimYfine,
     &                ProcGrid, NProcI, NProcJ, MPICOMM 
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Create pointers into the integer work space describing the storage 
C     of MSG communication setup.  Note that the pointer shift, pSI_MSG, 
C     is changed in every call.
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
      INTEGER NMSGi, NOG, NOGm, MPICOMM, NProcI, NPRocJ
      INTEGER iMSG_Geom(NMSGi), IGRD(NOGm,NBMG_pIGRD),
     &        pMSG(NBMG_pMSG,NOG), pMSGSO(NBMG_pMSG,NOG)
      INTEGER DimX(NProcI,NOGm), DimY(NProcJ,NOGm)
      INTEGER DimXfine(NProcI), DimYfine(NProcJ)
      INTEGER ProcGrid(NProcI,NProcJ)
      INTEGER NProc, MyProc, N, ierror, pSI_MSG

C --------------------------
C     Local Declarations:
C
      INTEGER I, J, iGs, jGs

C ==========================================================================

C     Note that the pointer shift, pSI_MSG, is changed in every call

      CALL BMG2_SymStd_SETUP_PtrMSG(
     &     IGRD(NOG,idL_BMG_NLx), IGRD(NOG,idL_BMG_NLy), pSI_MSG,
     &     NProc, NOG, pMSG)

      CALL BMG2_SymStd_SETUP_PtrMSGSO(
     &     IGRD(NOG,idL_BMG_NLx), IGRD(NOG,idL_BMG_NLy), pSI_MSG,
     &     NProc, NOG, pMSGSO)


C -------------------------------------------
C     Compute the dimensions for all local
C     grids on all processors
C -------------------------------------------

      DO I=1, NProcI
         DimX(I,NOG) = DimXfine(I)
      END DO

      DO J=1, NProcJ
         DimY(J,NOG) = DimYfine(J)
      END DO

      DO N=NOG-1, 1, -1
         
         iGs = 1
         
         DO I=1, NProcI

            if (mod(iGs,2).eq.1) then 
               DimX(I,N) = (DimX(I,N+1)+1)/2
            else
               IF (mod(DimX(I,N+1),2).eq.1) THEN
                  DimX(I,N) = (DimX(I,N+1)-1)/2
               ELSE
                  DimX(I,N) = (DimX(I,N+1)+1)/2
               ENDIF
            endif
            
            iGs = iGs + DimX(I,N+1)

         END DO

         jGs = 1

         DO J=1, NProcJ

            if (mod(jGs,2).eq.1) then 
               DimY(J,N) = (DimY(J,N+1)+1)/2
            else
               IF (mod(DimY(J,N+1),2).eq.1) THEN
                  DimY(J,N) = (DimY(J,N+1)-1)/2
               ELSE
                  DimY(J,N) = (DimY(J,N+1)+1)/2
               ENDIF
            endif
            
            jGs = jGs + DimY(J,N+1)

         END DO
         
      END DO


C -------------------------------------------
C     Create MSG grid information for U
C -------------------------------------------    

      DO N=NOG, 1, -1
         !
         CALL BMG2_SymStd_SETUP_MSGGrid(
     &        IGRD(N,idL_BMG_NGx), IGRD(N,idL_BMG_NGy),
     &        iMSG_Geom(pMSG(ipL_MSG_LocalArraySize,N)),
     &        iMSG_Geom(pMSG(ipL_MSG_GlobalCoordLocalData,N)),
     &        iMSG_Geom(pMSG(ipL_MSG_GlobalCoordActData,N)),
     &        iMSG_Geom(pMSG(ipL_MSG_ActDataStart,N)),
     &        DimX, DimY, ProcGrid, NProc, NProcI, NProcJ, NOGm, N,
     &        MPICOMM )
         !
      ENDDO

C --------------------------------------------
C     Create MSG grid information for SO
C --------------------------------------------

      DO N=NOG, 1, -1
         !
         CALL BMG2_SymStd_SETUP_MSGGridSO(
     &        IGRD(N,idL_BMG_NGx), IGRD(N,idL_BMG_NGy),
     &        iMSG_Geom(pMSGSO(ipL_MSG_LocalArraySize,N)),
     &        iMSG_Geom(pMSGSO(ipL_MSG_GlobalCoordLocalData,N)),
     &        iMSG_Geom(pMSGSO(ipL_MSG_GlobalCoordActData,N)),
     &        iMSG_Geom(pMSGSO(ipL_MSG_ActDataStart,N)),
     &        DimX, DimY, ProcGrid, NProc, NProcI, NProcJ, NOGm, N,
     &        MPICOMM )
         !
      ENDDO


C
C     Setup MSG communications for all grids
C

      DO N=NOG, 1, -1

c         write(*,*) N, NProc, MyProc

         iMSG_Geom(pMSG(ipL_MSG_NumAdjProc, N)) = 0

c         write(*,*)MyProc,':', (iWork(pMSG(ipL_MSG_ActDataStart,N)+i),
c     &        i=0,3*NProc-1)

         CALL MSG_tp_setup (
     &        iMSG_Geom(pMSG(ipL_MSG_LocalArraySize,N)),
     &        iMSG_Geom(pMSG(ipL_MSG_ActDataStart,N)),
     &        iMSG_Geom(pMSG(ipL_MSG_GlobalCoordLocalData,N)),
     &        iMSG_Geom(pMSG(ipL_MSG_GlobalCoordActData,N)),
     &        NProc, MyProc,
     &        iMSG_Geom(pMSG(ipL_MSG_NumAdjProc,N)),
     &        iMSG_Geom(pMSG(ipL_MSG_Proc,N)),
     &        iMSG_Geom(pMSG(ipL_MSG_Ipr,N)),
     &        iMSG_Geom(pMSG(ipL_MSG_Index,N)),
     &        1,1,ierror)

C         write(*,*) 'after MSG_tp_setup'


C         write(*,*) N,IGRD(N,idL_BMG_ICOORD),IGRD(N,idL_BMG_JCOORD)
         
         
         iMSG_Geom(pMSGSO(ipL_MSG_NumAdjProc, N)) = 0

         CALL MSG_tp_setup (
     &        iMSG_Geom(pMSGSO(ipL_MSG_LocalArraySize,N)),
     &        iMSG_Geom(pMSGSO(ipL_MSG_ActDataStart,N)),
     &        iMSG_Geom(pMSGSO(ipL_MSG_GlobalCoordLocalData,N)),
     &        iMSG_Geom(pMSGSO(ipL_MSG_GlobalCoordActData,N)),
     &        NProc, MyProc,
     &        iMSG_Geom(pMSGSO(ipL_MSG_NumAdjProc,N)),
     &        iMSG_Geom(pMSGSO(ipL_MSG_Proc,N)),
     &        iMSG_Geom(pMSGSO(ipL_MSG_Ipr,N)),
     &        iMSG_Geom(pMSGSO(ipL_MSG_Index,N)),
     &        1,1,ierror)


      ENDDO

      return
      end



