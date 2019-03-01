      SUBROUTINE BMG2_SymStd_restrict( 
     &                KF, KC, NOG,
     &                Q, QC, CI,
     &                Nx, Ny, Nxc, Nyc, iGs, jGs,
     &                iWorkMSG, NMSGi, pMSG, 
     &                MSG_Buffer, NMSGr, MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_restrict computes the restriction of a vector on the
C     fine grid, Q, to a vector on the coarse grid, QC.  The weights
C     involve the transpose of the interpolation operator from the
C     coarse grid to the fine grid.
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

#include      "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER  NMSGi, NMSGr, NOG, Nx, Nxc, Ny, Nyc

      INTEGER  iGs, jGs, KC, KF, iWorkMSG(NMSGi), 
     &         pMSG(NBMG_pMSG,NOG), MPICOMM
      REAL*RKIND   CI(Nxc,Nyc,8), Q(Nx,Ny), QC(Nxc,Nyc),
     &         MSG_Buffer(NMSGr)

C ----------------------------
C     Local Declarations
C
      INTEGER  i, ic, j, jc, ISTART, JSTART

C ==========================================================================
      
C ---------------------------------------------
C     Restrict the vector Q -> QC:
C ---------------------------------------------


      IF (mod(iGs,2).eq.1) THEN
         ISTART = 0
      ELSE
         ISTART = 1
      ENDIF

      IF (mod(jGs,2).eq.1) THEN
         JSTART = 0
      ELSE
         JSTART = 1
      ENDIF


      j=JSTART
      DO jc=2, Nyc-1
         j=j+2
         i=ISTART
         DO ic=2, Nxc-1
            i=i+2
            QC(ic,jc) = Ci(ic,jc,LNE)*Q(i-1,j-1)
     &                + Ci(ic,jc,LA)*Q(i,j-1)
     &                + Ci(ic+1,jc,LNW)*Q(i+1,j-1)
     &                + Ci(ic,jc,LR)*Q(i-1,j)
     &                + Q(i,j)
     &                + Ci(ic+1,jc,LL)*Q(i+1,j)
     &                + Ci(ic,jc+1,LSE)*Q(i-1,j+1)
     &                + Ci(ic,jc+1,LB)*Q(i,j+1)
     &                + Ci(ic+1,jc+1,LSW)*Q(i+1,j+1)
          ENDDO 
       ENDDO

       ! Note: no update of the ghost bdry of QC is necessary

c       ptrn = 1
       
c       call MSG_tbdx_send(QC, MSG_Buffer, 
c     &      iWorkMSG(pMSG(ipL_MSG_NumAdjProc,KC)),
c     &      iWorkMSG(pMSG(ipL_MSG_Proc,KC)),
c     &      iWorkMSG(pMSG(ipL_MSG_Ipr,KC)),
c     &      iWorkMSG(pMSG(ipL_MSG_Index,KC)),
c     &      ptrn, ierror)

c       call MSG_tbdx_receive(QC, MSG_Buffer,
c     &      iWorkMSG(pMSG(ipL_MSG_NumAdjProc,KC)),
c     &      iWorkMSG(pMSG(ipL_MSG_Proc,KC)),
c     &      iWorkMSG(pMSG(ipL_MSG_Ipr,KC)),
c     &      iWorkMSG(pMSG(ipL_MSG_Index,KC)),
c     &      ptrn, ierror)
       
c       call MSG_tbdx_close(QC, MSG_Buffer,
c     &      iWorkMSG(pMSG(ipL_MSG_NumAdjProc,KC)),
c     &      iWorkMSG(pMSG(ipL_MSG_Proc,KC)),
c     &      iWorkMSG(pMSG(ipL_MSG_Ipr,KC)),
c     &      iWorkMSG(pMSG(ipL_MSG_Index,KC)),
c     &      ptrn, ierror)


C ==========================================================================

       RETURN
       END
