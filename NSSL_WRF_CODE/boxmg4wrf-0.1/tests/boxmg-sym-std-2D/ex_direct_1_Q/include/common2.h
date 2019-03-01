C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Region data for creating problems with discontinuous coefficients.
C     
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   VARIABLES:
C  --------------------
C
C
C =======================================================================

      INTEGER    iREGmax
      PARAMETER  (iREGmax=20)

      INTEGER    iBL, iBR, iBT, iBB, iREG 
      REAL*8     Dix(iREGmax), Diy(iREGmax),
     &           Fi(iREGmax), Si(iREGmax),
     &           X1(iREGmax), X2(iREGmax), 
     &           Y1(iREGmax), Y2(iREGmax)

      REAL*8     xGf, xGs, yGf, yGs

      COMMON     /iDC2/ iREG, iBL, iBR, iBB, iBT
      COMMON     /rDC2/ xGs, xGf, yGs, yGf,
     &                  X1, X2, Y1, Y2,
     &                  Dix, Diy, Si, Fi
      
C =======================================================================

