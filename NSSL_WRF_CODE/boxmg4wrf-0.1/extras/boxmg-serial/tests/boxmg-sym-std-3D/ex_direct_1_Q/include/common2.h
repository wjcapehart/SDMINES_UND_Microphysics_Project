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
      PARAMETER  (iREGmax=1000)

      INTEGER    iBL, iBR, iByB, iByT, iBzB, iBzT, iREG 
      REAL*8     Dix(iREGmax), Diy(iREGmax), Diz(iREGmax), 
     &           Fi(iREGmax), Si(iREGmax),
     &           X1(iREGmax), X2(iREGmax), 
     &           Y1(iREGmax), Y2(iREGmax), 
     &           Z1(iREGmax), Z2(iREGmax)

      REAL*8      xGf, xGs, yGf, yGs, zGf, zGs

      COMMON     /iDC2/ iREG, iBL, iBR, iByB, iByT, iBzB, iBzT
      COMMON     /rDC2/ xGf, xGs, yGf, yGs, zGf, zGs, X1, X2, Y1, 
     &                  Y2, Z1, Z2, Dix, Diy, Diz, Si, Fi

C =======================================================================
