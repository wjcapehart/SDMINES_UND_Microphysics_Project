C     /DC2/ is the common block for putf.
C
C
      INTEGER    iREGmax
      PARAMETER  (iREGmax=20)

      INTEGER    iBL, iBR, iByB, iByT, iBzB, iBzT, iREG 
      REAL*8     Dix(iREGmax), Diy(iREGmax), Diz(iREGmax), 
     &           Fi(iREGmax), Si(iREGmax),
     &           X1(iREGmax), X2(iREGmax), 
     &           Y1(iREGmax), Y2(iREGmax), 
     &           Z1(iREGmax), Z2(iREGmax)

      COMMON     /iDC2/ iREG, iBL, iBR, iByB, iByT, iBzB, iBzT
      COMMON     /rDC2/ X1, X2, Y1, Y2, Z1, Z2, Dix, Diy, Diz, Si, Fi


