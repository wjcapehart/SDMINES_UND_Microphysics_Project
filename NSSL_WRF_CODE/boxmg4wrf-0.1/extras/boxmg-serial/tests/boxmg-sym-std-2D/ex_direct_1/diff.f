      subroutine diff (l, m, deltax, deltay, soln, mxnorm, OUT)
c
c--- This routine compares the exact solution to the computed solution
c--- and computes the max-norm of there difference.
c
      IMPLICIT NONE

      integer  l, m, i, j

      real*8   deltax, deltay, soln(0:l+1,0:m+1), mxnorm, P5, x, y,
     &         es, df
      LOGICAL  OUT
c
      P5=0.5D0
      mxnorm = 0

      IF (OUT) THEN
         write(*,*) ' i   j     x       y      exact soln           ',
     &              'soln         difference'
         write(*,*) '-----------------------------------------------',
     &              '------------------------'
      ENDIF

      do 20 j = 1, m
         y = 1 + deltay*(j-P5)
         do 10 i = 1, l
            x = deltax*(i-P5)
            es = x*x + y*y
            df = abs(soln(i,j)-es)
            if (mxnorm .lt. df) mxnorm = df
            IF (OUT) THEN
               write(*,100) i, j, x, y, es, soln(i,j), df
            ENDIF
 10      continue
 20   continue
      
 100  format(' ',i2,2x,i2,2x,f6.4,2x,f6.4,3x,f10.7,3x,f15.7,3x,e15.7)

      return
      end

