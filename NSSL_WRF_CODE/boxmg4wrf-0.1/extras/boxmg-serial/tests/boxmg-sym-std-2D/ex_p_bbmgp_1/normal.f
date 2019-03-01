
      subroutine normal (q,dphi,dtheta,l,m)
c
c------------------------------------------------------------------------
c
c   This subroutine normalizes the input vector by subtracting off the
c   null space component.
c
c------------------------------------------------------------------------
c
      implicit none

      integer l, m, i, j
      real*8  dphi, dtheta, q(0:l+1,0:m+1), sum 
c
c--  normalize the right hand side if the problem is singular.
c
c
      sum = 0
      do 100 j = 1, m
        do 100 i = 1, l
          sum = sum + q(i,j)
 100  continue
      sum = sum/(l*m)
      write(*,120) 'The null space component =',sum
      if (sum .ne. 0) then
        do 110 j = 1, m
          do 110 i = 1, l
            q(i,j) = q(i,j) - sum
 110    continue
      endif
c
 120  FORMAT(/,2X,A,1X,1P,E21.14,/)
c
      return
      end
