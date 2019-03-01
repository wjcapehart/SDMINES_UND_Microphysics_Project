      subroutine checker (soln, tsoln, l, m, dphi, dtheta, ipsbc)
c
c---------------------------------------------------------------------------
c
c     This subroutine checks the computed solution soln against the actual
c     solution, and then print the true error of the computed solutioon.
c
c--------------------------------------------------------------------------
c
c
      implicit none

      integer l, m, i, j, jmax, ipsbc
      real*8  dphi, dtheta, soln(0:l+1,0:m+1), tsoln(0:l+1,0:m+1),
     &        err, errmax, exsoln, theta, phi, P5

      P5=0.5D0
c
c
        write(*,100)
 100    format(' ',/,3x,'i',4x,'j',8x,'true soln',7x,'computed soln',
     &         9x,'error',/,' --------------',
     &         '------------------------------------------------')
c
      do 10 j = 1, m
        theta = (j-P5)*dtheta
        do 10 i = 1, l
          phi = (i-P5)*dphi
          tsoln(i,j) = exsoln(theta,phi)
  10  continue
c
      errmax = 0
      jmax=1
      do 20 j = 1, m
        do 20 i = 1, l
          err = tsoln(i,j) - soln(i,j)
          if (abs(err) .gt. errmax) then 
              errmax = abs(err)
              jmax=j
          endif
          write(*,110) i, j, tsoln(i,j), soln(i,j), err 
  20  continue
c
      write(*,120) 'Maximum error between the computed and true',
     &             'solution =', errmax
c
 110  FORMAT(' ',i3,2x,i3,5x,e13.6,5x,e13.6,5x,e13.6)
 120  FORMAT(/,2X,A,1X,A,1X,1P,E20.14,/)
c
      return
      end
