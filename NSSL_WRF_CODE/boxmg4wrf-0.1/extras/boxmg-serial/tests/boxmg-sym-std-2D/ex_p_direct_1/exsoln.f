      real*8 function exsoln (theta, phi)
c
c--------------------------------------------------------------------------
c
c   This function computes the exact solution to the problem in spherical
c   coordinates,
c
c    d[ sin(theta) du ]       1       d   [ du ]
c  - ------------------ - ---------- ---- ------  = - sin(theta) f
c         dtheta          sin(theta) dphi  dphi
c
c   where the solution is given by
c
c         u = [ sin(theta) ]**2 cos(phi)
c
c   and the right hand side f was obtained by using the discrete solution
c   in the above discretized PDE.
c   The boundary conditions used are periodic in phi and zero Neumann
c   in theta at -pi/2 and pi/2.
c
c--------------------------------------------------------------------------
c
      implicit none

      real*8   theta, phi
c
      exsoln = sin(theta)**2*cos(phi)
c
      return
      end
