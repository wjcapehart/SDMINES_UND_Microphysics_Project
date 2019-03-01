      subroutine putf(coef, rhs, l, m, dphi, dtheta, k, ipsbc)
c
c-------------------------------------------------------------------------
c
c   This subroutine computes the matrix coefficients and right hand
c   side for the system of equations obtained by discretizing the PDE
c   using central finite difference approximations.
c   The domain is the unit sphere where the grid used is staggered
c   (cell-centered) and logically rectangular and is border by a layer
c   of fictitious points needed by the solver routine.
c   The equation in spherical coordinates being discretized is 
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
c--INPUT
c
c  l      number of unkowns in the phi-direction, integer variable.
c
c  m      number of unknowns in the theta-direction, integer variable.
c
c  dphi   the grid spacing in the phi-direction, real variable.
c
c  dtheta the grid spacing in the theta-direction, real variable.
c
c  k      the grid level that is to be discretized, generating the
c         coefficient matrix and right hand side for that level, integer.
c
c  ipsbc  integer variable that is set to -2 to indicate that the problem
c         is periodic in the phi-direction, and also that the system of
c         equations is singular.
c
c--OUTPUT
c
c  coef   real array that contains the coefficients of the finite 
c         difference matrix. Dimensioned as coef(0:l+1,0:m+1,3)
c         where the 3 stands for the three non-zero entries of the
c         symmetric matrix being stored corresponding to the finite
c         difference 5-point stencil: center, west, and south points.
c         The grid has been bordered by a single layer of fictitious
c         points that are referred to by 0, l+1, and m+1.
c
c  rhs    real array that contains the right hand side for the discretized
c         system of equations.
c
c-------------------------------------------------------------------------
c
      implicit none

      integer   ko, kw, ks
      parameter (ko=1, kw=2, ks=3)
c
      integer   l, m, i, ipsbc, j, k
      real*8    coef(0:l+1,0:m+1,3), dpbydt, dphi, dtbydp, dtheta,
     &          exsoln, P5, rhs(0:l+1,0:m+1)
c
c
c-- Initialize right hand side and difference operator on grid k.
c
       P5=0.5D0
       dpbydt = dphi/dtheta
       dtbydp = dtheta/dphi
c
c---Compute the south stencil coefficients
c
      do 10 j = 1, m
        do 10 i = 1, l
          coef(i,j,ks) = (sin((j-1)*dtheta))*dpbydt
   10 continue
c
c---Compute the west stencil coefficients
c
      do 20 j = 1, m
        do 20 i = 1, l
          coef(i,j,kw) =(1/sin((j-P5)*dtheta))*dtbydp
   20 continue
c
c---Compute the central stencil coefficients
c
      do 30 j = 1, m
        do 30 i = 1, l
          coef(i,j,ko) = (sin((j-1)*dtheta) + sin(j*dtheta))*dpbydt
     *                 + (2/sin((j-P5)*dtheta))*dtbydp

          rhs(i,j) = dphi*(sin((j-1)*dtheta)
     *             *(exsoln((j-P5)*dtheta,(i-P5)*dphi)
     *             -exsoln((j-1-P5)*dtheta,(i-P5)*dphi))
     *             - sin(j*dtheta)*(exsoln((j+P5)*dtheta,(i-P5)*dphi)
     *             -exsoln((j-P5)*dtheta,(i-P5)*dphi)))/dtheta
     *             - dtheta*(exsoln((j-P5)*dtheta,(i-1-P5)*dphi)
     *             -2.*exsoln((j-P5)*dtheta,(i-P5)*dphi)
     *             +exsoln((j-P5)*dtheta,(i+P5)*dphi))
     *             /(sin((j-P5)*dtheta)*dphi)
   30 continue
c
c---copy periodic values from east side of grid to the west side.
c
      do 40 j=1,m
        coef(0,j,ko)=coef(l,j,ko)
        coef(0,j,ks)=coef(l,j,ks)
        coef(0,j,kw)=coef(l,j,kw)
        rhs(0,j)=rhs(l,j)
        
        coef(l+1,j,ko)=coef(1,j,ko)
        coef(l+1,j,ks)=coef(1,j,ks)
        coef(l+1,j,kw)=coef(1,j,kw)
        rhs(l+1,j)=rhs(1,j)
   40 continue
c
c---Set the fictitious points to match the boundary conditions
c
      do 50 i=0,l+1
        coef(i,0,ks)=0
        coef(i,m+1,ks)=0
        coef(i,0,ko)=1
        rhs(i,0)=0
        rhs(i,m+1)=0
        coef(i,m+1,ko)=1
   50 continue
c
c--  normalize the right hand side to make sure the problem is consistent.
c
        write (*,*) ' normalizing right side'
        call normal (rhs, dphi, dtheta, l, m)
c
      return
      end
