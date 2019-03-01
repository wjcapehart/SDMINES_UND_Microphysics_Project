      subroutine putf(so,qf,l,m,deltax,deltay,k)

c
c-----------------------------------------------------------------------
c
c***Author  Bandy, Victor A.
c             Computational Mathematics Group
c             University of Colorado at Denver
c
c***DESCRIPTION
c   PUTF is a user supplied subroutine that defines the difference operator
c   coefficients and the right hand side on grid k in two dimensions. 
c
c   This subroutine discretizes the PDE given below using central finite
c   difference approximations.
c
c      d [ P du ]   d [ P du ]            2    2
c    - --    --   - --    --   + u = -7 (x  + y  ) - 4
c      dx    dx     dy    dy
c
c     where  P = x**2 + y**2 + 1
c     and divergence flux boundary conditions on all four boundaries,
c
c       normal dot [ P grad u ] = g
c
c     and the exact solution of     u = x**2 + y**2
c     with domain X=[0,1], Y=[1,3]
c
c
c--INPUT
c
c  l      number of unkowns in the x-direction, integer variable.
c
c  m      number of unknowns in the y-direction, integer variable.
c
c  deltax the grid spacing in the x-direction, real variable.
c
c  deltax the grid spacing in the y-direction, real variable.
c
c  k      the grid level that is to be discretized, generating the
c         coefficient matrix and right hand side for that level, integer.
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
      implicit  none

      integer   ko, kw, ks
      parameter (ko=1, kw=2, ks=3)
c
      integer   l, m, k, i, j
      real*8    so(0:l+1,0:m+1,3), qf(0:l+1,0:m+1), deltax, deltay,
     &          P5, x(160), y(320), dx2, dy2
c
c
c-- Compute all south, west, and center stencil coefficients, and
c   compute the right hand side for each equation.
c   Note: the south coefficients for the points (i,1) i=1,l are set to
c        zero, and the coefficients for the points (1,j) j=1,m are set to
c        zero because they refer to ficticious points.
c
      P5 =  0.5D0

      do 10 i = 1, l
        x(i) = deltax*(i-P5)
 10   continue

      do 20 j = 1, m
        y(j) = 1 + deltay*(j-P5)
 20   continue

      dx2 = P5*deltax
      dy2 = P5*deltay

      do 30 j = 2, m
        do 30 i = 1, l
          so(i,j,ks) = x(i)**2 + (y(j) - dy2)**2 + 1
 30   continue

      do 40 j = 1, m
        do 40 i = 2, l
          so(i,j,kw) = (x(i) - dx2)**2 + y(j)**2 + 1
 40   continue

      do 50 j = 1, m
        do 50 i = 1, l
          so(i,j,ko) = (x(i) - dx2)**2 + 2*x(i)**2 + (x(i) + dx2)**2
     &               + (y(j) - dy2)**2 + 2*y(j)**2 + (y(j) + dy2)**2
     &               + 4 + deltax*deltay
          qf(i,j) = -deltax*deltay*(7*x(i)**2 + 7*y(j)**2 + 4)
 50   continue
c
c--- adjust the center coefficient and the right hand side for the
c--- boundary conditions.
c
      do 60 j = 1, m
        so(1,j,ko) = so(1,j,ko) - y(j)**2 - 1
        qf(1,j) = qf(1,j) - 0
c
        so(l,j,ko) = so(l,j,ko) - y(j)**2 - 2
        qf(l,j) = qf(l,j) + 2*deltax*(y(j)**2 + 2)
 60   continue

      do 70 i = 1, l
        so(i,1,ko) = so(i,1,ko) - x(i)**2 - 2
        qf(i,1) = qf(i,1) - 2*deltay*(x(i)**2 + 2)
c
        so(i,m,ko) = so(i,m,ko) - x(i)**2 - 10
        qf(i,m) = qf(i,m) + 6*deltay*(x(i)**2 + 10)
 70   continue
c
      return
      end
