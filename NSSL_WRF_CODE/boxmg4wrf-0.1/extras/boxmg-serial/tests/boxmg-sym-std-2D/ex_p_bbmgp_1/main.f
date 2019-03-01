      program main
c
c-------------------------------------------------------------------------
c***Author  Bandy, Victor A.
c             Computational Mathematics Group
c             University of Colorado at Denver
c
c   This program solves the PDE given below using the BBMGP interface to
c   the BOXMGP, black box multigrid, solver where the system of equations
c   is obtained by using central finite difference approximations.
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
c--variables
c
c  l,m   the number of unknowns in the phi and theta directions respectively.
c
c  dphi,dtheta  the phi and theta direction grid spacings respectively.
c
c  For other variables see the documentation header in BBMGP.
c
c-------------------------------------------------------------------------
c

      IMPLICIT NONE 

      INCLUDE 'BMG_SER_parameters.h'

      integer  idimsr, idimcf, idmiwk, idmrwk
      parameter ( idimsr = 218,
     *            idimcf = 802,
     *            idmrwk = 1336,    !!! Note the checks in bbmgp are broken
     *            idmiwk = 30 )
c
      integer  l, m, i,
     *         istncl, ipsbc, mgparm(10), iparm(4), iwork(idmiwk), ierr
      real*8   pi, dphi, dtheta,
     *         erparm(3), soln(idimsr), rhs(idimsr), coef(idimcf),
     *         rwork(idmrwk), tsoln(idimsr)
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)
c
c
c---initialize input arguments to define the problem and the multigrid
c   algorithm used to solve it.
c
      l = 10
      m = 10
      pi = 4.d0*datan(1.d0)
      dphi = 2.*pi/float(l)
      dtheta = pi/float(m)
c
      istncl = 5
      ipsbc = -2

      mgparm(1) = 1
      mgparm(2) = 1
      mgparm(3) = 1
      mgparm(4) = 2
      mgparm(5) = 0
      mgparm(6) = 10
      mgparm(7) = 1
      mgparm(8) = 1
      mgparm(9) = 3
      mgparm(10) = 0

      iparm(1) = 0
      iparm(2) = 1
      iparm(3) = 0
      iparm(4) = 0

      erparm(1) = 0.
      erparm(2) = 0.
      erparm(3) = 1.e-6
c
      write(*,*) 'dump mgparm = ', mgparm
      write(*,*) 'dump iparm = ', iparm
      write(*,*) 'dump erparm = ', erparm

      DO i=1, NBMG_SER_IOFLAG
         BMG_IOFLAG(i)=.FALSE.
      ENDDO
      
      BMG_IOFLAG(iBMG2_SER_OUT_WSPACE_SIZE)  = .TRUE.
      BMG_IOFLAG(iBMG2_SER_OUT_ITERATIONS)   = .TRUE.

      BMG_IOFLAG(iBMG2_SER_OUT_TIME_SETUP)   = .TRUE.
      BMG_IOFLAG(iBMG2_SER_OUT_TIME_CYCLING) = .TRUE.
      BMG_IOFLAG(iBMG2_SER_OUT_TIME_TOTAL)   = .FALSE.

      BMG_IOFLAG(iBMG2_SER_BUG_RES_RELAX)    = .TRUE.
      BMG_IOFLAG(iBMG2_SER_BUG_RES_CG_SOLVE) = .TRUE.

c
c---If core is not zeroed, then coef should be zeroed.
c
c---note: since we are doing FMG it is absolutely necessary to 
c         zero the right hand side. 

      do 10 i = 1, idimcf
         coef(i) = 0.
 10   continue

      do 20 i = 1, idimsr
        soln(i) = 0.
        rhs(i)=0.
 20   continue

c
c---define the coefficient matrix and the right hand side.
c

      call putf (coef, rhs, l, m, dphi, dtheta, 0, ipsbc)

c
c---end of initial setup.
c
c---Solve the system of equations by calling the BBMGP solver.
c
      write(*,*) ' --- calling boxmgp ---'
c
      call bbmgp( l, m, dphi, dtheta, istncl, ipsbc, 
     &            mgparm, iparm, erparm, BMG_IOFLAG,
     &            soln, rhs, idimsr, coef, idimcf,
     &            iwork, idmiwk, rwork, idmrwk, ierr )
c
      if ( ierr .ne. 0 ) then
         write(*,*) ' ***** ERROR: ierr = ',ierr
         do 30 i = 1, ierr
            write(*,100) i, iwork(i)
 100         format(' Error: iwork(',i2,') = ',i3)
 30      continue
         write(*,*) 'Error: idimsr should be > or = ',iwork(ierr+1)
         write(*,*) 'Error: idimcf should be > or = ',iwork(ierr+2)
         write(*,*) 'Error: idmrwk should be > or = ',iwork(ierr+3)
         write(*,*) 'Error: idmiwk should be > or = ',iwork(ierr+4)
c
         stop
      endif
c
      write(*,110) erparm(3)
 110  format(' l2 error of the residual = ',e15.7)
c
c---normalize the computed solution
c
      write (*,*) ' --- normalizing computed solution ---'
      call normal (soln, dphi, dtheta, l, m)
c
c---compare the computed solution to the exact solution.
c
      call checker (soln, tsoln, l, m, dphi, dtheta, ipsbc)
c
      end


      

