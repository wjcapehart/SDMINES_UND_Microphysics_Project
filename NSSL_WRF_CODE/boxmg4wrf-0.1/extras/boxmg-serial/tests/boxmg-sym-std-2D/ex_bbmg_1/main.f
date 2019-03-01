      program main
c
c-----------------------------------------------------------------------
c
c***Author  Bandy, Victor A.
c             Computational Mathematics Group
c             University of Colorado at Denver
c
c   This program solves the PDE given below using the BBMG interface to
c   the BOXMG, black box multigrid, solver where the system of equations
c   is obtained by using central finite difference approximations.
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
c--variables
c
c  l,m   the number of unknowns in the x and y directions respectively.
c
c  deltax,deltay  the x and y direction grid spacings respectively.
c
c  For the other variables see the documentation header in BBMG.
c
c
c----------------------------------------------------------------------
c
      IMPLICIT NONE

      INCLUDE 'BMG_SER_parameters.h'

      integer  idimsr, idimcf, idmiwk, idmrwk
      parameter (idimsr = 383,
     &           idimcf = 1915,
     &           idmiwk = 29,
     &           idmrwk = 2197)
c
      integer  l, m, istncl, mgparm(10), iparm(4), iwork(idmiwk),
     &         ierr, i
c
      REAL*8   deltax, deltay, erparm(3), soln(idimsr), rhs(idimsr),
     &         coef(idimcf), rwork(idmrwk), mxnorm
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)
      
c
c
c--- initialize input parameters
c
      l = 10
      m = 20
      deltax = 1D0/l
      deltay = 2D0/m
      istncl = 5

      mgparm(1) = 1
      mgparm(2) = 1
      mgparm(3) = 2
      mgparm(4) = 3
      mgparm(5) = 0
      mgparm(6) = 10
      mgparm(7) = 1
      mgparm(8) = 1
      mgparm(9) = 3
      mgparm(10) = 0

      iparm(1) = 0
      iparm(2) = 2
      iparm(3) = 0
      iparm(4) = 0

      erparm(1) = 0.
      erparm(2) = 0.
      erparm(3) = 1D-6
c
      write(*,*) ' dump mgparm = ',mgparm
      write(*,*) ' dump iparm = ',iparm
      write(*,*) ' dump erparm = ',erparm

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

      BMG_IOFLAG(iBMG2_SER_OUT_STOP_ERROR)   = .TRUE.
c
c--- initialize the solution and the coefficient matrix, as required
c    by BOXMG.
c
      do 10 i = 1, idimsr
        soln(i) = 0
        rhs(i) = 0
 10   continue

      do 20 i = 1, idimcf
        coef(i) = 0
 20   continue
c
c--- call PUTF to setup the system of equations
c
      call putf (coef, rhs, l, m, deltax, deltay, 0)
c
c--- call BBMG to solve the system of equations using multigrid
c
      write(*,*) ' The following output is from BOXMG:'

      call bbmg( l, m, deltax, deltay, istncl, 
     &           mgparm, iparm, erparm, BMG_IOFLAG,
     &           soln, rhs, idimsr, coef, idimcf,
     &           iwork, idmiwk, rwork, idmrwk, ierr )
c
      write(*,*) ' The following output is from the MAIN program:'

      if ( ierr .eq. 0 ) then
        write(*,*) ' The discrete l2 norm of the residual = ',erparm(3)
        call diff (l, m, deltax, deltay, soln, mxnorm)
        write(*,100) mxnorm
 100    format (' ','The max-norm of the difference between the exact ',
     &          /,' and the computed solution is ',e15.7)
      else
        write(*,*) ' ERROR: ierr = ',ierr
        do 30 i = 1, ierr
          write(*,40) i, iwork(i)
 40       format (' ERROR: iwork(',i2,') = ',i3)
 30     continue
        write(*,*) 'Error: idimsr should be > or = ',iwork(ierr+1)
        write(*,*) 'Error: idimcf should be > or = ',iwork(ierr+2)
        write(*,*) 'Error: idmrwk should be > or = ',iwork(ierr+3)
        write(*,*) 'Error: idmiwk should be > or = ',iwork(ierr+4)
      endif
 
      end



