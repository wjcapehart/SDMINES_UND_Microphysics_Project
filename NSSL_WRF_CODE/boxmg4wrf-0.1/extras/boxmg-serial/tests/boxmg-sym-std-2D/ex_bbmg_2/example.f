	program EXAMPLE
c
c       Sample main program for illustrating the use of BBMG.
c       The Partial differential equation solved for in this example is :
c
c       -D/DX [(X+Y)*DU/DX] - D/DY [(X-Y)*DU/DY] + X*Y*U = X-Y+(X**2)*(Y**2).
c
c       This equation has an exact solution u = x*y. 
c
c       The rectangular domain used for this problem is :
c       ----------------------------------------------
c
c       XS = 3, XF = 4, YS = 1, and YF = 2.  
c
c       The following boundary conditions were used :
c       --------------------------------------------
c
c                 U(3,Y) = 3*Y    	            (Derichlet)
c	          DU/DX(4,Y) = Y 		    (Neumann)
c	          DU/DY(X,1) - X*U(x,1) = X - X**2  (Robin)
c                 U(X,2) = 2*X  		    (Derichlet)
c
c       The number of unknowns in the x-direction is : NXPTS = 12
c       The number of unknowns in the y-direction is : NYPTS = 10
c
c       Computation of IDIMSR, IDIMCF, IDMIWK, IDMRWK :
c       -----------------------------------------------
c       First compute NGRIDS :
c	NGRIDS .ge. 1 + LOG((MIN(NXPTS,NYPTS)-1/(NXYO-1))/LOG(2).
c       So NGRIDS .ge. 3.  Set NGRIDS = 3.
c        
c       IDIMSR .le. 9*NGRIDS + (4./3.) * (NXPTS+4.5) * (NYPTS+4.5) - 40.
c       So IDIMSR .le. 306.  Set IDIMSR = 306.
c
c       IDIMCF = 3*IDIMSR + 2*(IDIMSR-(NXPTS+2)*(NYPTS+2)).
c       So IDIMSR = 1194.
c
c       IDMIWK = MAX(9*NGRIDS,29) = 29
c
c       IDMRWK = 2*IDIMSR + 8*NCMAX + NXCG*NYCG*(NXCG+3) + 2*NGRIDS,
c           where NXCG = max(1+(NXPTS-1)/(2**(NGRIDS-1)), NXYO), 
c                 NYCG = max(1+(NYPTS-1)/(2**(NGRIDS-1)), NXYO), 
c                 NCMAX = NFMAX - (NXPTS+2)*(NYPTS+2).
c       So IDMRWK = 1776
c
	
      IMPLICIT NONE

      INCLUDE 'BMG_SER_parameters.h'

      INTEGER IDIMSR, IDIMCF, IDMIWK, IDMRWK      
      PARAMETER ( IDIMSR = 249,
     &            IDIMCF = 1194,
     &            IDMIWK = 29,
     &            IDMRWK = 1455)

      INTEGER  I, J, NXPTS, NYPTS, ISTNCL, MGPARM(10), IPARM(4),
     &         IWORK(IDMIWK), IERR, DEBUG
      REAL*8   DELTAX, DELTAY, DIFF, ERPARM(3), ERR, ES,
     &         SOLN(IDIMSR), RHS(IDIMSR), SO(IDIMCF), RWORK(IDMRWK),
     &         TEMP(14,12), X, Y
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)

      REAL*8   ESOL
      EXTERNAL ESOL, PUTF

      REAL*8   XF, XS, YF, YS
      COMMON   /BBMGDM/ XS,XF,YS,YF  /DBUG/ DEBUG

      DEBUG=1

c     domain parameters for x and y variables
	
      XS=3
      XF=4
      YS=1
      YF=2

c     Finite difference grid parameters 

      NXPTS=12
      NYPTS=10
      DELTAX=(XF-XS)/NXPTS
      DELTAY=(YF-YS)/NYPTS
      ISTNCL=5

c     Specify multigrid parameters 

      MGPARM(1)=1
      MGPARM(2)=1
      MGPARM(3)=2
      MGPARM(4)=3
      MGPARM(5)=0
      MGPARM(6)=20
      MGPARM(7)=1
      MGPARM(8)=1
      MGPARM(9)=3
      MGPARM(10)=0

c     Specify truncation error, printing, and initialization options     

      IPARM(1)=0
      IPARM(2)=0
      IPARM(3)=0
      IPARM(4)=0

c     Specify error/stopping criteria for the multigrid cycling

      ERPARM(1)=0.
      ERPARM(2)=0.
      ERPARM(3)=1D-3

C     Specify IO flags.

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

      BMG_IOFLAG(iBMG_SER_OUT_SOLUTION)     = .TRUE.

      BMG_IOFLAG(iBMG2_SER_OUT_STOP_ERROR)  = .TRUE.
C     Must zero everything (very important on any machine except SUNs)

      DO 4 i = 1, IDIMSR
	 SOLN(i)=0
	 RHS(i) = 0
 4    CONTINUE

      DO 5 i=1, IDIMCF
	 SO(i)=0
 5    CONTINUE

      CALL PUTF(SO,RHS,NXPTS,NYPTS,DELTAX,DELTAY,0)

      CALL BBMG( NXPTS, NYPTS, DELTAX, DELTAY, ISTNCL,
     &           MGPARM, IPARM, ERPARM, BMG_IOFLAG,
     &           SOLN, RHS, IDIMSR, SO, IDIMCF,
     &           IWORK, IDMIWK, RWORK, IDMRWK, IERR )
      
c     check error parameters

      IF(IERR.NE.0) THEN
	 WRITE(*,1) IERR, (IWORK(I),I=1,IERR+4)
 1	 FORMAT(' IERR FROM BBMG :',I2,/,'ERROR PARAMTERS =',
     1   5I5)         
	 STOP
      ENDIF

c     Compute error in solution    
c     First copy SOLN array into TEMP matrix 

      DO 20 J=1,NYPTS+2
      DO 20 I=1,NXPTS+2
        TEMP(I,J) = SOLN((J-1)*(NXPTS+2) + I)
 20   CONTINUE
      WRITE(*,*) ' i   j     x       y      exact soln           ',
     &           'soln         difference'
      WRITE(*,*) '---------------------------------------------------',
     &           '--------------------'

      ERR = 0.
      DO 21 J=2,NYPTS+1
      DO 21 I=2,NXPTS+1
        X = XS + (I-1)*DELTAX
        Y = YS + (J-2)*DELTAY
        ES = ESOL(X,Y)
        DIFF = ABS(TEMP(I,J) - ES)
        ERR=MAX(DIFF,ERR)
        WRITE(*,100) I-1, J-1, X, Y, ES, TEMP(I,J), DIFF
 100    FORMAT(' ',I2,2X,I2,2X,F6.4,2X,F6.4,3X,F10.7,3X,F15.7,3X,E15.7)

 21   CONTINUE

      WRITE(*,2) NXPTS,NYPTS,IERR,ERR
 2    FORMAT('  ON THE GRID WITH NXPTS, AND NYPTS = ',2I4,/,
     1 'IERR =',I4,' AND THE ERROR WAS = ',E13.5)

c     ON THE SUN 4 THE WRITE STATEMENT ABOVE GENERATED THE FOLLOWING 
c     TWO LINES OF OUTPUT :
c
c     ON THE GRID WITH NXPTS, AND NYPTS =   12  10
c     IERR =   0 AND THE ERROR WAS =   0.23842E-05


      STOP
      END


      REAL*8 FUNCTION ESOL(X,Y)
      
      IMPLICIT NONE

      REAL*8  X, Y

      ESOL=X*Y

      RETURN
      END
