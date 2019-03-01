      SUBROUTINE BMG3_SER_SymStd_SOLVE_cg( 
     &         q, qf, ii, jj, kk, abd, bbd, nabd1, nabd2,
     &         BMG_IOFLAG, BMG_iPARMS 
     &         )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_SOLVE_cg does a direct solve on the coarsest grid. it
C     uses the LAPACK routine DPBTRS.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     QF        Refer to BMG3_SER_SymStd_SOLVE_boxmg
C     
C     II        Number of grid points in x direction, including
C               two fictitious points.
C     JJ        Number of grid points in y direction, including
C               two fictitious points.
C     KK        Number of grid points in z direction, including
C               two fictitious points.
C
C     ABD       Refer to BMG3_SER_SymStd_SOLVE_boxmg
C
C     NABD1     Refer to BMG3_SER_SymStd_SOLVE_boxmg
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     BBD       Workspace, refer to BMG3_SER_SymStd_SOLVE_boxmg
C
C  --------------------
C   OUTPUT:
C  --------------------
C
C     Q         Refer to BMG3_SER_SymStd_SOLVE_boxmg
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE
      
C -----------------------------
C     Includes
C
      INCLUDE 'BMG_SER_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER ii, jj, kk, nabd1, nabd2
      REAL*RKIND abd(nabd1,nabd2), bbd(nabd2), q(ii,jj,kk), qf(ii,jj,kk)

      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)
      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS)

C ----------------------------
C     Local Declarations
C
      integer i, i1, i2, ibw, j, j1, k, kt, k1, n, info

C =========================================================================

c     
c     direct solve on coarsest grid
c     

      i1=ii-1
      j1=jj-1
      k1=kk-1

      i2=i1-1
      n=i2*(j1-1)*(kk-2)
      ibw=i2*j1+1

      kt=0
      do k=2,k1
         do j=2,j1
            do i=2,i1
               kt=kt+1
               bbd(kt)=qf(i,j,k)
            enddo
         enddo
      enddo

C -------------------------------------------------------
C     Solve using the LAPACK routine DPBTRS
C -------------------------------------------------------

#if RKIND == 4
      CALL SPBTRS ('U', KT, IBW, 1, ABD, NABD1, BBD, NABD2, INFO) 
#else
      CALL DPBTRS ('U', KT, IBW, 1, ABD, NABD1, BBD, NABD2, INFO) 
#endif
      IF (INFO .NE. 0) THEN

         IF (BMG_IOFLAG(iBMG3_SER_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'Coarse grid solve failed!'
            WRITE(*,510) 'INFO = ', INFO
         END IF

         BMG_iPARMS(id_BMG3_SER_Ext_Err_Code) = INFO
         CALL BMG3_SER_SymStd_ErrTrap(BMG_iPARMS,20)
         RETURN

      ENDIF


      kt=0
      do k=2,k1
         do j=2,j1
            do i=2,i1
               kt=kt+1
               q(i,j,k)=bbd(kt)
            enddo
         enddo
      enddo

C =========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SER_SymStd_SOLVE_cg.f',/,5X,A)
 510  FORMAT (5X,A,1X,I3)

C ===========================================

      RETURN
      END
