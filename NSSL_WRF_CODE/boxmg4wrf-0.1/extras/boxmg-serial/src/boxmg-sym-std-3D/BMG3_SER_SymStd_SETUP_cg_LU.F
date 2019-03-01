      SUBROUTINE BMG3_SER_SymStd_SETUP_cg_LU( 
     &                SO, ii, jj, kk, NStncl, abd, nabd1, nabd2,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SER_SymStd_SETUP_cg_LU sets up the matrix on the coarsest grid,
C     and using the LAPACK routine DPBTRF, it forms the LU decomposition
C     of the matrix.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     SO        Refer to BMG3_SER_SymStd_SOLVE_boxmg
C
C     II        Number of grid points in x direction, including
C               two fictitious points.
C     JJ        Number of grid points in y direction, including
C               two fictitious points.
C     KK        Number of grid points in z direction, including
C               two fictitious points.
C
C     NABD1     Leading dimension of ABD.
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C     ABD       Refer to BMG3_SER_SymStd_SOLVE_boxmg
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================

      IMPLICIT NONE

C ----------------------------
C     Includes
C 
#include       "BMG_SER_constants.h"
      INCLUDE 'BMG_SER_stencils.h'
      INCLUDE 'BMG_SER_parameters.h'

C ----------------------------
C     Argument Declarations
C 
      INTEGER ii, info, jj, kk, nabd1, nabd2, NStncl
      REAL*RKIND  abd(nabd1,nabd2), so(ii,jj,kk,NStncl)
      LOGICAL BMG_IOFLAG(NBMG_SER_IOFLAG)
      INTEGER BMG_iPARMS(NBMG_SER_iPARMS)

C ----------------------------
C     Local Declarations
C
      INTEGER i, i1, i2, ibw, j, j1, k, k1, kl, n

C =========================================================================

C -------------------------------------------------------
C     Copy the operator on the coarsest grid into ABD 
C -------------------------------------------------------

      info = 0

      i1=ii-1
      j1=jj-1
      k1=kk-1

      i2=i1-1
C      n=i2*(j1-1)*(kk-2)
      ibw=i2*j1+1


      IF ( NStncl.EQ.14 ) THEN 

         kl=0
         do k=2,k1
            do j=2,j1
               do i=2,i1
                  !
                  kl=kl+1
                  !
                  abd(ibw+1,kl)    = so(i,j,k,kp)
                  abd(ibw,kl)      = -so(i,j,k,kpw)
                  abd(ibw-i1+3,kl) = -so(i+1,j,k,kpnw)
                  abd(ibw-i1+2,kl) = -so(i,j,k,kps)
                  abd(ibw-i1+1,kl) = -so(i,j,k,kpsw)
                  !
                  abd(ibw-(j1-2)*i2+2,kl) = -so(i+1,j+1,k,kbne)
                  abd(ibw-(j1-2)*i2+1,kl) = -so(i,j+1,k,kbn)
                  abd(ibw-(j1-2)*i2,kl)   = -so(i,j+1,k,kbnw)
                  abd(ibw-(j1-1)*i2+2,kl) = -so(i+1,j,k,kbe)
                  abd(ibw-(j1-1)*i2+1,kl) = -so(i,j,k,kb)
                  abd(ibw-(j1-1)*i2,kl)   = -so(i,j,k,kbw)
                  !
                  abd(3,kl) = -so(i+1,j,k,kbse)
                  abd(2,kl) = -so(i,j,k,kbs)
                  abd(1,kl) = -so(i,j,k,kbsw)
                  !
               enddo
            enddo
         enddo

      ELSE IF ( NStncl.EQ.4 ) THEN 
         
         kl=0
         do k=2,k1
            do j=2,j1
               do i=2,i1
                  !
                  kl=kl+1
                  !
                  abd(ibw+1,kl)    = so(i,j,k,kp)
                  abd(ibw,kl)      = -so(i,j,k,kpw)
                  abd(ibw-i1+3,kl) = rZERO
                  abd(ibw-i1+2,kl) = -so(i,j,k,kps)
                  abd(ibw-i1+1,kl) = rZERO
                  !
                  abd(ibw-(j1-2)*i2+2,kl) = rZERO
                  abd(ibw-(j1-2)*i2+1,kl) = rZERO
                  abd(ibw-(j1-2)*i2,kl)   = rZERO
                  abd(ibw-(j1-1)*i2+2,kl) = rZERO
                  abd(ibw-(j1-1)*i2+1,kl) = -so(i,j,k,kb)
                  abd(ibw-(j1-1)*i2,kl)   = rZERO
                  !
                  abd(3,kl) = rZERO
                  abd(2,kl) = rZERO
                  abd(1,kl) = rZERO
                  !
               enddo
            enddo
         enddo

      ELSE
         IF (BMG_IOFLAG(iBMG3_SER_OUT_STOP_ERROR)) THEN
            WRITE(*,500)
            WRITE(*,510) 'NEED: NStncl = 4 or 14 '
            WRITE(*,510) 'HAVE: NStncl = ', NStncl
         END IF
         
         CALL BMG3_SER_SymStd_ErrTrap(BMG_iPARMS,4)
         RETURN
         
      ENDIF


C -------------------------------------------------------
C     Factor using the LAPACK routine DPBTRF
C -------------------------------------------------------

#if RKIND == 4
      CALL SPBTRF('U', Kl, IBW, ABD, NABD1, INFO) 
#else
      CALL DPBTRF('U', Kl, IBW, ABD, NABD1, INFO) 
#endif
      IF (INFO .NE. 0) THEN

         IF (BMG_IOFLAG(iBMG3_SER_OUT_STOP_ERROR)) THEN
            WRITE(*,500) 'Coarse grid Cholesky decomposition failed!'
            WRITE(*,510) 'INFO = ', INFO
         END IF

         BMG_iPARMS(id_BMG3_SER_Ext_Err_Code) = INFO
         CALL BMG3_SER_SymStd_ErrTrap(BMG_iPARMS,5)
         RETURN

      ENDIF


C ========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SER_SymStd_SETUP_cg_LU.f',/,5X,A)
 510  FORMAT (5X,A,1X,I3)

C ===========================================

      RETURN
      END
