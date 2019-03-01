      SUBROUTINE BMG2_SymStd_SETUP_cg_LU(
     &                SO, II, JJ, NStncl, ABD, NABD1, NABD2,
     &                WS, NMSGr, NProcI, NProcJ, NProc, MyProc,
     &                ProcGrid, ProcCoord, LocArrSize, MPICOMM,
     &                BMG_IOFLAG, BMG_iPARMS
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SETUP_cg_LU sets up the matrix on the coarsest grid,
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
C     SO        Refer to BMG2_SymStd_SOLVE_boxmg
C     II        Number of grid points in x direction, including
C               two fictitious points.
C     JJ        Number of grid points in y direction, including
C               two fictitious points.
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
C     ABD       Refer to BMG2_SymStd_SOLVE_boxmg
C
C =======================================================================
C  --------------------
C   LOCAL:
C  --------------------
C
C
C ==========================================================================
#if RKIND == 4
#define FLOAT_PRECISION MPI_REAL
#else
#define FLOAT_PRECISION MPI_DOUBLE_PRECISION
#endif

      IMPLICIT NONE

C -----------------------------
C     Includes
C
      INCLUDE 'mpif.h'

      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'
#include      "BMG_constants.h"

C ----------------------------
C     Argument Declarations
C 
      INTEGER  II, JJ, NABD1, NABD2, NMSGr, 
     &         NProc, NProcI, NProcJ, NStncl
      INTEGER  ProcGrid(NProcI, NProcJ), ProcCoord(2,NProc),
     &         MPICOMM, MyProc, LocArrSize(3,*)
      REAL*RKIND   ABD(NABD1,NABD2), SO(II+1,JJ+1,NStncl), WS(NMSGr)
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C ----------------------------
C     Local Declarations
C
      INTEGER  I, INFO, I1, I2, J, J1, KK, N, INT_TEMP, IERR, K, 
     &         IIG, JJG, proc, p_WS
      INTEGER  LARGESTNODES, KKMAX
      INTEGER  INT_TEMP1, INT_TEMP2, P1, P2, P1SUM, P2SUM, PXP1, PYP2

C ==========================================================================

      I1=II-1
      J1=JJ-1

      I2=I1-1
      N=I2*(J1-1)

C ------------------------------------------------
C     Calculate the global number of points
C ------------------------------------------------

      !
      ! Global number in x
      !
      IIG=2
      DO i=1,NProcI
         IIG = IIG + LocArrSize(1,ProcGrid(i,1)) - 2
      ENDDO

      !
      ! Global number in y
      !
      JJG=2
      DO j=1,NProcJ
         JJG = JJG + LocArrSize(2,ProcGrid(1,j)) - 2
      ENDDO

C ------------------------------------------------
C     Find the largest local array
C ------------------------------------------------

      !
      ! Loop over local x dimensions
      !
      INT_TEMP = LocArrSize(1,1) - 2
      INT_TEMP1 = INT_TEMP
      DO proc=2, NProcI*NProcJ
         INT_TEMP = LocArrSize(1,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP1) THEN
            INT_TEMP1 = INT_TEMP
         END IF
      END DO

      !
      ! Loop over local y dimensions
      !
      INT_TEMP = LocArrSize(2,1) - 2
      INT_TEMP2 = INT_TEMP
      DO proc=2, NProcI*NProcJ
         INT_TEMP = LocArrSize(2,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP2) THEN
            INT_TEMP2 = INT_TEMP
         END IF
      END DO

      !
      ! Conservative: take largest from each
      !
      LARGESTNODES = INT_TEMP1 * INT_TEMP2 * 5 + 3

C ------------------------------------------------
C     Copy all information into the buffer
C ------------------------------------------------
      
      WS(1) = MyProc
      WS(2) = II
      WS(3) = JJ
      INT_TEMP = 4

      IF ( NStncl.EQ.5 ) THEN
         !
         DO J = 2, JJ - 1
            DO I = 2, II - 1
               WS(INT_TEMP)   =  SO(I,J,ko)
               WS(INT_TEMP+1) = -SO(I,J,kw)
               WS(INT_TEMP+2) = -SO(I+1,J,knw)
               WS(INT_TEMP+3) = -SO(I,J,ks)
               WS(INT_TEMP+4) = -SO(I,J,ksw)
               INT_TEMP = INT_TEMP + 5
            END DO
         END DO
         !
      ELSE IF ( NStncl.EQ.3 ) THEN
         !
         DO J = 2, JJ - 1
            DO I = 2, II - 1
               WS(INT_TEMP)   =  SO(I,J,ko)
               WS(INT_TEMP+1) = -SO(I,J,kw)
               WS(INT_TEMP+2) =  0.0D0
               WS(INT_TEMP+3) = -SO(I,J,ks)
               WS(INT_TEMP+4) =  0.0D0
               INT_TEMP = INT_TEMP + 5
            END DO
         END DO
         !
      ELSE
         IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
            WRITE(*,*) 
            WRITE(*,*) 'NEED: NStncl = 3 or 5 '
            WRITE(*,*) 'HAVE: NStncl = ', NStncl
         END IF
         
         CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,18)
         RETURN

      ENDIF
      
C ------------------------------------------------
C     Send/Receive information to/from everybody
C ------------------------------------------------

      IF (BMG_iPARMS(id_BMG2_CG_COMM) .eq. BMG_CG_ALLGATHER) THEN

         CALL MPI_ALLGATHER ( 
     &        WS(1), LARGESTNODES, FLOAT_PRECISION,
     &        WS(LARGESTNODES+1), LARGESTNODES, FLOAT_PRECISION,
     &        MPICOMM, IERR
     &        ) 
         
      ELSE

         CALL MPI_GATHER ( 
     &        WS(1), LARGESTNODES, FLOAT_PRECISION,
     &        WS(LARGESTNODES+1), LARGESTNODES, FLOAT_PRECISION,
     &        iZERO, MPICOMM, IERR
     &        ) 

      END IF

C ------------------------------------------------
C     Assemble the global coarse grid matrix
C ------------------------------------------------


      IF ( BMG_iPARMS(id_BMG2_CG_COMM) .eq. BMG_CG_ALLGATHER .OR.
     &     MyProc .eq. iONE ) THEN


         KKMAX=0
         
         DO proc=1, NProcI*NProcJ
            
            INT_TEMP = proc * LARGESTNODES + 1
            INT_TEMP = WS(INT_TEMP)
            
            P1 = ProcCoord(1,INT_TEMP)
            P2 = ProcCoord(2,INT_TEMP)
            
            P1SUM = 0
            DO i=1, P1-1
               P1SUM = P1SUM + LocArrSize(1,ProcGrid(i,1)) - 2
            END DO
            
            P2SUM = 0
            DO j=1, P2-1
               P2SUM = P2SUM + LocArrSize(2,ProcGrid(1,j)) - 2
            END DO
            
            PXP1 = LocArrSize(1,INT_TEMP) - 2
            PYP2 = LocArrSize(2,INT_TEMP) - 2
            
            DO j=1, PYP2
               DO i=1, PXP1
                  
                  KK = P1SUM + (i-1) + (IIG-2) * (P2SUM + j - 1) + 1
                  p_WS =  proc*LARGESTNODES+5*( (i-1)+(j-1)*PXP1) + 4
                  
                  ABD(IIG,KK)   = WS( p_WS )
                  ABD(IIG-1,KK) = WS( p_WS + 1 )
                  
                  ABD(3,KK) =  WS( p_WS + 2 )
                  ABD(2,KK) =  WS( p_WS + 3 )
                  ABD(1,KK) =  WS( p_WS + 4 )
                  
                  IF (KK.gt.KKMAX) KKMAX = KK
                  
               END DO
            END DO
            
         END DO



C ------------------------------------------------
C     DPBTRF is a LAPACK routine
C ------------------------------------------------
#if RKIND == 4
         CALL SPBTRF('U', KKMAX, IIG-1, ABD, NABD1, INFO) 
#else
         CALL DPBTRF('U', KKMAX, IIG-1, ABD, NABD1, INFO) 
#endif
         IF (INFO .NE. 0) THEN

            IF (BMG_IOFLAG(iBMG2_OUT_STOP_ERROR) .AND. MyProc.EQ.1) THEN
               WRITE(*,*) 'Coarse grid Cholesky decomposition failed!'
               WRITE(*,*) 'INFO = ', INFO
            END IF
            
            BMG_iPARMS(id_BMG2_Ext_Err_Code) = INFO
            CALL BMG2_SymStd_ErrTrap(BMG_iPARMS,17)
            RETURN   

         ENDIF
         
      END IF

C ==========================================================================


C =============================

      RETURN
      END
