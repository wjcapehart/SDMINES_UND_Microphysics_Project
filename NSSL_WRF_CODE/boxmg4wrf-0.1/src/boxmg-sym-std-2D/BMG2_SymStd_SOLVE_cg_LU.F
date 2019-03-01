      SUBROUTINE BMG2_SymStd_SOLVE_cg_LU(
     &                Q, QF, II, JJ, 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                ABD, BBD, NABD1, NABD2, NOG, 
     &                NProcI, NProcJ, NProc, MyProc,
     &                ProcGrid, ProcCoord, LocArrSize,
     &                iWork, NMSGi, pMSG, WS, NMSGr, MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SOLVE_cg does a direct solve on the coarsest grid. it
C     uses the LAPACK routine DPBTRS.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     QF        Refer to BMG2_SymStd_SOLVE_boxmg
C     
C     II        Number of grid points in x direction, including
C               two fictitious points.
C     JJ        Number of grid points in y direction, including
C               two fictitious points.
C
C     ABD       Refer to BMG2_SymStd_SOLVE_boxmg
C
C     NABD1     Refer to BMG2_SymStd_SOLVE_boxmg
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     BBD       Workspace, refer to BMG2_SymStd_SOLVE_boxmg
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C     Q         Refer to BMG2_SymStd_SOLVE_boxmg
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

C ----------------------------
C     Argument Declarations
C 
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'
#include      "BMG_constants.h"

      INCLUDE 'mpif.h'


      !
      !  BoxMG parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  II, JJ, NABD1, NABD2, NMSGi, NMSGr, NOG
      INTEGER  MPICOMM, MyProc, NProc, NProcI, NProcJ
      INTEGER  ProcGrid(NProcI, NProcJ), ProcCoord(2,NProc)
      INTEGER  LocArrSize(3,*), iWork(NMSGi), pMSG(NBMG_pMSG,NOG)
      REAL*RKIND   ABD(NABD1,NABD2), BBD(NABD2), Q(II,JJ), QF(II,JJ),
     &         WS(NMSGr)

C ----------------------------
C     Local Declarations
C
      INTEGER  I, I1, I2, J, J1, KK, N, INFO, INT_TEMP, IERR,
     &         IIG, JJG, proc, p_WS
      INTEGER  LARGESTNODES, KKMAX
      INTEGER  INT_TEMP1, INT_TEMP2, P1, P2, P1SUM, P2SUM, IIL, JJL
      INTEGER  TMP_RANK
      
      INTEGER  iBEG, iEND, jBEG, jEND

C =========================================================================

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
      do i=1,NProcI
         IIG = IIG + LocArrSize(1,ProcGrid(i,1)) - 2
      enddo

      !
      ! Global number in y
      !
      JJG=2
      do j=1,NProcJ
         JJG = JJG + LocArrSize(2,ProcGrid(1,j)) - 2
      enddo

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
      LARGESTNODES = INT_TEMP1 * INT_TEMP2  + 3

C ------------------------------------------------
C     Copy all information into the buffer
C ------------------------------------------------

      WS(1) = MyProc
      WS(2) = II
      WS(3) = JJ
      INT_TEMP = 4

      DO J = 2, JJ - 1
        DO I = 2, II - 1
           WS(INT_TEMP) = QF(I,J)
           INT_TEMP = INT_TEMP + 1
        END DO
      END DO

C ------------------------------------------------
C     Send/Receive information to/from everybody
C ------------------------------------------------

      IF (BMG_iPARMS(id_BMG2_CG_COMM) .eq. BMG_CG_ALLGATHER) THEN

         CALL MPI_ALLGATHER( 
     &        WS(1), LARGESTNODES, FLOAT_PRECISION,
     &        WS(LARGESTNODES+1), LARGESTNODES, FLOAT_PRECISION,
     &        MPICOMM, IERR
     &        ) 
         
      ELSE
         
         CALL MPI_GATHER( 
     &        WS(1), LARGESTNODES, FLOAT_PRECISION,
     &        WS(LARGESTNODES+1), LARGESTNODES, FLOAT_PRECISION,
     &        iZERO, MPICOMM, IERR
     &        )
         
      END IF

C ------------------------------------------------
C     Assemble the global right hand side
C ------------------------------------------------

      IF ( BMG_iPARMS(id_BMG2_CG_COMM) .eq. BMG_CG_ALLGATHER .OR. 
     &     MyProc .eq. iONE ) THEN

         KKMAX=0
         
         DO proc = 1, NProcI * NProcJ
            
            TMP_RANK = WS(proc*LARGESTNODES+1)
            
            P1 = ProcCoord(1,TMP_RANK)
            P2 = ProcCoord(2,TMP_RANK)
            
            P1SUM = 0
            P2SUM = 0
            
            DO i=1, P1-1
               P1SUM = P1SUM + LocArrSize(1,ProcGrid(i,1)) - 2
            END DO
            
            DO j=1, P2-1
               P2SUM = P2SUM + LocArrSize(2,ProcGrid(1,j)) - 2
            END DO
            
            IIL = LocArrSize(1,TMP_RANK) - 2
            JJL = LocArrSize(2,TMP_RANK) - 2
            
            DO j = 1, JJL
               DO i = 1, IIL
                  
                  KK = P1SUM + (i-1) + (IIG-2) * (P2SUM + j - 1) + 1
                  p_WS =  proc*LARGESTNODES + (i-1) + (j-1) * IIL + 4 
                  
                  BBD(KK) =  WS(p_WS)
                  
                  IF (KK.gt.KKMAX) KKMAX = KK
                  
               END DO
            END DO
            
         END DO
         
C ------------------------------------------------
C     DPBTRS is a LAPACK routine
C ------------------------------------------------
#if RKIND
         CALL SPBTRS ('U', KKMAX, IIG-1, 1, ABD, NABD1, BBD, KKMAX,INFO) 
#else
         CALL DPBTRS ('U', KKMAX, IIG-1, 1, ABD, NABD1, BBD, KKMAX,INFO) 
#endif
         IF (INFO .NE. 0) THEN
            write(*,*) 'Coarse grid solve failed!'
            STOP
         ENDIF

      END IF


      IF (BMG_iPARMS(id_BMG2_CG_COMM) .eq. BMG_CG_GATHER_SCATTER) THEN

         !
         ! Broadcast the solution vector
         !

         CALL MPI_Bcast(BBD(1), NABD2, FLOAT_PRECISION, 
     &        iZERO, MPICOMM, IERR)
         
      END IF

C ------------------------------------------------
C     Extract the local data
C ------------------------------------------------

      KK=0
 
      IIL = LocArrSize(1,MyProc) - 2
      JJL = LocArrSize(2,MyProc) - 2
      
      P1 = ProcCoord(1,MyProc)
      P2 = ProcCoord(2,MyProc)

      P1SUM = 0
      DO i=1, P1-1
         P1SUM = P1SUM + LocArrSize(1,ProcGrid(i,1)) - 2
      END DO

      P2SUM = 0
      DO j=1, P2-1
         P2SUM = P2SUM + LocArrSize(2,ProcGrid(1,j)) - 2
      END DO
      
      !
      !  Setup loop boundaries in x
      !      
      IF ( P1.EQ.1 ) THEN
         iBEG = 1
      ELSE
         iBEG = 0
      END IF
      
      IF (P1.EQ.NProcI) THEN
         iEND = IIL
      ELSE
         iEND = IIL+1
      END IF

      !
      !  Setup loop boundaries in y
      !
      IF ( P2.EQ.1 ) THEN
         jBEG = 1
      ELSE
         jBEG = 0
      END IF

      IF ( P2.EQ.NProcJ) THEN
         jEND = JJL
      ELSE
         jEND = JJL + 1
      ENDIF

      DO j=jBEG, jEND
         DO i=iBEG, iEND

            KK = P1SUM + (i-1) + (IIG-2) * (P2SUM + j - 1) + 1

            Q(i+1,j+1)=BBD(KK)           
            
         ENDDO
      ENDDO

C =========================================================================


C =======================================

      RETURN
      END
