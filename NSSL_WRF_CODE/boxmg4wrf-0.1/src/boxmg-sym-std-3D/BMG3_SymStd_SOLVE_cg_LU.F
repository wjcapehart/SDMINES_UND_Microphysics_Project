      SUBROUTINE BMG3_SymStd_SOLVE_cg_LU(
     &                Q, QF, II, JJ, KK, ABD, BBD, NABD1, NABD2, NOGm,
     &                NProcI, NProcJ, NProcK, NProc, MyProc,
     &                ProcGrid, ProcCoord, LocArrSize,
     &                WS, NMSGr, MPICOMM,
     &                BMG_iPARMS
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_SOLVE_cg_LU does a direct solve on the coarsest grid. 
C     It uses the LAPACK routine DPBTRS.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C     QF        Refer to BMG3_SymStd_SOLVE_boxmg
C     
C     II        Number of grid points in x direction, including
C               two fictitious points.
C     JJ        Number of grid points in y direction, including
C               two fictitious points.
C
C     ABD       Refer to BMG3_SymStd_SOLVE_boxmg
C
C     NABD1     Refer to BMG3_SymStd_SOLVE_boxmg
C
C =======================================================================
C  --------------------
C   INPUT/OUTPUT:
C  --------------------
C
C     BBD       Workspace, refer to BMG3_SymStd_SOLVE_boxmg
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C     Q         Refer to BMG3_SymStd_SOLVE_boxmg
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
      
      INCLUDE 'mpif.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'
#include    "BMG_constants.h"
      

C ----------------------------
C     Argument Declarations
C 
      INTEGER ii, jj, kk, nabd1, nabd2, NOGm

      INTEGER MPICOMM, MyProc, NMSGr, NProc, NProcI, NProcJ, NProcK
      INTEGER ProcGrid(NProcI, NProcJ, NProcK), ProcCoord(3,NProc),
     &        LocArrSize(3,*)
      REAL*RKIND  abd(nabd1,nabd2), bbd(nabd2), q(ii,jj,kk), 
     &        qf(ii,jj,kk), WS(NMSGr)
      INTEGER BMG_iPARMS(NBMG_iPARMS)

C ----------------------------
C     Local Declarations
C
      integer i, ibw, j, k, info
      INTEGER IIG, JJG, KKG
      INTEGER INT_TEMP, INT_TEMP1, INT_TEMP2, INT_TEMP3
      INTEGER KKMAX, iBEG, jBEG, kBEG, iEND, jEND, kEND
      INTEGER proc, P1, P2, P3, P1SUM, P2SUM, P3SUM
      INTEGER IERR, TMP_RANK, IIL, JJL, KKL, KK_t
      INTEGER p_WS, LARGESTNODES, MPI_IERROR
      INTEGER iMyProc, ictr

C =========================================================================

C ------------------------------------------------
C     Calculate the global number of points
C ------------------------------------------------

      !
      ! Global number in x
      !
      IIG=2
      do i=1,NProcI
         IIG = IIG + LocArrSize(1,ProcGrid(i,1,1)) - 2
      enddo
      
      !
      ! Global number in y
      !
      JJG=2
      do j=1,NProcJ
         JJG = JJG + LocArrSize(2,ProcGrid(1,j,1)) - 2
      enddo
      
      !
      ! Global number in z
      !
      KKG=2
      do k=1,NProcK
         KKG = KKG + LocArrSize(3,ProcGrid(1,1,k)) - 2
      enddo
      
C ------------------------------------------------
C     Find the largest local array
C ------------------------------------------------

      !
      ! Loop over local x dimensions
      !
      INT_TEMP = LocArrSize(1,1) - 2
      INT_TEMP1 = INT_TEMP
      DO proc=2, NProcI*NProcJ*NProcK
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
      DO proc=2, NProcI*NProcJ*NProcK
         INT_TEMP = LocArrSize(2,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP2) THEN
            INT_TEMP2 = INT_TEMP
         END IF
      END DO
      
      !
      ! Loop over local y dimensions
      !
      INT_TEMP = LocArrSize(3,1) - 2
      INT_TEMP3 = INT_TEMP
      DO proc=2, NProcI*NProcJ*NProcK
         INT_TEMP = LocArrSize(3,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP3) THEN
            INT_TEMP3 = INT_TEMP
         END IF
      END DO
      
          
      !
      ! Conservative: take largest from each
      !
      LARGESTNODES = INT_TEMP1 * INT_TEMP2 * INT_TEMP3  + 4
          
C ------------------------------------------------
C     Copy all information into the buffer
C ------------------------------------------------

      WS(1) = MyProc
      WS(2) = II
      WS(3) = JJ
      WS(4) = KK
      INT_TEMP = 5
      
      DO K = 2, KK - 1
         DO J = 2, JJ - 1
            DO I = 2, II - 1
               WS(INT_TEMP) = QF(I,J,K)
               INT_TEMP = INT_TEMP + 1
            END DO
         END DO
      END DO
      
C ------------------------------------------------
C     Send/Receive information to/from everybody
C ------------------------------------------------


      IF (BMG_iPARMS(id_BMG3_CG_COMM) .eq. BMG_CG_ALLGATHER) THEN
         
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
         
      ENDIF

C ------------------------------------------------
C     Assemble the global right hand side
C ------------------------------------------------

      IF ( BMG_iPARMS(id_BMG3_CG_COMM) .eq. BMG_CG_ALLGATHER .OR. 
     &     MyProc .eq. iONE ) THEN

         KKMAX=0
         
         DO proc = 1, NProcI * NProcJ * NProcK
            
            TMP_RANK = WS(proc*LARGESTNODES+1)
            
            P1 = ProcCoord(1,TMP_RANK)
            P2 = ProcCoord(2,TMP_RANK)
            P3 = ProcCoord(3,TMP_RANK)
            
            P1SUM = 0
 6          P2SUM = 0
            P3SUM = 0
            
            DO i=1, P1-1
               P1SUM = P1SUM + LocArrSize(1,ProcGrid(i,1,1)) - 2
            END DO
            
            DO j=1, P2-1
               P2SUM = P2SUM + LocArrSize(2,ProcGrid(1,j,1)) - 2
            END DO
            
            DO k=1, P3-1
               P3SUM = P3SUM + LocArrSize(3,ProcGrid(1,1,k)) - 2
            END DO
            
            IIL = LocArrSize(1,TMP_RANK) - 2
            JJL = LocArrSize(2,TMP_RANK) - 2
            KKL = LocArrSize(3,TMP_RANK) - 2
            
            DO k = 1, KKL
               DO j = 1, JJL
                  DO i = 1, IIL
                     
                     KK_t = (P1SUM + i-1) + (IIG-2) * (P2SUM + j - 1) 
     &                    + (IIG-2)*(JJG-2)*(P3SUM + k - 1) + 1
                     
                     p_WS =  proc*LARGESTNODES + 5 
     &                    + ((i-1) + (j-1)*IIL + (k-1)*IIL*JJL)
                     
                     BBD(KK_t) =  WS(p_WS)
                     
                     IF (KK_t.gt.KKMAX) KKMAX = KK_t
                     
                  END DO
               END DO
            END DO
            
            
         END DO      
         
C     -------------------------------------------
C     solve the linear system
C     -------------------------------------------
         
         ibw=(IIG-2)*(JJG-1)+1
         
#if RKIND == 4
         CALL SPBTRS ('U', KKMAX, IBW, 1, ABD, NABD1, BBD, KKMAX, INFO) 
#else
         CALL DPBTRS ('U', KKMAX, IBW, 1, ABD, NABD1, BBD, KKMAX, INFO) 
#endif
         IF (INFO .NE. 0) THEN

            WRITE(*,500) 'Coarse grid solve failed!'
            WRITE(*,510) 'INFO = ', INFO
           
            BMG_iPARMS(id_BMG3_Ext_Err_Code) = INFO
            CALL BMG3_SymStd_ErrTrap(BMG_iPARMS,20)
            RETURN

         ENDIF

      ENDIF
      
C ------------------------------------------------
C     Extract the local data 
C ------------------------------------------------
      
      IF (BMG_iPARMS(id_BMG3_CG_COMM) .eq. BMG_CG_GATHER_SCATTER) THEN

         !
         ! Broadcast the solution vector
         !

         CALL MPI_Bcast(BBD(1), NABD2, FLOAT_PRECISION, 
     &        iZERO, MPICOMM, IERR)
         
      END IF


      KK_t=0
      
      IIL = LocArrSize(1,MyProc) - 2
      JJL = LocArrSize(2,MyProc) - 2
      KKL = LocArrSize(3,MyProc) - 2
      
      P1 = ProcCoord(1,MyProc)
      P2 = ProcCoord(2,MyProc)
      P3 = ProcCoord(3,MyProc)
      
      
      P1SUM = 0
      DO i=1, P1-1
         P1SUM = P1SUM + LocArrSize(1,ProcGrid(i,1,1)) - 2
      END DO
      
      P2SUM = 0
      DO j=1, P2-1
         P2SUM = P2SUM + LocArrSize(2,ProcGrid(1,j,1)) - 2
      END DO
      
      P3SUM = 0
      DO k=1, P3-1
         P3SUM = P3SUM + LocArrSize(3,ProcGrid(1,1,k)) - 2
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
         
      !
      !  Setup loop boundaries in z
      !
      IF ( P3.EQ.1 ) THEN
         kBEG = 1
      ELSE
         kBEG = 0
      END IF

      IF ( P3.EQ.NProcK) THEN
         kEND = KKL
      ELSE
         kEND = KKL + 1
      ENDIF
         


      DO k=kBEG, kEND
         DO j=jBEG, jEND
            DO i=iBEG, iEND
               
               KK_t = (P1SUM + i-1) + (IIG-2) * (P2SUM + j - 1) +
     &              (IIG-2)*(JJG-2) * (P3SUM + k - 1) + 1
               
               Q(i+1,j+1,k+1)=BBD(KK_t)           
               
            END DO
         END DO
      END DO
         
C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SOLVE_cg.f',/,5X,A)
 510  FORMAT (5X,A,1X,I3)

C ===========================================

      RETURN
      END
