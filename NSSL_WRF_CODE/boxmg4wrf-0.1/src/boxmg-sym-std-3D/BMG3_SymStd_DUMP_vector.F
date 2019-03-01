      SUBROUTINE BMG3_SymStd_DUMP_vector( 
     &                BMG_IOFLAG, Q, kg, NOG,
     &                NLx, NLy, NLz, NGx, NGy, NGz,
     &                iGs, jGs, kGs, FILE_base, GHOST,
     &                ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &                MyProcI, MyProcJ, MyProcK, BMG_MPI_COMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_DUMP_vector outputs a component of a vector
C     interactively to the screen or dumps the entire vector into
C     an ASCII text file. Thus it is intended for debugging on
C     reasonably sized grids.
C
C     - assumes that Q(p_U) has been used in the call, where
C       p_U is returned from BMG3_SymStd_GET_pointers.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   INPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
C
C
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
      INCLUDE 'mpif.h'

#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  iGs, jGs, kGs, kg, NGx, NGy, NGz, 
     &         NOG, NStncl, NLx, NLy, NLz

      REAL*RKIND    Q(NLx,NLy,NLz)
      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG), GHOST
      CHARACTER FILE_base*(*)

      INTEGER  BMG_MPI_COMM, MyProcI, MyProcJ, MyProcK,
     &         NProc, NProcI, NProcJ, NProcK
      INTEGER  ProcGrid(NProcI,NProcJ,NProcK)

C --------------------------
C     Local Declarations:
C

      INTEGER   NTMP_iBUFFER
      PARAMETER ( NTMP_iBUFFER=7 )

      INTEGER   TMP_iBUFFER(NTMP_iBUFFER)
      INTEGER   IO_Proc, p_iL, p_jL, p_kL,
     &          p_Proc, p_iProc, p_jProc, p_kProc
      
      INTEGER   MPI_IERROR, MPI_MyProc, MyProc, iProc, jProc, kProc

      INTEGER   IEND, ISTART, JEND, JSTART, KEND, KSTART

      INTEGER   i, iL, j, jL, k, kL, IO_DONE, IO_CONT, IO_n
      CHARACTER ANS*1, FILEo*50

C ==========================================================================

      CALL MPI_COMM_RANK( BMG_MPI_COMM, MPI_MyProc, MPI_IERROR )

      MyProc = MPI_MyProc+1

      p_iL = 1
      p_jL = 2
      p_kL = 3
      
      p_iProc = 4
      p_jProc = 5
      p_kProc = 6
      p_Proc  = 7

      IF ( .FALSE. ) THEN

         IF ( MyProc.EQ.1 ) THEN 

            WRITE(*,125) '**** INTERACTIVE VECTOR EVALUATION *****'
            WRITE(*,126) 'Finest grid index = ', NOG
            WRITE(*,128) 'Current grid index = ', kg
            WRITE(*,128) 'Global index ranges including ghost points '
            WRITE(*,129) '[1,NGx_g] = [1,',NGx+2,']'
            WRITE(*,129) '[1,NGy_g] = [1,',NGy+2,']'
            WRITE(*,129) '[1,NGz_g] = [1,',NGz+2,']'
            WRITE(*,128) 'Local index ranges including ghost points '
            WRITE(*,129) '[1,NLx_g] = [1,',NLx,']'
            WRITE(*,129) '[1,NLy_g] = [1,',NLy,']'
            WRITE(*,129) '[1,NLz_g] = [1,',NLz,']'
            WRITE(*,*) 
            
         ENDIF

 10      CONTINUE

            IF ( MyProc.EQ.1 ) THEN

 20            CONTINUE  ! <<<<<<< BEGIN: Data entry

                  WRITE(*,140) 'Enter coordinates ProcGrid (i,j,k) ... '
                  READ(*,*) iProc, jProc, kProc

                  WRITE(*,140) 'Enter local array coordinates ... '
                  READ(*,*) iL, jL, kL
               
                  IO_Proc = ProcGrid(iProc,jProc,kProc)

                  TMP_iBUFFER(p_Proc)  = IO_Proc
                  
                  TMP_iBUFFER(p_iProc) = iProc
                  TMP_iBUFFER(p_jProc) = jProc
                  TMP_iBUFFER(p_kProc) = kProc
               
                  TMP_iBUFFER(p_iL) = iL
                  TMP_iBUFFER(p_jL) = jL
                  TMP_iBUFFER(p_kL) = kL
               
                  IF ( iL.LT.1 .OR. iL.GT.NLx
     &               .OR. jL.LT.1 .OR. jL.GT.NLy
     &               .OR. kL.LT.1 .OR. kL.GT.NLz ) THEN

                     !
                     ! Need local arrays instead!!!
                     !
                     WRITE(*,*) 'error index out of bounds'
                     GOTO 10

                  ENDIF

 30            CONTINUE      ! <<<<<<< END: Data entry

            ENDIF


            CALL MPI_Bcast( TMP_iBUFFER, NTMP_iBUFFER, MPI_INTEGER,
     &                      0, BMG_MPI_COMM, MPI_IERROR )



            IF ( MyProc.NE.1 ) THEN

               IO_Proc=TMP_iBUFFER(p_Proc)
               
               iProc = TMP_iBUFFER(p_iProc)
               jProc = TMP_iBUFFER(p_jProc)
               kProc = TMP_iBUFFER(p_kProc)
               
               iL = TMP_iBUFFER(p_iL)
               jL = TMP_iBUFFER(p_jL)
               kL = TMP_iBUFFER(p_kL)
               
            ENDIF



            IF ( MyProc.EQ.IO_Proc ) THEN 

               i=iL
               j=jL
               k=kL
               
               WRITE(*,150) Q(i,j,k)
            
               IO_DONE=1

            ENDIF

            CALL MPI_Bcast( IO_DONE, iONE, MPI_INTEGER,
     &                      IO_Proc-1, BMG_MPI_COMM, MPI_IERROR )

            IF ( MyProc.EQ.1 ) THEN
            
 40            CONTINUE

                  WRITE(*,145) 'Another stencil (Y/N) .... '
                  READ(*,'(A)') ANS

                  IF (INDEX('Yy',ANS).NE.0) THEN
                     IO_CONT = iONE
                  ELSE IF (INDEX('Nn',ANS).EQ.0) THEN
                     WRITE(*,*) 'Error .... '
                     GOTO 40
                  ELSE
                     IO_CONT = iZERO
                  ENDIF

 50            CONTINUE

            ENDIF

            CALL MPI_Bcast( IO_CONT, iONE, MPI_INTEGER,
     &                      0, BMG_MPI_COMM, MPI_IERROR )
         
            IF ( IO_CONT.EQ.iONE ) THEN
               GOTO 10
            ENDIF

 60      CONTINUE

      ELSE

         IF ( GHOST ) THEN 
            ISTART = 1
            JSTART = 1
            KSTART = 1
            IEND   = NLx
            JEND   = NLy
            KEND   = NLz
         ELSE
            ISTART = 2
            JSTART = 2
            KSTART = 2
            IEND   = NLx-1
            JEND   = NLy-1
            KEND   = NLz-1
         ENDIF

         IO_n = 9 + MyProc

         FILEo = 'output/'
     &           //FILE_base//'.'//CHAR(ICHAR('0')+kg)
     &                      //'.'//CHAR(ICHAR('0')+MyProcI)
     &                      //'-'//CHAR(ICHAR('0')+MyProcJ)
     &                      //'-'//CHAR(ICHAR('0')+MyProcK)

         OPEN ( IO_n, FILE=FILEo, STATUS='UNKNOWN' ) 
         
         DO k=KSTART, KEND
            DO j=JSTART, JEND
               DO i=ISTART, IEND
                  WRITE(IO_n,700) iGs+i-1, jGs+j-1, kGs+k-1,
     &                            Q(i,j,k)
               ENDDO
            ENDDO
         ENDDO
         
         CLOSE( IO_n )
         

      ENDIF

C ==========================================================================

 125  FORMAT (/,/,6X,A,/)
 126  FORMAT (/,8X,A,I2)
 128  FORMAT (8X,A,I2,/)
 129  FORMAT (10X,A,I5,A)

 130  FORMAT (/,4X,A)
 140  FORMAT (4X,A,$)
 145  FORMAT (/,4X,A,$)

 150  FORMAT (/,24X,1P,E12.5)
 152  FORMAT (/,8X,1P,3(E12.5,4X))
 155  FORMAT (8X,1P,3(E12.5,4X))
 160  FORMAT (24X,1P,E12.5,/)
 162  FORMAT (8X,1P,3(E12.5,4X),/)


 170  FORMAT (/,24X,1P,E12.5)
 175  FORMAT (16X,1P,2(E12.5,4X))
 180  FORMAT (24X,1P,E12.5)
 185  FORMAT (24X,1P,E12.5,/)

 700  FORMAT (I2,2(',',1X,I2),1P,',',1X,E16.8)

C ===========================================

      RETURN
      END
         
