      SUBROUTINE BMG2_SymStd_DUMP_vector( 
     &                BMG_IOFLAG, Q, kg, NOG,
     &                NLx, NLy, NGx, NGy, 
     &                iGs, jGs, FILE_base, GHOST,
     &                ProcGrid, NProcI, NProcJ, NProc,
     &                MyProcI, MyProcJ, BMG_MPI_COMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_DUMP_vector.f outputs a component of a vector
C     interactively to the screen or dumps the entire vector into
C     an ASCII text file. Thus it is intended for debugging on
C     reasonably sized grids.
C
C     - assumes that Q(p_U) has been used in the call, where
C       p_U is returned from BMG2_SymStd_GET_pointers.
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

#include      "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_parameters.h'

C ---------------------------
C     Argument Declarations:
C
      INTEGER  iGs, jGs, kg, NGx, NGy,  
     &         NOG, NLx, NLy

      REAL*RKIND    Q(NLx,NLy)
      LOGICAL   BMG_IOFLAG(NBMG_IOFLAG), GHOST
      CHARACTER FILE_base*(*)

      INTEGER  BMG_MPI_COMM, MyProcI, MyProcJ, 
     &         NProc, NProcI, NProcJ
      INTEGER  ProcGrid(NProcI,NProcJ)

C --------------------------
C     Local Declarations:
C

      INTEGER   NTMP_iBUFFER
      PARAMETER ( NTMP_iBUFFER=7 )

      INTEGER   TMP_iBUFFER(NTMP_iBUFFER)
      INTEGER   IO_Proc, p_iL, p_jL, 
     &          p_Proc, p_iProc, p_jProc
      
      INTEGER   MPI_IERROR, MPI_MyProc, MyProc, iProc, jProc

      INTEGER   IEND, ISTART, JEND, JSTART

      INTEGER   i, iL, j, jL, IO_DONE, IO_CONT, IO_n
      CHARACTER ANS*1, FILEo*50

C ==========================================================================

      CALL MPI_COMM_RANK( BMG_MPI_COMM, MPI_MyProc, MPI_IERROR )

      MyProc = MPI_MyProc+1

      p_iL = 1
      p_jL = 2
      
      p_iProc = 3
      p_jProc = 4
      p_Proc  = 5

      IF ( .FALSE. ) THEN

         IF ( MyProc.EQ.1 ) THEN 

            WRITE(*,125) '**** INTERACTIVE VECTOR EVALUATION *****'
            WRITE(*,126) 'Finest grid index = ', NOG
            WRITE(*,128) 'Current grid index = ', kg
            WRITE(*,128) 'Global index ranges including ghost points '
            WRITE(*,129) '[1,NGx_g] = [1,',NGx+2,']'
            WRITE(*,129) '[1,NGy_g] = [1,',NGy+2,']'
            WRITE(*,128) 'Local index ranges including ghost points '
            WRITE(*,129) '[1,NLx_g] = [1,',NLx,']'
            WRITE(*,129) '[1,NLy_g] = [1,',NLy,']'
            WRITE(*,*) 
            
         ENDIF

 10      CONTINUE

            IF ( MyProc.EQ.1 ) THEN

 20            CONTINUE  ! <<<<<<< BEGIN: Data entry

                  WRITE(*,140) 'Enter coordinates ProcGrid (i,j,k) ... '
                  READ(*,*) iProc, jProc

                  WRITE(*,140) 'Enter local array coordinates ... '
                  READ(*,*) iL, jL
               
                  IO_Proc = ProcGrid(iProc,jProc)

                  TMP_iBUFFER(p_Proc)  = IO_Proc
                  
                  TMP_iBUFFER(p_iProc) = iProc
                  TMP_iBUFFER(p_jProc) = jProc
               
                  TMP_iBUFFER(p_iL) = iL
                  TMP_iBUFFER(p_jL) = jL
               
                  IF ( iL.LT.1 .OR. iL.GT.NLx
     &               .OR. jL.LT.1 .OR. jL.GT.NLy ) THEN

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
               
               iL = TMP_iBUFFER(p_iL)
               jL = TMP_iBUFFER(p_jL)
               
            ENDIF



            IF ( MyProc.EQ.IO_Proc ) THEN 

               i=iL
               j=jL
               
               WRITE(*,150) Q(i,j)
            
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
            IEND   = NLx
            JEND   = NLy
         ELSE
            ISTART = 2
            JSTART = 2
            IEND   = NLx-1
            JEND   = NLy-1
         ENDIF

         IO_n = 6 + MyProc

         FILEo = 'output/'
     &           //FILE_base//'.'//CHAR(ICHAR('0')+kg)
     &                      //'.'//CHAR(ICHAR('0')+MyProcI)
     &                      //'-'//CHAR(ICHAR('0')+MyProcJ)

         OPEN ( IO_n, FILE=FILEo, STATUS='UNKNOWN' ) 
         
         DO j=JSTART, JEND
            DO i=ISTART, IEND
               WRITE(IO_n,700) iGs+i-1, jGs+j-1, Q(i,j)
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

 700  FORMAT (I3,1(',',1X,I3),1P,',',1X,E16.8)

C ===========================================

      RETURN
      END
         
