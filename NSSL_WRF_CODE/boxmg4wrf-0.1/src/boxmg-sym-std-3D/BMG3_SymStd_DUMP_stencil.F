      SUBROUTINE BMG3_SymStd_DUMP_stencil( 
     &                BMG_IOFLAG, SO, NStncl, kg, NOG,
     &                NLx, NLy, NLz, NGx, NGy, NGz,
     &                iGs, jGs, kGs,
     &                ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &                MyProcI, MyProcJ, MyProcK, BMG_MPI_COMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_DUMP_stencil outputs the stencil interactively
C     to the screen or dumps it to an ASCII text file.  Thus, it
C     it is intended for debugging the stencil generation on 
C     modest sizes grids. 
C     
C     - assumes that SO(p_SO) has been used in the call, where
C       p_SO is returned by BMG3_SymStd_GET_pointers.
C     - assumes that you have provided the value of NStncl
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
C   INPUT/OUTPUT:
C  --------------------
C
C
C =======================================================================
C  --------------------
C   OUTPUT:
C  --------------------
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

      REAL*RKIND   SO(NLx+1,NLy+1,NLz+1,NStncl)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

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

      INTEGER   i, iL, j, jL, k, kL, IO_DONE, IO_CONT, IO_n
      CHARACTER ANS*1, FILEo*30

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

      IF ( BMG_IOFLAG(iBMG3_OUT_STENCIL_TTY) ) THEN

         IF ( MyProc.EQ.1 ) THEN 

            WRITE(*,125) '**** INTERACTIVE STENCIL EVALUATION *****'
            WRITE(*,126) 'Finest grid index = ', NOG
            WRITE(*,128) 'Current grid index = ', kg
            WRITE(*,128) 'Global index ranges including ghost points '
            WRITE(*,129) '[1,NGx_g] = [1,',NGx+2,']'
            WRITE(*,129) '[1,NGy_g] = [1,',NGy+2,']'
            WRITE(*,129) '[1,NGz_g] = [1,',NGz+2,']'
            WRITE(*,128) 'Local index ranges including ghost points '
            WRITE(*,129) '[1,NLx_g] = [1,',NLx+1,']'
            WRITE(*,129) '[1,NLy_g] = [1,',NLy+1,']'
            WRITE(*,129) '[1,NLz_g] = [1,',NLz+1,']'
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

               IF (NStncl.EQ.4) THEN

                  WRITE(*,*)
                  WRITE(*,150) -SO(i,j+1,k,kps)
                  WRITE(*,*)
                  
                  WRITE(*,150) -SO(i,j,k+1,kb)
                  WRITE(*,155) -SO(i,j,k,kpw),
     &                          SO(i,j,k,kp),
     &                         -SO(i+1,j,k,kpw)
                  WRITE(*,160) -SO(i,j,k,kb)
               
                  WRITE(*,*)
                  WRITE(*,150) -SO(i,j,k,kps)
                  WRITE(*,*)
               
               ELSEIF (NStncl.EQ.14) THEN
                  
                  WRITE(*,152) -SO(i,j+1,k+1,kbse),
     &                         -SO(i,j+1,k+1,kbs),
     &                         -SO(i+1,j+1,k+1,kbsw)
                  WRITE(*,155) -SO(i,j+1,k,kpnw),
     &                         -SO(i,j+1,k,kps),
     &                         -SO(i+1,j+1,k,kpsw)
                  WRITE(*,162) -SO(i,j+1,k,kbnw),
     &                         -SO(i,j+1,k,kbn),
     &                         -SO(i+1,j+1,k,kbne)
         
                  WRITE(*,152) -SO(i,j,k+1,kbe),
     &                         -SO(i,j,k+1,kb),
     &                         -SO(i+1,j,k+1,kbw)
                  WRITE(*,155) -SO(i,j,k,kpw),
     &                          SO(i,j,k,kp),
     &                         -SO(i+1,j,k,kpw)
                  WRITE(*,162) -SO(i,j,k,kbw),
     &                         -SO(i,j,k,kb),
     &                         -SO(i+1,j,k,kbe)
               
                  WRITE(*,152) -SO(i,j,k+1,kbne),
     &                         -SO(i,j,k+1,kbn),
     &                         -SO(i+1,j,k+1,kbnw)
                  WRITE(*,155) -SO(i,j,k,kpsw),
     &                         -SO(i,j,k,kps),
     &                         -SO(i+1,j,k,kpnw)
                  WRITE(*,162) -SO(i,j,k,kbsw),
     &                         -SO(i,j,k,kbs),
     &                         -SO(i+1,j,k,kbse)

               ENDIF
            
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

         IO_n = 9 + MyProc

         FILEo = 'output/'
     &           //'Stencil.'//CHAR(ICHAR('0')+kg) 
     &                  //'.'//CHAR(ICHAR('0')+MyProcI)
     &                  //'-'//CHAR(ICHAR('0')+MyProcJ)
     &                  //'-'//CHAR(ICHAR('0')+MyProcK)

         OPEN ( IO_n, FILE=FILEo, STATUS='UNKNOWN' ) 
         
         IF (NStncl.EQ.4) THEN
            
            DO k=2, NLz-1
               DO j=2, NLy-1
                  DO i=2, NLx-1
                     WRITE(IO_n,700) (kGs+k-3)*(NGx*NGy) 
     &                               + (jGs+j-3)*NGx + iGs+i-2,
     &                                iGs+i-1, jGs+j-1, kGs+k-1,
     &                               -SO(i,j+1,k,kps),
     &                               -SO(i,j,k+1,kb),
     &                               -SO(i,j,k,kpw),
     &                                SO(i,j,k,kp),
     &                               -SO(i+1,j,k,kpw),
     &                               -SO(i,j,k,kb),
     &                               -SO(i,j,k,kps)
                  ENDDO
               ENDDO
            ENDDO

         ELSEIF (NStncl.EQ.14) THEN

            DO k=2, NLz-1
               DO j=2, NLy-1
                  DO i=2, NLx-1
                     WRITE(IO_n,800) (kGs+k-3)*(NGx*NGy) 
     &                               + (jGs+j-3)*NGx + iGs+i-2,
     &                                iGs+i-1, jGs+j-1, kGs+k-1, 'N',
     &                               -SO(i,j+1,k+1,kbse),
     &                               -SO(i,j+1,k+1,kbs),
     &                               -SO(i+1,j+1,k+1,kbsw),
     &                               -SO(i,j+1,k,kpnw),
     &                               -SO(i,j+1,k,kps),
     &                               -SO(i+1,j+1,k,kpsw),
     &                               -SO(i,j+1,k,kbnw),
     &                               -SO(i,j+1,k,kbn),
     &                               -SO(i+1,j+1,k,kbne)
         
                     WRITE(IO_n,800) (kGs+k-3)*(NGx*NGy) 
     &                               + (jGs+j-3)*NGx + iGs+i-2,
     &                                iGs+i-1, jGs+j-1, kGs+k-1, 'O',
     &                               -SO(i,j,k+1,kbe),
     &                               -SO(i,j,k+1,kb),
     &                               -SO(i+1,j,k+1,kbw),
     &                               -SO(i,j,k,kpw),
     &                               SO(i,j,k,kp),
     &                               -SO(i+1,j,k,kpw),
     &                               -SO(i,j,k,kbw),
     &                               -SO(i,j,k,kb),
     &                               -SO(i+1,j,k,kbe)
               
                     WRITE(IO_n,800) (kGs+k-3)*(NGx*NGy) 
     &                               + (jGs+j-3)*NGx + iGs+i-2,
     &                                iGs+i-1, jGs+j-1, kGs+k-1, 'S',
     &                               -SO(i,j,k+1,kbne),
     &                               -SO(i,j,k+1,kbn),
     &                               -SO(i+1,j,k+1,kbnw),
     &                               -SO(i,j,k,kpsw),
     &                               -SO(i,j,k,kps),
     &                               -SO(i+1,j,k,kpnw),
     &                               -SO(i,j,k,kbsw),
     &                               -SO(i,j,k,kbs),
     &                               -SO(i+1,j,k,kbse)

                  ENDDO
               ENDDO
            ENDDO
               
         ENDIF
            
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

 700  FORMAT (I6,3(',',1X,I2),1P,7(',',1X,E14.7))
 800  FORMAT (I6,3(',',1X,I2),',',1X,A,1P,9(',',2X,E14.7))

C ===========================================

      RETURN
      END
         
