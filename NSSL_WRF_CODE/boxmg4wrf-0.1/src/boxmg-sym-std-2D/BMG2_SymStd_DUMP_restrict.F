      SUBROUTINE BMG2_SymStd_DUMP_restrict( 
     &                BMG_IOFLAG, CI, KFG, KCG, NOG,
     &                NLx, NLy, NGx, NGy, 
     &                NLx_c, NLy_c, NGx_c, NGy_c,
     &                iGs, jGs, iGs_c, jGs_c, 
     &                ProcGrid, NProcI, NProcJ, NProc,
     &                MyProcI, MyProcJ, BMG_MPI_COMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_DUMP_restrict outputs the interpolation weights
C     interactively to the screen or dumps the entire array
C     of interpolation weights into an ASCII text file. Thus it is 
C     intended for debugging the stencil generation on reasonable size 
C     grids.
C     
C     - assumes that CI(p_CIC) has been used in the call, where
C       p_CIC is returned from BMG2_SymStd_GET_pointers.
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
      INTEGER  iGs, iGs_c, jGs, jGs_c, 
     &         KCG, KFG,
     &         NGx, NGy, NGx_c, NGy_c, 
     &         NOG, NStncl,
     &         NLx, NLy, NLx_c, NLy_c

      REAL*RKIND   CI(NLx_c,NLy_c,26)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER  BMG_MPI_COMM, MyProcI, MyProcJ,
     &         NProc, NProcI, NProcJ
      INTEGER  ProcGrid(NProcI,NProcJ)

C --------------------------
C     Local Declarations:
C

      INTEGER   NTMP_iBUFFER
      PARAMETER ( NTMP_iBUFFER=7 )

      INTEGER   TMP_iBUFFER(NTMP_iBUFFER)
      INTEGER   IO_Proc, p_iLc, p_jLc, 
     &          p_Proc, p_iProc, p_jProc
      
      INTEGER   MPI_IERROR, MPI_MyProc, MyProc, iProc, jProc

      INTEGER   ic, iLc, jc, jLc, IO_DONE, IO_CONT, IO_n
      CHARACTER ANS*1, FILEo*30

C ==========================================================================

      CALL MPI_COMM_RANK( BMG_MPI_COMM, MPI_MyProc, MPI_IERROR )

      MyProc = MPI_MyProc+1

      p_iLc = 1
      p_jLc = 2
      
      p_iProc = 3
      p_jProc = 4
      p_Proc  = 5

      IF ( BMG_IOFLAG(iBMG2_OUT_RESTRICT_TTY) ) THEN

         IF ( MyProc.EQ.1 ) THEN 

            WRITE(*,125) '**** INTERACTIVE RESTRICTION EVALUATION *****'
            WRITE(*,126) 'Finest grid index = ', NOG
            WRITE(*,126) 'Fine grid index   = ', KFG
            WRITE(*,128) 'Coarse grid index = ', KCG
            WRITE(*,127) 'Global coarse-grid index ranges ' 
     &                   //'(including ghost points) ... '
            WRITE(*,129) '[1,NGx_c] = [1,',NGx_c,']'
            WRITE(*,129) '[1,NGy_c] = [1,',NGy_c,']'
            WRITE(*,127) 'Local fine-grid index ranges '
     &                   //'(including ghost points) ... '
            WRITE(*,129) '[1,NLx_c] = [1,',NLx_c+1,']'
            WRITE(*,129) '[1,NLy_c] = [1,',NLy_c+1,']'
            WRITE(*,*) 
            
         ENDIF

 10      CONTINUE

            IF ( MyProc.EQ.1 ) THEN

 20            CONTINUE  ! <<<<<<< BEGIN: Data entry

                  WRITE(*,140) 'Enter coordinates ProcGrid (i,j,k) ... '
                  READ(*,*) iProc, jProc

                  WRITE(*,140) 'Enter the local coarse-grid '
     &                         //'coordinates ... '
                  READ(*,*) iLc, jLc
               
                  IO_Proc = ProcGrid(iProc,jProc)
               
                  TMP_iBUFFER(p_Proc)  = IO_Proc
                  
                  TMP_iBUFFER(p_iProc) = iProc
                  TMP_iBUFFER(p_jProc) = jProc
               
                  TMP_iBUFFER(p_iLc) = iLc
                  TMP_iBUFFER(p_jLc) = jLc
               
                  IF ( iLc.LT.1 .OR. iLc.GT.NLx_c
     &               .OR. jLc.LT.1 .OR. jLc.GT.NLy_c ) THEN

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
               
               iLc = TMP_iBUFFER(p_iLc)
               jLc = TMP_iBUFFER(p_jLc)
               
            ENDIF



            IF ( MyProc.EQ.IO_Proc ) THEN 

               ic=iLc
               jc=jLc
               
               WRITE(*,152) CI(ic,jc+1,LSE),
     &                      CI(ic,jc+1,LB),
     &                      CI(ic+1,jc+1,LSW)

               WRITE(*,155) CI(ic,jc,LR),
     &                      rONE,
     &                      CI(ic+1,jc,LL)

               WRITE(*,162) CI(ic,jc,LNE),
     &                      CI(ic,jc,LA),
     &                      CI(ic+1,jc,LNW)
 
               IO_DONE=1

            ENDIF

            CALL MPI_Bcast( IO_DONE, iONE, MPI_INTEGER,
     &                      IO_Proc-1, BMG_MPI_COMM, MPI_IERROR )

            IF ( MyProc.EQ.1 ) THEN
            
 40            CONTINUE

                  WRITE(*,145) 'Another restriction (Y/N) .... '
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

         IO_n = 6 + MyProc

         FILEo = 'output/'
     &           //'Restrict.'//CHAR(ICHAR('0')+KCG) 
     &                  //'.'//CHAR(ICHAR('0')+MyProcI)
     &                  //'-'//CHAR(ICHAR('0')+MyProcJ)

         OPEN ( IO_n, FILE=FILEo, STATUS='UNKNOWN' ) 
         
         DO jc=2, NLy_c-1
            DO ic=2, NLx_c-1

               WRITE(IO_n,800) iGs_c+ic-1, jGs_c+jc-1,
     &                         CI(ic,jc+1,LSE),
     &                         CI(ic,jc+1,LB),
     &                         CI(ic+1,jc+1,LSW),
     &                         CI(ic,jc,LR),
     &                         rONE,
     &                         CI(ic+1,jc,LL),
     &                         CI(ic,jc,LNE),
     &                         CI(ic,jc,LA),
     &                         CI(ic+1,jc,LNW)

            ENDDO
         ENDDO

         CLOSE( IO_n )

      ENDIF

C ==========================================================================

 125  FORMAT (/,/,6X,A,/)
 126  FORMAT (/,8X,A,I2)
 127  FORMAT (8X,A,/)
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

 800  FORMAT (I2,',',1X,I2,1X,1P,9(',',1X,E14.7))

C ===========================================

      RETURN
      END
         
