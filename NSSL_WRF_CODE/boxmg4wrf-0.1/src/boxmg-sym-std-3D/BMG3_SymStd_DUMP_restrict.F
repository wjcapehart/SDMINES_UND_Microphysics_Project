      SUBROUTINE BMG3_SymStd_DUMP_restrict( 
     &                BMG_IOFLAG, CI, KFG, KCG, NOG,
     &                NLx, NLy, NLz, NGx, NGy, NGz,
     &                NLx_c, NLy_c, NLz_c, NGx_c, NGy_c, NGz_c,
     &                iGs, jGs, kGs, iGs_c, jGs_c, kGs_c,
     &                ProcGrid, NProcI, NProcJ, NProcK, NProc,
     &                MyProcI, MyProcJ, MyProcK, BMG_MPI_COMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_DUMP_restrict outputs the interpolation weights
C     interactively to the screen or dumps the entire array
C     of interpolation weights into an ASCII text file. Thus it is 
C     intended for debugging the stencil generation on reasonable size 
C     grids.
C     
C     - assumes that CI(p_CIC) has been used in the call, where
C       p_CIC is returned from BMG3_SymStd_GET_pointers.
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
      INTEGER  iGs, iGs_c, jGs, jGs_c, kGs, kGs_c,
     &         KCG, KFG,
     &         NGx, NGy, NGz, NGx_c, NGy_c, NGz_c,
     &         NOG, NStncl,
     &         NLx, NLy, NLz, NLx_c, NLy_c, NLz_c

      REAL*RKIND   CI(NLx_c,NLy_c,NLz_c,26)
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
      INTEGER   IO_Proc, p_iLc, p_jLc, p_kLc,
     &          p_Proc, p_iProc, p_jProc, p_kProc
      
      INTEGER   MPI_IERROR, MPI_MyProc, MyProc, iProc, jProc, kProc

      INTEGER   ic, iLc, jc, jLc, kc, kLc, IO_DONE, IO_CONT, IO_n
      CHARACTER ANS*1, FILEo*30

C ==========================================================================

      CALL MPI_COMM_RANK( BMG_MPI_COMM, MPI_MyProc, MPI_IERROR )

      MyProc = MPI_MyProc+1

      p_iLc = 1
      p_jLc = 2
      p_kLc = 3
      
      p_iProc = 4
      p_jProc = 5
      p_kProc = 6
      p_Proc  = 7

      IF ( BMG_IOFLAG(iBMG3_OUT_RESTRICT_TTY) ) THEN

         IF ( MyProc.EQ.1 ) THEN 

            WRITE(*,125) '**** INTERACTIVE RESTRICTION EVALUATION *****'
            WRITE(*,126) 'Finest grid index = ', NOG
            WRITE(*,126) 'Fine grid index   = ', KFG
            WRITE(*,128) 'Coarse grid index = ', KCG
            WRITE(*,127) 'Global coarse-grid index ranges ' 
     &                   //'(including ghost points) ... '
            WRITE(*,129) '[1,NGx_c] = [1,',NGx_c,']'
            WRITE(*,129) '[1,NGy_c] = [1,',NGy_c,']'
            WRITE(*,129) '[1,NGz_c] = [1,',NGz_c,']'
            WRITE(*,127) 'Local fine-grid index ranges '
     &                   //'(including ghost points) ... '
            WRITE(*,129) '[1,NLx_c] = [1,',NLx_c+1,']'
            WRITE(*,129) '[1,NLy_c] = [1,',NLy_c+1,']'
            WRITE(*,129) '[1,NLz_c] = [1,',NLz_c+1,']'
            WRITE(*,*) 
            
         ENDIF

 10      CONTINUE

            IF ( MyProc.EQ.1 ) THEN

 20            CONTINUE  ! <<<<<<< BEGIN: Data entry

                  WRITE(*,140) 'Enter coordinates ProcGrid (i,j,k) ... '
                  READ(*,*) iProc, jProc, kProc

                  WRITE(*,140) 'Enter the local coarse-grid '
     &                         //'coordinates ... '
                  READ(*,*) iLc, jLc, kLc
               
                  IO_Proc = ProcGrid(iProc,jProc,kProc)
               
                  TMP_iBUFFER(p_Proc)  = IO_Proc
                  
                  TMP_iBUFFER(p_iProc) = iProc
                  TMP_iBUFFER(p_jProc) = jProc
                  TMP_iBUFFER(p_kProc) = kProc
               
                  TMP_iBUFFER(p_iLc) = iLc
                  TMP_iBUFFER(p_jLc) = jLc
                  TMP_iBUFFER(p_kLc) = kLc
               
                  IF ( iLc.LT.1 .OR. iLc.GT.NLx_c
     &               .OR. jLc.LT.1 .OR. jLc.GT.NLy_c
     &               .OR. kLc.LT.1 .OR. kLc.GT.NLz_c ) THEN

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
               
               iLc = TMP_iBUFFER(p_iLc)
               jLc = TMP_iBUFFER(p_jLc)
               kLc = TMP_iBUFFER(p_kLc)
               
            ENDIF



            IF ( MyProc.EQ.IO_Proc ) THEN 

               ic=iLc
               jc=jLc
               kc=kLc

               WRITE(*,152)  CI(ic,jc,kc+1,lbne),
     &                       CI(ic,jc,kc+1,lxzse),
     &                       CI(ic,jc+1,kc+1,lbse)
               WRITE(*,155)  CI(ic,jc,kc,lxyne),
     &                       CI(ic,jc,kc,lxyr),
     &                       CI(ic,jc+1,kc,lxyse)
               WRITE(*,162)  CI(ic,jc,kc,ltne),
     &                       CI(ic,jc,kc,lxzne),
     &                       CI(ic,jc+1,kc,ltse)

               WRITE(*,152)  CI(ic,jc,kc+1,lyzsw),
     &                       CI(ic,jc,kc+1,lxzb),
     &                       CI(ic,jc+1,kc+1,lyzse)
               WRITE(*,155)  CI(ic,jc,kc,lxya),
     &                       rONE,
     &                       CI(ic,jc+1,kc,lxyb)
               WRITE(*,162)  CI(ic,jc,kc,lyznw),
     &                       CI(ic,jc,kc,lxza),
     &                       CI(ic,jc+1,kc,lyzne)

               WRITE(*,152) CI(ic+1,jc,kc+1,lbnw),
     &                      CI(ic+1,jc,kc+1,lxzsw),
     &                      CI(ic+1,jc+1,kc+1,lbsw)
               WRITE(*,155) CI(ic+1,jc,kc,lxynw),
     &                      CI(ic+1,jc,kc,lxyl),
     &                      CI(ic+1,jc+1,kc,lxysw)
               WRITE(*,162) CI(ic+1,jc,kc,ltnw),
     &                      CI(ic+1,jc,kc,lxznw),
     &                      CI(ic+1,jc+1,kc,ltsw)
 
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

         IO_n = 9 + MyProc

         FILEo = 'output/'
     &           //'Restrict.'//CHAR(ICHAR('0')+KCG) 
     &                  //'.'//CHAR(ICHAR('0')+MyProcI)
     &                  //'-'//CHAR(ICHAR('0')+MyProcJ)
     &                  //'-'//CHAR(ICHAR('0')+MyProcK)

         OPEN ( IO_n, FILE=FILEo, STATUS='UNKNOWN' ) 
         
         DO kc=2, NLz_c-1
            DO jc=2, NLy_c-1
               DO ic=2, NLx_c-1

                  WRITE(IO_n,800) iGs_c+ic-1, jGs_c+jc-1, kGs_c+kc-1,
     &                            'N',
     &                            CI(ic,jc,kc+1,lbne),
     &                            CI(ic,jc,kc+1,lxzse),
     &                            CI(ic,jc+1,kc+1,lbse),
     &                            CI(ic,jc,kc,lxyne),
     &                            CI(ic,jc,kc,lxyr),
     &                            CI(ic,jc+1,kc,lxyse),
     &                            CI(ic,jc,kc,ltne),
     &                            CI(ic,jc,kc,lxzne),
     &                            CI(ic,jc+1,kc,ltse)

                  WRITE(IO_n,800) iGs_c+ic-1, jGs_c+jc-1, kGs_c+kc-1,
     &                            'O',
     &                            CI(ic,jc,kc+1,lyzsw),
     &                            CI(ic,jc,kc+1,lxzb),
     &                            CI(ic,jc+1,kc+1,lyzse),
     &                            CI(ic,jc,kc,lxya),
     &                            rONE,
     &                            CI(ic,jc+1,kc,lxyb),
     &                            CI(ic,jc,kc,lyznw),
     &                            CI(ic,jc,kc,lxza),
     &                            CI(ic,jc+1,kc,lyzne)

                  WRITE(IO_n,800) iGs_c+ic-1, jGs_c+jc-1, kGs_c+kc-1,
     &                            'S',
     &                            CI(ic+1,jc,kc+1,lbnw),
     &                            CI(ic+1,jc,kc+1,lxzsw),
     &                            CI(ic+1,jc+1,kc+1,lbsw),
     &                            CI(ic+1,jc,kc,lxynw),
     &                            CI(ic+1,jc,kc,lxyl),
     &                            CI(ic+1,jc+1,kc,lxysw),
     &                            CI(ic+1,jc,kc,ltnw),
     &                            CI(ic+1,jc,kc,lxznw),
     &                            CI(ic+1,jc+1,kc,ltsw)


               ENDDO
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

 700  FORMAT (I2,2(',',1X,I2),1P,7(',',1X,E14.7))
 800  FORMAT (I2,2(',',1X,I2),1X,A,1P,9(',',1X,E14.7))

C ===========================================

      RETURN
      END
         
