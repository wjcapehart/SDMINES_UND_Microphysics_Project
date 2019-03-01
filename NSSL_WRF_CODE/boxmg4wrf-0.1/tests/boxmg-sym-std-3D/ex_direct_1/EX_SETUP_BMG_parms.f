      SUBROUTINE EX_SETUP_BMG_parms( 
     &              CYCLEFILEi, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &              MPI_MyProc, NProc, EX_MPI_COMM, EX_MPI_IERROR
     &              )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     EX_SETUP_BMG_parms on the master ( MPI_MyProc=0 ) first sets the
C     default cycle parameters and then reads customization parameters
C     from a data file and then broadcasts the arrays.
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
C -----------------------------

      INCLUDE 'mpif.h'

      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_parameters.h'

C ----------------------------
C     Argument Declarations
C ----------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      REAL*8   BMG_rPARMS(NBMG_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

      INTEGER   EX_MPI_COMM, EX_MPI_IERROR, MPI_MyProc, NProc
      CHARACTER CYCLEFILEi*(*)

C ----------------------------
C     Local Declarations
C ----------------------------

      INTEGER   CE, i, IO_in, MyProc, READ_STATUS
      LOGICAL   exist

C ==========================================================================

C ------------------------------------------- 
C     MPI indexing:
C ------------------------------------------- 

      MyProc = MPI_MyProc

C -------------------------------------------
C     Set default parameter values
C -------------------------------------------
      
      IF ( MyProc.EQ.0 ) THEN
      
         CALL BMG3_SymStd_SETUP_parms( 
     &             BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &             ) 
         
      ENDIF

C ------------------------------------------- 
C     Read parameter data
C ------------------------------------------- 

      IF ( MyProc.EQ.0 ) THEN 

         !
         ! Eliminate trailing blanks
         !
         CE=LEN(CYCLEFILEi)
         DO i=LEN(CYCLEFILEi), 1, -1
            IF (CYCLEFILEi(i:i).EQ.' ') THEN
               CE=i
            ENDIF
         ENDDO
         CE=CE-1
         
         !
         ! Check that the file exists 
         !
         INQUIRE( FILE=CYCLEFILEi(1:CE), EXIST=exist )

         IF ( exist ) THEN

            !
            !  Open data file
            !
            IO_in = 10
            OPEN( IO_in, FILE=CYCLEFILEi(1:CE), STATUS='OLD' ) 

            !
            !  Read the parameters
            !
            READ(IO_in,*) BMG_iPARMS(id_BMG3_SETUP)
            READ(IO_in,*) BMG_iPARMS(id_BMG3_STENCIL)
            READ(IO_in,*) BMG_iPARMS(id_BMG3_RELAX),
     &                    BMG_iPARMS(id_BMG3_NRELAX_DOWN),
     &                    BMG_iPARMS(id_BMG3_NRELAX_UP),
     &                    BMG_iPARMS(id_BMG3_RELAX_SYM )

            READ(IO_in,*) BMG_iPARMS(id_BMG3_CYCLE_CLASS),
     &                    BMG_iPARMS(id_BMG3_NCYCLE_TYPE)

            READ(IO_in,*) BMG_iPARMS(id_BMG3_CG_SOLVER),
     &                    BMG_iPARMS(id_BMG3_CG_MIN_DIM)
         
            READ(IO_in,*) BMG_iPARMS(id_BMG2_CYCLE_CLASS), 
     &                    BMG_iPARMS(id_BMG2_NCYCLE_TYPE)

            READ(IO_in,*) BMG_iPARMS(id_BMG2_CG_SOLVER),
     &                    BMG_iPARMS(id_BMG2_CG_MIN_DIM)

            READ(IO_in,*) BMG_iPARMS(id_BMG3_MAX_ITERS),
     &                    BMG_rPARMS(id_BMG3_STOP_TOL)

            !
            !  Close the data file
            !
            CLOSE(IO_in)

            !
            !  Override default stopping criteria
            !
            BMG_iPARMS(id_BMG3_STOP_TEST) = BMG_STOP_REL_RES_L2

            !
            !  Override debbugging flags
            !
            BMG_IOFLAG(iBMG3_BUG_RES_RELAX)    = .TRUE.
            BMG_IOFLAG(iBMG3_BUG_RES_CG_SOLVE) = .TRUE.

            !
            !  Set the Coarse Grid Comm Pattern
            !
            BMG_iPARMS(id_BMG3_CG_COMM) = BMG_CG_GATHER_SCATTER

            BMG_iPARMS(id_BMG3_SYNC_INITIAL_GUESS) 
     &           = BMG_SYNC_INITIAL_GUESS

            READ_STATUS = 1

         ELSE

            READ_STATUS = -1

         ENDIF

      ENDIF

C ------------------------------------------- 
C     Broadcast read status
C ------------------------------------------- 

      CALL MPI_Bcast( READ_STATUS, 1, MPI_INTEGER,
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

      IF ( READ_STATUS.NE.1 ) THEN
         IF ( MyProc.EQ.0 ) THEN
            WRITE(*,*) 'ERROR: File ',CYCLEFILEi(1:CE),
     &                 ' does not exist!'
         ENDIF
         CALL MPI_FINALIZE(EX_MPI_IERROR)
         STOP
      ENDIF

C ------------------------------------------- 
C     Broadcast parameters arrays
C ------------------------------------------- 

      CALL MPI_Bcast( BMG_iPARMS, NBMG_iPARMS, MPI_INTEGER,
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

      CALL MPI_Bcast( BMG_rPARMS, NBMG_rPARMS, MPI_REAL8,
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

      CALL MPI_Bcast( BMG_IOFLAG, NBMG_IOFLAG, MPI_LOGICAL,
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

C ------------------------------------------- 
C     Override some I/O on the master
C ------------------------------------------- 

      IF ( MyProc.EQ.0 ) THEN

         !
         !  Override I/O Parameters
         !
         BMG_IOFLAG(iBMG3_OUT_WSPACE_POINT) = .FALSE.
         BMG_IOFLAG(iBMG3_OUT_WSPACE_SIZE)  = .FALSE.

         BMG_IOFLAG(iBMG3_OUT_ITERATIONS)   = .TRUE.
         BMG_IOFLAG(iBMG3_OUT_TIME_SETUP)   = .TRUE.
         BMG_IOFLAG(iBMG3_OUT_TIME_CYCLING) = .TRUE.
         BMG_IOFLAG(iBMG3_OUT_TIME_TOTAL)   = .TRUE.

         BMG_IOFLAG(iBMG3_BUG_PARAMETERS)   = .FALSE.
         !
         !  Debugging parameters?
         !
         IF ( BMG_IOFLAG(iBMG3_BUG_PARAMETERS) ) THEN
            CALL BMG3_SymStd_DUMP_parms( 
     &                       BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &                       ) 
         ENDIF

      ENDIF

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: EX_SETUP_BMG_parms.f',
     &        //,5X,A,/)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================


      RETURN
      END
