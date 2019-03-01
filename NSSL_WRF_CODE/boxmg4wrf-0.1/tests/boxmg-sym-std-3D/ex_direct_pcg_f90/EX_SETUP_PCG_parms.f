      SUBROUTINE EX_SETUP_PCG_parms( 
     &              PCGFILEi, BMG_PCG_iPARMS, BMG_PCG_rPARMS,
     &              MPI_MyProc, NProc, EX_MPI_COMM, EX_MPI_IERROR
     &              )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     EX_SETUP_PCG_parms.f on the master ( MyProc=1 ) first sets 
C     the default PCG parameters, then reads customization
C     parameters from a data file, and finally broadcasts the arrays.
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
      INCLUDE 'BMG_PCG_parameters.h'

C ----------------------------
C     Argument Declarations
C ----------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER  BMG_PCG_iPARMS(NBMG_PCG_iPARMS)
      REAL*8   BMG_PCG_rPARMS(NBMG_PCG_rPARMS)

      INTEGER   EX_MPI_COMM, EX_MPI_IERROR, MPI_MyProc, NProc
      CHARACTER PCGFILEi*(*)

C ----------------------------
C     Local Declarations
C ----------------------------

      INTEGER i, IO_in, MyProc

C ==========================================================================

C ------------------------------------------- 
C     MPI indexing:
C ------------------------------------------- 

      MyProc = MPI_MyProc
      

C ------------------------------------------- 
C     Read parameter data
C ------------------------------------------- 

      IF ( MyProc.EQ.0 ) THEN 
 
         !
         !  Open data file
         !
         IO_in = 10
         OPEN( IO_in, FILE=PCGFILEi, STATUS='OLD' ) 

         !
         !  Read the parameters
         !
         READ(IO_in,*) BMG_PCG_iPARMS(id_BMG_PCG_NMG_CYCLES)
         READ(IO_in,*) BMG_PCG_iPARMS(id_BMG_PCG_STOP_TEST)
         READ(IO_in,*) BMG_PCG_iPARMS(id_BMG_PCG_PRECON)
         READ(IO_in,*) BMG_PCG_iPARMS(id_BMG_PCG_MAX_ITERS)
         READ(IO_in,*) BMG_PCG_iPARMS(id_BMG_PCG_BMG_SETUP)
         READ(IO_in,*) BMG_PCG_rPARMS(id_BMG_PCG_STOP_TOL)
         READ(IO_in,*) BMG_PCG_iPARMS(id_BMG_PCG_OUT_FREQ)

         !  Close the data file
         !
         CLOSE(IO_in)

          !
      ENDIF

C ------------------------------------------- 
C     Broadcast parameters arrays
C ------------------------------------------- 

      CALL MPI_Bcast( BMG_PCG_iPARMS, NBMG_PCG_iPARMS, MPI_INTEGER,
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

      CALL MPI_Bcast( BMG_PCG_rPARMS, NBMG_PCG_rPARMS, MPI_REAL8,
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

C ------------------------------------------- 
C     Override some I/O on the master
C ------------------------------------------- 


C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: EX_SETUP_BMG_parms.f',
     &        //,5X,A,/)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================


      RETURN
      END
