      SUBROUTINE EX_SETUP_BMG_SER_parms( 
     &              CYCLEFILEi, BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG
     &              )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     EX_SETUP_BMG_SER_parms sets the default cycle parameters and then
C     reads customization parameters from a data file.
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

      INCLUDE 'BMG_SER_constants.h'
      INCLUDE 'BMG_SER_parameters.h'

C ----------------------------
C     Argument Declarations
C ----------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER  BMG_iPARMS(NBMG_SER_iPARMS)
      REAL*8   BMG_rPARMS(NBMG_SER_rPARMS)
      LOGICAL  BMG_IOFLAG(NBMG_SER_IOFLAG)

      CHARACTER CYCLEFILEi*(*)

C ----------------------------
C     Local Declarations
C ----------------------------

      INTEGER   CE, i, IO_in, READ_STATUS
      LOGICAL   exist

C ==========================================================================

C -------------------------------------------
C     Set default parameter values
C -------------------------------------------
      
      
      CALL BMG2_SER_SymStd_SETUP_parms( 
     &     BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
     &     ) 
      
C ------------------------------------------- 
C     Read parameter data
C ------------------------------------------- 

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
         READ(IO_in,*) BMG_iPARMS(id_BMG2_SER_SETUP)
         READ(IO_in,*) BMG_iPARMS(id_BMG2_SER_STENCIL)
         READ(IO_in,*) BMG_iPARMS(id_BMG2_SER_RELAX),
     &                 BMG_iPARMS(id_BMG2_SER_NRELAX_DOWN),
     &                 BMG_iPARMS(id_BMG2_SER_NRELAX_UP),
     &                 BMG_iPARMS(id_BMG2_SER_RELAX_SYM )

         READ(IO_in,*) BMG_iPARMS(id_BMG2_SER_CYCLE_CLASS),
     &                 BMG_iPARMS(id_BMG2_SER_NCYCLE_TYPE)

         READ(IO_in,*) BMG_iPARMS(id_BMG2_SER_CG_MIN_DIM)
         
         READ(IO_in,*) BMG_iPARMS(id_BMG2_SER_MAX_ITERS),
     &                 BMG_rPARMS(id_BMG2_SER_STOP_TOL)

         READ(IO_in,*) BMG_iPARMS(id_BMG2_SER_CG_TYPE),
     &                 BMG_iPARMS(id_BMG2_SER_CG_CONSTRUCT)

         !
         !  Close the data file
         !
         CLOSE(IO_in)

         !
         !  Override default stopping criteria
         !
         BMG_iPARMS(id_BMG2_SER_STOP_TEST) = BMG_SER_STOP_REL_RES_L2

         !
         !  Override debbugging flags
         !
         BMG_IOFLAG(iBMG2_SER_BUG_RES_RELAX)    = .TRUE.
         BMG_IOFLAG(iBMG2_SER_BUG_RES_CG_SOLVE) = .TRUE.

         READ_STATUS = 1

      ELSE

         READ_STATUS = -1

      ENDIF

C ------------------------------------------- 
C     Override some I/O 
C ------------------------------------------- 

      !
      !  Override I/O Parameters
      !
      BMG_IOFLAG(iBMG2_SER_OUT_WSPACE_POINT) = .FALSE.
      BMG_IOFLAG(iBMG2_SER_OUT_WSPACE_SIZE)  = .FALSE.

      BMG_IOFLAG(iBMG2_SER_OUT_ITERATIONS)   = .TRUE.
      BMG_IOFLAG(iBMG2_SER_OUT_TIME_SETUP)   = .TRUE.
      BMG_IOFLAG(iBMG2_SER_OUT_TIME_CYCLING) = .TRUE.
      BMG_IOFLAG(iBMG2_SER_OUT_TIME_TOTAL)   = .TRUE.
      
      BMG_IOFLAG(iBMG2_SER_BUG_PARAMETERS)   = .FALSE.
      !
      !  Debugging parameters?
      !
c$$$         IF ( BMG_IOFLAG(iBMG2_SER_BUG_PARAMETERS) ) THEN
c$$$            CALL BMG2_SER_SymStd_DUMP_parms( 
c$$$     &                       BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG 
c$$$     &                       ) 
c$$$         ENDIF

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: EX_SETUP_BMG_SER_parms.f',
     &        //,5X,A,/)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================


      RETURN
      END
