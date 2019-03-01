      SUBROUTINE EX_SETUP_GRID_size( 
     &              GRIDFILEi, Nx, Ny, Nz
     &              )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     EX_SETUP_GRID_size.f reads x, y, and z dimensions for the problem.
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

C ----------------------------
C     Argument Declarations
C ----------------------------

      !
      ! BoxMG Cycle and I/O Parameters
      !
      INTEGER  Nx, Ny, Nz

      CHARACTER GRIDFILEi*(*)

C ----------------------------
C     Local Declarations
C ----------------------------

      INTEGER   CE, i, IO_in, READ_STATUS
      LOGICAL   exist

C ==========================================================================

      !
      ! Eliminate trailing blanks
      !
      CE=LEN(GRIDFILEi)
      DO i=LEN(GRIDFILEi), 1, -1
         IF (GRIDFILEi(i:i).EQ.' ') THEN
            CE=i
         ENDIF
      ENDDO
      CE=CE-1
         
      !
      ! Check that the file exists 
      !
      INQUIRE( FILE=GRIDFILEi(1:CE), EXIST=exist )
      
      IF ( exist ) THEN

         !
         !  Open data file
         !
         IO_in = 10
         OPEN( IO_in, FILE=GRIDFILEi(1:CE), STATUS='OLD' ) 

         !
         !  Read the parameters
         !
         READ(IO_in,*) Nx, Ny, Nz

         !
         !  Close the data file
         !
         CLOSE(IO_in)
         
         READ_STATUS = 1

      ELSE

         READ_STATUS = -1

      ENDIF

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: EX_SETUP_GRID_size.f',
     &        //,5X,A,/)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================


      RETURN
      END
