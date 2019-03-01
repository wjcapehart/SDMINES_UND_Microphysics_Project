      SUBROUTINE EX_SETUP_PDE_parms(
     &                    PDEFILEi
     &                    )
      
C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     EX_SETUP_PDE_parms reads the PDE paramters from a data file.
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

      INCLUDE 'common2.h'

C ----------------------------
C     Argument Declarations
C ----------------------------

      CHARACTER  PDEFILEi*(*)

C ----------------------------
C     Local Declarations
C ----------------------------

      INTEGER   NTMP_iBUFFER, NTMP_rBUFFER
      PARAMETER ( NTMP_iBUFFER = 5, NTMP_rBUFFER = 8*iREGmax + 4 )

      INTEGER   CE, i, IO_in, NCHK_iBUFFER, NCHK_rBUFFER,
     &          READ_STATUS
      LOGICAL   exist
      
      !
      !  pointer to the number of regions
      !
      INTEGER   p_iREG

      !
      !  Region boundary pointers
      !
      INTEGER   p_X1, p_X2, p_Y1, p_Y2

      !
      !  Region coefficient pointers
      ! 
      INTEGER   p_Dix, p_Diy, p_Si, p_Fi

      !
      !  Global domain boundary 
      ! 
      INTEGER   p_xGf, p_xGs, p_yGf, p_yGs

      !
      !  Boundary condition index pointers
      !
      INTEGER   p_ibl, p_ibr, p_ibb, p_ibt

C ==========================================================================

C -----------------------------------
C     Initialize local pointers
C -----------------------------------
      
      !
      !  Number of regions
      !
      p_iREG = 1

      !
      ! Boundary condition indices
      !
      p_ibl  = p_iREG + 1
      p_ibr  = p_ibl + 1
      p_ibb  = p_ibr + 1
      p_ibt  = p_ibb + 1

      NCHK_iBUFFER = p_ibt

      !
      !  PDE coefficients
      !
      p_Dix  = 1
      p_Diy  = p_Dix + iREGmax 

      p_Si   = p_Diy + iREGmax
      p_Fi   = p_Si  + iREGmax

      !
      !  Region boundaries
      !
      p_X1   = p_Fi  + iREGmax
      p_X2   = p_X1  + iREGmax 
      p_Y1   = p_X2  + iREGmax
      p_Y2   = p_Y1  + iREGmax

      !
      !  Global domain boundaries
      !
      p_xGs   = p_Y2 + iREGmax
      p_xGf   = p_xGs + 1
      p_yGs   = p_xGf + 1
      p_yGf   = p_yGs + 1

      NCHK_rBUFFER = p_yGf

      IF ( NCHK_iBUFFER.NE.NTMP_iBUFFER 
     &     .OR. NCHK_rBUFFER.NE.NTMP_rBUFFER ) THEN

         WRITE(*,*) 'ERROR: EX_SETUP_PDE_parms.f ... '
         WRITE(*,*) '  HAVE: NCHK_iBUFFER  = ', NCHK_iBUFFER
         WRITE(*,*) '  HAVE: NTMP_iBUFFER  = ', NTMP_iBUFFER
         WRITE(*,*) '  HAVE: NCHK_rBUFFER  = ', NCHK_rBUFFER
         WRITE(*,*) '  HAVE: NTMP_rBUFFER  = ', NTMP_rBUFFER

         STOP

      ENDIF
      
C -----------------------------------
C     Initialize arrays in common
C -----------------------------------

      DO i=1, iREGmax
         !
         ! Diffusion coefficient (diagonal tensor)
         !
         Dix(i) = rZERO
         Diy(i) = rZERO
         !
         ! Removal/Absorption coefficent
         !
         Si(i)  = rZERO
         !
         ! Source/RHS 
         !
         Fi(i)  = rZERO
         !
         ! Region coordinates (tensor product)
         !
         X1(i)  = rZERO
         X2(i)  = rZERO
         Y1(i)  = rZERO
         Y2(i)  = rZERO
      ENDDO

C -----------------------------------
C     Read parameter data
C -----------------------------------

      !
      ! Eliminate trailing blanks
      !
      CE=LEN(PDEFILEi)
      DO i=LEN(PDEFILEi), 1, -1
         IF (PDEFILEi(i:i).EQ.' ') THEN
            CE=i
         ENDIF
      ENDDO
      CE=CE-1
         
      !
      ! Check that the file exists 
      !
      INQUIRE( FILE=PDEFILEi(1:CE), EXIST=exist )

      IF ( exist ) THEN

         !
         !  Open data file
         !
         IO_in = 10
         OPEN( IO_in, FILE=PDEFILEi(1:CE), STATUS='OLD' ) 

         !
         !  Global physical domain
         !
         READ(IO_in,*) xGs, xGf, yGs, yGf
         !
         !  Global boundary condition indices
         !
         READ(IO_in,*) iBL,iBR,iBB,iBT

         !
         !  The number of regions:
         !
         READ(IO_in,*) iREG
            
         !
         !  Read data for each region:
         !      
         DO i=1, iREG

            READ(IO_in,*) x1(i),x2(i),y1(i),y2(i)
            READ(IO_in,*) dix(i),diy(i),si(i),fi(i)
            
         ENDDO

         !
         !  Close the data file
         !
         CLOSE(IO_in)

         !
         !  READ was successful
         !
         READ_STATUS = 1

      ELSE

         !
         !  READ failed because the file didn't exist
         !
         READ_STATUS = -1

      ENDIF


      !
      !  Output the PDE parameters
      !
         
      WRITE(*,1000) 'The number of regions is ', iREG
      WRITE(*,1020) 'Global Physical Domain:',
     &     '[',xGs,',',xGf,'] x [',
     &     yGs,',',yGf,']' 
      DO i=1, iREG

         WRITE(*,1010) 'Parameters for region ', i
         WRITE(*,1020) 'Physical Domain:',
     &        '[',x1(i),',',x2(i),'] x [',
     &        y1(i),',',y2(i),']' 
         WRITE(*,1030) '(Diagonal) Diffusion Tensor:',
     &        Dix(i), Diy(i)
         WRITE(*,1040) 'Absorption Coefficient:', si(i)
         WRITE(*,1040) 'Source Strength:', fi(i)

      ENDDO

      WRITE(*,1050) 'Boundary Condition Indices'
      WRITE(*,1060) 'x=Xs', ibl
      WRITE(*,1060) 'x=Xf', ibr
      WRITE(*,1060) 'y=Ys', ibb
      WRITE(*,1060) 'y=Yf', ibt

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: EX_SETUP_PDE_parms.f',
     &        //,5X,A,/)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C --------------------------------
C     Region Parameters
C
 1000 FORMAT(/,2X,A,I3,/)
 1005 FORMAT(/,2X,A,$)
 1010 FORMAT(/,2X,30("="),/,4X,A,I3,/,2X,30("="),/)
 1020 FORMAT(2X,A,T35,A,4(F5.2,A),/)
 1030 FORMAT(2X,A,T35,1P,2(E14.7))
 1040 FORMAT(2X,A,T35,1P,E14.7)

C --------------------------------
C     Boundary Condition Indices
C
 1050 FORMAT(/,/,2X,30("="),/,4X,A,/,2X,30("="),/)
 1060 FORMAT(6X,A,T20,I2)
 1070 FORMAT(6X,A,T20,I2,/)

C ===========================================

      RETURN
      END
