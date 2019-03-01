      SUBROUTINE EX_SETUP_PDE_parms(
     &                    PDEFILEi, 
     &                    MPI_MyProc, NProc, EX_MPI_COMM, EX_MPI_IERROR
     &                    )
      
C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     EX_SETUP_PDE_parms.f on the master ( MyProc=1 ) reads the PDE
C     paramters from a data file and then broadcasts the arrays.
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

      INCLUDE 'common2.h'

C ----------------------------
C     Argument Declarations
C ----------------------------

      INTEGER    EX_MPI_COMM, EX_MPI_IERROR, MPI_MyProc, NProc
      CHARACTER  PDEFILEi*(*)

C ----------------------------
C     Local Declarations
C ----------------------------

      INTEGER   NTMP_iBUFFER, NTMP_rBUFFER
      PARAMETER ( NTMP_iBUFFER = 5, NTMP_rBUFFER = 8*iREGmax + 4 )

      INTEGER   TMP_iBUFFER(NTMP_iBUFFER)
      REAL*8    TMP_rBUFFER(NTMP_rBUFFER)

      INTEGER   CE, i, IO_in, MyProc, NCHK_iBUFFER, NCHK_rBUFFER,
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

C ------------------------------------------- 
C     MPI indexing:
C ------------------------------------------- 

      MyProc = MPI_MyProc

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

         IF ( MyProc.EQ.0 ) THEN
            WRITE(*,*) 'ERROR: EX_SETUP_PDE_parms.f ... '
            WRITE(*,*) '  HAVE: NCHK_iBUFFER  = ', NCHK_iBUFFER
            WRITE(*,*) '  HAVE: NTMP_iBUFFER  = ', NTMP_iBUFFER
            WRITE(*,*) '  HAVE: NCHK_rBUFFER  = ', NCHK_rBUFFER
            WRITE(*,*) '  HAVE: NTMP_rBUFFER  = ', NTMP_rBUFFER
         ENDIF

         CALL MPI_FINALIZE(EX_MPI_IERROR)
         STOP

      ENDIF
      
C -----------------------------------
C     Initialize arrays in common
C -----------------------------------

      IF ( MyProc.EQ.0 ) THEN

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

      ENDIF

C -----------------------------------
C     Read parameter data
C -----------------------------------

      IF ( MyProc.EQ.0 ) THEN 

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

      ENDIF

C ------------------------------------------- 
C     Broadcast read status
C ------------------------------------------- 

      CALL MPI_Bcast( READ_STATUS, 1, MPI_INTEGER,
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

      IF ( READ_STATUS.NE.1 ) THEN
         IF ( MyProc.EQ.0 ) THEN
            WRITE(*,*) 'ERROR: File ', PDEFILEi(1:CE),
     &                 ' does not exist!'
         ENDIF
         CALL MPI_FINALIZE(EX_MPI_IERROR)
         STOP
      ENDIF

C ---------------------------------------
C     Pack Buffer:
C ---------------------------------------
      
      IF ( MyProc.EQ.0 ) THEN

         DO i=1, iREGmax
            !
            ! Diffusion coefficient (diagonal tensor)
            !
            TMP_rBUFFER(p_Dix+i-1) = Dix(i)
            TMP_rBUFFER(p_Diy+i-1) = Diy(i)
            !
            ! Removal/Absorption coefficent
            !
            TMP_rBUFFER(p_Si+i-1)  = Si(i)
            !
            ! Source/RHS 
            !
            TMP_rBUFFER(p_Fi+i-1)  = Fi(i)
            !
            ! Region coordinates (tensor product?)
            !
            TMP_rBUFFER(p_X1+i-1)  = X1(i)
            TMP_rBUFFER(p_X2+i-1)  = X2(i)
            TMP_rBUFFER(p_Y1+i-1)  = Y1(i)
            TMP_rBUFFER(p_Y2+i-1)  = Y2(i)
            
         ENDDO

         !
         ! Dimensions of the global domain
         !
         TMP_rBUFFER(p_xGs)  = xGs
         TMP_rBUFFER(p_xGf)  = xGf
         TMP_rBUFFER(p_yGs)  = yGs
         TMP_rBUFFER(p_yGf)  = yGf

         !
         ! Number of regions
         !
         TMP_iBUFFER(p_iREG) = iREG

         !
         ! Boundary condition indices
         !
         TMP_iBUFFER(p_ibl) = ibl
         TMP_iBUFFER(p_ibr) = ibr
         TMP_iBUFFER(p_ibb) = ibb
         TMP_iBUFFER(p_ibt) = ibt

      ENDIF

C ---------------------------------------
C     Broadcast:
C ---------------------------------------

      CALL MPI_Bcast( TMP_rBUFFER, NTMP_rBUFFER, MPI_REAL8,
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

      CALL MPI_Bcast( TMP_iBUFFER, NTMP_iBUFFER, MPI_INTEGER,
     &                0, EX_MPI_COMM, EX_MPI_IERROR )

C ---------------------------------------
C     Unpack BUFFER:
C ---------------------------------------

      IF ( MyProc.NE.0 ) THEN

         DO i=1, iREGmax
            !
            ! Diffusion coefficient (diagonal tensor)
            !
            Dix(i) = TMP_rBUFFER(p_Dix+i-1)
            Diy(i) = TMP_rBUFFER(p_Diy+i-1)
            !
            ! Removal/Absorption coefficent
            !
            Si(i)  = TMP_rBUFFER(p_Si+i-1)
            !
            ! Source/RHS 
            !
            Fi(i)  = TMP_rBUFFER(p_Fi+i-1)
            !
            ! Region coordinates (tensor product)
            !
            X1(i)  = TMP_rBUFFER(p_X1+i-1)
            X2(i)  = TMP_rBUFFER(p_X2+i-1)
            Y1(i)  = TMP_rBUFFER(p_Y1+i-1)
            Y2(i)  = TMP_rBUFFER(p_Y2+i-1)

         ENDDO

         !
         ! Dimensions of the global domain
         !
         xGs = TMP_rBUFFER(p_xGs)
         xGf = TMP_rBUFFER(p_xGf)
         yGs = TMP_rBUFFER(p_yGs)
         yGf = TMP_rBUFFER(p_yGf)

         !
         ! Number of regions
         !
         iREG = TMP_iBUFFER(p_iREG)

         !
         ! Boundary condition indices
         !
         iBL = TMP_iBUFFER(p_ibl)
         iBR = TMP_iBUFFER(p_ibr)
         iBB = TMP_iBUFFER(p_ibb)
         iBT = TMP_iBUFFER(p_ibt)

      ENDIF

      !
      !  Output the PDE parameters
      !
      IF ( MyProc.EQ.0 ) THEN 
         
         WRITE(*,1000) 'The number of regions is ', iREG
         WRITE(*,1020) 'Global Physical Domain:',
     &                    '[',xGs,',',xGf,'] x [',
     &                        yGs,',',yGf,']' 
         DO i=1, iREG

            WRITE(*,1010) 'Parameters for region ', i
            WRITE(*,1020) 'Physical Domain:',
     &                    '[',x1(i),',',x2(i),'] x [',
     &                        y1(i),',',y2(i),']' 
            WRITE(*,1030) '(Diagonal) Diffusion Tensor:',
     &                    Dix(i), Diy(i)
            WRITE(*,1040) 'Absorption Coefficient:', si(i)
            WRITE(*,1040) 'Source Strength:', fi(i)

         ENDDO

         WRITE(*,1050) 'Boundary Condition Indices'
         WRITE(*,1060) 'x=Xs', ibl
         WRITE(*,1060) 'x=Xf', ibr
         WRITE(*,1060) 'y=Ys', ibb
         WRITE(*,1060) 'y=Yf', ibt

      ENDIF
         
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
