      SUBROUTINE EX_SETUP_PDE_parms(
     &                    PDEFILEi, 
     &                    MPI_MyProc, NProc, EX_MPI_COMM, EX_MPI_IERROR
     &                    )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C     
C     EX_SETUP_PDE_parms on the master ( MyProc=1 ) reads the PDE
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
      PARAMETER ( NTMP_iBUFFER = 7, NTMP_rBUFFER = 11*iREGmax + 6 )

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
      INTEGER   p_X1, p_X2, p_Y1, p_Y2, p_Z1, p_Z2

      !
      !  Region coefficient pointers
      ! 
      INTEGER   p_Dix, p_Diy, p_Diz, p_Si, p_Fi

      !
      !  Global domain boundary 
      ! 
      INTEGER   p_xGf, p_xGs, p_yGf, p_yGs, p_zGf, p_zGs

      !
      !  Boundary condition index pointers
      !
      INTEGER   p_ibl, p_ibr, p_ibyb, p_ibyt, p_ibzb, p_ibzt

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
      p_ibyb = p_ibr + 1
      p_ibyt = p_ibyb + 1
      p_ibzb = p_ibyt + 1
      p_ibzt = p_ibzb + 1

      NCHK_iBUFFER = p_ibzt

      !
      !  PDE coefficients
      !
      p_Dix  = 1
      p_Diy  = p_Dix + iREGmax 
      p_Diz  = p_Diy + iREGmax

      p_Si   = p_Diz + iREGmax
      p_Fi   = p_Si  + iREGmax

      !
      !  Region boundaries
      !
      p_X1   = p_Fi  + iREGmax
      p_X2   = p_X1  + iREGmax 
      p_Y1   = p_X2  + iREGmax
      p_Y2   = p_Y1  + iREGmax
      p_Z1   = p_Y2  + iREGmax
      p_Z2   = p_Z1  + iREGmax

      !
      !  Global domain boundaries
      !
      p_xGs   = p_Z2 + iREGmax
      p_xGf   = p_xGs + 1
      p_yGs   = p_xGf + 1
      p_yGf   = p_yGs + 1
      p_zGs   = p_yGf + 1
      p_zGf   = p_zGs + 1

      NCHK_rBUFFER = p_zGf

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
            Diz(i) = rZERO
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
            Z1(i)  = rZERO
            Z2(i)  = rZERO
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
            READ(IO_in,*) xGs, xGf, yGs, yGf, zGs, zGf
            !
            !  Global boundary condition indices
            !
            READ(IO_in,*) ibl,ibr,ibyb,ibyt,ibzb,ibzt

            !
            !  The number of regions:
            !
            READ(IO_in,*) iREG
            
            !
            !  Read data for each region:
            !      
            DO i=1, iREG

               READ(IO_in,*) x1(i),x2(i),y1(i),y2(i),z1(i),z2(i)
               READ(IO_in,*) dix(i),diy(i),diz(i),si(i),fi(i)

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
            TMP_rBUFFER(p_Diz+i-1) = Diz(i)
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
            TMP_rBUFFER(p_Z1+i-1)  = Z1(i)
            TMP_rBUFFER(p_Z2+i-1)  = Z2(i)
            
         ENDDO

         !
         ! Dimensions of the global domain
         !
         TMP_rBUFFER(p_xGs)  = xGs
         TMP_rBUFFER(p_xGf)  = xGf
         TMP_rBUFFER(p_yGs)  = yGs
         TMP_rBUFFER(p_yGf)  = yGf
         TMP_rBUFFER(p_zGs)  = zGs
         TMP_rBUFFER(p_zGf)  = zGf

         !
         ! Number of regions
         !
         TMP_iBUFFER(p_iREG) = iREG

         !
         ! Boundary condition indices
         !
         TMP_iBUFFER(p_ibl)  = ibl
         TMP_iBUFFER(p_ibr)  = ibr
         TMP_iBUFFER(p_ibyb) = ibyb
         TMP_iBUFFER(p_ibyt) = ibyt
         TMP_iBUFFER(p_ibzb) = ibzb
         TMP_iBUFFER(p_ibzt) = ibzt

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
            Diz(i) = TMP_rBUFFER(p_Diz+i-1)
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
            Z1(i)  = TMP_rBUFFER(p_Z1+i-1)
            Z2(i)  = TMP_rBUFFER(p_Z2+i-1)

         ENDDO

         !
         ! Dimensions of the global domain
         !
         xGs = TMP_rBUFFER(p_xGs)
         xGf = TMP_rBUFFER(p_xGf)
         yGs = TMP_rBUFFER(p_yGs)
         yGf = TMP_rBUFFER(p_yGf)
         zGs = TMP_rBUFFER(p_zGs)
         zGf = TMP_rBUFFER(p_zGf)

         !
         ! Number of regions
         !
         iREG = TMP_iBUFFER(p_iREG)

         !
         ! Boundary condition indices
         !
         ibl  = TMP_iBUFFER(p_ibl)
         ibr  = TMP_iBUFFER(p_ibr)
         ibyb = TMP_iBUFFER(p_ibyb)
         ibyt = TMP_iBUFFER(p_ibyt)
         ibzb = TMP_iBUFFER(p_ibzb)
         ibzt = TMP_iBUFFER(p_ibzt)

      ENDIF

      !
      !  Output the PDE parameters
      !
      IF ( MyProc.EQ.0 ) THEN 
         
         WRITE(*,1000) 'The number of regions is ', iREG
         WRITE(*,1020) 'Global Physical Domain:',
     &                    '[',xGs,',',xGf,'] x [',
     &                        yGs,',',yGf,'] x [',
     &                        zGs,',',zGf,']' 
         DO i=1, iREG

            WRITE(*,1010) 'Parameters for region ', i
            WRITE(*,1020) 'Physical Domain:',
     &                    '[',x1(i),',',x2(i),'] x [',
     &                        y1(i),',',y2(i),'] x [',
     &                        z1(i),',',z2(i),']' 
            WRITE(*,1030) '(Diagonal) Diffusion Tensor:',
     &                    Dix(i), Diy(i), Diz(i)
            WRITE(*,1040) 'Absorption Coefficient:', si(i)
            WRITE(*,1040) 'Source Strength:', fi(i)

         ENDDO

         WRITE(*,1050) 'Boundary Condition Indices'
         WRITE(*,1060) 'x=Xs', ibl
         WRITE(*,1060) 'x=Xf', ibr
         WRITE(*,1060) 'y=Ys', ibyb
         WRITE(*,1060) 'y=Yf', ibyt
         WRITE(*,1060) 'z=Zs', ibzb
         WRITE(*,1070) 'z=Zf', ibzt

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
 1020 FORMAT(2X,A,T35,A,6(F5.2,A),/)
 1030 FORMAT(2X,A,T35,1P,3(E14.7))
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
