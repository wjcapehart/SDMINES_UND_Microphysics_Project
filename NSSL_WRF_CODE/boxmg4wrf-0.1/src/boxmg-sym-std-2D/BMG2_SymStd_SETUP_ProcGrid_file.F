      SUBROUTINE BMG2_SymStd_SETUP_ProcGrid_file(
     &                GRIDFILEi, BMG_MSG_iGRID, BMG_MSG_pGRID,
     &                BMG_iPARMS, NBMG_MSG_iGRIDm, NBMG_MSG_iGRID,
     &                MPI_MyProc, NProc, BMG_MPI_COMM, MPI_IERROR
     &                ) 

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG2_SymStd_SETUP_ProcGrid_file reads and communicates minimal
C     data from a file that describes the processor grid and creates
C     related information for use locally.
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

      INCLUDE 'BMG_parameters.h'
      INCLUDE 'BMG_workspace.h'

C -----------------------------
C     Argument Declarations
C
      INTEGER  BMG_iPARMS(NBMG_iPARMS)
      
      INTEGER  NBMG_MSG_iGRIDm, NBMG_MSG_iGRID
      INTEGER  BMG_MSG_iGRID(NBMG_MSG_iGRIDm),
     &         BMG_MSG_pGRID(NBMG_MSG_pGRID)

      INTEGER   BMG_MPI_COMM, MPI_MyProc, MPI_IERROR, NProc
      CHARACTER GRIDFILEi*(*)

C -----------------------------
C     Local Declarations
C 
      INTEGER  GRIDDIST, i, iGs, j, jGs, MyProc, MyProcI, MyProcJ,
     &         NGx, NGy, NLx, NLy, NProcI, NProcJ, pSI,
     &         p_iGs, p_jGs, p_NLx, p_NLy, p_PC, p_PG,
     &         p_ProcCoord, p_ProcGrid
      
C ==========================================================================

C ------------------------------------------- 
C     MPI -> MSG indexing:
C ------------------------------------------- 

      MyProc = MPI_MyProc+1

C ------------------------------------------- 
C     Open the data file, read the header
C -------------------------------------------     

      !
      !  Master:
      !
      IF ( MyProc .EQ. 1 ) THEN

         OPEN(10,FILE=GRIDFILEi, STATUS='OLD' )
         READ(10,*) GRIDDIST

      ENDIF

      !
      !  All:
      !
      CALL MPI_Bcast( GRIDDIST, 1, MPI_INTEGER,
     &                0, BMG_MPI_COMM, MPI_IERROR  )


C --------------------------------------------
C     Process and distribute the grid data
C --------------------------------------------

      IF ( GRIDDIST.EQ.BMG_MSG_GRIDDIST_AUTOMATIC ) THEN

         
      
      ELSEIF ( GRIDDIST.EQ.BMG_MSG_GRIDDIST_CUSTOM
     &        .OR. GRIDDIST.EQ.BMG_MSG_GRIDDIST_UNIFORM ) THEN
         !
         !  Master:
         !
         IF ( MyProc.EQ.1 ) THEN
            ! 
            ! Read the number of processors in i and j
            !
            READ(10,*) NProcI, NProcJ
            !
         ENDIF

         !
         !  All:
         !
         CALL MPI_Bcast( NProcI, 1, MPI_INTEGER, 
     &                   0, BMG_MPI_COMM, MPI_IERROR )
         
         CALL MPI_BCast( NProcJ, 1, MPI_INTEGER, 
     &                   0, BMG_MPI_COMM, MPI_IERROR )
         
         IF ( NProcI*NProcJ .NE. NProc ) THEN
            WRITE(*,*) 'Procgrid Error : ',NProcI,' x ',NProcJ
     &           ,' != ',NProc
            STOP
         ENDIF

         !
         ! Compute storage requirements and pointers
         !
         pSI = NBMG_MSG_iGRID_data + 1
         !
         BMG_MSG_pGRID(ip_BMG_MSG_NLx_Grid)  = pSI
         pSI = pSI + NProcI
         !
         BMG_MSG_pGRID(ip_BMG_MSG_NLy_Grid)  = pSI
         pSI = pSI + NProcJ
         !
         BMG_MSG_pGRID(ip_BMG_MSG_iGs_Grid)  = pSI
         pSI = pSI + NProcI
         !
         BMG_MSG_pGRID(ip_BMG_MSG_jGs_Grid)  = pSI
         pSI = pSI + NProcJ
         !
         BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)  = pSI
         pSI = pSI + NProcI*NProcJ
         !
         BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord) = pSI
         pSI = pSI + 2*NProcI*NProcJ
         !
         NBMG_MSG_iGRID = pSI-1

         IF ( BMG_iPARMS(id_BMG2_POINTERS).EQ.BMG_USE_pointers ) THEN
            
            NBMG_MSG_iGRIDm = NBMG_MSG_iGRID

            CLOSE(10)
            RETURN

         ELSEIF ( BMG_iPARMS(id_BMG2_POINTERS).EQ.BMG_NO_pointers ) THEN
            
            IF ( NBMG_MSG_iGRID.GT.NBMG_MSG_iGRIDm ) THEN 

               IF (MyProc.EQ.1) THEN
                  WRITE(*,500) 'Space for BMG_MSG_iGRID is too small!'
                  WRITE(*,510) 'HAVE: NBMG_MSG_iGRIDm = ',
     &                         NBMG_MSG_iGRIDm
                  WRITE(*,520) 'NEED: NBMG_MSG_iGRIDm = ',
     &                         NBMG_MSG_iGRID
               ENDIF

               CALL MPI_FINALIZE(MPI_IERROR)
               STOP

            ENDIF

         ELSE

            IF (MyProc.EQ.1) THEN
               WRITE(*,500) 'Value of BMG_iPARMS(id_BMG2_pointers) is '
     &                      //'out of range ...'
               WRITE(*,510) 'HAVE: BMG_iPARMS(id_BMG2_POINTERS) = ',
     &                      BMG_iPARMS(id_BMG2_POINTERS)
               WRITE(*,505) 'NEED: BMG_iPARMS(id_BMG2_POINTERS) = '
     &                      //'0 or 1'  
            ENDIF

            CALL MPI_FINALIZE(MPI_IERROR)
            STOP

         ENDIF

         !
         !  Make local copies of pointers
         !
         p_NLx = BMG_MSG_pGRID(ip_BMG_MSG_NLx_Grid)
         p_NLy = BMG_MSG_pGRID(ip_BMG_MSG_NLy_Grid)
         p_iGs = BMG_MSG_pGRID(ip_BMG_MSG_iGs_Grid)
         p_jGs = BMG_MSG_pGRID(ip_BMG_MSG_jGs_Grid)
         !
         p_ProcGrid  = BMG_MSG_pGRID(ip_BMG_MSG_ProcGrid)
         p_ProcCoord = BMG_MSG_pGRID(ip_BMG_MSG_ProcCoord)

         !
         !  Master:
         !
         IF ( MyProc.EQ.1 ) THEN
            !
            IF ( GRIDDIST.EQ.BMG_MSG_GRIDDIST_CUSTOM ) THEN
               !
               ! Read arrays of fine-grid local dimensions 
               !
               DO i=1, NProcI
                  READ(10,*) BMG_MSG_iGRID(p_NLx+i-1)
               ENDDO
               !
               DO j=1, NProcJ
                  READ(10,*) BMG_MSG_iGRID(p_NLy+j-1)
               ENDDO
               !
               ! Close the grid file.
               !
               CLOSE (10)
            ELSE
               !
               ! Read fine-grid local dimensions
               !
               READ(10,*) NLx, NLy
               !
               ! Close the grid file.
               !
               CLOSE (10)
               !
               ! Fill arrays of fine-grid local dimensions
               !
               DO i=1, NProcI
                  BMG_MSG_iGRID(p_NLx+i-1) = NLx
               ENDDO
               !
               DO j=1, NProcJ
                  BMG_MSG_iGRID(p_NLy+j-1) = NLy
               ENDDO
               !
            ENDIF

         ENDIF
         
         !
         !  All:
         !
         CALL MPI_Bcast( BMG_MSG_iGRID(p_NLx), NProcI, MPI_INTEGER,
     &                   0, BMG_MPI_COMM, MPI_IERROR )
         
         CALL MPI_BCast( BMG_MSG_iGRID(p_NLy), NProcJ, MPI_INTEGER,
     &                   0, BMG_MPI_COMM, MPI_IERROR ) 
         
         !
         ! Setup processor grid: 
         !
         !   MyProc  range is [1,NProcI*NProcJ]
         !   MyProcI range is [1,NProcI]
         !   MyProcJ range is [1,NProcJ] 
         !
         MyProcJ = (MyProc-1)/NProcI + 1
         MyProcI = MyProc - (MyProcJ-1) * NProcI
         
         p_PG = p_ProcGrid
         DO j=1,NProcJ
            DO i=1,NProcI
               BMG_MSG_iGRID(p_PG) = (j-1)*NProcI + i 
               p_PG = p_PG + 1
            ENDDO
         ENDDO
         
         !
         ! Store coordinates of each processor
         ! 
         p_PC = p_ProcCoord
         DO j=1,NProcJ
            DO i=1,NProcI
               BMG_MSG_iGRID(p_PC)   = i
               BMG_MSG_iGRID(p_PC+1) = j
               p_PC = p_PC + 2
            ENDDO
         ENDDO

         !
         !  Global problem size
         ! 
         NGx = 0
         do i=1, NProcI
            NGx = NGx + BMG_MSG_iGRID(p_NLx+i-1)
         enddo
         
         NGy = 0
         do i=1, NProcJ
            NGy = NGy + BMG_MSG_iGRID(p_NLy+i-1)
         enddo
         
         !
         ! Global index (i,j) of the lower left corner
         !
         iGs = 1
         DO i=1, NprocI
            BMG_MSG_iGRID(p_iGs+i-1) = iGs
            iGs = iGs + BMG_MSG_iGRID(p_NLx+i-1)
         ENDDO

         jGs = 1
         DO i=1, NProcJ
            BMG_MSG_iGRID(p_jGs+i-1) = jGs
            jGs = jGs + BMG_MSG_iGRID(p_NLy+i-1)
         ENDDO
         
         !
         !  Sanity check
         !
         IF ( NGx.NE.iGs-1 .OR. NGy.NE.jGs-1 ) THEN

            IF ( MyProc.EQ.1 ) THEN
               WRITE(*,*) ' THIS IS STRANGE: NGx = ', NGx
               WRITE(*,*) '                  iGs = ', iGs
               WRITE(*,*) '                  NGy = ', NGy
               WRITE(*,*) '                  jGs = ', jGs
            ENDIF
               
            CALL MPI_FINALIZE(MPI_IERROR)
            STOP

         ENDIF

      ENDIF

C --------------------------------------------
C     Packup the local grid data
C --------------------------------------------

      ! 
      !  Store: NGx, NGy, NLx, NLy, iGs, jGs
      ! 
      BMG_MSG_iGRID(id_BMG_MSG_NGx) = NGx
      BMG_MSG_iGRID(id_BMG_MSG_NGy) = NGy
      BMG_MSG_iGRID(id_BMG_MSG_NLx) = BMG_MSG_iGRID(p_NLx+MyProcI-1)
      BMG_MSG_iGRID(id_BMG_MSG_NLy) = BMG_MSG_iGRID(p_NLy+MyProcJ-1)
      BMG_MSG_iGRID(id_BMG_MSG_iGs) = BMG_MSG_iGRID(p_iGs+MyProcI-1)
      BMG_MSG_iGRID(id_BMG_MSG_jGs) = BMG_MSG_iGRID(p_jGs+MyProcJ-1)
      !
      !  Store: Nproc, NProcI, NProcJ
      !
      BMG_MSG_iGRID(id_BMG_MSG_NProc)   = NProc
      BMG_MSG_iGRID(id_BMG_MSG_NProcI)  = NProcI
      BMG_MSG_iGRID(id_BMG_MSG_NProcJ)  = NProcJ
      !
      !  Store: MyProc, MyProcI, MyProcJ
      !
      BMG_MSG_iGRID(id_BMG_MSG_MyProc)  = MyProc
      BMG_MSG_iGRID(id_BMG_MSG_MyProcI) = MyProcI
      BMG_MSG_iGRID(id_BMG_MSG_MyProcJ) = MyProcJ
      !
      !  Store/set communicators:
      !
      BMG_MSG_iGRID(id_BMG_MSG_COMM)   = BMG_MPI_COMM
      BMG_MSG_iGRID(id_BMG_MSG_COMM_x) = 0
      BMG_MSG_iGRID(id_BMG_MSG_COMM_y) = 0

C ==========================================================================

 500  FORMAT (/,'FATAL ERROR: BMG2_SymStd_SETUP_ProcGrid_file.f',
     &        //,5X,A,/)
 505  FORMAT (5X,A,/)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

C ===========================================

      RETURN
      END
