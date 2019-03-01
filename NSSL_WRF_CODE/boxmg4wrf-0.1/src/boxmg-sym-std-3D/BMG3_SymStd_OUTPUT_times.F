      SUBROUTINE BMG3_SymStd_OUTPUT_times(
     &                BMG_rPARMS, MyProc, NProc, MPICOMM
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     Compute and output timing information.
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
#if RKIND == 4
#define FLOAT_PRECISION MPI_REAL
#else
#define FLOAT_PRECISION MPI_DOUBLE_PRECISION
#endif

      IMPLICIT NONE
      
      INCLUDE 'mpif.h'
#include    "BMG_constants.h"
      INCLUDE 'BMG_parameters.h'

      INTEGER MyProc, NProc, MPICOMM
      REAL*RKIND  BMG_rPARMS(NBMG_rPARMS)

      INTEGER I, ierr
      REAL*RKIND  AVG_rPARMS(NBMG_rPARMS)
      REAL*RKIND  MAX_rPARMS(NBMG_rPARMS)
      REAL*RKIND  MIN_rPARMS(NBMG_rPARMS)

C ==========================================================================

      DO I=NBMG2_rPARMS+2, NBMG_rPARMS
         CALL MPI_AllReduce(BMG_rPARMS(i),AVG_rPARMS(i),
     &        iONE, FLOAT_PRECISION, MPI_SUM, 
     &        MPICOMM, ierr)
         AVG_rPARMS(i) = AVG_rPARMS(i)/NProc
      ENDDO

      DO I=NBMG2_rPARMS+2, NBMG_rPARMS
         CALL MPI_AllReduce(BMG_rPARMS(i),MAX_rPARMS(i),
     &        iONE, FLOAT_PRECISION, MPI_MAX,
     &        MPICOMM, ierr)
      ENDDO

      DO I=NBMG2_rPARMS+2, NBMG_rPARMS
         CALL MPI_AllReduce(BMG_rPARMS(i),MIN_rPARMS(i),
     &        iONE, FLOAT_PRECISION, MPI_MIN,
     &        MPICOMM, ierr)
      ENDDO

      IF ( MyProc.EQ.iZERO ) THEN

         WRITE(*,*)
         WRITE(*,*) 'Average Times:'
         WRITE(*,*) '=============='
         
         WRITE(*,110)  AVG_rPARMS(id_BMG3_TIME_SETUP_TOTAL),
     &        AVG_rPARMS(id_BMG3_TIME_SOLVE_TOTAL)
         
         WRITE(*,111) AVG_rPARMS(id_BMG3_TIME_SETUP_FINE_STENCIL),
     &         AVG_rPARMS(id_BMG3_TIME_SETUP_CG_ITLI)

         WRITE(*,112) AVG_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI),
     &        AVG_rPARMS(id_BMG3_TIME_SETUP_RELAX)

         WRITE(*,113) AVG_rPARMS(id_BMG3_TIME_SETUP_CG_LU),
     &        AVG_rPARMS(id_BMG3_TIME_SOLVE_CG)

         WRITE(*,114) AVG_rPARMS(id_BMG3_TIME_relax),
     &        AVG_rPARMS(id_BMG3_TIME_restrict)

         WRITE(*,115) AVG_rPARMS(id_BMG3_TIME_interp_add),
     &        AVG_rPARMS(id_BMG3_TIME_SETUP_MSG)
         
         WRITE(*,116) AVG_rPARMS(id_BMG3_TIME_SETUP_PTR_GRID),
     &        AVG_rPARMS(id_BMG3_TIME_SETUP_PARTS)

         WRITE(*,117) AVG_rPARMS(id_BMG3_TIME_PCG_TOTAL),
     &        AVG_rPARMS(id_BMG3_TIME_PCG_PRECON)

         WRITE(*,*)
         WRITE(*,*) 'Max Times:'
         WRITE(*,*) '=========='
         
         WRITE(*,110)  MAX_rPARMS(id_BMG3_TIME_SETUP_TOTAL),
     &        MAX_rPARMS(id_BMG3_TIME_SOLVE_TOTAL)
         
         WRITE(*,111) MAX_rPARMS(id_BMG3_TIME_SETUP_FINE_STENCIL),
     &         MAX_rPARMS(id_BMG3_TIME_SETUP_CG_ITLI)

         WRITE(*,112) MAX_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI),
     &        MAX_rPARMS(id_BMG3_TIME_SETUP_RELAX)

         WRITE(*,113) MAX_rPARMS(id_BMG3_TIME_SETUP_CG_LU),
     &        MAX_rPARMS(id_BMG3_TIME_SOLVE_CG)

         WRITE(*,114) MAX_rPARMS(id_BMG3_TIME_relax),
     &        MAX_rPARMS(id_BMG3_TIME_restrict)

         WRITE(*,115) MAX_rPARMS(id_BMG3_TIME_interp_add),
     &        MAX_rPARMS(id_BMG3_TIME_SETUP_MSG)
         
         WRITE(*,116) MAX_rPARMS(id_BMG3_TIME_SETUP_PTR_GRID),
     &        MAX_rPARMS(id_BMG3_TIME_SETUP_PARTS)

         WRITE(*,117) MAX_rPARMS(id_BMG3_TIME_PCG_TOTAL),
     &        MAX_rPARMS(id_BMG3_TIME_PCG_PRECON)

         WRITE(*,*)
         WRITE(*,*) 'Min Times:'
         WRITE(*,*) '=========='
         
         WRITE(*,110)  MIN_rPARMS(id_BMG3_TIME_SETUP_TOTAL),
     &        MIN_rPARMS(id_BMG3_TIME_SOLVE_TOTAL)
         
         WRITE(*,111) MIN_rPARMS(id_BMG3_TIME_SETUP_FINE_STENCIL),
     &         MIN_rPARMS(id_BMG3_TIME_SETUP_CG_ITLI)

         WRITE(*,112) MIN_rPARMS(id_BMG3_TIME_SETUP_INTERP_OI),
     &        MIN_rPARMS(id_BMG3_TIME_SETUP_RELAX)

         WRITE(*,113) MIN_rPARMS(id_BMG3_TIME_SETUP_CG_LU),
     &        MIN_rPARMS(id_BMG3_TIME_SOLVE_CG)

         WRITE(*,114) MIN_rPARMS(id_BMG3_TIME_relax),
     &        MIN_rPARMS(id_BMG3_TIME_restrict)

         WRITE(*,115) MIN_rPARMS(id_BMG3_TIME_interp_add),
     &        MIN_rPARMS(id_BMG3_TIME_SETUP_MSG)
         
         WRITE(*,116) MIN_rPARMS(id_BMG3_TIME_SETUP_PTR_GRID),
     &        MIN_rPARMS(id_BMG3_TIME_SETUP_PARTS)

         WRITE(*,117) MIN_rPARMS(id_BMG3_TIME_PCG_TOTAL),
     &        MIN_rPARMS(id_BMG3_TIME_PCG_PRECON)

         WRITE(*,*)

      ENDIF
      
C ==========================================================================

 110  FORMAT(4X,'Setup (total)       = ',F16.4,3X,
     &     'Solve (total) = ',F16.4)
 111  FORMAT(4X,'SETUP_fine_stencil  = ',F16.4,3X,
     &     'SETUP_cg_ITLI = ',F16.4)
 112  FORMAT(4X,'SETUP_interp_IO     = ',F16.4,3X,
     &     'SETUP_relax   = ',F16.4)
 113  FORMAT(4X,'SETUP_cg_LU         = ',F16.4,3X,
     &     'SOLVE_cg      = ',F16.4)
 114  FORMAT(4X,'relax               = ',F16.4,3X,
     &     'restrict      = ',F16.4)
 115  FORMAT(4X,'interp_add          = ',F16.4,3X,
     &     'SETUP_MSG     = ',F16.4)
 116  FORMAT(4X,'SETUP_PtrGrid       = ',F16.4,3X,
     &     'SETUP_parts   = ',F16.4)
 117  FORMAT(4X,'PCG (total)         = ',F16.4,3X,
     &     'PCG_precon    = ',F16.4)
     
C ===========================================

      RETURN
      END


