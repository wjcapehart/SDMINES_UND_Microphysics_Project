      SUBROUTINE BMG3_SymStd_relax_planes_yz( 
     &                BMG_iPARMS, BMG_rPARMS, BMG_IOFLAG,
     &                Nxf, Nyf, Nzf, kg, IGRD, Q, QF, RES, NF,
     &                SO, NSO, SOR, NSOR, CI, NCI, IFD,
     &                NOGm, NOG, UPDOWN, RES_L2,
     &                BMG_iPARMS_PL_yz, BMG_pWORK_PL_yz,
     &                BMG_iWORK_PL, NBMG_iWORK_PL,
     &                BMG_rWORK_PL, NBMG_rWORK_PL
     &                )

C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     BMG3_SymStd_relax_planes_yz performs one relaxation step of 
C     red-black yz-plane relaxation (in the correct order)
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

C ------------------------------------------------
C     Includes
C
#include    "BMG_constants.h"
      INCLUDE 'BMG_stencils.h'
      INCLUDE 'BMG_workspace.h'
      INCLUDE 'BMG_parameters.h'

C ------------------------------------------------
C     Argument Declarations
C 
      INTEGER  kg, NCI, NBMG_iWORK_PL, NBMG_rWORK_PL,
     &         NF, NOG, NOGm, NSO, NSOR, Nxf, Nyf, Nzf

      INTEGER  BMG_iPARMS(NBMG_iPARMS), 
     &         BMG_iPARMS_PL_yz(NBMG_iPARMS,NOG),
     &         BMG_iWORK_PL(NBMG_iWORK_PL),
     &         BMG_pWORK_PL_yz(NBMG_pWORK,Nxf,NOG),
     &         IFD, IGRD(NOGm,24), IRELAX_SYM, UPDOWN

      REAL*RKIND   BMG_rPARMS(NBMG_rPARMS), BMG_rWORK_PL(NBMG_rWORK_PL), 
     &         CI(NCI), Q(NF), QF(NF), RES(NF),   
     &         SO(NSO), SOR(NSOR)
      
      REAL*8   RES_L2
      
      LOGICAL  BMG_IOFLAG(NBMG_IOFLAG)

C -------------------------------------------------
C    Local Declarations
C
      INTEGER  iPL, iPL_beg, iPL_end, Nx, Ny, Nz,
     &         LSTART, LEND, LSTRIDE,
     &         p_CIC, p_SO, p_SOR, p_U

      INTEGER  ngx,ngy,ngz

      INTEGER  NFm_yz, NOGm_yz, NSOm_yz

      INTEGER  NC_yz, NCBW_yz, NCI_yz, NCU_yz, NF_yz, NOG_yz,
     &         NSO_yz, NSOR_yz, NStncl_2D, NStncl_3D 

      INTEGER  p_CI_yz, p_CSO_yz, p_CU_yz, p_IG_yz, p_SO_yz, p_SOR_yz,
     &         p_Q_yz, p_RES_yz, p_U_yz

C ==========================================================================
      
      IRELAX_SYM = BMG_iPARMS_PL_yz(id_BMG2_RELAX_SYM,kg)

      IF ( ( IRELAX_SYM.EQ.BMG_RELAX_NONSYM )
     &   .OR. ( UPDOWN.EQ.BMG_DOWN )         ) THEN
         LSTART  =  2
         LEND    =  3
         LSTRIDE =  1
      ELSE
         LSTART  =  3
         LEND    =  2
         LSTRIDE = -1
      ENDIF
      
      !
      !  Set stencil size
      !
      IF (kg.NE.NOG .OR. IFD.NE.1 ) THEN
         NStncl_2D = 5
         NStncl_3D = 14
      ELSE
         NStncl_2D = 3
         NStncl_3D = 4
      ENDIF

      !
      ! Store maximum dimensions for (x,y) - planes
      ! 
      NOGm_yz = BMG_iPARMS_PL_yz(id_BMG2_DIM_NOG,NOG)
      NFm_yz  = BMG_iPARMS_PL_yz(id_BMG2_DIM_NF,NOG)
      NSOm_yz = BMG_iPARMS_PL_yz(id_BMG2_DIM_NSO,NOG)

c$$$      WRITE(*,*) 'Nxf, Nyf  = ', Nxf-2, Nyf-2
c$$$      WRITE(*,*) 'NOGm_yz   = ', NOGm_yz
c$$$      WRITE(*,*) 'NFm_yz    = ', NFm_yz
c$$$      WRITE(*,*) 'NSOm_yz   = ', NSOm_yz

C -------------------------------------------------------------------------
C    >>>>>>>>>>>>>>>>>>>>>> BEGIN: yz sweep <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C -------------------------------------------------------------------------

      CALL BMG3_SymStd_GET_pointers( kg, IGRD, NOGm,
     &            p_U, p_SO, p_SOR, p_CIC, Nx, Ny, Nz,ngx,ngy,ngz)

c$$$      WRITE(*,*) 'Relaxing on Grid kg = ', kg
c$$$      WRITE(*,*) 'Nx, Nz, Nz  = ', Nx, Ny, Nz

      NOG_yz  = BMG_iPARMS_PL_yz(id_BMG2_DIM_NOG,kg)
      NF_yz   = BMG_iPARMS_PL_yz(id_BMG2_DIM_NF,kg)
      NC_yz   = BMG_iPARMS_PL_yz(id_BMG2_DIM_NC,kg)
      NSO_yz  = BMG_iPARMS_PL_yz(id_BMG2_DIM_NSO,kg)
      NSOR_yz = BMG_iPARMS_PL_yz(id_BMG2_DIM_NSOR,kg)
      NCI_yz  = BMG_iPARMS_PL_yz(id_BMG2_DIM_NCI,kg)
      NCBW_yz = BMG_iPARMS_PL_yz(id_BMG2_DIM_NCBW,kg)
      NCU_yz  = BMG_iPARMS_PL_yz(id_BMG2_DIM_NCU,kg)
                                !

c$$$      WRITE(*,*) '**** Entering the loop with with ... '
c$$$      WRITE(*,*) 
c$$$      WRITE(*,*) ' NOG_yz  = ', BMG_iPARMS_PL_yz(id_BMG2_DIM_NOG,kg)
c$$$      WRITE(*,*) ' NF_yz   = ', BMG_iPARMS_PL_yz(id_BMG2_DIM_NF,kg)
c$$$      WRITE(*,*) ' NC_yz   = ', BMG_iPARMS_PL_yz(id_BMG2_DIM_NC,kg)
c$$$      WRITE(*,*) ' NSO_yz  = ', BMG_iPARMS_PL_yz(id_BMG2_DIM_NSO,kg)
c$$$      WRITE(*,*) ' NSOR_yz = ', BMG_iPARMS_PL_yz(id_BMG2_DIM_NSOR,kg)
c$$$      WRITE(*,*) ' NCI_yz  = ', BMG_iPARMS_PL_yz(id_BMG2_DIM_NCI,kg)
c$$$      WRITE(*,*) ' NCBW_yz = ', BMG_iPARMS_PL_yz(id_BMG2_DIM_NCBW,kg)
c$$$      WRITE(*,*) ' NCU_yz  = ', BMG_iPARMS_PL_yz(id_BMG2_DIM_NCU,kg)


      DO iPL_beg=LSTART, LEND, LSTRIDE
         !
         iPL_end=2*((Nx-1-iPL_beg)/2)+iPL_beg
         !
         DO iPL=iPL_beg, iPL_end, 2

            p_SO_yz  = BMG_pWORK_PL_yz(ip_SO,iPL,kg)
            p_U_yz   = BMG_pWORK_PL_yz(ip_U,iPL,kg)
            p_Q_yz   = BMG_pWORK_PL_yz(ip_Q,iPL,kg)
            p_RES_yz = BMG_pWORK_PL_yz(ip_RES,iPL,kg)
            p_SOR_yz = BMG_pWORK_PL_yz(ip_SOR,iPL,kg)
            p_CI_yz  = BMG_pWORK_PL_yz(ip_CI,iPL,kg)
            p_CSO_yz = BMG_pWORK_PL_yz(ip_CSO,iPL,kg)
            p_CU_yz  = BMG_pWORK_PL_yz(ip_CU,iPL,kg)
            p_IG_yz  = BMG_pWORK_PL_yz(ip_IG,iPL,kg)

             ! solve the plane
           
c$$$            WRITE(*,*) 
c$$$            WRITE(*,*) 'iPL = ', iPL
c$$$            WRITE(*,*) 
c$$$            WRITE(*,*) '**** Calling Boxmg with ... '
c$$$            WRITE(*,*) 
c$$$            WRITE(*,*) ' NOG_yz  = ',
c$$$     &                 BMG_iPARMS_PL_yz(id_BMG2_DIM_NOG,kg)
c$$$            WRITE(*,*) ' NF_yz   = ',
c$$$     &                 BMG_iPARMS_PL_yz(id_BMG2_DIM_NF,kg)
c$$$            WRITE(*,*) ' NC_yz   = ', 
c$$$     &                 BMG_iPARMS_PL_yz(id_BMG2_DIM_NC,kg)
c$$$            WRITE(*,*) ' NSO_yz  = ', 
c$$$     &                 BMG_iPARMS_PL_yz(id_BMG2_DIM_NSO,kg)
c$$$            WRITE(*,*) ' NSOR_yz = ', 
c$$$     &                 BMG_iPARMS_PL_yz(id_BMG2_DIM_NSOR,kg)
c$$$            WRITE(*,*) ' NCI_yz  = ',
c$$$     &                 BMG_iPARMS_PL_yz(id_BMG2_DIM_NCI,kg)
c$$$            WRITE(*,*) ' NCBW_yz = ', 
c$$$     &                 BMG_iPARMS_PL_yz(id_BMG2_DIM_NCBW,kg)
c$$$            WRITE(*,*) ' NCU_yz  = ', 
c$$$     &                 BMG_iPARMS_PL_yz(id_BMG2_DIM_NCU,kg)
c$$$            WRITE(*,*) 
c$$$            WRITE(*,*) ' p_SO_yz  = ', BMG_pWORK_PL_yz(ip_SO,iPL,kg)
c$$$            WRITE(*,*) ' p_U_yz   = ', BMG_pWORK_PL_yz(ip_U,iPL,kg)
c$$$            WRITE(*,*) ' p_Q_yz   = ', BMG_pWORK_PL_yz(ip_Q,iPL,kg)
c$$$            WRITE(*,*) ' p_RES_yz = ', BMG_pWORK_PL_yz(ip_RES,iPL,kg)
c$$$            WRITE(*,*) ' p_SOR_yz = ', BMG_pWORK_PL_yz(ip_SOR,iPL,kg)
c$$$            WRITE(*,*) ' p_CI_yz  = ', BMG_pWORK_PL_yz(ip_CI,iPL,kg)
c$$$            WRITE(*,*) ' p_CSO_yz = ', BMG_pWORK_PL_yz(ip_CSO,iPL,kg)
c$$$            WRITE(*,*) ' p_CU_yz  = ', BMG_pWORK_PL_yz(ip_CU,iPL,kg)
c$$$            WRITE(*,*) ' p_IG_yz  = ', BMG_pWORK_PL_yz(ip_IG,iPL,kg)

            !
            ! Copy the iPL{th} plane of the 3D solution into 2D
            !
            CALL BMG3_SymStd_COPY_rV_32_yz( 
     &                Q(p_U), BMG_rWORK_PL(p_U_yz), iPL, Nx, Ny, Nz 
     &                )

            
            CALL BMG3_SymStd_COPY_RHS_yz( 
     &                SO(p_SO), Q(p_U), QF(p_U), BMG_rWORK_PL(p_Q_yz), 
     &                iPL, Nx, Ny, Nz, NStncl_3D
     &                )



            BMG_rPARMS(id_BMG2_STOP_TOL) = rZERO

            CALL BMG2_SymStd_SOLVE_boxmg(
     &                Ny-2, Nz-2, 
     &                BMG_iPARMS_PL_yz(1,kg), BMG_rPARMS, BMG_IOFLAG,
     &                BMG_rWORK_PL(p_U_yz), BMG_rWORK_PL(p_Q_yz),
     &                BMG_rWORK_PL(p_RES_yz), NF_yz, NC_yz,
     &                BMG_rWORK_PL(p_SO_yz), NSO_yz,
     &                BMG_rWORK_PL(p_SOR_yz), NSOR_yz,
     &                BMG_rWORK_PL(p_CI_yz), NCI_yz, 
     &                BMG_rWORK_PL(p_CSO_yz), BMG_rWORK_PL(p_CU_yz),
     &                NCBW_yz, NCU_yz, BMG_iWORK_PL(p_iG_yz),
     &                NOG_yz, NOG_yz
     &                )
            !
            ! Copy the 2D solution back in the iPL{th} plane of 3D
            !
            CALL BMG3_SymStd_COPY_rV_23_yz( 
     &                Q(p_U), BMG_rWORK_PL(p_U_yz), iPL, Nx, Ny, Nz 
     &                )
            !
         ENDDO
         !
      ENDDO

      IF ( BMG_IOFLAG(iBMG3_BUG_RES_RELAX) ) THEN
         CALL BMG3_SymStd_residual( kg, NOG, IFD,
     &             Q(p_U), QF(p_U), SO(p_SO), RES(p_U), Nx, Ny, Nz
     &             )
         CALL BMG3_SymStd_UTILS_norm_l2( 
     &             RES(p_U), Nx, Ny, Nz, RES_L2
     &             )
         WRITE(*,100) 'Residual (l2-norm) after yz sweep = ', RES_L2
      ENDIF

C -------------------------------------------------------------------------
C    >>>>>>>>>>>>>>>>>>>>>>> END: yz sweep <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C -------------------------------------------------------------------------

C ==========================================================================

 100  FORMAT(1X,A,1X,1P,E15.7)

C ===========================================

      RETURN
      END









