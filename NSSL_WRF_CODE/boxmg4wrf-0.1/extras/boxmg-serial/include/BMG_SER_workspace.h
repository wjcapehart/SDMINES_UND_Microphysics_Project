C ==========================================================================
C  --------------------
C   DESCRIPTION:
C  --------------------
C
C     This include file provides a single resource for defining commonly
C     used parameters in the BOXMG family of codes.
C
C =======================================================================
C $license_flag$
C =======================================================================
C  --------------------
C   VARIABLES:
C  --------------------
C
C
C ==========================================================================

C ==========================================================================
C ----------------------------------------------
C     Workspace Pointer Array Indeces
C ----------------------------------------------

      INTEGER   ip_CI,
     &          ip_CSO,
     &          ip_CU,
     &          ip_hG,
     &          ip_iG,
     &          ip_Q,
     &          ip_RES,
     &          ip_SO,
     &          ip_SOR,
     &          ip_U,
     &          ip_BMG_SER_PCG_P,
     &          ip_BMG_SER_PCG_R,
     &          ip_BMG_SER_PCG_Z
 
      PARAMETER( ip_U         = 1,
     &           ip_Q         = 2,
     &           ip_RES       = 3,
     &           ip_SO        = 4,
     &           ip_SOR       = 5,
     &           ip_CI        = 6,
     &           ip_CU        = 7,
     &           ip_CSO       = 8,
     &           ip_hG        = 9,
     &           ip_iG        = 10,
     &           ip_BMG_SER_PCG_P = 11,
     &           ip_BMG_SER_PCG_R = 12,
     &           ip_BMG_SER_PCG_Z = 13  )

C -------------------------------------------------------
C     Number of workspace pointers
C -------------------------------------------------------

      INTEGER   NBMG_SER_pWORK
      PARAMETER ( NBMG_SER_pWORK = 13 )

C ==========================================================================
C -----------------------------------------------------
C     Plane Solve Workspace Pointer Array Indeces
C -----------------------------------------------------

      INTEGER   ip_BMG_iPARMS_PL_xy,
     &          ip_BMG_iPARMS_PL_yz,
     &          ip_BMG_iPARMS_PL_xz,
     &          ip_BMG_pWORK_PL_xy, 
     &          ip_BMG_pWORK_PL_yz, 
     &          ip_BMG_pWORK_PL_xz, 
     &          ip_BMG_iWORK_pSI_xy,
     &          ip_BMG_iWORK_pSI_yz,
     &          ip_BMG_iWORK_pSI_xz,
     &          ip_BMG_iWORK_pSR_xy,
     &          ip_BMG_iWORK_pSR_yz,
     &          ip_BMG_iWORK_pSR_xz,
     &          ip_BMG_iWORK_PL_SF,          ! index shift to integer data
     &          ip_BMG_rWORK_PL_SF,          ! index shift to real data
     &          ip_BMG_SER_CGTEMP_yo,
     &          ip_BMG_SER_CGTEMP_zo

      PARAMETER ( ip_BMG_iPARMS_PL_xy = 1,
     &            ip_BMG_iPARMS_PL_yz = 2,
     &            ip_BMG_iPARMS_PL_xz = 3,
     &            ip_BMG_pWORK_PL_xy  = 4, 
     &            ip_BMG_pWORK_PL_yz  = 5, 
     &            ip_BMG_pWORK_PL_xz  = 6, 
     &            ip_BMG_iWORK_pSI_xy = 7,
     &            ip_BMG_iWORK_pSI_yz = 8,
     &            ip_BMG_iWORK_pSI_xz = 9,
     &            ip_BMG_iWORK_pSR_xy = 10,
     &            ip_BMG_iWORK_pSR_yz = 11,
     &            ip_BMG_iWORK_pSR_xz = 12,
     &            ip_BMG_iWORK_PL_SF  = 13,
     &            ip_BMG_rWORK_PL_SF  = 14,
     &            ip_BMG_SER_CGTEMP_yo    = 15,
     &            ip_BMG_SER_CGTEMP_zo    = 16 )

      
      INTEGER   NBMG_iWORK_PL_ptrs
      PARAMETER ( NBMG_iWORK_PL_ptrs = 16 )
      
C ==========================================================================
C ----------------------------------------------
C    Include in workspace
C ----------------------------------------------

      INTEGER   i_InWORK_SO,
     &          i_InWORK_U,
     &          i_InWORK_Q,
     &          i_InWORK_RES,
     &          i_InWORK_PCG_P,
     &          i_InWORK_PCG_R,
     &          i_InWORK_PCG_Z

      PARAMETER( i_InWORK_SO    = 1, 
     &           i_InWORK_U     = 2, 
     &           i_InWORK_Q     = 3,
     &           i_InWORK_RES   = 4,
     &           i_InWORK_PCG_P = 5,
     &           i_InWORK_PCG_R = 6,
     &           i_InWORK_PCG_Z = 7  )

C -------------------------------------------------------
C     Number of InWORK logical flags
C -------------------------------------------------------

      INTEGER   NBMG_SER_InWORK
      PARAMETER( NBMG_SER_InWORK = 7 )

C ==========================================================================
C  -----------------------------------------------
C     IGRD pointer indeces
C  -----------------------------------------------

      INTEGER  ipL_BMG_SER_U, 
     &         ipL_BMG_SER_SO,
     &         ipL_BMG_SER_SOR,
     &         ipL_BMG_SER_CI

      PARAMETER ( ipL_BMG_SER_U       = 1, 
     &            ipL_BMG_SER_SO      = 2,
     &            ipL_BMG_SER_SOR     = 3,
     &            ipL_BMG_SER_CI      = 4  )

C ------------------------------------------------
C    IGRD data indeces
C ------------------------------------------------

      INTEGER  idL_BMG_SER_IVW,
     &         idL_BMG_SER_Nx,
     &         idL_BMG_SER_Ny,
     &         idL_BMG_SER_Nz

      PARAMETER ( idL_BMG_SER_IVW     = 5,
     &            idL_BMG_SER_Nx      = 6,
     &            idL_BMG_SER_Ny      = 7,
     &            idL_BMG_SER_Nz      = 8  )

C ----------------------------------------------
C     Number of workspace pointers
C ----------------------------------------------

      INTEGER   NBMG_SER_pIGRD
      PARAMETER ( NBMG_SER_pIGRD = 8 )

C ==========================================================================





