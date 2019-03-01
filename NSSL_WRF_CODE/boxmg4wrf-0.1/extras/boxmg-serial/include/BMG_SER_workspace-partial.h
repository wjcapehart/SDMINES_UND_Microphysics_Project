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

C =======================================================================
C ----------------------------------------------
C     Workspace Pointer Array Indeces
C ----------------------------------------------

c$$$      INTEGER   ip_CI,
c$$$     &          ip_CSO,
c$$$     &          ip_CU,
c$$$     &          ip_hG,
c$$$     &          ip_iG,
c$$$     &          ip_Q,
c$$$     &          ip_RES,
c$$$     &          ip_SO,
c$$$     &          ip_SOR,
c$$$     &          ip_U,
c$$$     &          ip_BMG_SER_PCG_P,
c$$$     &          ip_BMG_SER_PCG_R,
c$$$     &          ip_BMG_SER_PCG_Z
c$$$ 
c$$$      PARAMETER( ip_U         = 1,
c$$$     &           ip_Q         = 2,
c$$$     &           ip_RES       = 3,
c$$$     &           ip_SO        = 4,
c$$$     &           ip_SOR       = 5,
c$$$     &           ip_CI        = 6,
c$$$     &           ip_CU        = 7,
c$$$     &           ip_CSO       = 8,
c$$$     &           ip_hG        = 9,
c$$$     &           ip_iG        = 10,
c$$$     &           ip_BMG_SER_PCG_P = 11,
c$$$     &           ip_BMG_SER_PCG_R = 12,
c$$$     &           ip_BMG_SER_PCG_Z = 13  )

C -------------------------------------------------------
C     Number of workspace pointers
C -------------------------------------------------------

      INTEGER   NBMG_SER_pWORK
      PARAMETER ( NBMG_SER_pWORK = 13 )

C =======================================================================
C  ----------------------------------------------------
C     Plane Solve Workspace Pointer Array Indeces
C -----------------------------------------------------

c$$$      INTEGER   ip_BMG_iPARMS_PL_xy,
c$$$     &          ip_BMG_iPARMS_PL_yz,
c$$$     &          ip_BMG_iPARMS_PL_xz,
c$$$     &          ip_BMG_pWORK_PL_xy, 
c$$$     &          ip_BMG_pWORK_PL_yz, 
c$$$     &          ip_BMG_pWORK_PL_xz, 
c$$$     &          ip_BMG_iWORK_pSI_xy,
c$$$     &          ip_BMG_iWORK_pSI_yz,
c$$$     &          ip_BMG_iWORK_pSI_xz,
c$$$     &          ip_BMG_iWORK_pSR_xy,
c$$$     &          ip_BMG_iWORK_pSR_yz,
c$$$     &          ip_BMG_iWORK_pSR_xz,
c$$$     &          ip_BMG_iWORK_PL_SF,          ! index shift to integer data
c$$$     &          ip_BMG_rWORK_PL_SF,          ! index shift to real data
c$$$     &          ip_BMG_SER_CGTEMP_yo,
c$$$     &          ip_BMG_SER_CGTEMP_zo
c$$$
c$$$      PARAMETER ( ip_BMG_iPARMS_PL_xy = 1,
c$$$     &            ip_BMG_iPARMS_PL_yz = 2,
c$$$     &            ip_BMG_iPARMS_PL_xz = 3,
c$$$     &            ip_BMG_pWORK_PL_xy  = 4, 
c$$$     &            ip_BMG_pWORK_PL_yz  = 5, 
c$$$     &            ip_BMG_pWORK_PL_xz  = 6, 
c$$$     &            ip_BMG_iWORK_pSI_xy = 7,
c$$$     &            ip_BMG_iWORK_pSI_yz = 8,
c$$$     &            ip_BMG_iWORK_pSI_xz = 9,
c$$$     &            ip_BMG_iWORK_pSR_xy = 10,
c$$$     &            ip_BMG_iWORK_pSR_yz = 11,
c$$$     &            ip_BMG_iWORK_pSR_xz = 12,
c$$$     &            ip_BMG_iWORK_PL_SF  = 13,
c$$$     &            ip_BMG_rWORK_PL_SF  = 14,
c$$$     &            ip_BMG_SER_CGTEMP_yo    = 15,
c$$$     &            ip_BMG_SER_CGTEMP_zo    = 16 )
c$$$
c$$$      
c$$$      INTEGER   NBMG_iWORK_PL_ptrs
c$$$      PARAMETER ( NBMG_iWORK_PL_ptrs = 16 )
      
C =======================================================================
C ----------------------------------------------
C    Include in workspace
C ----------------------------------------------

c$$$      INTEGER   i_InWORK_SO,
c$$$     &          i_InWORK_U,
c$$$     &          i_InWORK_Q,
c$$$     &          i_InWORK_RES,
c$$$     &          i_InWORK_PCG_P,
c$$$     &          i_InWORK_PCG_R,
c$$$     &          i_InWORK_PCG_Z
c$$$
c$$$      PARAMETER( i_InWORK_SO    = 1, 
c$$$     &           i_InWORK_U     = 2, 
c$$$     &           i_InWORK_Q     = 3,
c$$$     &           i_InWORK_RES   = 4,
c$$$     &           i_InWORK_PCG_P = 5,
c$$$     &           i_InWORK_PCG_R = 6,
c$$$     &           i_InWORK_PCG_Z = 7  )

C -------------------------------------------------------
C     Number of InWORK logical flags
C -------------------------------------------------------

      INTEGER   NBMG_SER_InWORK
      PARAMETER( NBMG_SER_InWORK = 7 )

C =======================================================================
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

C =======================================================================





