! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     This include file provides a single resource for defining commonly
!     used workspace parameters in the BOXMG family of codes.
!
! =======================================================================
! $license_flag$
! =======================================================================
!  --------------------
!   VARIABLES:
!  --------------------
!
!
! ==========================================================================

! -------------------------------------------------------
!     Workspace Pointer Array Indeces: serial
! -------------------------------------------------------

      INTEGER   ip_CI,                                                  &
     &          ip_CSO,                                                 &
     &          ip_CU,                                                  &
     &          ip_hG,                                                  &
     &          ip_iG,                                                  &
     &          ip_Q,                                                   &
     &          ip_RES,                                                 &
     &          ip_SO,                                                  &
     &          ip_SOR,                                                 &
     &          ip_U,                                                   &
     &          ip_BMG_PCG_P,                                           &
     &          ip_BMG_PCG_R,                                           &
     &          ip_BMG_PCG_Z
 
      PARAMETER( ip_U         = 1,                                      &
     &           ip_Q         = 2,                                      &
     &           ip_RES       = 3,                                      &
     &           ip_SO        = 4,                                      &
     &           ip_SOR       = 5,                                      &
     &           ip_CI        = 6,                                      &
     &           ip_CU        = 7,                                      &
     &           ip_CSO       = 8,                                      &
     &           ip_hG        = 9,                                      &
     &           ip_iG        = 10,                                     &
     &           ip_BMG_PCG_P = 11,                                     &
     &           ip_BMG_PCG_R = 12,                                     &
     &           ip_BMG_PCG_Z = 13  )

! -------------------------------------------------------
!     Workspace Pointer Array Indeces: parallel
! -------------------------------------------------------

      INTEGER    ip_LS,                                                 &
     &           ip_MSG,                                                &
     &           ip_MSGSO,                                              &
     &           ip_MSG_BUF,                                            &
     &           ip_ProcGRID

      PARAMETER ( ip_MSG       = 14,                                    &
     &            ip_MSGSO     = 15,                                    &
     &            ip_MSG_BUF   = 16,                                    &
     &            ip_ProcGrid  = 17,                                    &
     &            ip_LS        = 18  )

! -------------------------------------------------------
!     Number of workspace pointers
! -------------------------------------------------------

      INTEGER   NBMG_pWORK
      PARAMETER ( NBMG_pWORK = 18 )

! =======================================================================
! -------------------------------------------------------
!     Plane Solve Workspace Pointer Array Indeces
! -------------------------------------------------------

      INTEGER   ip_BMG_iPARMS_PL_xy,                                    &
     &          ip_BMG_iPARMS_PL_yz,                                    &
     &          ip_BMG_iPARMS_PL_xz,                                    &
     &          ip_BMG_pWORK_PL_xy,                                     &
     &          ip_BMG_pWORK_PL_yz,                                     &
     &          ip_BMG_pWORK_PL_xz,                                     &
     &          ip_BMG_iWORK_pSI_xy,                                    &
     &          ip_BMG_iWORK_pSI_yz,                                    &
     &          ip_BMG_iWORK_pSI_xz,                                    &
     &          ip_BMG_iWORK_pSR_xy,                                    &
     &          ip_BMG_iWORK_pSR_yz,                                    &
     &          ip_BMG_iWORK_pSR_xz,                                    &
     &          ip_BMG_iWORK_PL_SF,                                     & 
     &          ip_BMG_rWORK_PL_SF,                                     &
     &          ip_BMG_CGTEMP_yo,                                       &
     &          ip_BMG_CGTEMP_zo

      PARAMETER ( ip_BMG_iPARMS_PL_xy = 1,                              &
     &            ip_BMG_iPARMS_PL_yz = 2,                              &
     &            ip_BMG_iPARMS_PL_xz = 3,                              &
     &            ip_BMG_pWORK_PL_xy  = 4,                              &
     &            ip_BMG_pWORK_PL_yz  = 5,                              & 
     &            ip_BMG_pWORK_PL_xz  = 6,                              & 
     &            ip_BMG_iWORK_pSI_xy = 7,                              &
     &            ip_BMG_iWORK_pSI_yz = 8,                              &
     &            ip_BMG_iWORK_pSI_xz = 9,                              &
     &            ip_BMG_iWORK_pSR_xy = 10,                             &
     &            ip_BMG_iWORK_pSR_yz = 11,                             &
     &            ip_BMG_iWORK_pSR_xz = 12,                             &
     &            ip_BMG_iWORK_PL_SF  = 13,                             &
     &            ip_BMG_rWORK_PL_SF  = 14,                             &
     &            ip_BMG_CGTEMP_yo    = 15,                             &
     &            ip_BMG_CGTEMP_zo    = 16 )

      
      INTEGER   NBMG_iWORK_PL_ptrs
      PARAMETER ( NBMG_iWORK_PL_ptrs = 16 )

! =======================================================================
! -------------------------------------------------------
!     Plane Solve Workspace Pointer Array Indeces
! -------------------------------------------------------

      INTEGER   id_BMG_iWORK_CS_NCBW,                                   & 
     &          id_BMG_iWORK_CS_NCU,                                    &
     &          id_BMG_NOGm_CS,                                         &
     &          id_BMG_iWORK_CS,                                        &
     &          id_BMG_rWORK_CS,                                        &
     &          id_BMG_iWORK_PL_CS,                                     &
     &          id_BMG_rWORK_PL_CS,                                     &
     &          ip_BMG_iWORK_CS_CSO,                                    &
     &          ip_BMG_iWORK_CS_CU,                                     &
     &          ip_BMG_iPARMS_CS,                                       &
     &          ip_BMG_rPARMS_CS,                                       &
     &          ip_BMG_pWORK_CS,                                        &
     &          ip_BMG_iWORK_CS,                                        &
     &          ip_BMG_rWORK_CS,                                        &
     &          ip_BMG_iWORK_PL_CS,                                     &
     &          ip_BMG_rWORK_PL_CS,                                     &
     &          ip_BMG_iWORK_CS_SF,                                     & 
     &          ip_BMG_rWORK_CS_SF                                      

      PARAMETER ( id_BMG_iWORK_CS_NCBW    = 1,                          & 
     &            id_BMG_iWORK_CS_NCU     = 2,                          &
     &            id_BMG_NOGm_CS          = 3,                          &
     &            id_BMG_iWORK_CS         = 4,                          &
     &            id_BMG_rWORK_CS         = 5,                          &
     &            id_BMG_iWORK_PL_CS      = 6,                          &
     &            id_BMG_rWORK_PL_CS      = 7,                          &
     &            ip_BMG_iWORK_CS_CSO     = 8,                          &
     &            ip_BMG_iWORK_CS_CU      = 9,                          &
     &            ip_BMG_iPARMS_CS        = 10,                         &
     &            ip_BMG_rPARMS_CS        = 11,                         &
     &            ip_BMG_pWORK_CS         = 12,                         &
     &            ip_BMG_iWORK_CS         = 13,                         &
     &            ip_BMG_rWORK_CS         = 14,                         &
     &            ip_BMG_iWORK_PL_CS      = 15,                         &
     &            ip_BMG_rWORK_PL_CS      = 16,                         &
     &            ip_BMG_iWORK_CS_SF      = 17,                         & 
     &            ip_BMG_rWORK_CS_SF      = 18  )                           

      INTEGER   NBMG_iWORK_CS_ptrs
      PARAMETER ( NBMG_iWORK_CS_ptrs = 18 )

! =======================================================================
! -------------------------------------------------------
!     Include in workspace
! -------------------------------------------------------

      INTEGER   i_InWORK_SO,                                            &
     &          i_InWORK_U,                                             &
     &          i_InWORK_Q,                                             &
     &          i_InWORK_RES,                                           &
     &          i_InWORK_PCG_P,                                         & 
     &          i_InWORK_PCG_R,                                         & 
     &          i_InWORK_PCG_Z

      PARAMETER( i_InWORK_SO    = 1,                                    & 
     &           i_InWORK_U     = 2,                                    & 
     &           i_InWORK_Q     = 3,                                    &
     &           i_InWORK_RES   = 4,                                    &
     &           i_InWORK_PCG_P = 5,                                    &
     &           i_InWORK_PCG_R = 6,                                    &
     &           i_InWORK_PCG_Z = 7  )

! -------------------------------------------------------
!     Number of InWORK logical flags
! -------------------------------------------------------

      INTEGER   NBMG_InWORK
      PARAMETER( NBMG_InWORK = 7 )
     
! =======================================================================
! -------------------------------------------------------
!     IGRD pointer indeces: serial
! -------------------------------------------------------

      INTEGER  ipL_BMG_U,                                               &
     &         ipL_BMG_SO,                                              &
     &         ipL_BMG_SOR,                                             &
     &         ipL_BMG_CI

      PARAMETER ( ipL_BMG_U       = 1,                                  & 
     &            ipL_BMG_SO      = 2,                                  &
     &            ipL_BMG_SOR     = 3,                                  &
     &            ipL_BMG_CI      = 4  )

! -------------------------------------------------------
!    IGRD data indeces: serial/parallel
! -------------------------------------------------------

      INTEGER  idL_BMG_IVW,                                             &
     &         idL_BMG_NLx,                                             &
     &         idL_BMG_NLy,                                             &
     &         idL_BMG_NLz,                                             &
     &         idL_BMG_NGx,                                             &
     &         idL_BMG_NGy,                                             &
     &         idL_BMG_NGz,                                             &
     &         idL_BMG_Icoord,                                          &
     &         idL_BMG_Jcoord,                                          &
     &         idL_BMG_Kcoord

      PARAMETER ( idL_BMG_IVW     =  5,                                 &
     &            idL_BMG_NLx     =  6,                                 &
     &            idL_BMG_NLy     =  7,                                 &
     &            idL_BMG_NLz     =  8,                                 &
     &            idL_BMG_NGx     =  9,                                 &
     &            idL_BMG_NGy     = 10,                                 &
     &            idL_BMG_NGz     = 11,                                 & 
     &            idL_BMG_Icoord  = 12,                                 &
     &            idL_BMG_Jcoord  = 13,                                 &
     &            idL_BMG_Kcoord  = 14  )

! -------------------------------------------------------
!     Number of workspace pointers
! -------------------------------------------------------

      INTEGER   NBMG_pIGRD
      PARAMETER ( NBMG_pIGRD = 14 )

! =======================================================================
! -------------------------------------------------------
!     MSG workspace pointer indices
! -------------------------------------------------------

      INTEGER ipL_MSG_Index,                                            &
     &        ipL_MSG_LocalArraySize,                                   &
     &        ipL_MSG_Proc,                                             &
     &        ipL_MSG_Ipr,                                              &
     &        ipL_MSG_NumAdjProc,                                       &
     &        ipL_MSG_ActDataStart,                                     &
     &        ipL_MSG_GlobalCoordLocalData,                             &
     &        ipL_MSG_GlobalCoordActData,                               &
     &        ipL_MSG_ProcGrid,                                         &
     &        ipL_MSG_ProcGridCoord_x,                                  &
     &        ipL_MSG_ProcGridCoord_y,                                  &
     &        ipL_MSG_ProcGridCoord_z

      PARAMETER ( ipL_MSG_Index                = 1,                     &
     &            ipL_MSG_LocalArraySize       = 2,                     &
     &            ipL_MSG_Proc                 = 3,                     &
     &            ipL_MSG_Ipr                  = 4,                     &
     &            ipL_MSG_NumAdjProc           = 5,                     &
     &            ipL_MSG_ActDataStart         = 6,                     &
     &            ipL_MSG_GlobalCoordLocalData = 7,                     &
     &            ipL_MSG_GlobalCoordActData   = 8,                     &
     &            ipL_MSG_ProcGrid             = 9,                     &
     &            ipL_MSG_ProcGridCoord_x      = 10,                    &
     &            ipL_MSG_ProcGridCoord_y      = 11,                    &
     &            ipL_MSG_ProcGridCoord_z      = 12  )

! -------------------------------------------------------
!     Number of MSG workspace pointers
! -------------------------------------------------------

      INTEGER   NBMG_pMSG
      PARAMETER ( NBMG_pMSG = 12 )
      
! ==================================================================
! -------------------------------------------------------
!     1D Line Solve workspace pointer indices
! -------------------------------------------------------

      INTEGER ipL_LS_XDataDist,                                         &
     &        ipL_LS_YDataDist,                                         &
     &        ipL_LS_ZDataDist
      
      PARAMETER ( ipL_LS_XDataDist = 1,                                 & 
     &            ipL_LS_YDataDist = 2,                                 &
     &            ipL_LS_ZDataDist = 3  )

! -------------------------------------------------------
!     Number of Line Solve workspace pointers
! -------------------------------------------------------

      INTEGER  NBMG_pLS
      PARAMETER ( NBMG_pLS = 3 )

! ==================================================================
! -------------------------------------------------------
!     Indeces of pointers for MSG workspace:
! -------------------------------------------------------

      INTEGER  ip_BMG_MSG_iWORK_pLS,                                    &
     &         ip_BMG_MSG_iWORK_pU,                                     &
     &         ip_BMG_MSG_iWORK_pSO,                                    & 
     &         ip_BMG_MSG_iWORK_pyo,                                    &
     &         ip_BMG_MSG_iWORK_pzo,                                    &
     &         ip_BMG_MSG_iWORK_LS,                                     &
     &         ip_BMG_MSG_iWORK_U,                                      &
     &         ip_BMG_MSG_iWORK_SO,                                     &
     &         ip_BMG_MSG_iWORK_yo,                                     &
     &         ip_BMG_MSG_iWORK_zo

      PARAMETER ( ip_BMG_MSG_iWORK_pLS = 1,                             &
     &            ip_BMG_MSG_iWORK_pU  = 2,                             &
     &            ip_BMG_MSG_iWORK_pSO = 3,                             & 
     &            ip_BMG_MSG_iWORK_pyo = 4,                             &
     &            ip_BMG_MSG_iWORK_pzo = 5,                             &
     &            ip_BMG_MSG_iWORK_LS  = 6,                             &
     &            ip_BMG_MSG_iWORK_U   = 7,                             &
     &            ip_BMG_MSG_iWORK_SO  = 8,                             &
     &            ip_BMG_MSG_iWORK_yo  = 9,                             &
     &            ip_BMG_MSG_iWORK_zo  = 10 )

      INTEGER   NBMG_MSG_pWORK
      PARAMETER ( NBMG_MSG_pWORK = 10 )

! ==================================================================
! -------------------------------------------------------
!     Data pointers for MPI Processor Grid
! -------------------------------------------------------

      INTEGER  id_BMG_MSG_NGx,                                          &   
     &         id_BMG_MSG_NGy,                                          &    
     &         id_BMG_MSG_NGz,                                          &
     &         id_BMG_MSG_NLx,                                          &
     &         id_BMG_MSG_NLy,                                          &
     &         id_BMG_MSG_NLz,                                          &
     &         id_BMG_MSG_iGs,                                          &
     &         id_BMG_MSG_jGs,                                          &
     &         id_BMG_MSG_kGs,                                          &
     &         id_BMG_MSG_NProc,                                        &  
     &         id_BMG_MSG_NProcI,                                       &
     &         id_BMG_MSG_NProcJ,                                       &
     &         id_BMG_MSG_NProcK,                                       &
     &         id_BMG_MSG_MyProc,                                       &
     &         id_BMG_MSG_MyProcI,                                      &
     &         id_BMG_MSG_MyProcJ,                                      &
     &         id_BMG_MSG_MyProcK,                                      &
     &         id_BMG_MSG_COMM,                                         &
     &         id_BMG_MSG_COMM_x,                                       &
     &         id_BMG_MSG_COMM_y,                                       &
     &         id_BMG_MSG_COMM_xy,                                      &
     &         id_BMG_MSG_COMM_yz,                                      &
     &         id_BMG_MSG_COMM_xz

      PARAMETER( id_BMG_MSG_NGx      = 1,                               &
     &           id_BMG_MSG_NGy      = 2,                               &
     &           id_BMG_MSG_NGz      = 3,                               &
     &           id_BMG_MSG_NLx      = 4,                               &
     &           id_BMG_MSG_NLy      = 5,                               &
     &           id_BMG_MSG_NLz      = 6,                               &
     &           id_BMG_MSG_iGs      = 7,                               &
     &           id_BMG_MSG_jGs      = 8,                               &
     &           id_BMG_MSG_kGs      = 9,                               &
     &           id_BMG_MSG_NProc    = 10,                              &
     &           id_BMG_MSG_NProcI   = 11,                              &
     &           id_BMG_MSG_NProcJ   = 12,                              &
     &           id_BMG_MSG_NProcK   = 13,                              & 
     &           id_BMG_MSG_MyProc   = 14,                              &
     &           id_BMG_MSG_MyProcI  = 15,                              & 
     &           id_BMG_MSG_MyProcJ  = 16,                              &
     &           id_BMG_MSG_MyProcK  = 17,                              &
     &           id_BMG_MSG_COMM     = 18,                              &
     &           id_BMG_MSG_COMM_x   = 19,                              &
     &           id_BMG_MSG_COMM_y   = 20,                              &
     &           id_BMG_MSG_COMM_xy  = 21,                              &
     &           id_BMG_MSG_COMM_yz  = 22,                              &
     &           id_BMG_MSG_COMM_xz  = 23  )

      INTEGER   NBMG_MSG_iGRID_data
      PARAMETER ( NBMG_MSG_iGRID_data = 23 ) 

! -------------------------------------------------------
!     Indeces to pointers for MPI Processor Grid
! -------------------------------------------------------

      INTEGER  ip_BMG_MSG_NLx_Grid,                                     &
     &         ip_BMG_MSG_NLy_Grid,                                     &
     &         ip_BMG_MSG_NLz_Grid,                                     &
     &         ip_BMG_MSG_iGs_Grid,                                     &
     &         ip_BMG_MSG_jGs_Grid,                                     &
     &         ip_BMG_MSG_kGs_Grid,                                     &
     &         ip_BMG_MSG_ProcGrid,                                     &
     &         ip_BMG_MSG_ProcCoord

      PARAMETER ( ip_BMG_MSG_NLx_Grid = 1,                              &
     &            ip_BMG_MSG_NLy_Grid = 2,                              &
     &            ip_BMG_MSG_NLz_Grid = 3,                              &
     &            ip_BMG_MSG_iGs_Grid = 4,                              &
     &            ip_BMG_MSG_jGs_Grid = 5,                              &
     &            ip_BMG_MSG_kGs_Grid = 6,                              &
     &            ip_BMG_MSG_ProcGrid = 7,                              &
     &            ip_BMG_MSG_ProcCoord = 8  )

      INTEGER   NBMG_MSG_pGRID
      PARAMETER ( NBMG_MSG_pGRID = 8 )

! -------------------------------------------------------
!     Grid Distribution methods
! -------------------------------------------------------

      INTEGER   BMG_MSG_GRIDDIST_UNIFORM,                               &
     &          BMG_MSG_GRIDDIST_CUSTOM,                                &
     &          BMG_MSG_GRIDDIST_AUTOMATIC

      PARAMETER ( BMG_MSG_GRIDDIST_UNIFORM   = 1,                       &
     &            BMG_MSG_GRIDDIST_CUSTOM    = 2,                       &
     &            BMG_MSG_GRIDDIST_AUTOMATIC = 3  )

! ==================================================================


! >>>>>>>>>> BELOW HERE IS STUFF THAT IS GOING AWAY !!!!

! ==================================================================

      INTEGER MaxNumProc
      PARAMETER ( MaxNumProc = 100 )

! ==================================================================
! -------------------------------------------------------
!     Workspace pointers for processor grid information
! -------------------------------------------------------

      INTEGER ipL_PGrid,                                                &
     &        ipL_NProcI,                                               &
     &        ipL_NProcJ,                                               &
     &        ipL_NProcK,                                               &
     &        ipL_MyProcI,                                              &
     &        ipL_MyProcJ,                                              &
     &        ipL_MyProcK

      PARAMETER ( ipL_PGrid   = 1,                                      &
     &            ipL_NProcI  = 2,                                      &
     &            ipL_NProcJ  = 3,                                      &
     &            ipL_NProcK  = 4,                                      &
     &            ipL_MyProcI = 5,                                      &
     &            ipL_MyProcJ = 6,                                      &
     &            ipL_MyProcK = 7  )


! ---------------------------------------------------
!     Number of Processor Grid workspace pointers:
! ---------------------------------------------------

      INTEGER   NBMG_ProcGridInfo
      PARAMETER ( NBMG_ProcGridInfo = 7 )

! ==================================================================


