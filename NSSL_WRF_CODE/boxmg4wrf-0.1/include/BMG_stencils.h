! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     This include file provides a single resource for defining the
!     stencil compass references.
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

! -----------------------
!    2D symmetric:
! -----------------------
      
      INTEGER   ko, knw, ks, ksw, kw
      PARAMETER ( ko=1, kw=2, ks=3, ksw=4, knw=5 )

! -----------------------
!    2D interpolation:
! -----------------------

      INTEGER   LA, LB, LL, LNE, LNW, LR, LSE, LSW      
      PARAMETER ( LL=1, LR=2, LA=3, LB=4, LSW=5, LNW=6, LNE=7, LSE=8 )

! -----------------------
!    3D symmetric:
! -----------------------

      INTEGER   kb, kbe, kbn, kbne, kbnw, kbw, kbs, kbse, kbsw,         & 
     &          kp, kpnw, kps, kpsw, kpw

      PARAMETER( kp=1, kpw=2, kps=3, kb=4, kpsw=5, kpnw=6,              &
     &           kbw=7, kbnw=8, kbn=9, kbne=10, kbe=11,                 &
     &           kbse=12, kbs=13, kbsw=14 )

! ------------------------
!    3d Interpolation:
! ------------------------

      INTEGER  lbne, lbnw, lbse, lbsw, ltne, ltnw, ltse, ltsw,          &
     &         lxya, lxyb, lxyl, lxyr, lxyne, lxynw, lxyse,             &
     &         lxysw, lxza, lxzb, lxzne, lxznw, lxzse, lxzsw,           &
     &         lyzne, lyznw, lyzse, lyzsw

      PARAMETER( lxyl=1, lxyr=2, lxya=3, lxyb=4, lxza=5, lxzb=6,        &
     &           lxyne=7, lxyse=8, lxysw=9, lxynw=10, lxzsw=11,         &
     &           lxznw=12, lxzne=13, lxzse=14, lyzsw=15, lyznw=16,      &
     &           lyzne=17, lyzse=18, lbsw=19, lbnw=20, lbne=21,         &
     &           lbse=22, ltsw=23, ltnw=24, ltne=25, ltse=26 )

! ------------------------
!     SOR workspace:
! ------------------------

      INTEGER   msor, msos, msot, mtot
      PARAMETER( mtot=1, msor=2, msos=3, msot=4 )

      INTEGER    ml1, md1, ml2, md2, ml3, md3, md, ml
      PARAMETER( ml1=2, md1=3, ml2=4, md2=5,                            &
     &           ml3=6, md3=7, md=1, ml=2 )
      
      INTEGER    ipL_BMG_LUD1, ipL_BMG_LUD2, ipL_BMG_LUD3,              &
     &           ipL_BMG_LUL1, ipL_BMG_LUL2, ipL_BMG_LUL3

      PARAMETER ( ipL_BMG_LUL1=1, ipL_BMG_LUD1=2,                       &
     &            ipL_BMG_LUL2=3, ipL_BMG_LUD2=4,                       &
     &            ipL_BMG_LUL3=5, ipL_BMG_LUD3=6 )

! ==========================================================================






