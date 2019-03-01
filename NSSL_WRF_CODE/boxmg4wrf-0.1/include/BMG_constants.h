! ==========================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     This include file provides a single resource for defining
!     commonly used constants of the appropriate precision (TYPE).
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


      INTEGER    iZERO, iONE, iTWO, iTHREE, iFOUR, iFIVE,               &
     &           iSIX, iSEVEN, iEIGHT, iNINE, iTEN

      PARAMETER( iZERO = 0, iONE = 1, iTWO  = 2, iTHREE = 3,            &
     &            iFOUR = 4, iFIVE = 5, iSIX = 6, iSEVEN = 7,           &
     &            iEIGHT = 8, iNINE = 9, iTEN = 10             )

! split the lines to avoid problem with line changing length after substitution
      REAL(kind=RKIND) rZERO, rONE, rTWO, rTHREE, rFOUR
      REAL(kind=RKIND) rFIVE, rSIX, rSEVEN, rEIGHT, rNINE, rTEN
      
      PARAMETER( rZERO = 0, rONE = 1, rTWO = 2, rTHREE = 3,             &
     &           rFOUR = 4, rFIVE = 5, rSIX = 6, rSEVEN = 7,            &
     &           rEIGHT = 8, rNINE = 9, rTEN = 10              )


      REAL*8     dZERO

      PARAMETER( dZERO = 0)


! =======================================================================
