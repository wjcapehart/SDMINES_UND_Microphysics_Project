!_____________________________________________
!
!    include file to describe the MSG functions 
!
!                 version 2.0 
!
!---------------------------------------------

      real*8 MSG_timer
      external MSG_timer 
      integer MSG_myproc
      external MSG_myproc
      integer MSG_nproc
      external MSG_nproc
      integer MSG_VERSION
      common /MSG_GLOBAL_DATA/ MSG_VERSION 
!
