!_______________________________________________________
!
!      MSG 2.0 include file 
!   
!      version May 14, 1997
!_______________________________________________________
!
!
!     include the standard MPI definitions
!
      include 'mpif.h'
!
!   status tables for MPI non blocking routines
!
      integer SendStatus(MPI_STATUS_SIZE),             &
     &        RecvStatus(MPI_STATUS_SIZE)
!
!   arrays sendid and recvid contain the IDs
!   used for channel communications
!
      integer MSG_sendid(MAX_PROCS, MAX_PATTERNS),     &
     &        MSG_recvid(MAX_PROCS, MAX_PATTERNS)
!
!   MSGSegment is the pointer to the receive buffer within
!   the array allocated for a buffer
!
      integer MSGSegment(MAX_PATTERNS)
!
!   MSG_COMM identifies the communicator to be used for the MSG
!
      integer MSG_COMM
!
!   MSG_COMM_PARENT identifies the communicator given by the
!   caller (MSG_COMM_WORLD by default)
! 
      integer MSG_COMM_PARENT
!
!   MSG_COMM_PARENT_FLAG tells whether the parent communicator has
!   been changed from the default MSG_COMM_WORLD
!
      integer MSG_COMM_PARENT_FLAG
      integer MSG_COMM_PARENT_MODIFIED
      parameter (MSG_COMM_PARENT_MODIFIED = 12345678)   
!
!   MSG_TRANSFER_TYPE indicates algorithm of data transfer:
!   1 stands for "all to all" transfer
!   0 stands for "series of shifts" transfer
!
      integer MSG_TRANSFER_TYPE(MAX_PATTERNS)
!
!   MSG_COMM_TYPE indicates wheter blocking or non-blocking
!   communication should be used
!   0 .. stands for non-blocking
!   1 .. stands for blocking

      integer MSG_BLOCKING

!-----------------------------------------------------------
      common /MSG_sendrec/ MSG_sendid, MSG_recvid,               &
     &                     SendStatus, RecvStatus,               &
     &                     MSGSegment, MSG_COMM,                 &
     &                     MSG_COMM_PARENT, MSG_COMM_PARENT_FLAG,&
     &                     MSG_TRANSFER_TYPE, MSG_BLOCKING
!
