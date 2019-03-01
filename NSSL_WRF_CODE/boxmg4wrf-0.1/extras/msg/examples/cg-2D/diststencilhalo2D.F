      subroutine DistributeStencilHalo2D (s, GlobalCoordLocalData, 
     &     GlobalCoordActData, ActualDataStart, LocalArraySize,
     &     StencilSize, MSGMyProc, NumProc, ierror)
      implicit none
      include 'mpif.h'
      include 'MSG.h'
      include 'cg.h'
      real*8  s(maxndimx+2,maxndimy+2,MaxStencilSize)
      integer GlobalCoordLocalData(2,3,MaxNumProc)
      integer GlobalCoordActData(2,3,MaxNumProc) 
      integer ActualDataStart(3,MaxNumProc)
      integer LocalArraySize(3,MaxNumProc)
      integer StencilSize,MSGMyProc,NumProc,ierror
c---------------------------------------------------------------------
c
c     Distributes the stencil for boundary points to their neighbors
c
c     arguments:
c
c     s = array containing the stencil (input/output)
c
c     GlobalCoordLocalData = array containing the global coordinates of
c                           the local data array (input, see MSG2.0 
c                           manual)
c
c     GlobalCoordActData = array conatining global coordinates of the 
c                         actual data, i.e. local plus halo (input, see
c                         MSG2.0 manual)
c
c     ActualDatastart = array containing the indices of the start of 
c                       actual data within the local data array (input,
c                       see MSG2.0 manual)
c
c     StencilSize = number of non zeros in the stencil, this is the 
c                   leading dimension of the array s (input)
c
c     MSGMyProc = processor number in MSG (>=1) (input)
c
c     NumProc = number of processors (input)
c
c     Halo = number of ghost lines around the local data (input)
c
c     ierror = error codes from MPI/MSG calls, if .ne. 0, something
c              went wrong (output)
c
c---------------------------------------------------------------------
      real*8 Buffer(10*2*(maxndimx+maxndimy))
      integer NumAdjProc
      integer Index(10*2*(maxndimx+2+maxndimy+2))
      integer Ipr(27)
      integer Proc(10)
      
      integer i
c---------------------------------------------------------------------
c
c     local variables
c
c     MatGcld = array conatining the global coordinates of the stencil
c               data within s
c
c     MatGcad = array containing the global coordinates of the actual
c               stencil data within s (halo and local)
c
c     Las = dimensions of the local arrays s
c
c     Ads = coodrinates where the actual data starts within s
c
c     Buffer = buffer for MSG communication routines
c
c     NumAdjProc = number of adjacent processors (see MSG2.0 manual)
c
c     prtn = number of MSG communication pattern 
c
c     Index = index array for MSG communication, set by MSG_tp_setup
c             (see MSG2.0 manual)
c
c     Ipr = table of pointers, set by MSG_tp_setup (see MSG2.0 manual 
c
c     Proc = list of adjacent processors, set by MSG_tp_setup (see MSG2.0
c            manual)
c
c     i = index variable
c
c---------------------------------------------------------------------
c
c     we need to set up a communication pattern for the matrix that is 
c     based on the given communication pattern for the data
c
c
c     set up the communication pattern for distribution of the stencil
c

      do i=1,MaxStencilSize

         NumAdjProc = 0
         
         call MSG_tp_setup(LocalArraySize,ActualDataStart,
     &        GlobalCoordLocalData,GlobalCoordActData,NumProc,MSGMyProc,
     &        NumAdjProc,Proc,Ipr,Index,1,1,ierror)
         if (ierror .ne. 0) goto 999
         
c     write(*,*) MSGMyProc,(Proc(i),i=1,NumAdjProc)
         
         
c     
c     distribute the boundary stencils
c     
         call MSG_tbdx_send(s(1,1,i),Buffer,NumAdjProc,Proc,Ipr,Index,
     &        MSGptrn2,ierror)
         if (ierror .ne. 0) goto 999
c     
c     receive the boundary stencils
c     
         call MSG_tbdx_receive(s(1,1,i),Buffer,NumAdjProc,Proc,Ipr,
     &        Index,MSGptrn2,ierror)
         if (ierror .ne. 0) goto 999
c     
c     close all communication channels associated with this pattern
c     
         call MSG_tbdx_close(s(1,1,i),Buffer,NumAdjProc,Proc,Ipr,Index,
     &        MSGptrn2,ierror)
      enddo
c     
 999  return
      end











