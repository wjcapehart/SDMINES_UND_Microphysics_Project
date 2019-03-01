      subroutine DistributeStencilHalo3D (s, GlobalGridLocalData, 
     &     GlobalGridActData, ActualDataStart, LocalArraySize, 
     &     StencilSize, MSGMyProc, NumProc, Halo, ierror)
      implicit none
      include 'mpif.h'
      include 'MSG.h'
      include 'cg.h'
      real*8 s(maxndimx+2,maxndimy+2,maxndimz+2,MaxStencilSize)
      integer GlobalGridLocalData(2,3,MaxNumProc)
      integer GlobalGridActData(2,3,MaxNumProc) 
      integer ActualDataStart(3,MaxNumProc)
      integer LocalArraySize(3,MaxNumProc)
      integer StencilSize,MSGMyProc,NumProc,Halo,ierror
c---------------------------------------------------------------------
c
c     Distributes the stencil for boundary points to their neighbors
c
c     arguments:
c
c     s = array containing the stencil (input/output)
c
c     GlobalGridLocalData = array containing the global coordinates of
c                           the local data array (input, see MSG2.0 
c                           manual)
c
c     GlobalGridActData = array conatining global coordinates of the 
c                         actual data, i.e. local plus halo (input, see
c                         MSG2.0 manual)
c
c     ActualDatastart = array containing the indices of the start of 
c                       actual data within the local data array (input,
c                       see MSG2.0 manual)
c
c     LocalArraySize = array containing the size the local data array
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
      real*8 Buffer(2*MaxHalo*MaxStencilSize*(maxndimx+2)**2
     &     + 2*MaxHalo*(maxndimy+2)**2 + 2*MaxHalo*(maxndimz+2)**2)
      integer NumAdjProc
      integer Index(2*MaxHalo*StencilSize*(maxndimx+2)**2 
     &     + 2*MaxHalo*(maxndimy+2)**2+ 2*MaxHalo*(maxndimz+2)**2)
      integer Ipr(55)
      integer Proc(27)
      
      integer i
c---------------------------------------------------------------------
c
c     local variables
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
c     we use 2 as the pattern number for the communication of the
c     stencil data
c


      do i=1,StencilSize
         
c     NumAdjProc must be set to zero beofre the initial call to
c     MSG_tp_setup
c     
         NumAdjProc = 0
c     
c     set up the communication pattern for distribution of the stencil
c     
         call MSG_tp_setup(LocalArraySize,ActualDataStart,
     &        GlobalGridLocalData,GlobalGridActData,NumProc,MSGMyProc,
     &        NumAdjProc,Proc,Ipr,Index,1,1,ierror)
         if (ierror .ne. 0) goto 999
c     
c     distribute the boundary stencils
c     
         call MSG_tbdx_send(s(1,1,1,i),Buffer,NumAdjProc,Proc,Ipr,Index,
     &        MsGptrn2,ierror)
         if (ierror .ne. 0) goto 999
c     
c     receive the boundary stencils
c     
         call MSG_tbdx_receive(s(1,1,1,i),Buffer,NumAdjProc,Proc,Ipr,
     &        Index,MSGptrn2,ierror)
         if (ierror .ne. 0) goto 999

c     
c     close all communication channels associated with this pattern
c
         call MSG_tbdx_close(s(1,1,1,i),Buffer,NumAdjProc,Proc,Ipr,
     &        Index,MSGptrn2,ierror)
         
      enddo
c
 999  return
      end









