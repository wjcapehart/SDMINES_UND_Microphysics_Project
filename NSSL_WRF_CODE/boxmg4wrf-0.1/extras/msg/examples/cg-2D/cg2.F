      program conjgrad2D
      implicit none
      include 'mpif.h'
      include 'MSG.h'
      include 'cg.h'
      integer ndim, i, j, k, ierror
      integer locndimx,locndimy
      integer MSGMyProc, NumProc, MyProc
      integer MyProci, MyProcj
      integer Proc(8), Ipr(17), NumAdjProc,ldim
      integer SqrtNumProc

      integer Index(4*(maxndimx+maxndimz+4))
      integer LocalArraySize(3,MaxNumProc)
      integer ActDatastart(3,MaxNumProc)
      integer GlobalCoordLocalData(2,3,MaxNumProc)
      integer GlobalCoordActData(2,3,MaxNumProc)


      real*8 r(maxndimx+2,maxndimy+2)
      real*8 x(maxndimx+2,maxndimy+2)
      real*8 b(maxndimx+2,maxndimy+2)
      real*8 p(maxndimx+2,maxndimy+2)
      real*8 w(maxndimx+2,maxndimy+2)
      real*8 z(maxndimx+2,maxndimy+2)

      real*8 buf(8*maxndimx+4)

      real*8 h, beta, alpha
      real*8 rzprod2, rzprod1, rzprod
      real*8 rbuf,bnorm

      real*8 s(maxndimx+2,maxndimy+2,MaxStencilSize)


      call MPI_INIT(ierror)
      call MSG_ENABLE(MSGMyProc, NumProc)
      MyProc = MSGMyProc-1

      SqrtNumProc = sqrt(real(NumProc))
      if (SqrtNumProc**2 .ne. NumProc) then
         write(*,*)NumProc,' is not square'
         goto 999
      endif

      MyProcj = MyProc/SqrtNumProc 
      MyProci = MyProc - MyProcj*SqrtNumProc

      MyProci = MyProci+1
      MyProcj = MyProcj+1

c     calculate the size of the local arrays
      if (MyProc .eq. 0) then
         write(*,*)'enter number of gridlines '
         read(*,*) ndim
      endif
      call MPI_BCAST(ndim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      if (ierror.ne.0) write(*,*)'ierror = ',ierror

      locndimx = ndim / SqrtNumProc
      if (MyProci .eq. SqrtNumProc) then
         locndimx = locndimx+(ndim-locndimx*SqrtNumProc)
      endif

      locndimy = ndim / SqrtNumProc
      if (MyProcj .eq. SqrtNumProc) then
         locndimy = locndimy+(ndim-locndimy*SqrtNumProc)
      endif

c      write(*,*) MSGMyProc, MyProci, MyProcj, locndimx, locndimy

c
c     initialize the stencil matrix
c
c     in this example we are using a
c     nine point stencil for the Laplacian
c     in 2D
c
      do i=2,locndimx+1
         do j=2,locndimy+1
            s(i,j,k0)=8.0
            s(i,j,kw)=-1.0
            s(i,j,ks)=-1.0
            s(i,j,knw)=-1.0
            s(i,j,ksw)=-1.0
         enddo
      enddo
c
c     processors on the upper boundary of the 
c     processor grid need to initially have
c     all the information for their points on 
c     the physical boundary
c
      if (MyProci.eq. SqrtNumProc) then
         do j=2,locndimy+1
            s(locndimx+2,j,k0)=8.0
            s(locndimx+2,j,kw)=-1.0
            s(locndimx+2,j,ks)=-1.0
            s(locndimx+2,j,knw)=-1.0
            s(locndimx+2,j,ksw)=-1.0
         enddo
      endif

      if (MyProcj.eq. SqrtNumProc) then
         do i=2,locndimx+1
            s(i,locndimy+2,k0)=8.0
            s(i,locndimy+2,kw)=-1.0
            s(i,locndimy+2,ks)=-1.0
            s(i,locndimy+2,knw)=-1.0
            s(i,locndimy+2,ksw)=-1.0
         enddo
      endif

      if ( (MyProci .eq. SqrtNumProc) .and.
     &     (MyProcj .eq. SqrtNumProc) ) then
         s(locndimx+2,locndimy+2,k0)=8.0
         s(locndimx+2,locndimy+2,kw)=-1.0
         s(locndimx+2,locndimy+2,ks)=-1.0
         s(locndimx+2,locndimy+2,knw)=-1.0
         s(locndimx+2,locndimy+2,ksw)=-1.0
      endif


      h=1.0/ndim

c     initialize right hand side 
      do j=2, locndimy+1
         do i=2, locndimx+1
            b(i,j) = 3*h**2
         enddo
      enddo
         
      

c     set up variables for CG iteration
      k=0
      
      do j=1,locndimy+2
         do i=1, locndimx+2
            x(i,j) = 0.0
         enddo
      enddo
      
c     calculate the initial residual 
c     (it is the right hand side because we
c     assume a zero initial guess)
      do j=2,locndimy+1
         do i=2,locndimx+1
            r(i,j) = b(i,j)
         enddo
      enddo
      
      bnorm = 0.0
      do j=2,locndimy+1
         do i=2,locndimx+1
            bnorm = bnorm + b(i,j)**2
         enddo
      enddo
      call MPI_AllReduce(bnorm, rbuf, 1, MPI_REAL8, MPI_SUM, 
     &                   MPI_COMM_WORLD,ierror)
      bnorm = rbuf
      if (MyProc .eq. 0) write(*,*)'b-norm = ',bnorm


c     initialize variables for call to MSG
      NumAdjProc = 0
      ldim = ndim/SqrtNumProc

      do i=1,SqrtNumProc
         do j=1,SqrtNumProc
            LocalArraySize(1,(j-1)*SqrtNumProc+i)=maxndimx+2
            LocalArraySize(2,(j-1)*SqrtNumProc+i)=maxndimy+2
            LocalArraySize(3,(j-1)*SqrtNumProc+i)=1
            
            ActDataStart(1,(j-1)*SqrtNumProc+i)=1
            ActDataStart(2,(j-1)*SqrtNumProc+i)=1
            ActDataStart(3,(j-1)*SqrtNumProc+i)=1

            GlobalCoordLocalData(1,1,(j-1)*SqrtNumProc+i)
     &           =(i-1)*ldim+2
            GlobalCoordLocalData(1,2,(j-1)*SqrtNumProc+i)
     &           =(j-1)*ldim+2
            GlobalCoordLocalData(1,3,(j-1)*SqrtNumProc+i)=1

            if (i.ne.SqrtNumProc) then
               GlobalCoordLocalData(2,1,(j-1)*SqrtNumProc+i)=i*ldim+1
            else
               GlobalCoordLocalData(2,1,(j-1)*SqrtNumProc+i)=ndim+1
            endif
            if (j.ne.SqrtNumProc) then
               GlobalCoordLocalData(2,2,(j-1)*SqrtNumProc+i)=j*ldim+1
            else
               GlobalCoordLocalData(2,2,(j-1)*SqrtNumProc+i)=ndim+1
            endif
            GlobalCoordLocalData(2,3,(j-1)*SqrtNumProc+i)=1
            
            GlobalCoordActData(1,1,(j-1)*SqrtNumProc+i)=(i-1)*ldim+1
            GlobalCoordActData(1,2,(j-1)*SqrtNumProc+i)=(j-1)*ldim+1
            GlobalCoordActData(1,3,(j-1)*SqrtNumProc+i)=1
 
            if (i.eq.SqrtNumProc) then
               GlobalCoordActData(2,1,(j-1)*SqrtNumProc+i)=ndim+2
            else
               GlobalCoordActData(2,1,(j-1)*SqrtNumProc+i)=i*ldim+2
            endif
            if (j.eq.SqrtNumProc) then
               GlobalCoordActData(2,2,(j-1)*SqrtNumProc+i)=ndim+2
            else
               GlobalCoordActData(2,2,(j-1)*SqrtNumProc+i)=j*ldim+2
            endif
            GlobalCoordActData(2,3,(j-1)*SqrtNumProc+i)=1
         enddo
      enddo

      call MSG_tp_setup(LocalArraySize, ActDataStart,
     &                  GlobalCoordLocalData, GlobalCoordActData,
     &                  NumProc, MSGMyProc, NumAdjProc, Proc, Ipr,
     &                  Index, 1, 1, ierror)  
c
c     distribute the stencil information to the neighbor
c     processors that need it
      call DistributeStencilHalo2D (s,GlobalCoordLocalData, 
     &     GlobalCoordActData, ActDataStart, LocalArraySize, 
     &     5,MSGMyProc,NumProc,ierror)

c      do i=locndimx,1,-1
c      write(*,*) MyProci, MyProcj, (s(i,j,knw),j=1,locndimy)
c      enddo
c      write(*,*) '-------'

c
c     the CG iteration starts here
c
c     preconditioner ... no preconditioner here, yet
 111  do j=2,locndimy+1
         do i=2,locndimx+1
            z(i,j) = r(i,j)
         enddo
      enddo
      
      k = k+1

      if (k .eq. 1) then
c     calculate (r,z)
         rzprod1 = 0.
         do j=2,locndimy+1
            do i=2,locndimx+1
               rzprod1 = rzprod1 + r(i,j)*z(i,j)
            enddo
         enddo
         call MPI_AllReduce(rzprod1, rbuf, 1, MPI_REAL8, MPI_SUM, 
     &        MPI_COMM_WORLD,ierror)
         rzprod1 = rbuf

c     set p equal to z
         do j=2,locndimy+1
            do i=2,locndimx+1
               p(i,j) = z(i,j)
            enddo
         enddo
      else
         rzprod2 = rzprod1

c     calculate beta = (r,z) / (r_old,z_old)
         rzprod1 = 0.
         do j=2,locndimy+1
            do i=2,locndimx+1
               rzprod1 = rzprod1 + r(i,j)*z(i,j)
            enddo
         enddo
         call MPI_AllReduce(rzprod1, rbuf, 1, MPI_REAL8, MPI_SUM, 
     &        MPI_COMM_WORLD,ierror)
         rzprod1 = rbuf
         beta = rzprod1/rzprod2

c     set p = z + beta*p
         do j=2,locndimy+1
            do i=2,locndimx+1
               p(i,j) = z(i,j) + beta*p(i,j)
            enddo
         enddo
      endif


c     update the ghost boundaries of p
      call MSG_tbdx_send(p,buf,NumAdjProc,Proc,Ipr,Index,MSGptrn1
     &     ,ierror)
      if (ierror.ne.0) write(*,*) MyProc,'-a:',ierror
      
      call MSG_tbdx_receive(p,buf,NumAdjProc,Proc,Ipr,Index,MSGptrn1
     &     ,ierror)
      if (ierror.ne.0) write(*,*) MyProc,'-b:',ierror
      

c     calculate w = A p
      do j=2,locndimy+1
         do i=2,locndimx+1
            w(i,j) = s(i,j,k0)*p(i,j) +s(i,j,kw)     *p(i-1,j)
     &                                +s(i+1,j,kw)   *p(i+1,j)
     &                                +s(i,j,ks)     *p(i,j-1)
     &                                +s(i,j+1,ks)   *p(i,j+1)
     &                                +s(i,j,ksw)    *p(i-1,j-1)
     &                                +s(i,j+1,knw)  *p(i-1,j+1)
     &                                +s(i+1,j+1,ksw)*p(i+1,j+1)
     &                                +s(i+1,j,knw)  *p(i+1,j-1)
         enddo
      enddo

c     calculate alpha = (r,z)/(p,q)
      alpha = 0.0
      do j=2,locndimy+1
         do i=2,locndimx+1
            alpha = alpha+p(i,j)*w(i,j)
         enddo
      enddo
      call MPI_AllReduce(alpha, rbuf, 1, MPI_REAL8, MPI_SUM, 
     &                   MPI_COMM_WORLD,ierror)
      alpha = rbuf
      alpha = rzprod1/alpha
      

c     calculate x = x + alpha*p
c     and       r = r - alpha*w
      do j=2,locndimy+1
         do i=2,locndimx+1
            x(i,j) = x(i,j) + alpha*p(i,j)
            r(i,j) = r(i,j) - alpha*w(i,j)
         enddo
      enddo

c     some stopping criterion
      if ((rzprod1 .le. 1.0e-3*bnorm).or.(k.ge.500)) goto 999
      goto 111


c     report the results
 999  if (MyProc .eq. 0) then 
         write(*,*)'number of iterations = ',k
         write(*,*)'norm of residual     = ',rzprod1
      endif

c     finalize things and exit
      call MSG_DISABLE(ierror)
      call MPI_FINALIZE(ierror)

      end





