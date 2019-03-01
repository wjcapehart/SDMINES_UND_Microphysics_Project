      program conjgrad3D
      implicit none
      include 'mpif.h'
      include 'MSG.h'
      include 'cg.h'
      integer ndim, i, j, k, kk, ierror
      integer locndimx,locndimy,locndimz
      integer MSGMyProc, NumProc, MyProc
      integer MyProci, MyProcj, MyProck
      integer Proc(27), Ipr(55), NumAdjProc,ldim,ptrn
      integer DirNumProc

      integer Index
     &     (2*6*((maxndimx+2)**2+(maxndimy+2)**2+(maxndimz+2)**2))
      integer LocalArraySize(3,MaxNumProc)
      integer ActDatastart(3,MaxNumProc)
      integer GlobalCoordLocalData(2,3,MaxNumProc)
      integer GlobalCoordActData(2,3,MaxNumProc)
      real*8  itime,etime


      real*8 r(maxndimx+2,maxndimy+2,maxndimz+2)
      real*8 x(maxndimx+2,maxndimy+2,maxndimz+2)
      real*8 b(maxndimx+2,maxndimy+2,maxndimz+2)
      real*8 p(maxndimx+2,maxndimy+2,maxndimz+2)
      real*8 w(maxndimx+2,maxndimy+2,maxndimz+2)
      real*8 z(maxndimx+2,maxndimy+2,maxndimz+2)

      real*8 buf
     &     (6*2*((maxndimx+2)**2+(maxndimy+2)**2+(maxndimz+2)**2))

      real*8 h, beta, alpha
      real*8 rzprod2, rzprod1, rzprod
      real*8 rbuf,bnorm

      real*8 s(maxndimx+2,maxndimy+2,maxndimz+2,MaxStencilSize)
      real*8 work(maxndimx+2,maxndimy+2,maxndimz+2,2)

      real*8 eps

      call MPI_INIT(ierror)

      call MSG_ENABLE(MSGMyProc, NumProc)
      MyProc = MSGMyProc-1
      ptrn=1
c      write(*,*) MyProc,' is alive'
      

      DirNumProc = (real(NumProc))**(1.0/3)
      if (DirNumProc**3 .ne. NumProc) then
         write(*,*)NumProc,' is not cube'
         goto 999
      endif


      MyProck = MyProc/DirNumProc**2 
      MyProcj = (MyProc- MyProck*DirNumProc**2)/DirNumProc 
      MyProci = MyProc - MyProcj*DirNumProc - MyProck*(DirNumProc**2)

      MyProci = MyProci+1
      MyProcj = MyProcj+1
      MyProck = MyProck+1

c      write(*,*) MyProc,': ',DirNumProc, MyProci, MyProcj,MyProck
           

c     calculate the geometry of the variables
      if (MyProc .eq. 0) then
         write(*,*)'enter number of gridlines '
         read(*,*) ndim
         write(*,*)'eps '
         read(*,*) eps
      endif
      call MPI_BCAST(ndim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(eps,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

c      write(*,*) ndim, eps
      call MPI_barrier(MPI_COMM_WORLD,ierror)
      itime = MPI_Wtime()


      locndimx = ndim / DirNumProc
      if (MyProci .eq. DirNumProc) then
         locndimx = locndimx+(ndim-locndimx*DirNumProc)
      endif

      locndimy = ndim / DirNumProc
      if (MyProcj .eq. DirNumProc) then
         locndimy = locndimy+(ndim-locndimy*DirNumProc)
      endif

      locndimz = ndim / DirNumProc
      if (MyProck .eq. DirNumProc) then
         locndimz = locndimz+(ndim-locndimz*DirNumProc)
      endif
      
c      write(*,*) MSGMyProc,locndimx,locndimy,locndimz
      
c      write(*,*) MyProc,': ', MyProci,MyProcj,MyProck 

c     initialize the stencil matrix

      do k=2,locndimz+1
         do j=2,locndimy+1
            do i=2,locndimx+1
               s(i,j,k,kp)=2.0+4.0*eps
               s(i,j,k,kpw)=-eps
               s(i,j,k,kps)=-eps
               s(i,j,k,kb)=-1.0
            enddo
         enddo
      enddo

      if (MyProci.eq. DirNumProc) then
         do k=2,locndimz+1
            do j=2,locndimy+1
               s(locndimx+2,j,k,kp)=2.0+4.0*eps
               s(locndimx+2,j,k,kpw)=-eps
               s(locndimx+2,j,k,kps)=-eps
               s(locndimx+2,j,k,kb)=-1.0
            enddo
         enddo
      endif

      if (MyProcj.eq. DirNumProc) then
         do k=2,locndimz+1
            do i=2,locndimx+1
               s(i,locndimy+2,k,kp)=2.0+4.0*eps
               s(i,locndimy+2,k,kpw)=-eps
               s(i,locndimy+2,k,kps)=-eps
               s(i,locndimy+2,k,kb)=-1.0
            enddo
         enddo
      endif

      if (MyProck.eq. DirNumProc) then
         do j=2,locndimy+1
            do i=2,locndimx+1
               s(i,j,locndimz+2,kp)=2.0+4.0*eps
               s(i,j,locndimz+2,kpw)=-eps
               s(i,j,locndimz+2,kps)=-eps
               s(i,j,locndimz+2,kb)=-1.0
            enddo
         enddo
      endif

      if ((MyProci.eq.DirNumProc).and.(MyProcj.eq.DirNumProc)) then
         do k=2,locndimz+1
            s(locndimx+2,locndimy+2,k,kp)=2.0+4.0*eps
            s(locndimx+2,locndimy+2,k,kpw)=-eps
            s(locndimx+2,locndimy+2,k,kps)=-eps
            s(locndimx+2,locndimy+2,k,kb)=-1.0
         enddo
      endif
      
      if ((MyProcj.eq.DirNumProc).and.(MyProck.eq.DirNumProc)) then
         do i=2,locndimx+1
            s(i,locndimy+2,locndimz+2,kp)=2.0+4.0*eps
            s(i,locndimy+2,locndimz+2,kpw)=-eps
            s(i,locndimy+2,locndimz+2,kps)=-eps
            s(i,locndimy+2,locndimz+2,kb)=-1.0
         enddo
      endif
       
      if ((MyProci.eq.DirNumProc).and.(MyProck.eq.DirNumProc)) then
         do j=2,locndimy+1
            s(locndimx+2,j,locndimz+2,kp)=2.0+4.0*eps
            s(locndimx+2,j,locndimz+2,kpw)=-eps
            s(locndimx+2,j,locndimz+2,kps)=-eps
            s(locndimx+2,j,locndimz+2,kb)=-1.0
         enddo
      endif
 
      if ( (MyProci .eq. DirNumProc) .and.
     &     (MyProcj .eq. DirNumProc)  .and.
     &     (MyProck .eq. DirNumProc) ) then
         s(locndimx+2,locndimy+2,locndimz+2,kp)=2.0+4.0*eps
         s(locndimx+2,locndimy+2,locndimz+2,kpw)=-eps
         s(locndimx+2,locndimy+2,locndimz+2,kps)=-eps
         s(locndimx+2,locndimy+2,locndimz+2,kb)=-1.0
      endif

      h=1.0/ndim

c     initialize right hand side

      do k=2, locndimz+1
         do j=2, locndimy+1
            do i=2, locndimx+1
               b(i,j,k) = h**2
            enddo
         enddo
      enddo
      

c     set up variables for CG iteration

      kk=0
      
      do k=1, locndimz+2
         do j=1,locndimy+2
            do i=1, locndimx+2
               x(i,j,k) = 0.0
            enddo
         enddo
      enddo
      
c     zero the variable r
      do k=1,locndimz+2
         do j=1,locndimy+2
            do i=1,locndimx+2
               r(i,j,k) = 0.0
            enddo
         enddo
      enddo

c     calculate the initial residual
      do k=2,locndimz+1
         do j=2,locndimy+1
            do i=2,locndimx+1
               r(i,j,k) = b(i,j,k)
            enddo
         enddo
      enddo

      bnorm = 0.0
      do k=2,locndimz+1
         do j=2,locndimy+1
            do i=2,locndimx+1
               bnorm = bnorm + b(i,j,k)**2
            enddo
         enddo
      enddo

      call MPI_AllReduce(bnorm, rbuf, 1, MPI_REAL8, MPI_SUM, 
     &     MPI_COMM_WORLD,ierror)
      bnorm = rbuf


c
c     set the all the neccessary data for MSG
c
      NumAdjProc = 0      
      ldim = ndim/DirNumProc

      do k=1,DirNumProc
      do j=1,DirNumProc
      do i=1,DirNumProc
         LocalArraySize(1,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)
     &        =maxndimx+2
         LocalArraySize(2,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)
     &        =maxndimy+2
         LocalArraySize(3,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)
     &        =maxndimz+2
         
         ActDataStart(1,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=1
         ActDataStart(2,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=1
         ActDataStart(3,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=1
         
         GlobalCoordLocalData
     &        (1,1,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)
     &        =(i-1)*ldim+2
         GlobalCoordLocalData
     &        (1,2,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)
     &        =(j-1)*ldim+2
         GlobalCoordLocalData
     &        (1,3,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)
     &        =(k-1)*ldim+2
         
         if (i.ne.DirNumProc) then
            GlobalCoordLocalData
     &           (2,1,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=i*ldim+1
         else
            GlobalCoordLocalData
     &           (2,1,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=ndim+1
         endif
         if (j.ne.DirNumProc) then
            GlobalCoordLocalData
     &           (2,2,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=j*ldim+1
         else
            GlobalCoordLocalData
     &           (2,2,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=ndim+1
         endif
         if (k.ne.DirNumProc) then
            GlobalCoordLocalData
     &           (2,3,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=k*ldim+1
         else
            GlobalCoordLocalData
     &           (2,3,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=ndim+1
         endif
            
         GlobalCoordActData(1,1,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)
     &        =(i-1)*ldim+1
         GlobalCoordActData(1,2,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)
     &        =(j-1)*ldim+1
         GlobalCoordActData(1,3,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)
     &        =(k-1)*ldim+1
         
         if (i.eq.DirNumProc) then
            GlobalCoordActData
     &           (2,1,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=ndim+2
         else
            GlobalCoordActData
     &           (2,1,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=i*ldim+2
         endif
         if (j.eq.DirNumProc) then
            GlobalCoordActData
     &           (2,2,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=ndim+2
         else
            GlobalCoordActData
     &           (2,2,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=j*ldim+2
         endif
         if (k.eq.DirNumProc) then
            GlobalCoordActData
     &           (2,3,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=ndim+2
         else
            GlobalCoordActData
     &           (2,3,(k-1)*DirNumProc**2+(j-1)*DirNumProc+i)=k*ldim+2
         endif
      enddo
      enddo
      enddo

      call MSG_tp_setup(LocalArraySize, ActDataStart,
     &                  GlobalCoordLocalData, GlobalCoordActData,
     &                  NumProc, MSGMyProc, NumAdjProc, Proc, Ipr,
     &                  Index, 1, 1, ierror)  

 
c
c     distribute the stencil to all neighbor processors
c
      call DistributeStencilHalo3D (s,GlobalCoordLocalData, 
     &     GlobalCoordActData, ActDataStart, LocalArraySize, 
     &     MaxStencilSize, MSGMyProc, NumProc, 1,ierror)




c     this is the CG iteration
 111  continue

      call LineZPre(s,r,z,work,locndimx,locndimy,locndimz)
c      call JacobiPre(s,r,z,locndimx,locndimy,locndimz)
c      call IdPre(s,r,z,locndimx,locndimy,locndimz)

      kk = kk+1

      if (kk .eq. 1) then
         rzprod1 = 0.
         do k=2,locndimz+1
            do j=2,locndimy+1
               do i=2,locndimx+1
                  rzprod1 = rzprod1 + r(i,j,k)*z(i,j,k)
               enddo
            enddo
         enddo
         call MPI_AllReduce(rzprod1, rbuf, 1, MPI_REAL8, MPI_SUM, 
     &        MPI_COMM_WORLD,ierror)
         rzprod1 = rbuf
         
         do k=2,locndimz+1
            do j=2,locndimy+1
               do i=2,locndimx+1
                  p(i,j,k) = z(i,j,k)
               enddo
            enddo
         enddo
      else
         rzprod2 = rzprod1
         rzprod1 = 0.
         do k=2,locndimz+1
            do j=2,locndimy+1
               do i=2,locndimx+1
                  rzprod1 = rzprod1 + r(i,j,k)*z(i,j,k)
               enddo
            enddo
         enddo
         call MPI_AllReduce(rzprod1, rbuf, 1, MPI_REAL8, MPI_SUM, 
     &        MPI_COMM_WORLD,ierror)
         rzprod1 = rbuf

         beta = rzprod1/rzprod2
         do k=2,locndimz+1
            do j=2,locndimy+1
               do i=2,locndimx+1
                  p(i,j,k) = z(i,j,k) + beta*p(i,j,k)
               enddo
            enddo
         enddo
      endif

      ptrn=1
      call MSG_tbdx_send(p,buf,NumAdjProc,Proc,Ipr,Index,ptrn,ierror)
      if (ierror.ne.0) then
         write(*,*) MyProc,'-a:',ierror
         goto 999
      endif
      
      ptrn=1
      call MSG_tbdx_receive(p,buf,NumAdjProc,Proc,Ipr,Index,ptrn,ierror)
      if (ierror.ne.0) then
         write(*,*) MyProc,'-b:',ierror
         goto 999
      endif
      
      do k=2,locndimz+1
         do j=2,locndimy+1
            do i=2,locndimx+1
               w(i,j,k) = s(i,j,k,kp)*p(i,j,k)
     &              +s(i,j,k,kpw)  *p(i-1,j,k)
     &              +s(i+1,j,k,kpw)*p(i+1,j,k)
     &              +s(i,j,k,kps)  *p(i,j-1,k)
     &              +s(i,j+1,k,kps)*p(i,j+1,k)
     &              +s(i,j,k,kb)   *p(i,j,k-1)
     &              +s(i,j,k+1,kb) *p(i,j,k+1)
            enddo
         enddo
      enddo


      alpha = 0.0
      do k=2,locndimz+1
         do j=2,locndimy+1
            do i=2,locndimx+1
               alpha = alpha+p(i,j,k)*w(i,j,k)
            enddo
         enddo
      enddo
      call MPI_AllReduce(alpha, rbuf, 1, MPI_REAL8, MPI_SUM, 
     &     MPI_COMM_WORLD,ierror)
      alpha = rbuf
      alpha = rzprod1/alpha
      
      do k=2,locndimz+1
         do j=2,locndimy+1
            do i=2,locndimx+1
               x(i,j,k) = x(i,j,k) + alpha*p(i,j,k)
               r(i,j,k) = r(i,j,k) - alpha*w(i,j,k)
            enddo
         enddo
      enddo

      if ((rzprod1 .le. 1.0e-3*bnorm).or.(kk.ge.1000)) goto 999
      goto 111

      
 999  call MPI_barrier(MPI_COMM_WORLD,ierror)
      etime = MPI_Wtime()


      if (MyProc .eq. 0) then 
         write(*,*)'number of iterations = ',kk
         write(*,*)'norm of residual     = ',rzprod1
         write(*,*)'time elapsed         = ',etime-itime
      endif
      call MSG_DISABLE(ierror)
      call MPI_FINALIZE(ierror)

      end


