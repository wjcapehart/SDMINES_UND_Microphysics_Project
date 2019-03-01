      subroutine LineZPre (s,r,z,work,locndimx,locndimy,locndimz)
      implicit none
      
      include 'cg.h'

      real*8 s(maxndimx+2, maxndimy+2, maxndimz+2, MaxStencilSize)
      real*8 r(maxndimx+2,maxndimy+2,maxndimz+2)
      real*8 z(maxndimx+2,maxndimy+2,maxndimz+2)
      real*8 work(maxndimx+2,maxndimy+2,maxndimz+2,2)
      real*8 y(maxndimz+2)
      integer locndimx,locndimy,locndimz
      
      integer i,j,k
      logical firstpass 
      save firstpass
      data firstpass /.true./


c      write(*,*) ' in LinZPre ', locndimx, locndimy, locndimz

      if (firstpass .eqv. .true.) then
 
c     LU decomposition of the lines
c     notation is according to Atkinson, Numerical Analysis, pp 527
c           work(i,j,k,1) ... alpha(k)
c           work(i,j,k,2) ... gamma(k)
c           s(i,j,k,kb)   ... b(k)
c           s(i,j,k+1,kb) ... c(k)
c           s(i,j,k,kp)   ... a(k)

         do j=2,locndimy+1
            do i=2,locndimx+1
c              interior points start at k = 2
               work(i,j,2,1) = s(i,j,2,kp)
               work(i,j,2,2) = s(i,j,2+1,kb)/work(i,j,2,1)
               do k=3,locndimz
                  work(i,j,k,1) = s(i,j,k,kp) 
     &                 - s(i,j,k,kb)*work(i,j,k-1,2)
                  work(i,j,k,2) = s(i,j,k+1,kb)/work(i,j,k,1)
               enddo
               work(i,j,locndimz+1,1) = s(i,j,locndimz+1,kp)
     &              - s(i,j,locndimz+1,kb)*work(i,j,locndimz,2)
            enddo
         enddo
         firstpass = .false.
      endif

c     forward elimination
      do j=2,locndimy+1
         do i=2,locndimx+1
            y(2) = r(i,j,2)/work(i,j,2,1)
            do k=3,locndimz+1
               y(k) = (r(i,j,k) - s(i,j,k,kb)*y(k-1))
     &               / work(i,j,k,1)
            enddo

c     backsubstitution
            z(i,j,locndimz+1) = y(locndimz+1)
            do k=locndimz,2,-1
               z(i,j,k) = y(k) - work(i,j,k,2)*z(i,j,k+1)
            enddo
         enddo
      enddo


      return
      end
