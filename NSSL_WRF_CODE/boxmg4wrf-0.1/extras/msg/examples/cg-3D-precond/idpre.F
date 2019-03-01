      subroutine IdPre (s,r,z,locndimx,locndimy,locndimz)
      implicit none
      
      include 'cg.h'

      real*8 s(maxndimx+2, maxndimy+2, maxndimz+2, MaxStencilSize)
      real*8 r(maxndimx+2,maxndimy+2,maxndimz+2)
      real*8 z(maxndimx+2,maxndimy+2,maxndimz+2)
      integer locndimx,locndimy,locndimz

      integer i,j,k

      do k=2,locndimz+1
         do j=2,locndimz+1
            do i=2,locndimx+1
               z(i,j,k) = r(i,j,k)
            enddo
         enddo
      enddo

      return
      end
