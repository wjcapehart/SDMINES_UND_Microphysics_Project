      subroutine JacobiPre (s,r,z,locndimx,locndimy,locndimz)
      implicit none
      
      include 'cg.h'

      real*8 s(maxndimx+2, maxndimy+2, maxndimz+2, MaxStencilSize)
      real*8 r(maxndimx+2,maxndimy+2,maxndimz+2)
      real*8 z(maxndimx+2,maxndimy+2,maxndimz+2)
      integer locndimx,locndimy,locndimz

      integer i,j,k

c     this is Jacobi preconditioning

      do k=2,locndimz+1
         do j=2,locndimz+1
            do i=2,locndimx+1
               z(i,j,k) = r(i,j,k)/s(i,j,k,kp)
c     &              +(s(i,j,k+1,kb) *r(i,j,k+1)
c     &              + s(i,j,k,kb)   *r(i,j,k-1)
c     &              + s(i,j+1,k,kps)*r(i,j+1,k)
c     &              + s(i,j,k,kps)  *r(i,j-1,k)
c     &              + s(i+1,j,k,kpw)*r(i+1,j,k)
c     &              + s(i,j,k,kpw)  *r(i-1,j,k))/s(i,j,k,kp)**2

               
            enddo
         enddo
      enddo

      return
      end
