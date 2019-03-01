      subroutine vecprod (prod,u,v,L,M)
      
      implicit none

      integer L,M
      real*8 u(0:L+1,0:M+1), v(0:L+1,0:M+1), prod

      integer i,j
   

      prod = 0.0
      do j=1,M
         do i=1,L
            prod = prod+u(i,j)*v(i,j)
         enddo
      enddo

      return
      end
