      subroutine vecprod (prod,u,v,ii,jj,kk)
      
      implicit none

      integer ii,jj,kk
      real*8 u(0:ii-1,0:jj-1,0:kk-1), v(0:ii-1,0:jj-1,0:kk-1),
     &     prod

      integer i,j,k
   

      prod = 0.0
      do j=1,jj-2
         do i=1,ii-2
            do k=1,kk-2
               prod = prod+u(i,j,k)*v(i,j,k)
            enddo
         enddo
      enddo
      
      return
      end
