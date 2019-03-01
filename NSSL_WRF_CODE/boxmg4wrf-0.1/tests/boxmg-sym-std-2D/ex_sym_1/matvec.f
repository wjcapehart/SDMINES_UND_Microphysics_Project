      subroutine matvec (SO, u, bu, L, M)

      implicit none

      integer   ko, kw, ks
      parameter (ko=1, kw=2, ks=3)
      
      integer L,M
      real*8 u(0:L+1,0:M+1), bu(0:L+1,0:M+1), SO(0:L+1,0:M+1,3)

      integer i,j


      do j=1,M
         do i=1,L
            bu(i,j) = so(i,j,ko)*u(i,j)
     &               -so(i,j,kw)*u(i-1,j)
     &               -so(i,j,ks)*u(i,j-1)
     &               -so(i+1,j,kw)*u(i+1,j)
     &               -so(i,j+1,ks)*u(i,j+1)
         enddo
      enddo

      return
      end
