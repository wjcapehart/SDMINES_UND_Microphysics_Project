      subroutine putf(so,qf,l,m,deltax,deltay,k)

      implicit  none

      integer   ko, kw, ks
      parameter (ko=1, kw=2, ks=3)
c
      integer   l, m, k, i, j
      real*8    so(0:l+1,0:m+1,3), qf(0:l+1,0:m+1), deltax, deltay

c

      do 30 j = 1, m
        do 30 i = 1, l
          so(i,j,ks) = 1
 30   continue

      do 40 j = 1, m
        do 40 i = 1, l
          so(i,j,kw) = 1
 40   continue

      do 50 j = 1, m
        do 50 i = 1, l
          so(i,j,ko) =  4 
          qf(i,j) = deltax*deltay
 50   continue
c
      return
      end
