c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
       
      subroutine dmaket(m,n,dim,des,lddes,npoly,t,ldt,
     * wptr,info,ptab,ldptab)
      integer m,n,dim,lddes,npoly,ldt,wptr(dim),info,ptab(ldptab,dim)
      double precision des(lddes,dim),t(ldt,npoly)
c
c Purpose: create t matrix and append s1 to it.
c
c On Entry:
c   m			order of the derivatives in the penalty
c   n			number of rows in des
c   dim			number of columns in des
c   des(lddes,dim)	variables to be splined
c   lddes		leading dimension of des as declared in the
c			calling program
c   ldt			leading dimension of t as declared in the
c			calling program
c
c   npoly		dimension of polynomial part of spline
c On Exit:
c   t(ldt,npoly+ncov1)	[t:s1]
c   info 		error indication
c   			   0 : successful completion
c		 	   1 : error in creation of t
c Work Arrays:
c   wptr(dim)		integer work vector
c
c Subprograms Called Directly:
c	Other - mkpoly
c
c
      integer i,j,k,kk,tt,nt,bptr,eptr
c
      info = 0
c      npoly = mkpoly(m,dim)
      do 5 j=1,n
        t(j,1)=1.0
 5       continue
      nt = 1
      if (npoly .gt. 1) then
          do 10 j=1,dim
             nt = j + 1
             wptr(j) = nt
             ptab(nt,j)= ptab(nt,j) +1
             do 15 kk = 1, n
                t(kk,nt)= des(kk,j)
   15        continue     
c             call dcopy(n,des(1,j),1,t(1,nt),1)
   10     continue
c
c     get cross products of x's in null space for m>2
c
c     WARNING: do NOT change next do loop unless you fully understand:
c              This first gets x1*x1, x1*x2, x1*x3, then
c              x2*x2, x2*x3, and finally x3*x3 for dim=3,n=3
c              wptr(1) is always at the beginning of the current
c	       level of cross products, hence the end of the
c	       previous level which is used for the next.
c	       wptr(j) is at the start of xj * (previous level)
c
          do 50 k=2,m-1
             do 40 j=1,dim
                bptr = wptr(j)
                wptr(j) = nt + 1
                eptr = wptr(1) - 1
                do 30 tt=bptr,eptr
                   nt = nt + 1
                        do 21 jj= 1,dim 
                        ptab(nt,jj)= ptab(tt,jj)
 21                     continue
                   ptab( nt,j)= 1+ ptab( nt,j)
                   do 20 i=1,n
                      t(i,nt) = des(i,j) * t(i,tt)
   20              continue
   30           continue
   40        continue
   50     continue
          if (nt .ne. npoly) then
	      info = 1
	      return
          endif
      endif
c
      end
