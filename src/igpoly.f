c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 

      subroutine igpoly(nx, xg,ny,yg,np,xp,yp,ind)
!----------------------------------------------------------------------
!     This subroutine determines whether or not an 2-d point (xg(i),yg(j))
!     element is inside a (closed) set of points xp,yp that
!     are assumed to form polygon.
!     result is an  matrix indicating comparision for all 
!     combinations of xg and yg
!----------------------------------------------------------------------

      integer np      		! # of points in polygon
      integer nx,ny		! # points to check 
      real xg(nx) 		! x grid values to check
      real yg(ny)               ! y grid values to check 
      real    xp(np)		! 2d-locations of polygon
      real    yp(np)		
      real  x1, x2, y1,y2       ! min and max of x and y
      real  temp, xt, yt
      integer ind(nx,ny)      ! THE ANSWER : ind(i)=1 if point xd(i),yd(i) is 
                                !  in polygon 0 otherwise 
      integer in
      x1= xp(1)
      x2= xp(2)
      y1= yp(1)
      y2= yp(1)
!
!     find the minima and maxima of the polygon coordinates
!     i.e. the smallest rectangle containing the polygon.
      do j = 1,np

        temp= xp(j)
        if( temp.lt.x1)  then 
           x1 = temp
        endif
        if( temp.gt.x2) then
           x2 = temp
        endif

        temp= yp(j)
        if( temp.lt.y1)  then 
           y1 = temp
        endif
        if( temp.gt.y2) then
           y2 = temp
        endif
      enddo

! loop over all 
! grid points

      do i = 1, nx
       do j = 1,ny
        xt= xg(i)
        yt= yg(j)
! quick test that point is inside the bounding rectangle
! of the polygon

        if( (xt.le. x2).and. (xt.ge.x1).and.
     *       (yt.le.y2).and.(yt.ge.y1) ) then
         call inpoly2(xt,yt,np,xp,yp, in) 
         ind(i,j)=in
        else
         ind(i,j)= 0
        endif

       enddo
      enddo

      return 
      end
