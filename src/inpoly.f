c fields, Tools for spatial data
c Copyright (C) 2017, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 

      subroutine inpoly(nd, xd,yd,np,xp,yp,ind)
!----------------------------------------------------------------------
!     This subroutine determines whether or not an 2-d point (xd(j),yd(j))
!     element is inside a (closed) set of points xp,yp that
!     are assumed to form polygon.
!----------------------------------------------------------------------

      integer np
!     # of points in polygon
      integer nd
!     # points to check 
      real xd(nd)
!     2d-locations to check
      real yd(nd)
      real    xp(np)
!     2d-locations of polygon
      real  yp(np)
      real  x1, x2, y1,y2
!     min and max of x and y
      real  temp, xt, yt
      integer ind(nd)
!     THE ANSWER : ind(i)=1 if point xd(i),yd(i) is 
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
    
      do j = 1,nd
      xt= xd(j)
      yt= yd(j)
! quick test that point is inside the bounding rectangle
! if not it is not inside polygon
      if( (xt.le. x2).and. (xt.ge.x1).and.
     *       (yt.le.y2).and.(yt.ge.y1) ) then
        call inpoly2(xt,yt,np,xp,yp, in) 
        ind(j)=in
      else
        ind(j)= 0
      endif

      enddo

      return 
      end

      subroutine inpoly2(xpnt,ypnt,np,xp,yp,in)
C
      parameter (pi=3.14159265358979,ttpi=2.*pi)
C
      dimension xp(np),yp(np)
      real xpnt, ypnt
      integer in
C
C----------------------------------------------------------------------
C
C  THE VALUE OF THIS FUNCTION IS NON-ZERO IF AND ONLY IF (XPNT,YPNT) 
C  IS INSIDE OR *ON* THE POLYGON DEFINED BY THE POINTS (XP(I),YP(I)), 
C  FOR I FROM 1 TO NP.
C
C  THE INPUT POLYGON DOES NOT HAVE TO BE A CLOSED POLYGON, THAT IS
C  IT DOES NOT HAVE TO BE THAT (XP(1),YP(1) = (XP(NP,YP(NP)).
C
C----------------------------------------------------------------------
C
C  DETERMINE THE NUMBER OF POINTS TO LOOK AT (DEPENDING ON WHETHER THE
C  CALLER MADE THE LAST POINT A DUPLICATE OF THE FIRST OR NOT).
C
      if (xp(np).eq.xp(1) .and. yp(np).eq.yp(1)) then
        npts = np-1
      else
        npts = np
      end if

      in = 0
!     ASSUME POINT IS OUTSIDE

C --- ------------------------------------------------------------------
C --- CHECK TO SEE IF THE POINT IS ON THE POLYGON.
C --- ------------------------------------------------------------------

      do ipnt = 1,npts
         if (xpnt .eq. xp(ipnt) .and. ypnt .eq. yp(ipnt) ) then
         in = 1
      goto 999
! EARLY EXIT
         endif
      enddo

C --- ------------------------------------------------------------------
C --- COMPUTE THE TOTAL ANGULAR CHANGE DESCRIBED BY A RAY EMANATING 
C --- FROM THE POINT (XPNT,YPNT) AND PASSING THROUGH A POINT THAT 
C --- MOVES AROUND THE POLYGON.
C --- ------------------------------------------------------------------

      anch = 0.

      inxt = npts
      xnxt = xp(npts)
      ynxt = yp(npts)
      anxt = atan2(ynxt-ypnt, xnxt-xpnt)

      do 100 ipnt=1,npts

         ilst = inxt
         xlst = xnxt
         ylst = ynxt
         alst = anxt

         inxt = ipnt
         xnxt = xp(inxt)
         ynxt = yp(inxt)
         anxt = atan2(ynxt-ypnt, xnxt-xpnt)

         adif = anxt-alst

         if (abs(adif) .gt. pi) adif = adif - sign(ttpi,adif)

         anch = anch + adif

 100  continue

C --- ------------------------------------------------------------------
C --- IF THE POINT IS OUTSIDE THE POLYGON, THE TOTAL ANGULAR CHANGE 
C --- SHOULD BE EXACTLY ZERO, WHILE IF THE POINT IS INSIDE THE POLYGON,
C --- THE TOTAL ANGULAR CHANGE SHOULD BE EXACTLY PLUS OR MINUS TWO PI.
C --- WE JUST TEST FOR THE ABSOLUTE VALUE OF THE CHANGE BEING LESS 
C --- THAN OR EQUAL TO PI.
C --- ------------------------------------------------------------------

      if (abs(anch) .ge. pi) in = 1

 999  continue
      return
      end
