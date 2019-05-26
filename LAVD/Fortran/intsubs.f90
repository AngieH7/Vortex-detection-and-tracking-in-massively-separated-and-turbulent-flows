module intsubs 
   use grid
   use inout
   implicit none
   
contains
   
   subroutine setup_traj(tstart,traj)
      real, dimension(ox,oy,2) :: traj
      integer :: tstart
      integer :: i,j,k, status
      integer*4 :: hi, hi2, hi3, hi4
      character(40) :: fname
      integer :: fp=62

      if(freshstart) then

         startstep = -1

         do i=1,ox
            do j=1,oy
               
               traj(i,j,1) = x1(i)
               traj(i,j,2) = y1(j)
                             
            end do
         end do

      else  

         if(dir.lt.0) then
            write(fname,'(a,i6.6,a,i6.6,a)') trim(odir) // trim(oname) // 'traj_neg_',tstart,'_',startstep,'.dat'
         else
            write(fname,'(a,i6.6,a,i6.6,a)') trim(odir) // trim(oname) // 'traj_pos_',tstart,'_',startstep,'.dat'
         end if
         
         print *,"Reading trajectory file ",trim(fname)
         open(unit=fp, file=trim(fname), status="old", form="unformatted", &
              iostat=status)
         !  print *, 'file opened'
         if (status .ne. 0) then
            print *,"Error: could not open "//trim(fname)//" for reading.", status
            stop
         end if
         
         read(fp) hi, hi2, hi4
         read(fp) ((traj(i,j,1),i=1,ox),j=1,oy), &
              ((traj(i,j,2),i=1,ox),j=1,oy)
         close(fp)
         
      end if

  end subroutine setup_traj
  
  subroutine step_vels_dis(t,vel,traj)
    integer, intent(in) :: t
    real, dimension(ox,oy,2), intent(inout) :: vel
    real, dimension(ox,oy,2), intent(in) :: traj
    real, dimension(2) :: vect
    integer, dimension(2) :: spos
    integer :: i,j,m

    do i=1,ox
       do j=1,oy

             spos = getbound(traj(i,j,:))
             
             vel(i,j,1) = interp2d(spos, traj(i,j,:), getsquare(spos, velt(:,:,1)))
             vel(i,j,2) = interp2d(spos, traj(i,j,:), getsquare(spos, velt(:,:,2)))

             if(traj(i,j,1).gt.x(nx)) then
                vel(i,j,1) = u(nx,ny)
             end if

       end do
    end do
    
  end subroutine step_vels_dis
  
  subroutine step_vor_dis(t,vor,traj)
    integer, intent(in) :: t
    real, dimension(ox,oy), intent(inout) :: vor
    real, dimension(ox,oy,2), intent(in) :: traj
    integer, dimension(2) :: spos
    integer :: i,j,m

    do i=1,ox
       do j=1,oy

             spos = getbound(traj(i,j,:))
             
             vor(i,j) = interp2d(spos, traj(i,j,:), getsquare(spos, vort(:,:)))

             if(traj(i,j,1).gt.x(nx)) then
                vor(i,j) = vort(nx,ny)
             end if

       end do
    end do
    
  end subroutine step_vor_dis

  subroutine step_traj_dis(t,tjump,vel,traj)
    integer, intent(in) :: t,tjump
    real, dimension(ox,oy,2), intent(in) :: vel
    real, dimension(ox,oy,2), intent(inout) :: traj
    real*4, dimension(:,:,:), allocatable :: ev
    real*4, dimension(:,:,:), allocatable :: mv
    real, dimension(2) :: loc
    integer, dimension(2) :: spos
    integer :: i,j,m, igrd
    real :: k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y

    allocate(ev(nx,ny,2), mv(nx,ny,2))
    
    do i=1,nx
       do j=1,ny
          mv(i,j,1) = interp4(tjump*idata*delt, (tjump+dir)*idata*delt, & 
               u1(i,j), u(i,j), (t+.5)*dir*delta)
          mv(i,j,2) = interp4(tjump*idata*delt, (tjump+dir)*idata*delt, & 
               v1(i,j), v(i,j), (t+.5)*dir*delta)
          
          ev(i,j,1) = interp4(tjump*idata*delt, (tjump+dir)*idata*delt, & 
               u1(i,j), u(i,j), (t+1)*dir*delta)
          ev(i,j,2) = interp4(tjump*idata*delt, (tjump+dir)*idata*delt, & 
               v1(i,j), v(i,j), (t+1)*dir*delta)
       end do
    end do

    
    ! RK 4
    do i=1,ox
       do j=1,oy
          
          k1x = dir*vel(i,j,1)*delta
          k1y = dir*vel(i,j,2)*delta
          
          loc(1) = traj(i,j,1)+(k1x/2)
          loc(2) = traj(i,j,2)+(k1y/2)
          
          spos = getbound(loc)
          
          k2x = dir*delta*interp2d(spos, loc, getsquare(spos, mv(:,:,1)))
          k2y = dir*delta*interp2d(spos, loc, getsquare(spos, mv(:,:,2)))
          
          loc(1) = traj(i,j,1)+(k2x/2)
          loc(2) = traj(i,j,2)+(k2y/2)
          
          spos = getbound(loc)
          
          k3x = dir*delta*interp2d(spos, loc, getsquare(spos, mv(:,:,1)))
          k3y = dir*delta*interp2d(spos, loc, getsquare(spos, mv(:,:,2)))
          
          loc(1) = traj(i,j,1)+(k3x)
          loc(2) = traj(i,j,2)+(k3y)
          
          spos = getbound(loc)
          
          k4x = dir*delta*interp2d(spos, loc, getsquare(spos, ev(:,:,1)))
          k4y = dir*delta*interp2d(spos, loc, getsquare(spos, ev(:,:,2)))

          traj(i,j,1) = traj(i,j,1) + k1x/6 + k2x/3 + k3x/3 + k4x/6
          traj(i,j,2) = traj(i,j,2) + k1y/6 + k2y/3 + k3y/3 + k4y/6

          
       end do
    end do

  end subroutine step_traj_dis
  
  function interp(a1,a2,b1,b2,am) result(bm)
    real, intent(in) :: a1,a2, am
    real, intent(in) :: b1,b2
    real :: bm
    
    bm = b2 + ((b1-b2)*(am-a2)/(a1-a2))
    
  end function interp

  function interp4(a1,a2,b1,b2,am) result(bm)
    real, intent(in) :: a1,a2, am
    real*4, intent(in) :: b1,b2
    real*4 :: bm
    
    bm = b2 + ((b1-b2)*(am-a2)/(a1-a2))
    
  end function interp4

  function interp2d(spos,pos,square2d) result(yea)
    integer, dimension(2), intent(in) :: spos
    real, dimension(2), intent(in) :: pos
    real, dimension(2,2), intent(in) :: square2d
    real :: yea
    real :: xp, yp
    real :: fy1,fy2

    xp = pos(1)
    yp = pos(2)

    ! Need to make sure that the points are still in the domain, otherwise take care of them 
    ! (here, by sticking to walls.


    if(xp.gt.x(nx)) xp = x(nx)
    if(xp.lt.x(1))  xp = x(1)
    
    if(yp.gt.y(ny)) yp = y(ny)
    if(yp.lt.y(1))  yp = y(1)
       
    !interpolate the x = x1 and the x = x2 faces
    
    fy1 = interp(x(spos(1)), x(spos(1)+1), square2d(1,1), square2d(2,1), xp)
    fy2 = interp(x(spos(1)), x(spos(1)+1), square2d(1,2), square2d(2,2), xp)
    
    yea = interp(y(spos(2)), y(spos(2)+1), fy1, fy2, yp)
    
  end function interp2d

  function getsquare(spos,thing) result(square)
    integer, dimension(2), intent(in) :: spos
    real*4, dimension(nx,ny), intent(in) :: thing
    real, dimension(2,2) :: square
   
    !obtain cube of 8 data points to pass to interp3d
    square(1,1) = real(thing(spos(1),spos(2)),8)
    square(1,2) = real(thing(spos(1),spos(2)+1),8)
    square(2,1) = real(thing(spos(1)+1,spos(2)),8)
    square(2,2) = real(thing(spos(1)+1,spos(2)+1),8)

  end function getsquare
  
  function getbound(loc) result(spos)
    real, dimension(2), intent(in) :: loc
    integer, dimension(2) :: spos
    real :: xp, yp
    integer :: j,igrd
    real :: thetap, theta1, theta2

    xp = loc(1)
    yp = loc(2)

    ! Need to make sure that the points are still in the domain, otherwise take care of them 
    ! (here, by sticking to walls.)

    if(xp.gt.x(nx)) xp = x(nx)
    if(xp.lt.x(1))  xp = x(1)
    if(yp.gt.y(ny)) yp = y(ny)
    if(yp.lt.y(1))  yp = y(1)

    spos(1) = nx-1
    spos(2) = ny-1

    do j=1,nx
       if(xp.gt.x(j)) then
       else
          spos(1) = j-1
          exit
       end if
    end do

    do j=1,ny
       if(yp.gt.y(j)) then
       else
          spos(2) = j-1
          exit
       end if
    end do

    if(spos(1).eq.0) spos(1)=1
    if(spos(2).eq.0) spos(2)=1
    
    if(spos(1).eq.nx) spos(1) = nx-1
    if(spos(2).eq.ny) spos(2) = ny-1
    
  end function getbound

  subroutine velst(t,tjump)
    integer, intent(in) :: t, tjump
    integer :: i,j,igrd

    do i=1,nx
       do j=1,ny
             velt(i,j,1) = interp4(tjump*idata*delt, (tjump+dir)*idata*delt, & 
                  u1(i,j), u(i,j), t*delta*dir)
             velt(i,j,2) = interp4(tjump*idata*delt, (tjump+dir)*idata*delt, & 
                  v1(i,j), v(i,j), t*delta*dir)
       end do
    end do

  end subroutine velst

 subroutine vorst(t,tjump)
    integer, intent(in) :: t, tjump
    integer :: i,j,igrd

    do i=1,nx
       do j=1,ny
             vort(i,j) = interp4(tjump*idata*delt, (tjump+dir)*idata*delt, & 
                  om1(i,j), om(i,j), t*delta*dir)
       end do
    end do

  end subroutine vorst

 end module intsubs
