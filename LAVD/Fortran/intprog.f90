program intprog
  use grid
  use inout
  use analysis
  use intsubs
  implicit none
  integer :: n,m,i,j,t,tjump,filenum
  real, dimension(:,:,:), allocatable :: vel, traj, eulcrit
  real, dimension(:,:), allocatable :: ftle
  real, dimension(:,:), allocatable :: vor
  integer :: tstart,enc
  
  ! Read input variables
  call initialize
  ! Set up data and FTLE grids
  call setup_grid
  
  ! Allocate data velocity variables
  allocate(u(nx,ny))
  allocate(v(nx,ny))
  allocate(om(nx,ny))
  allocate(om1(nx,ny))
  allocate(psi(nx,ny))
  allocate(u1(nx,ny))
  allocate(v1(nx,ny))
  allocate(velt(nx,ny,2))
  allocate(vort(nx,ny))
  
  allocate(vel(ox,oy,2))
  allocate(vor(ox,oy))
  allocate(ivd(ox,oy))
  allocate(lavd(ox,oy))
  allocate(traj(ox,oy,2))
  allocate(ftle(ox,oy))
  allocate(eulcrit(ox,oy,5))
  
  call writemesh
  
  do tstart = nstart,tend,tinc
     
     ! Initialize FTLE variables on the FTLE grid
     
     ftle = 0.
     lavd = 0.
     
     call setup_traj(tstart, traj)
     
     tjump = dir*floor((startstep+1)*delta/idata/delt)
     !     print *, 'tjump = ', tjump
     
     call readdata(tstart+(idata*tjump))
     u1 = u
     v1 = v
     om1 = om
     
     call readdata(tstart+idata*(tjump+dir))
     
     ! Start FTLE integration
     !     enc = abs(49000 - tstart)/100
     !     print *,enc
     !     do t=startstep+1,enc
     do t = startstep+1,int((inttime/(delta)))
        
        print *, 't ', t
        
        if(abs(dir*t*delta).gt.abs(idata*delt*(tjump+dir))) then
           tjump = tjump + dir
           
           u1 = u
           v1 = v
           om1 = om
           
           call readdata(tstart+(idata*(tjump+dir)))
        end if
        
        call velst(t, tjump)
        call vorst(t, tjump)
          
        ! Interpolate velocity in space from data grid to current location 
        ! of FTLE trajectories
        
        call step_vels_dis(t,vel,traj)
        call step_vor_dis(t,vor,traj)
        
        ! Use updated velocity fields to step forward in time.
if (t.ne.0)then         
                call step_traj_dis(t, tjump, vel, traj) 
end if
       
       call get_ivd(ivd,vor)

        call lavd_calc(ivd,lavd)
print *,lavd(ox,oy)
        if(mod(t,ilavd).eq.0.and.t.ne.0) then
           call writelavd(t, tstart,lavd)
        end if
        
        
        ! if(mod(t,enc).eq.0.and.t.ne.0)then
        if(mod(t,iftle).eq.0.and.t.ne.0) then
           !ftle = ftle_calc(t, traj)
           !call writeftle(ftle, t, tstart)
        end if
        
        ! if(mod(t,itraj).eq.0) then
        if(t.eq.0)then
     !      call writetraj(vel,traj,t,tstart)
        end if
        
        if(t.eq.0) then
           !call eulerian(eulcrit, vel)
           !call writeeul(eulcrit, tstart)
           !call writevor(vor, tstart)
           !        call write_nx(tstart)
           !call writeivd(ivd, tstart)
        end if
 
     end do
     
  end do
  
  print *, 'program finished!'
  
  deallocate(vel,traj,u,v,x,y,x1,y1,u1,v1,ftle,eulcrit,vor,vort,om,om1,velt,ivd,lavd)
end program intprog
