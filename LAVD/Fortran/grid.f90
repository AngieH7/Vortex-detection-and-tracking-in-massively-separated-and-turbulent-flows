!Set up the grids for both the data and the ftle calculation

module grid
  real, dimension(:), allocatable :: x,y
  integer :: nx, ny
  real, dimension(:), allocatable :: x1,y1
  integer :: ox, oy
  integer :: oxs=1,oys=1
  integer :: tend, dir, tinc
  real :: xrat, yrat
  integer :: wherecut
  real :: delta, inttime, pi, delt
  real :: xscale, yscale, xstart, ystart
  character(64) :: oname, iname
  integer :: itraj, iftle, idata, eul, startstep,ilavd
  character(125) :: idir, odir
  logical :: freshstart
  
  namelist /inputdata/ nx, ny, freshstart, startstep, itraj, &
       iftle, ilavd, idata, eul, delt, idir, iname, odir, oname, &
       nstart, tend, tinc, delta, inttime, dir, xscale, &
       yscale, xrat, yrat, xstart, ystart

contains
  
  subroutine setup_grid
    implicit none
    integer :: i,j,m, status
    character(1) :: thing
    real :: dx1, dy1
    character (125) fname
    integer, parameter :: fp=67
    pi = 4*atan(1.0)

    allocate(x(nx), y(ny))

    write(fname, '(a,i7.7,a)') trim(idir) // 'part',nstart,'.dat'
    print *, 'Reading file (for grid): ', fname
    open(unit=10,file=trim(fname),form='formatted')
    read(10,'(a)') thing
    read(10,'(a)') thing
    read(10,'(a)') thing
    do j=1,ny
       do i=1,nx
          read(10,'(e16.5, 1X, e16.5)') x(i), y(j)
       end do
    end do
    close(10)
    
    ox = int(xscale*nx)
    oy = int(yscale*ny)
    
    print *, "nx = ", nx, "ox = ", ox, '[ ', x(1),' ',x(nx),' ]'
    print *, "ny = ", ny, "oy = ", oy, '[ ', y(1),' ',y(ny),' ]'

    allocate(x1(ox),y1(oy))
    
    dx1=abs((x(nx)-x(1)) * xrat / ((float(nx)*xscale)-1))
    dy1=abs((y(ny)-y(1)) * yrat / ((float(ny)*yscale)-1))
    
    ! trajectory grid -- base on factor of resolution, 
    ! startpoints, and fraction of domain
    
    do j=1,ox
       x1(j) = x(1) + ((j-1)*dx1) + xstart
    end do

    do j=1,oy
       y1(j) = y(1) + ((j-1)*dy1) + ystart
    end do

  end subroutine setup_grid

  subroutine initialize
    implicit none
    integer :: status, j
    integer :: fp=310
    
    open(unit=fp, file="strain.inp", status="old", iostat=status)
 
    print *, "Reading strain.inp..."
    read(unit=fp, nml=inputdata)
 
    close(fp)

  end subroutine initialize
  

end module grid
