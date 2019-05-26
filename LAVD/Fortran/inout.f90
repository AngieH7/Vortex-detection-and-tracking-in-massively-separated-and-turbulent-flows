!Routines and functions for outputting the M|Z criterion and FTLE.
module inout
  use grid
  implicit none
  real*4, dimension(:,:), allocatable :: u,v, om, psi
  real*4, dimension(:,:), allocatable :: u1,v1,om1,vort
  real*4, dimension(:,:,:), allocatable :: velt
real, dimension(:,:), allocatable :: ivd,lavd

contains
  
  subroutine readdata(t)
    integer :: i,j
    integer, intent(in) :: t
    character(150) :: fname
    character(1) :: thing
    integer :: time, n1, n2, n3, n4
    real*4 :: test1

    write(fname, '(a,i7.7,a)') trim(idir) // trim(iname), t,'.dat'
    
    open(28, file=trim(fname), form='formatted')
    
    write(6,*)'Reading file : ', fname
    read(28, '(a)') thing
    read(28, '(a)') thing
    read(28, '(a)') thing
    do j=1,ny
       do i=1,nx
          read(28,'(e16.5,e16.5,e16.5,e16.5,e16.5,e16.5)') test1, test1, u(i,j), v(i,j), om(i,j), psi(i,j)
       end do
    end do
print *,'om',om(100,100)
    close(28)

  end subroutine readdata

  subroutine writetraj(vel, traj, t, tstart)
    real, dimension(ox,oy,3), intent(in) :: vel, traj
    integer, intent(in) :: t, tstart
    character(50) :: fname
    integer :: i,j

    if(dir.lt.0) then
       write(fname,'(a,i6.6,a,i6.6,a)') trim(odir) // trim(oname) // 'traj_neg_',tstart,'_',t,'.dat'
    else
       write(fname,'(a,i6.6,a,i6.6,a)') trim(odir) // trim(oname) // 'traj_pos_',tstart,'_',t,'.dat'
    end if
    
   print *, 'Writing trajectory file ', fname
   open(24,file=trim(fname),form="unformatted") 
   write(24) ox,oy,2
   write(24) ((real(traj(i,j,1),4),i=1,ox),j=1,oy), &
        ((real(traj(i,j,2),4),i=1,ox),j=1,oy)

   close(24)    

   if(dir.lt.0) then
      write(fname,'(a,i6.6,a,i6.6,a)') trim(odir) // trim(oname) // 'vel_neg_',tstart,'_',t,'.dat'
   else
      write(fname,'(a,i6.6,a,i6.6,a)') trim(odir) // trim(oname) // 'vel_pos_',tstart,'_',t,'.dat'
   end if
   
   print *, 'Writing trajectory velocity file ', fname
   open(25,file=trim(fname),form="unformatted") 
   write(25) ox,oy,2
   write(25) ((real(vel(i,j,1),4),i=1,ox),j=1,oy), &
        &    ((real(vel(i,j,2),4),i=1,ox),j=1,oy)

   close(25)    
    
  end subroutine writetraj

  subroutine writeftle(ftle, t, tstart)
    integer :: i,j
    real, dimension(ox,oy), intent(in) :: ftle
    integer, intent(in) :: t, tstart
    character(50) :: fname

    if(dir.lt.0) then
       write(fname, '(a,i6.6,a,i6.6,a)') trim(odir) // trim(oname) // 'ftle_neg_',tstart, '_', t,'.dat'
    else
       write(fname, '(a,i6.6,a,i6.6,a)') trim(odir) // trim(oname) // 'ftle_pos_',tstart, '_', t,'.dat'
    end if

    print *, 'Writing ftle file', fname
    open(34, file=trim(fname), form='unformatted')
    write(34) ox,oy,1
    write(34) ((real(ftle(i,j),4),i=1,ox), j=1,oy)
    close(34)

  end subroutine writeftle

  subroutine writelavd(t, tstart,lavd)
    integer :: i,j
    real, dimension(ox,oy), intent(in) :: lavd
    integer, intent(in) :: t, tstart
    character(50) :: fname

       write(fname, '(a,i6.6,a,i6.6,a)') trim(odir) // trim(oname) // 'lavd_',tstart,'_',t,'.dat'

    print *, 'Writing lavd file', fname
    open(84, file=trim(fname), form='unformatted')
    write(84) ox,oy,1
    write(84) ((real(lavd(i,j),4),i=1,ox), j=1,oy)
    close(84)
print *,'lavd(ox,oy)=',lavd(ox,oy)
  end subroutine writelavd

  subroutine writeeul(eulcrit, tstart)
    integer :: i,j
    real, dimension(ox,oy,5), intent(in) :: eulcrit
    integer, intent(in) :: tstart
    character(40) :: fname
    integer :: nzpstart, nzp

    ! Will only be printing z-sections, so not going to do all the extra stuff
    ! in the x- and y- directions

    write(fname, '(a,i6.6,a,i6.6,a)') trim(odir) // trim(oname) // 'eul_zrat_',tstart,'.dat'

    print *, 'Writing eulerian file', fname
    open(34, file=trim(fname), form='unformatted')
    write(34) ox,oy,2
    write(34) ((real(eulcrit(i,j,1),4),i=1,ox), j=1,oy), &
         ((real(eulcrit(i,j,2),4),i=1,ox), j=1,oy)
!         ((real(eulcrit(i,j,3),4),i=1,ox), j=1,oy), &
!         ((real(eulcrit(i,j,4),4),i=1,ox), j=1,oy), &
!         ((real(eulcrit(i,j,5),4),i=1,ox), j=1,oy)
    close(34)

  end subroutine writeeul
  
  subroutine write_nx(tstart)
    integer :: i,j
    integer, intent(in) :: tstart
    character(40) :: fname

    write(fname, '(a,i6.6,a)') trim(odir) // trim(oname) // 'eul_nx_',tstart,'.dat'

    print *, 'Writing eulerian nx file', fname
    open(35, file=trim(fname), form='unformatted')
    write(35) nx,ny,2
    write(35) ((om(i,j),i=1,nx), j=1,ny), &
         ((psi(i,j),i=1,nx), j=1,ny)
    close(35)

    write(fname, '(a)') trim(odir) // trim(oname) // 'grid_nx.dat'

    print *, 'Writing eulerian nx grid file', fname
    open(35, file=trim(fname), form='unformatted')
    write(35) nx,ny
    write(35) ((real(x(i),4),i=1,nx), j=1,ny), &
         ((real(y(j),4),i=1,nx), j=1,ny)
    close(35)

  end subroutine write_nx
  
  subroutine writemesh
    integer :: i,j,k

    print *, 'Writing mesh file'
 
    open(31,file=trim(odir)//'/grid.dat',form='UNFORMATTED')
    write(31) ox,oy
    write(31) ((real(x1(i),4),i=1,ox),j=1,oy), &
         &    ((real(y1(j),4),i=1,ox),j=1,oy)
 
    close(31)

  end subroutine writemesh
  
  subroutine writeivd(ivd, tstart)
    integer :: i,j
    real, dimension(ox,oy), intent(in) :: ivd
    integer, intent(in) :: tstart
    character(50) :: fname

    write(fname, '(a,i6.6,a,i6.6,a)') trim(odir) // trim(oname) // 'ivd_',tstart,'.dat'

    print *, 'Writing ivd file', fname
    open(37, file=trim(fname), form='unformatted')
    write(37) ox,oy,1
    write(37) ((real(ivd(i,j),4),i=1,ox), j=1,oy)
    close(37)

  end subroutine writeivd

  subroutine writevor(vor, tstart)
    integer :: i,j
    real, dimension(ox,oy), intent(in) :: vor
    integer, intent(in) :: tstart
    character(50) :: fname

    write(fname, '(a,i6.6,a,i6.6,a)') trim(odir) // trim(oname) // 'vor_',tstart,'.dat'

    print *, 'Writing vor file', fname
    open(737, file=trim(fname), form='unformatted')
    write(737) ox,oy,1
    write(737) ((real(vor(i,j),4),i=1,ox), j=1,oy)
    close(737)

  end subroutine writevor

end module inout
