module analysis 
  use grid 
  implicit none
contains
  
  subroutine eulerian(eulcrit, vel)
    use inout
    integer :: i,j,k
    real, dimension(3,3) :: gradu, s, omega, ft1, ft2, mult, tensor
    real, dimension(ox,oy,5), intent(inout) :: eulcrit
    real, dimension(ox,oy,2), intent(in) :: vel
    real, dimension(3) :: eigs
    real :: orgp, orgq, orgr, tilr, tilq
    real :: theta, ptrm, mtrm, swmax, maxe
    
    swmax = 0.0

    
    do i=1,ox
       do j=1,oy
          
          if(i.eq.1) then
             gradu(1,1) = (vel(i+1,j,1) - vel(i,j,1)) / (x1(i+1) - x1(i))
             gradu(2,1) = (vel(i+1,j,2) - vel(i,j,2)) / (x1(i+1) - x1(i))
             gradu(3,1) = 0.0
          elseif(i.eq.ox) then
             gradu(1,1) = (vel(i,j,1) - vel(i-1,j,1)) / (x1(i) - x1(i-1))
             gradu(2,1) = (vel(i,j,2) - vel(i-1,j,2)) / (x1(i) - x1(i-1))
             gradu(3,1) = 0.0
          else
             gradu(1,1) = (vel(i+1,j,1) - vel(i-1,j,1)) / (x1(i+1) - x1(i-1))
             gradu(2,1) = (vel(i+1,j,2) - vel(i-1,j,2)) / (x1(i+1) - x1(i-1))
             gradu(3,1) = 0.0
          end if
          if(j.eq.1) then
             gradu(1,2) = (vel(i,j+1,1) - vel(i,j,1)) / (y1(j+1) - y1(j))
             gradu(2,2) = (vel(i,j+1,2) - vel(i,j,2)) / (y1(j+1) - y1(j))
             gradu(3,2) = 0.0
          elseif(j.eq.oy) then
             gradu(1,2) = (vel(i,j,1) - vel(i,j-1,1)) / (y1(j) - y1(j-1))
             gradu(2,2) = (vel(i,j,2) - vel(i,j-1,2)) / (y1(j) - y1(j-1))
             gradu(3,2) = 0.0
          else
             gradu(1,2) = (vel(i,j+1,1) - vel(i,j-1,1)) / (y1(j+1) - y1(j-1))
             gradu(2,2) = (vel(i,j+1,2) - vel(i,j-1,2)) / (y1(j+1) - y1(j-1))
             gradu(3,2) = 0.0
          end if
          
          gradu(1,3) = 0.0
          gradu(2,3) = 0.0
          gradu(3,3) = 0.0
          
          omega(1,1) = 0.
          omega(1,2) = 0.5*(gradu(1,2)-gradu(2,1))
          omega(1,3) = 0.5*(gradu(1,3)-gradu(3,1))
          omega(2,2) = 0.
          omega(2,3) = 0.5*(gradu(2,3)-gradu(3,2))
          omega(3,3) = 0.
          omega(2,1) = -omega(1,2)
          omega(3,1) = -omega(1,3)
          omega(3,2) = -omega(2,3)
    
          s(1,1) = gradu(1,1)
          s(1,2) = 0.5*(gradu(1,2)+gradu(2,1))
          s(1,3) = 0.5*(gradu(1,3)+gradu(3,1))
          s(2,2) = gradu(2,2)
          s(2,3) = 0.5*(gradu(2,3)+gradu(3,2))
          s(3,3) = gradu(3,3)
          s(2,1) = s(1,2)
          s(3,1) = s(1,3)
          s(3,2) = s(2,3)
          
          ft1 = matmul(omega, transpose(omega))
          ft2 = matmul(s,transpose(s))
    
          ! Q criterion
          eulcrit(i,j,1) = .5*(ft1(1,1) + ft1(2,2) + ft1(3,3) - ft2(1,1) - ft2(2,2) - ft2(3,3))

          !!! Vorticity in second spot now
          eulcrit(i,j,2) = gradu(2,1) - gradu(1,2)
          
!!$          orgp = - (gradu(1,1) + gradu(2,2) + gradu(3,3))
!!$          mult = matmul(gradu, gradu)
!!$          orgq = .5 * (orgp**2 - (mult(1,1) + mult(2,2) + mult(3,3)))
!!$          mult = matmul(mult, gradu)
!!$          orgr = (-orgp**3 + 3.*orgp*orgq - (mult(1,1) + mult(2,2) + mult(3,3)))/3.
!!$          tilr = orgr + 2.*orgp**3/27. - orgp*orgq/3.
!!$          theta = acos( (-.5*tilr) / ((-1./3.)*tilq)**(1.5) ) 
!!$
!!$          ! Delta criterion
!!$          eulcrit(i,j,2) = orgq**3/27. + .25*orgr**2
!!$
!!$          if(eulcrit(i,j,2).gt.0) then
!!$
!!$             ptrm = .5*orgr + sqrt(eulcrit(i,j,2))
!!$             ptrm = sign(abs(ptrm)**(1./3.), ptrm)
!!$             mtrm = .5*orgr - sqrt(eulcrit(i,j,2))
!!$             mtrm = sign(abs(mtrm)**(1./3.), mtrm)
!!$
!!$             ! Swirl criterion, lambda_ci^2 (have to remember to re-add P)
!!$             ! this didn't matter before because P=0 for incompressible
!!$             eulcrit(i,j,4) = (.75*(ptrm-mtrm) - (orgp/3.))**2
!!$
!!$             if(sqrt(eulcrit(i,j,4)).gt..1)then
!!$                ! Chakraborty criterion
!!$                eulcrit(i,j,5) = (.5*(ptrm+mtrm)-(orgp/3.))/sqrt(eulcrit(i,j,4))
!!$             else
!!$                eulcrit(i,j,5) = 0.
!!$             end if
!!$          else
!!$             eulcrit(i,j,4) = 0.
!!$             eulcrit(i,j,5) = 0.
!!$          end if
!!$
!!$          swmax = max(swmax, eulcrit(i,j,4))
!!$          
!!$          ! Lambda-2 criterion
!!$          tensor = matmul(s,s) + matmul(omega,omega)
!!$          
!!$          orgp = - (tensor(1,1) + tensor(2,2) + tensor(3,3))
!!$          mult = matmul(tensor, tensor)
!!$          orgq = .5 * (orgp**2 - (mult(1,1) + mult(2,2) + mult(3,3)))
!!$          mult = matmul(mult, tensor)
!!$          orgr = (-orgp**3 + 3.*orgp*orgq - (mult(1,1) + mult(2,2) + mult(3,3)))/3.
!!$          tilr = orgr + 2.*orgp**3/27. - orgp*orgq/3.
!!$          theta = acos( (-.5*tilr) / ((-1./3.)*tilq)**(1.5) ) 
!!$    
!!$          ! We know tensor is positive definite, so can just use other
!!$          ! set of eigenvalue equations for discriminant > 0.
!!$          eigs(1) = 2 * ((-1./3.)*tilq)**.5 * cos(theta/3.) 
!!$          eigs(2) = 2 * ((-1./3.)*tilq)**.5 * cos((theta+(2*pi))/3.) 
!!$          eigs(3) = 2 * ((-1./3.)*tilq)**.5 * cos((theta+(4*pi))/3.)
!!$          
!!$          eigs(1) = eigs(1) - (1./3.)*orgp 
!!$          eigs(2) = eigs(2) - (1./3.)*orgp 
!!$          eigs(3) = eigs(3) - (1./3.)*orgp
!!$          
!!$          maxe = max(eigs(1),eigs(2),eigs(3))
!!$          if(eigs(1).gt.eigs(2)) then
!!$             if(eigs(2).gt.eigs(3)) then
!!$                eulcrit(i,j,3) = eigs(2)
!!$             else
!!$                if(eigs(1).gt.eigs(3)) then
!!$                   eulcrit(i,j,3) = eigs(3)
!!$                else
!!$                   eulcrit(i,j,3) = eigs(1)
!!$                end if
!!$             end if
!!$          else
!!$             if(eigs(1).gt.eigs(3)) then
!!$                eulcrit(i,j,3) = eigs(1)
!!$             else
!!$                if(eigs(2).gt.eigs(3)) then
!!$                   eulcrit(i,j,3) = eigs(3)
!!$                else
!!$                   eulcrit(i,j,3) = eigs(2)
!!$                end if
!!$             end if
!!$          end if
          
       end do
    end do
    
    
!!$    do i=1,ox
!!$       do j=1,oy
!!$          ! Normalizing swirl by maximum value
!!$             eulcrit(i,j,4) = eulcrit(i,j,4) / swmax
!!$                    
!!$       end do
!!$    end do

  end subroutine eulerian
  
  function ftle_calc(num,traj) result(ftle) 
    integer :: num 
    real, dimension(ox,oy,2), intent(in) :: traj 
    real, dimension(ox,oy) :: ftle
    real, dimension(2) :: xa,xb,xc,ya,yb,yc, st9 
    real :: xoa,xob,xoc,yoa,yob,yoc
    integer :: i,j
    
    ! Setup finite difference scheme -- don’t actually use xb, yb, zb 
    ! anymore, since use a central different scheme. do i=1,ox
    
    do i=1,ox
       do j=1,oy

          if(i.eq.1)then 
             xa = traj(i,j,:) 
             xb = traj(i,j,:) 
             xc = traj(i+1,j,:) 
             xoa = x1(i) 
             xob = x1(i)
             xoc = x1(i+1)
          elseif(i.eq.ox)then 
             xa = traj(i-1,j,:) 
             xb = traj(i,j,:) 
             xc = traj(i,j,:) 
             xoa = x1(i-1) 
             xob = x1(i) 
             xoc = x1(i)
          else 
             xa = traj(i-1,j,:) 
             xb = traj(i,j,:) 
             xc = traj(i+1,j,:) 
             xoa = x1(i-1) 
             xob = x1(i) 
             xoc = x1(i+1)
          end if
          
          if(j.eq.1)then 
             ya = traj(i,j,:) 
             yb = traj(i,j,:) 
             yc = traj(i,j+1,:)
             yoa = y1(j) 
             yob = y1(j) 
             yoc = y1(j+1)
          elseif(j.eq.oy)then 
             ya = traj(i,j-1,:) 
             yb = traj(i,j,:) 
             yc = traj(i,j,:) 
             yoa = y1(j-1) 
             yob = y1(j) 
             yoc = y1(j)
          else 
             ya = traj(i,j-1,:) 
             yb = traj(i,j,:) 
             yc = traj(i,j+1,:) 
             yoa = y1(j-1) 
             yob = y1(j) 
             yoc = y1(j+1)
          end if

          st9 = stretcher(xc, xa, yc, ya, xoc, xoa, yoc, yoa) 
          
          ftle(i,j) = Log(max(st9(1), st9(2)))/(2*num*delta)
             
       end do
    end do
    
  end function ftle_calc

   subroutine  lavd_calc(ivd,lavd)
     integer :: i,j,k
     real, dimension(ox,oy), intent(in) :: ivd 
     real, dimension(ox,oy), intent(inout) :: lavd
     real total
     total = 0.0   
     do i=1,ox
        do j=1,oy
           total = total+lavd(i,j)
        end do
     end do
    ! print *,'before total=',total
     do i=1,ox
        do j=1,oy
          lavd(i,j) = lavd(i,j)+ivd(i,j)
       end do
    end do
    total = 0.0   
    do i=1,ox
       do j=1,oy
          total = total+lavd(i,j)
       end do
    end do
    !print *,'after total=',total
    !   print *,'ivd=',ivd(1,1)
      print *,'lavd=',lavd(ox,oy)
    
  
  end subroutine lavd_calc
  
  subroutine get_ivd(ivd, vor)
    real :: omavery,total
    integer :: i,j,k
    real, dimension(ox,oy), intent(inout) :: ivd
    real, dimension(ox,oy), intent(in) :: vor
    total = 0.0   
    do i=1,ox
       do j=1,oy
          total = total+vor(i,j)
       end do
    end do
   ! print *,'total=',total
    omavery = total/(ox*oy)
   ! print *,'omavery =',omavery
    
    do i=1,ox
       do j=1,oy
          ivd(i,j) = abs(vor(i,j)-omavery)
       end do
    end do
    print *,'ivd=',ivd(i,j)
    print *,'ivd=',ivd(ox,oy)
    
  end subroutine get_ivd
  
  
  function stretcher(xd1,xd2,yd1,yd2,x01,x02,y01,y02)result(eigs) 
    real, dimension(2), intent(in) :: xd1,xd2,yd1,yd2
    real, intent(in) :: x01,x02,y01,y02
    real, dimension(2,2) :: ftlemat, mult 
    real, dimension(2) :: eigs
    real :: a,b,c

    ! Set up Cauchy-Green strain tensor 

    ftlemat(1,1) = (xd1(1)-xd2(1))/(x01-x02) 
    ftlemat(2,1) = (xd1(2)-xd2(2))/(x01-x02) 
    ftlemat(1,2) = (yd1(1)-yd2(1))/(y01-y02) 
    ftlemat(2,2) = (yd1(2)-yd2(2))/(y01-y02) 
    
    ftlemat = matmul(transpose(ftlemat), ftlemat)

    ! Since ftlemat is symmetric, all eigenvalues will be real 
    ! Use quadratic equation to calculate eigenvalues. 

    a = 1.
    b = -(ftlemat(1,1)+ftlemat(2,2))
    c = (ftlemat(1,1)*ftlemat(2,2))-(ftlemat(2,1)*ftlemat(1,2))

    eigs(1) = -b + sqrt(b**2 - (4*a*c)) / 2 / a
    eigs(2) = -b - sqrt(b**2 - (4*a*c)) / 2 / a

  end function stretcher
end module analysis
