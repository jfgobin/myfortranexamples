! lorenzimplicit
! Uses implicit (backward) Euleur to compute the next (x,y,z)
! This involves solving a system of 3 nonlinear equations of 3 variables
! At each step
!
!

subroutine lorenzimplicit(x,sigma,rho,beta,dt, x1)
 real(8), dimension(3), intent(in) :: x				! Current position
 real(8), intent(in) :: sigma, rho, beta, dt			! Parameters
 real(8), dimension(3), intent(out) :: x1			! Next position
 real(8), dimension(3,3) :: J, Jinv				! Jacobian and inverse of Jacobian
 real(8), dimension(3) :: fx,x2, jf				! Value of the left sides of the three Lorenz equations
 real(8) :: distance,d 						! put in implicit form
 integer :: k,l,c
 distance=1
 x2=x
 c=1
 do
!  J=reshape((/1+sigma*dt, -dt*(rho-x2(3)), -dt*x2(2), &
!              -sigma*dt,  (1+dt),          -dt*x2(1), &
!              0.0_8,      x2(1),           (1+beta*dt)/), (/3,3/))
  J=reshape((/1+sigma*dt,     -sigma*dt, 0.0_8,        &
             -dt*(rho-x2(3)), (1+dt),    x2(1),        &
	     -dt*x2(2),       -dt*x2(1), (1+beta*dt)/), (/3,3/))
  call det(J,3,d)
  if (d == 0.0 ) then
   x1=x2							! The matrix can't be inverted - minima or solution
   exit
  end if
  call inverse(J,3,Jinv)
  fx(1)=(1+sigma*dt)*x2(1)-sigma*dt*x2(2)-x(1)
  fx(2)=(1+dt)*x2(2)-x2(1)*dt*(rho-x2(3))-x(2)
  fx(3)=(1+beta*dt)*x2(3)-dt*x2(1)*x2(2)-x(3)
  jf=0.0
  do k=1,3
   do l=1,3
    jf(k)= jf(k)+Jinv(l,k)*fx(l)
   end do
  end do
  do k=1,3
   x1(k)=x2(k)-jf(k)
  end do
  distance=sqrt( (x1(1)-x2(1))**2 + (x1(2)-x2(2))**2 +(x1(3)-x2(3))**2)
  x2=x1
  c=c+1
  if ((distance < 1E-4) .or. (c > 99)) then
   goto 100
  end if
 end do
100 x1=x2 							! Exit label
end subroutine lorenzimplicit 
