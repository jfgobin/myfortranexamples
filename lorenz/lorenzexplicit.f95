! Lorenz Explicit
! Uses Forward Euler to compute the next (x,y,z)
!

subroutine lorenzexplicit(x,sigma,rho,beta,dt, x1)
    ! Computes the vector dx using the Lorenz system
    real(8), dimension(3), intent(in) :: x
    real(8), intent(in) :: sigma,rho,beta,dt
    real(8), dimension(3), intent(out) :: x1
    x1(1) = x(1) + dt*(sigma * (x(2) - x(1)))
    x1(2) = x(2) + dt*(x(1)*(rho-x(3))-x(2))
    x1(3) = x(3) + dt*(x(1)*x(2)-beta*x(3))
    return
end subroutine lorenzexplicit

