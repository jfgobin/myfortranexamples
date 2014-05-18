! LORENZ
!
! Compute (x,y,z) values for of Lorenz's Strange Attractor
! using explicit (forward), implicit (backward) and implicit/explicit (forward-backward) Euler
!

program lorenz
    real(8), dimension(3) :: xe, xe1, xi, xi1, xei, xei1, xei2
    real(8) :: beta, sigma, rho, dt
    integer :: i
    character(len=80) :: fmt
    fmt="(F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5, F10.5)"
    dt=0.001
    xe =(/1,1,1/)
    xi=xe
    xei=xe
    rho=28
    sigma=10
    beta=2.6666666667
    do i = 1,100000
        call lorenzexplicit(xe,sigma,rho,beta,dt,xe1)
	call lorenzimplicit(xi,sigma,rho,beta,dt,xi1)
	call lorenzexplicit(xei,sigma,rho,beta,dt,xei1)
	call lorenzimplicit(xei,sigma,rho,beta,dt,xei2)
        write(*, fmt),  xe(1), xe(2), xe(3), xi(1), xi(2), xi(3), xei(1), xei(2), xei(3)
	xe=xe1
	xi=xi1
	xei=(xei1+xei2)/2.0
    end do
    stop
end program
