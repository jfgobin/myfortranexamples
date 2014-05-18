function gcd(a, b) result(c)
! Returns the GCD of a and b
 integer :: a,b,c,u,l,m
 if ( a > b ) then
   u = a
   l = b
 else
   u = b
   l = a
 end if
 do while (l > 0)
    m = modulo(u,l)
    u=l
    l=m
 end do
 c=u
 return
end function

subroutine simplify(a, b, c, d)
! Returns the fraction a/b in its simplified form 
! c/d where c and d are relatively prime
 integer, intent(in) :: a,b
 integer, intent(out) :: c,d
 integer :: n,gcd
 n=gcd(a,b)
 c=a/n
 d=b/n
 return
end subroutine

subroutine findfractionpref(pref, rept, mpref, a, b)
! Returns the fraction a/b such as its division gives the pattern prefreptreptrept ....
! With the necessary multiplier
 integer, intent(in) :: pref, rept, mpref
 integer, intent(out) :: a, b
 integer :: d1, d2, num, den
 d1=1+floor(log10(real(pref)))
 d2=1+floor(log10(real(rept)))
 num=((pref*10**(d2)+rept)-pref)
 den=10**(d1+d2)-10**(d1)
 if ( mpref > 0) then
  den=den*10**mpref
 else
  num=num*10**mpref
 end if
 call simplify(num, den,a , b)
 return
end subroutine

subroutine findfractionnopref(rept, a,b)
! Returns the fraction a/b such as its division gives the pattern reptreptrept...
 integer, intent(in) :: rept
 integer, intent(out) :: a,b
 integer :: d1, num, den
 d1=1+floor(log10(real(rept)))
 num=rept
 den=10**d1-1
 call simplify(num,den,a,b)
 return
end subroutine 


program fractionfinder
 implicit none
 integer :: pref,rept,a,b
 character :: c
 do
  print *, 'Does your fraction include a prefix (yY/nN) or Q to quit (qQ)?'
  read(*,'(A1)'), c
  if ((c == 'y').or.(c == 'Y')) then
   print *, 'Prefix part?'
   read(*, '(I12)'), pref
   print *, 'Repeated part?'
   read(*,'(I12)'), rept
   call findfractionpref(pref,rept,0, a,b)
   print *, 'The requested fraction is ', a, '/', b
  else if ((c == 'n').or.(c == 'N')) then
   print *, 'Repeated part?'
   read(*, '(I8)'), rept
   call findfractionnopref(rept,a,b)
   print *, 'The requested fraction is ', a, '/', b
  else if ((c == 'q').or.(c == 'Q')) then
   goto 100
  end if
 end do
100 print *, 'Bye bye!'
 stop
end program
