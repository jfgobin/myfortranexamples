! matrix - matrix operations
!
! det(X,n,d) - computes the determinant of matrix X
! inverse(X,Y) - computes the inverse of matrix X


! det(ox,n,d)
! 
! Computes the determinant of matrix on (size n by n) using Gauss Jordan decomposition
! The original matrix is transformed using row operations into a U (upper triangular) matrix
! whose determinant is the product of the diagonal elements.
!
subroutine det(ox,n,d)
 integer, intent(in) :: n
 real(8), dimension(n,n), intent(in) :: ox
 real(8), dimension(n,n) :: x
 real(8), intent(out) :: d
 real(8)  :: row, pivot
 integer :: sgn, k, l, m
 ! The determinant is computed by creating an upper matrix
 ! And simply multiplying the elements on the diagonal
 ! First, let's copy the matrix into a new one, to 
 ! honour the INTENT(IN) statement
 x=ox
 sgn=1					! If two rows are swapped, then sgn is multiplied by -1 
 do k = 1,n-1				! We will create an upper triangular matrix by 
  if (x(k,k) == 0) then			! Suppressing the (n-k-1) values of the rows below the kth
   do l=k,n				! If the (k,k) element is 0, let's swap the line with another one.
    if (x(k,l) /= 0) then 		! If it is not possible (i.e all the remaining element at the kth column are
     do m=1,n  				! zero, break and return 0
      row=x(m,k)
      x(m,k)=x(m,l)
      x(m,l)=row
     end do
     sgn=-sgn
     exit
     else
      d=0.0
      return
    end if
   end do
  end if				! At this point, the line intend to use as a pivot as a non zero value ... let's
  do l=k+1,n				! use that pivot to suppress all values on the kth column for the row below
   pivot=x(k,l)/x(k,k)
   do m=k,n
    x(m,l)=x(m,l)-pivot*x(m,k)
   end do
  end do
 end do
 d=1.0
 do k = 1,n
  d=d*x(k,k)
 end do
 d=d*sgn
end subroutine det

! Inverse(x,n,xi)
! 
! Computes the inverse of matrix ox (size n by n) using Gauss Jordan decomposition:
! The original matrix is transformed using row and column operations into the identity
! matrix, while the (originally) identity matrix receives the same operations (same order).
subroutine inverse(ox,n,xi) 
 integer, intent(in) :: n
 real(8), dimension(n,n), intent(in) :: ox
 real(8), dimension(n,n), intent(out) :: xi
 real(8), dimension(n,n) :: tempx, x
 integer :: k,l,m
 x=ox
 tempx=0				! First, let's init the temp matrix as the identity one
 do k = 1,n
  tempx(k,k)=1
 end do
 do k=1,n				! I start by creating a U matrix from the original one (same code as for the determinant)
  if (x(k,k) == 0) then			! If the current element is 0, let's exchange with a non zero value
   do l=k+1,n
    if(x(l,k) /= 0) then			
     call swaprows(x,n,k,l)
     call swaprows(tempx,n,k,l)
     exit
    else					! This matrix is not invertible
     xi=0
     return
    end if
   end do
  end if
  pivot=x(k,k)
  do l=1,n				! First, transform the element (k,k) into a 1 by dividing it (and the rest of the line as well
   x(l,k)=x(l,k)/pivot			! That way, I can avoid a later division.
   tempx(l,k)=tempx(l,k)/pivot		
  end do
  do l=1,n				! Let's suppress all the values on the column except myself
   if ( l /= k) then
    pivot=x(k,l)
    do m=1,n
     x(m,l)=x(m,l)-pivot*x(m,k)
     tempx(m,l)=tempx(m,l)-pivot*tempx(m,k)
    end do
   end if
  end do
 end do					! At this point, we have an upper triangular matrix
 xi=tempx
end subroutine inverse   

! Swap the k-th and l-th rows in a matrix that is n by n. 
! The transformation is done in place.  
subroutine swaprows(x,n,k,l)
 integer, intent(in) :: n,k,l
 real(8), dimension(n,n), intent(inout) :: x
 real(8) :: v
 integer :: i
 do i = 1,n
  v=x(i,k)
  x(i,k)=x(i,l)
  x(i,l)=v
 end do
end subroutine
