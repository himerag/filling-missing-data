module dineof_utils


contains

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dindiff(x,B,alpha,numit)
  use ufileformat
  implicit none

  real,intent(inout)         :: B(:)
  real,intent(in)            :: x(:)
  real, pointer              :: BF(:),xe(:)
  integer                    :: nsize,k,numit
  real                       :: alpha,valex


nsize = size(B,1)

allocate(BF(size(B,1)+1))

allocate(xe(size(x,1)+1))
BF=0

! extended x
xe(1) = 1.5*x(1) - .5 * x(2)
xe(2:nsize) = (x(1:nsize-1) + x(2:nsize))/2
xe(nsize+1) = 1.5*x(nsize) - .5 * x(nsize-1)


do k=1,numit
!  F(2:nsize) = alpha * (f(2:nsize) - f(1:nsize-1));
!  f = f +  (F(2:nsize+1) - F(1:nsize));  
  BF(2:nsize) = alpha * (B(2:nsize) - B(1:nsize-1))/(x(2:nsize) - x(1:nsize-1))
  B = B + (BF(2:nsize+1) - BF(1:nsize))/ (xe(2:nsize+1) - xe(1:nsize))
end do

deallocate(BF,xe)

end subroutine dindiff

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module
