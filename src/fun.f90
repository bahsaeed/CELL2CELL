module fun

implicit none

contains
 function reciprocal(a) result (b)
    real(8), intent(in) :: a(3,3)
    real(8)             :: b(3,3)
    real(8)  twopi_v
    b(:,1)=cross(a(:,2),a(:,3))
    b(:,2)=cross(a(:,3),a(:,1))
    b(:,3)=cross(a(:,1),a(:,2))
    twopi_v=8d0*atan(1d0)/dot_product(a(:,1),b(:,1))
    b=twopi_v*b
 end function


 function cross(a, b)
    real(8), dimension(3) :: cross
    real(8), dimension(3), intent(in) :: a, b

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
 end function cross

 function matinv3(a) result(b)
    real(8), intent(in) :: a(3,3)
    real(8)             :: b(3,3)
    real(8)             :: detinv

    detinv = 1/(a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2)&
              - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1)&
              + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1))

    b(1,1) = +detinv * (a(2,2)*a(3,3) - a(2,3)*a(3,2))
    b(2,1) = -detinv * (a(2,1)*a(3,3) - a(2,3)*a(3,1))
    b(3,1) = +detinv * (a(2,1)*a(3,2) - a(2,2)*a(3,1))
    b(1,2) = -detinv * (a(1,2)*a(3,3) - a(1,3)*a(3,2))
    b(2,2) = +detinv * (a(1,1)*a(3,3) - a(1,3)*a(3,1))
    b(3,2) = -detinv * (a(1,1)*a(3,2) - a(1,2)*a(3,1))
    b(1,3) = +detinv * (a(1,2)*a(2,3) - a(1,3)*a(2,2))
    b(2,3) = -detinv * (a(1,1)*a(2,3) - a(1,3)*a(2,1))
    b(3,3) = +detinv * (a(1,1)*a(2,2) - a(1,2)*a(2,1))
 end function

 function expect(n,x,a,y)
    integer*4,intent(in) :: n
    complex*16,intent(in):: x(n),a(n,n),y(n)
    complex*16           :: b(n),c(n),t,expect
    integer*4 i,j
    expect= dcmplx(0d0,0d0)
    b=dcmplx(0d0,0d0)
    c=conjg(x)
    do i=1,n
     do j=1,n
      b(i)=b(i)+a(j,i)*y(j)
     enddo
    enddo
    do i=1,n
     t=c(i)*b(i)
     expect=expect+t
    enddo
 end function expect

end module fun
