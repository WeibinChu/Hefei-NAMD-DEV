program test
    implicit none
    integer, dimension(3) :: x=[1,2,3],y=[4,5,6],z=[7,8,9]
    write(*,*) calc(x,y,z)
contains

    elemental function calc(a,b,c)
        integer, intent(in) :: a, b, c
        integer :: calc
        calc = a + b + c
    end function
end program test