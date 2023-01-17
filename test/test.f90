program test
    implicit none
    complex :: a = (2,-1), b = (3, -2)
    real, dimension(3) :: x = [1,2,3]
    character(len=48) :: out = '(2(F11.7,SP,F9.6,"i"))'

    write(*,out) a,b
end program test