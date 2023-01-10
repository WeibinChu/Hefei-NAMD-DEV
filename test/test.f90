program test
    implicit none
    complex :: a = (1,0)
    real :: x
    character(len=256) :: out_fmt
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f9.4))" )' )  10

    write(*,*) out_fmt
end program test