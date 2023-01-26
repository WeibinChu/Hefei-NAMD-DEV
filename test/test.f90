program test_sgesv
    implicit none
    real :: a(3,3),b(3)
    integer :: v(3),iflag
    external sgesv
    a=reshape([2.0,0.0,0.0,0.0,3.0,0.0,0.0,0.0,4.0],[3,3])
    b=[998.0,999.0,1000.0]
    print *,'a=',a
    print *,'b=',b
    call sgesv(3,1,a,3,v,b,3,iflag)
    print *,'solve=',b
end program test_sgesv