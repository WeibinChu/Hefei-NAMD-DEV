program test_sgesv
    use omp_lib
    implicit none
    integer :: i, j, x
    integer :: ct, cr, cm, t1, t2
    real :: t3, t4

    ! x = 999
    ! call omp_set_num_threads(8)
    ! x = omp_get_max_threads()
    ! write(*,*) x

    ! !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
    ! do i = 1, 20
    !     x = x+i
    !     write(*,*) i, x, omp_get_thread_num()
    ! end do
    ! !$OMP END PARALLEL DO

    call system_clock(t1, cr, cm)
    call cpu_time(t3)

    do i=1, 1000
        do j=1, 1000
            write(100,*) ''
        end do
    end do

    call cpu_time(t4)
    call system_clock(t2, cr, cm)

    write(*,*) t1, t2, t3, t4, cr,cm
    write(*,*) 'CPU', t4-t3
    write(*,*) 'SYS', MOD(t2-t1, cm) / REAL(cr)


end program test_sgesv