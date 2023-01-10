module utils
  use prec
  use constants
  implicit none

  private :: interpolate_c, interpolate_r

  interface interpolate
    procedure :: interpolate_c, interpolate_r
  end interface

  type, private :: mathLib
  contains
    ! private
    ! procedure, private :: interpolate_r, interpolate_c
    ! generic, public :: interpolate => interpolate_r, interpolate_c
  end type
  type(mathLib), protected :: maths

contains

  ! initialize the random seed from the system clock
  ! code from: http://fortranwiki.org/fortran/show/random_seed
  subroutine init_random_seed()
    implicit none
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)

    deallocate(seed)
  end subroutine

  elemental function interpolate_r(x, a0, a1)
    real(q), intent(in) :: x, a0, a1
    real(q)             :: interpolate_r

    interpolate_r = (1.0_q-x) * a0 + x * a1
  end function
    
  elemental function interpolate_c(x, a0, a1)
    real(q), intent(in)    :: x
    complex(q), intent(in) :: a0, a1
    complex(q)             :: interpolate_c

    interpolate_c = (con%uno-x) * a0 + x * a1
  end function

  elemental function cmplx_norm(x)
    complex(q), intent(in) :: x
    real(q)                :: cmplx_norm

    cmplx_norm = REAL(CONJG(x) * x, kind=q)
  end function

end module utils