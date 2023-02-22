module utils
  use prec
  use constants
#ifdef ENABLEMPI
  use mpi
#endif
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
  subroutine init_random_seed(salt)
    implicit none
    integer, intent(in) :: salt
    integer :: i, n, clock, ierr
    integer, dimension(:), allocatable :: seed
    real :: r

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    seed = seed + (modulo(clock, 65536)+300) * salt
    call random_seed(put = seed)

    deallocate(seed)

#ifdef ENABLEMPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
    call random_number(r)
    write(*,'(A,I6,A,F6.3)') '[I] Salt: ', salt, ', First random number: ', r
#ifdef ENABLEMPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

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

  pure function zDIAG(N, V)
    implicit none
    integer, intent(in)   :: N
    complex(q), intent(in), dimension(N) :: V
    complex(q), dimension(N,N) :: zDIAG

    logical, dimension(N,N) :: MASK
    integer :: i

    zDIAG = 0.0_q
    MASK = reshape([(MOD(i,N) == i/N, i=0,N*N-1)], shape=[N,N])
    zDIAG = unpack(V, MASK, zDIAG)
  end function

  subroutine printToScreen(msg, fh)
    character(len=*), intent(in) :: msg
    integer, intent(in) :: fh
    integer :: ierr
#ifdef ENABLEMPI
    call MPI_FILE_WRITE_SHARED(fh, msg, len(msg), MPI_CHAR, MPI_STATUS_IGNORE, ierr)
    ! call MPI_FILE_WRITE_ORDERED(fh, msg, len(msg), MPI_CHAR, MPI_STATUS_IGNORE, ierr)
#else
    write(fh,*) msg
#endif
  end subroutine

end module utils