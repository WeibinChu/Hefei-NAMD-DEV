module hamil
  use prec
  use fileio
  use utils
  use couplings
  use constants
  implicit none

  private :: calcDij_i, calcDij_r

  interface calcDij
    procedure :: calcDij_i, calcDij_r
  end interface

  type TDKS
    integer :: ndim
    ! _[c,p,n] means current, previous, next
    complex(kind=q), allocatable, dimension(:) :: psi_c, psi_p, psi_n
    complex(kind=q), allocatable, dimension(:,:) :: psi_a
    ! the result of hamiltonian acting on a vector
    ! complex(kind=q), allocatable, dimension(:) :: hpsi
    !! population
    real(kind=q), allocatable, dimension(:,:) :: pop_a, mppop_a
    !!!!real(kind=q), allocatable, dimension(:) :: pop_c, pop_n
    ! real(kind=q) :: norm_c !! unused

    complex(kind=q), allocatable, dimension(:,:) :: ham_c

    ! KS eigenvalues & Non-adiabatic couplings
    real(kind=q), allocatable, dimension(:,:) :: eigKs    !! These two are abandoned in the long version.
    real(kind=q), allocatable, dimension(:,:,:) :: NAcoup

    !! surface hopping related
    real(kind=q), allocatable, dimension(:,:) :: sh_pops, sh_mppops
    real(kind=q), allocatable, dimension(:,:,:) :: sh_prob !! P_ij (j,i,N-1)

    !! decoherence induced surface hopping
    real(kind=q), allocatable, dimension(:,:) :: dish_pops, dish_mppops !! dish_mppops & recom_pops are not used.
    real(kind=q), allocatable, dimension(:,:) :: recom_pops
    real(kind=q), allocatable, dimension(:)   :: dish_decmoment !! t_i(t)
    ! whether the memory has been allocated
    logical :: LALLO = .FALSE.

  end type

contains

  subroutine initTDKS(ks)
    implicit none

    type(TDKS), intent(inout)  :: ks

    integer :: N

    ! memory allocation
    ks%ndim = inp%NBASIS
    N = inp%NBASIS

    if (.NOT. ks%LALLO) then
      allocate(ks%psi_c(N))
      allocate(ks%psi_p(N))
      allocate(ks%psi_n(N))
      ! allocate(ks%hpsi(N))

      allocate(ks%ham_c(N,N))

      ! allocate(ks%eigKs(N, inp%NSW))
      ! allocate(ks%NAcoup(N, N, inp%NSW))

      select case (inp%ALGO)
      case ('FSSH')
        allocate(ks%psi_a(N, inp%NAMDTIME))
        allocate(ks%pop_a(N, inp%NAMDTIME))
        allocate(ks%sh_pops(N, inp%NAMDTIME))
        allocate(ks%sh_prob(N, N, inp%NAMDTIME-1))
        if (inp%LSPACE) then
          allocate(ks%mppop_a(inp%NBADNS, inp%NAMDTIME))
          allocate(ks%sh_mppops(inp%NBADNS, inp%NAMDTIME))
        end if
      case ('DISH')
        allocate(ks%dish_decmoment(N))
      case default
      end select

      ! Now copy olap%eig&Dij => ks%eig%Dij
      ! ks%eigKs = olap%Eig
      ! ks%NAcoup = olap%Dij

      ks%LALLO = .TRUE.
    end if

    ! init, put carrier on iniband
    ks%psi_c = con%cero
    ks%psi_c(inp%INIBAND) = con%uno

    select case (inp%ALGO)
    case ('FSSH')
      ks%psi_a = con%cero
      ks%psi_a(inp%INIBAND, 1) = con%uno
      ks%pop_a = 0.0_q 
      ks%pop_a(inp%INIBAND, 1) = 1.0_q
      ks%sh_pops = 0.0_q
      ks%sh_prob = 0.0_q
    case ('DISH')
    end select

  end subroutine

  ! constructing the hamiltonian by replicating NAC
  subroutine make_hamil(tion, tele, olap, ks)
    implicit none

    integer, intent(in) :: tion, tele
    type(overlap), intent(in) :: olap
    type(TDKS), intent(inout) :: ks

    integer :: RTIME,XTIME !! left & right
    integer :: i
    RTIME = tion
    XTIME = RTIME + 1

    ks%ham_c = (0.0_q, 0.0_q)
    if (tele <= (inp%NELM / 2)) then
      ks%ham_c(:,:) = interpolate((tele + inp%NELM/2.0_q - 0.5_q) / inp%NELM, &
                                          olap%Dij(:,:,RTIME-1), olap%Dij(:,:,RTIME))
    else 
      ks%ham_c(:,:) = interpolate((tele - inp%NELM/2.0_q - 0.5_q) / inp%NELM, &
                                          olap%Dij(:,:,RTIME), olap%Dij(:,:,XTIME))
    end if

    ! multiply by -i * hbar
    ks%ham_c = -con%I * con%hbar * ks%ham_c 
    
    ! the energy eigenvalue part
    do i=1, ks%ndim
      ks%ham_c(i,i) = interpolate((tele - 0.5_q) / inp%NELM, &
                                  olap%Eig(i,RTIME), olap%Eig(i,XTIME))
    end do

  end subroutine

  subroutine make_hamil2(tion, tele, olap, ks)
    implicit none

    integer, intent(in) :: tion, tele
    type(overlap), intent(in) :: olap
    type(TDKS), intent(inout) :: ks

    integer :: RTIME,XTIME !! left & right
    integer :: i
    RTIME = tion
    XTIME = RTIME + 1

    ks%ham_c = (0.0_q, 0.0_q)
    if (tele <= (inp%NELM / 2)) then
      ks%ham_c(:,:) = interpolate((tele + inp%NELM/2.0_q) / inp%NELM, &
                                          olap%Dij(:,:,RTIME-1), olap%Dij(:,:,RTIME))
    else 
      ks%ham_c(:,:) = interpolate((tele - inp%NELM/2.0_q) / inp%NELM, &
                                          olap%Dij(:,:,RTIME), olap%Dij(:,:,XTIME))
    end if

    ! multiply by -i * hbar
    ks%ham_c = -con%I * con%hbar * ks%ham_c 
    
    ! the energy eigenvalue part
    do i=1, ks%ndim
      ks%ham_c(i,i) = interpolate(REAL(tele, kind=q) / inp%NELM, &
                                  olap%Eig(i,RTIME), olap%Eig(i,XTIME))
    end do

  end subroutine

  subroutine make_hamil_wrong(tion, tele, olap, ks)
    implicit none

    integer, intent(in) :: tion, tele
    type(overlap), intent(in) :: olap
    type(TDKS), intent(inout) :: ks

    integer :: RTIME,XTIME
    integer :: i
    RTIME = tion
    XTIME = RTIME + 1

    ks%ham_c = (0.0_q, 0.0_q)
    ks%ham_c(:,:) = interpolate(REAL(tele, kind=q) / inp%NELM, &
                                olap%Dij(:,:,RTIME), olap%Dij(:,:,XTIME))

    ! multiply by -i * hbar
    ks%ham_c = -con%I * con%hbar * ks%ham_c 
    
    ! the energy eigenvalue part
    do i=1, ks%ndim
      ks%ham_c(i,i) = interpolate(REAL(tele, kind=q) / inp%NELM, &
                                  (olap%Eig(i,RTIME)+olap%Eig(i,XTIME))/2.0_q, &
                                  (olap%Eig(i,XTIME)+olap%Eig(i,XTIME+1))/2.0_q)
    end do

  end subroutine

  elemental function calcDij_r(dij1, dij2, dij3, tele)
    implicit none
    real(kind=q), intent(in) :: dij1, dij2, dij3
    real(kind=q), intent(in) :: tele
    real(kind=q) :: calcDij_r

    if (tele <= (inp%NELM / 2)) then
      calcDij_r = interpolate((tele + inp%NELM/2.0_q) / inp%NELM, &
                                          dij1, dij2)
    else 
      calcDij_r = interpolate((tele - inp%NELM/2.0_q) / inp%NELM, &
                                          dij2, dij3)
    end if
  end function

  elemental function calcDij_i(dij1, dij2, dij3, tele)
    implicit none
    real(kind=q), intent(in) :: dij1, dij2, dij3
    integer, intent(in) :: tele
    real(kind=q) :: calcDij_i

    calcDij_i = calcDij_r(dij1, dij2, dij3, REAL(tele, kind=q))
  end function

  end module
