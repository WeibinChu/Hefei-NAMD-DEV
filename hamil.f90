module hamil
  use prec
  use fileio
  use utils
  use couplings
  use constants
  implicit none

  type TDKS
    integer :: ndim
    ! _[c,p,n] means current, previous, next
    complex(kind=q), allocatable, dimension(:) :: psi_c, psi_p
    complex(kind=q), allocatable, dimension(:,:) :: psi_a
    ! the result of hamiltonian acting on a vector
    ! complex(kind=q), allocatable, dimension(:) :: hpsi
    !! population
    real(kind=q), allocatable, dimension(:,:) :: pop_a
    !!!!real(kind=q), allocatable, dimension(:) :: pop_c, pop_n
    ! real(kind=q) :: norm_c !! unused

    complex(kind=q), allocatable, dimension(:,:) :: ham_c
    ! complex(kind=q), allocatable, dimension(:) :: alpha_c
    ! real(kind=q), allocatable, dimension(:,:) :: beta_c

    ! KS eigenvalues & Non-adiabatic couplings
    real(kind=q), allocatable, dimension(:,:) :: eigKs
    real(kind=q), allocatable, dimension(:,:,:) :: NAcoup

    !! surface hopping related
    !! Bkm = REAL(DCONJG(Akm) * Ckm)
    real(kind=q), allocatable, dimension(:) :: Bkm
    real(kind=q), allocatable, dimension(:,:) :: sh_pops
    real(kind=q), allocatable, dimension(:,:) :: sh_prop
    real(kind=q), allocatable, dimension(:,:,:) :: sh_prob !! P_ij (j,i,N-1)

    !! decoherence induced surface hopping
    real(kind=q), allocatable, dimension(:,:) :: dish_pops
    real(kind=q), allocatable, dimension(:,:) :: recom_pops
    ! whether the memory has been allocated
    logical :: LALLO = .FALSE.

  end type

contains

  subroutine initTDKS(ks, olap)
    implicit none

    type(TDKS), intent(out)  :: ks
    type(overlap), intent(in)  :: olap

    integer :: N

    ! memory allocation
    ks%ndim = inp%NBASIS
    N = inp%NBASIS

    if (.NOT. ks%LALLO) then
      allocate(ks%psi_c(N))
      ! allocate(ks%hpsi(N))

      allocate(ks%ham_c(N,N))
      ! allocate(ks%alpha_c(N))
      ! allocate(ks%beta_c(N,N))

      allocate(ks%eigKs(N, inp%NSW))
      allocate(ks%NAcoup(N, N, inp%NSW))

      select case (inp%ALGO)
      case ('FSSH')
        allocate(ks%psi_p(N))
        allocate(ks%psi_a(N, inp%NAMDTIME))
        allocate(ks%pop_a(N, inp%NAMDTIME))
        allocate(ks%sh_pops(N, inp%NAMDTIME))
        allocate(ks%sh_prop(N, inp%NAMDTIME))
        allocate(ks%sh_prob(N, N, inp%NAMDTIME-1))
        allocate(ks%Bkm(N))
      case ('DISH')
        allocate(ks%dish_pops(N, inp%NAMDTIME))
        allocate(ks%recom_pops(N,inp%NAMDTIME))
      case default !! should be check in fileio.f90
      end select

      ! Now copy olap%eig&Dij => ks%eig%Dij
      ks%eigKs = olap%Eig
      ks%NAcoup = olap%Dij / (2*inp%POTIM)

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
      ks%sh_prop = 0.0_q
      ks%sh_prob = 0.0_q
    case ('DISH')
      ks%dish_pops = 0.0_q
      ks%recom_pops = 0.0_q
    end select

  end subroutine

  ! constructing the hamiltonian by replicating NAC
  subroutine make_hamil_rtime(tion, TELE, ks)
    implicit none

    type(TDKS), intent(inout) :: ks
    integer, intent(in) :: tion, TELE

    integer :: RTIME,XTIME !! left & right
    integer :: i
    RTIME = tion
    XTIME = RTIME + 1

    ks%ham_c = (0.0_q, 0.0_q)
    ! ks%alpha_c = (0.0_q, 0.0_q)
    ! ks%beta_c = 0.0_q
    if (TELE <= (inp%NELM / 2)) then
      ks%ham_c(:,:) = interpolate((TELE + inp%NELM/2 - 0.5_q) / inp%NELM, &
                                          ks%NAcoup(:,:,RTIME-1), ks%NAcoup(:,:,RTIME))
      ! ks%ham_c(:,:) = ks%NAcoup(:,:,RTIME)
      ! ks%ham_c(:,:) = ks%NAcoup(:,:,RTIME - 1) + (ks%NAcoup(:,:,RTIME) - &
      !               & ks%NAcoup(:,:,RTIME-1)) * (TELE+inp%NELM/2 - 0.5_q) / inp%NELM
    else 
      ks%ham_c(:,:) = interpolate((TELE - inp%NELM/2 - 0.5_q) / inp%NELM, &
                                          ks%NAcoup(:,:,RTIME), ks%NAcoup(:,:,XTIME))
      ! ks%ham_c(:,:) = ks%NAcoup(:,:,RTIME)
      ! ks%ham_c(:,:) = ks%NAcoup(:,:,RTIME ) + (ks%NAcoup(:,:,XTIME) - &
      !               & ks%NAcoup(:,:,RTIME)) * (TELE-inp%NELM/2 - 0.5_q ) / inp%NELM
    end if

    ! multiply by -i * hbar
    ks%ham_c = -con%I * con%hbar * ks%ham_c 
    
    ! the energy eigenvalue part
    do i=1, ks%ndim
      ks%ham_c(i,i) = interpolate((TELE - 0.5_q) / inp%NELM, &
                                  ks%eigKs(i,RTIME), ks%eigKs(i,XTIME))
      ! ks%ham_c(i,i) = ks%eigKs(i,RTIME)
      ! ks%ham_c(i,i) = ks%eigKs(i,RTIME) +  (ks%eigKs(i,XTIME) - ks%eigKs(i,RTIME)) * (TELE - 0.5_q) / inp%NELM
    end do

  end subroutine

end module
