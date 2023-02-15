module dish
  use prec
  use constants
  use fileio
  use utils
  use couplings
  use hamil
  use TimeProp
#ifdef ENABLEMPI
  use mpi
#endif
  implicit none

  private :: calcDect, whichToDec, projector

contains

  subroutine calcDect(ks, DECOTIME, COEFFISQ)
    implicit none
    type(TDKS), intent(in) :: ks

    real(kind=q), dimension(ks%ndim), intent(out) :: DECOTIME !! tau_i(t)
    real(kind=q), dimension(ks%ndim), intent(out) :: COEFFISQ !! |c_i(t)|^2
    integer :: i

    COEFFISQ(:) = REAL(CONJG(ks%psi_c(:)) * ks%psi_c(:), kind=q)

    ! Calculating Decoherence Time tau_a(t)
    ! 1/tau_a(t) = SUM_{i!=a}|c_i(t)|^2 * r_ai
    do i=1, ks%ndim
      ! DEPHMATR(j,i)=DEPHMATR(i,j)
      DECOTIME(i) = 1.0_q / SUM(inp%DEPHMATR(:,i) * COEFFISQ(:))
    end do
    !write(90,*),DECOTIME
  end subroutine

  subroutine whichToDec(ks, DECOTIME, which, shuffle)
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), dimension(ks%ndim), intent(in) :: DECOTIME  !! tau_i(t)
    integer, intent(out)                         :: which
    integer, dimension(ks%ndim), intent(inout)   :: shuffle

    integer :: tmp
    integer :: i, j
    !! integer :: tm(inp%NBASIS), decstate(inp%NBASIS)
    !! the colomn of decstate is the sum number of decoherence states (n)
    !! the element of decstate is the serial number of decoherence states (which)
    real(kind=q) :: r

    !! n = 0
    !! decstate = 0
    which = 0


    do j=ks%ndim,1,-1
      call random_number(r)
      i=INT(ks%ndim*r)+1
      tmp=shuffle(i)
      shuffle(i)=shuffle(j)
      shuffle(j)=tmp
    end do


    do j=1, ks%ndim
      i=shuffle(j)
      if ( DECOTIME(i) <= ks%dish_decmoment(i) ) then !! tau_i(t) <= t_i(t)
        which = i
        ks%dish_decmoment(which) = 0.0_q    !! update the decoherence moment
        exit
      end if
    end do
  end subroutine

  subroutine projector(ks, COEFFISQ, which, tion, indion, cstat, iend, fgend, ks0)
    implicit none
    type(TDKS), intent(inout) :: ks, ks0
    real(kind=q), dimension(ks%ndim), intent(inout) :: COEFFISQ !! |c_i(t)|^2

    integer, intent(in) :: which, tion, indion, iend
    integer, intent(inout) :: cstat
    logical, intent(inout) :: fgend

    integer :: i
    real(kind=q) :: r, kbT, lower, upper, norm_c, popOnWhich
    real(kind=q), dimension(ks%ndim) :: popBoltz, dE

    kbT = inp%TEMP * BOLKEV
    call random_number(r)

    !! hopping probability with Boltzmann factor
    popBoltz = COEFFISQ

    dE = ((ks0%eigKs(  :  ,tion) + ks0%eigKs(  :  , tion+1)) - &
          (ks0%eigKs(cstat,tion) + ks0%eigKs(cstat, tion+1))) /2.0_q
    if (inp%LHOLE) then
      dE = -dE
    end if
    where ( dE > 0.0_q )
      popBoltz = popBoltz * exp(-dE / kbT)
    end where
    popOnWhich = popBoltz(which)
    ! popBoltz = COEFFISQ(which)

    ! dE = ((ks0%eigKs(which,tion) + ks0%eigKs(which, tion+1)) - &
    !       (ks0%eigKs(cstat,tion) + ks0%eigKs(cstat, tion+1))) /2.0_q
    ! if (inp%LHOLE) then
    !   dE = -dE
    ! end if
    ! if ( dE > 0.0_q ) then
    !   popBoltz = popBoltz * exp(-dE / kbT)
    ! end if
    !write(*,*) which,popBoltz

    if (r <= popOnWhich) then
    !! project in, hop: cstat -> which
      ks%psi_c = con%cero
      ks%psi_c(which) = con%uno

      if ((.NOT. fgend) .AND. which == iend) then
        ks0%dish_pops(cstat, (indion+1):inp%NAMDTIME) = ks0%dish_pops(cstat, (indion+1):inp%NAMDTIME) + 1.0_q/inp%NTRAJ
        fgend = .TRUE.
      end if

      cstat = which

    else
    !! project out, hop: cstat -> j != which, with prob |c_j(t)|^2
    !! Here is a problem, when projecting out, hop upward hasn't been scaled by Boltzmann factor
    !! Assuming if r > MAX(upper) then reject this hop, this should obey the DBC.
      ks%psi_c(which) = con%cero
      COEFFISQ(which) = 0.0_q
      norm_c = SUM(COEFFISQ) !! 1-|c_a|^2
      ks%psi_c(:) = ks%psi_c(:) / SQRT(norm_c)

      !if (cstat==which) then
      !! renormalize COEFFISQ, so SUM(|c_j(t)|^2) = 1
      ! COEFFISQ = COEFFISQ / norm_c 
      ! call random_number(r)
      popBoltz(which) = 0.0_q
      do i=1, ks%ndim
        if (i == 1) then
          ! lower = 0.0_q
          ! upper = COEFFISQ(i)
          lower = popOnWhich 
          upper = popOnWhich + popBoltz(i)
        else
          ! lower = upper
          ! upper = upper + COEFFISQ(i)
          lower = upper
          upper = upper + popBoltz(i)
        end if

        if (lower <= r .AND. r < upper) then
          if ((.NOT. fgend) .AND. i == iend) then
            ! Here use decimal part of dish_pops to store the recom_pops matrix.
            ks0%dish_pops(cstat, (indion+1):inp%NAMDTIME) = ks0%dish_pops(cstat, (indion+1):inp%NAMDTIME) + 1.0_q/inp%NTRAJ
            fgend = .TRUE.
          end if
          cstat = i !! Hop to i
          exit
        end if
      end do
      !end if
    end if

  end subroutine


  subroutine runDISH(ks, olap)
    implicit none
    type(TDKS), intent(inout) :: ks
    type(overlap), intent(in) :: olap

    integer :: i, N, tion, indion, tele
    logical :: init
    real(kind=q) :: edt, norm
    complex(kind=q), dimension(ks%ndim,ks%ndim) :: ham
    type(TDKS) :: ksi

    integer :: istat, iend
    integer :: which
    type(TDKS), allocatable, dimension(:) :: ks_list
    logical,    allocatable, dimension(:) :: fgend !! recomb indicator, tells recomb from which state at which time
    integer,    allocatable, dimension(:) :: cstat
    integer, dimension(ks%ndim)                :: shuffle
    real(kind=q), dimension(ks%ndim)           :: COEFFISQ  !! |c_i(t)|^2
    real(kind=q), dimension(ks%ndim)           :: DECOTIME  !! tau_i(t)

    integer :: ierr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! A = Diag, B = Hpsi, C = SH, D = write
    ! #PROG  1  A...B..A...B..C....D
    !        2  A...B..A...B..C....
    !        3  A...B..A...B..C....
    !        4  A...B..A...B..C....
    !        5  A...B..A...B..C....
#ifdef ENABLEMPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

    N = ceiling(REAL(inp%NTRAJ)/inp%NPROG)
    allocate(ks_list(N))
    allocate(fgend(N))
    allocate(cstat(N))

    ! At the first step, current state always equal initial state
    istat = inp%INIBAND
    iend  = inp%LORB
    cstat = istat
    fgend = .FALSE.
    init  = .TRUE.

    deallocate(ks%eigKs)
    deallocate(ks%NAcoup)
    ks_list = ks
    allocate(ks%eigKs(ks%ndim, inp%NSW))
    allocate(ks%NAcoup(ks%ndim, ks%ndim, inp%NSW))
    allocate(ks%dish_pops(ks%ndim, inp%NAMDTIME))
    if (inp%IPROG == 0) allocate(ks%recom_pops(ks%ndim,inp%NAMDTIME))
    ks%eigKs = olap%Eig
    ks%NAcoup = olap%Dij
    ks%dish_pops = 0.0_q
    if (inp%IPROG == 0) ks%recom_pops = 0.0_q

    shuffle = [(i, i=1,ks%ndim)]

    edt = inp%POTIM / inp%NELM

    do indion = 1, inp%NAMDTIME - 1
      ! NAMDTINI >= 2, tion+inp%NAMDTINI-1 <= NSW(FSSH), dump 1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !          S                  E
      ! E1--D1--E2--D2--E3--D3------EN--DN--XX
      !                         E1--D1--E2--D2--E3--D3------EN--DN--XX
      !                              S                      E
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! for DISH, one cycle has NSW-2 intervals.
      tion = 2 + MOD(indion+inp%NAMDTINI-1-2, inp%NSW-2) ! 2 <= tion <= 2+(NSW-3) = NSW-1

      do tele = 1, inp%NELM
        select case (inp%ALGO_INT)
        case (0)
          call make_hamil(tion, tele, ks)
          call TrotterMat(ks, edt)
        case (11) ! this method will accumulate error! May renormalize wfc.
          call make_hamil(tion, tele, ks)
          call EulerMat(ks, edt)
        case (10)
          if (init) then
            ks%psi_p = ks%psi_c
            call make_hamil(tion, tele, ks)
            call EulerMat(ks, edt)
          else
            call make_hamil2(tion, tele-1, ks)
            call EulerModMat(ks, edt)
          end if
        case (2)
          call make_hamil(tion, tele, ks)
          call DiagonizeMat(ks, edt)
        end select
        ham = ks%ham_c

        ! TODO init conflict with OMP
        !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ksi) IF (inp%NPARDISH > 1)
        do i = 1, N
          select case (inp%ALGO_INT)
          case (10)
            if (.NOT. init) then
              ksi = ks_list(i)
              ksi%psi_n = ksi%psi_p + HamPsi(ham, ksi%psi_p, 'C')
              ksi%psi_p = ksi%psi_c
              ksi%psi_c = ksi%psi_n
              ks_list(i) = ksi
            else
              ks_list(i)%psi_c = ks_list(i)%psi_c + HamPsi(ham, ks_list(i)%psi_c, 'C')
              init = .FALSE.
            end if
          case (11)
            ks_list(i)%psi_c = ks_list(i)%psi_c + HamPsi(ham, ks_list(i)%psi_c, 'C')
          case default
            ks_list(i)%psi_c = HamPsi(ham, ks_list(i)%psi_c, 'N')
          end select
        end do
        !!$OMP END PARALLEL DO
      end do
     
      !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ksi,norm,which,COEFFISQ,DECOTIME) FIRSTPRIVATE(shuffle) IF (inp%NPARDISH > 1)
      do i = 1, N
        ksi = ks_list(i)
        norm = REAL(SUM(CONJG(ksi%psi_c) * ksi%psi_c), kind=q) 
        if (norm <= 0.99_q .OR. norm >= 1.01_q)  then
            write(*,*) "[E] Error in Electronic Propagation"
            stop
        end if

        call calcDect(ksi, DECOTIME, COEFFISQ)   

        ksi%dish_decmoment(:) = ksi%dish_decmoment(:) + inp%POTIM
        call whichToDec(ksi, DECOTIME, which, shuffle)  
        
        if ( which > 0 ) then
          call projector(ksi, COEFFISQ, which, tion, indion, cstat(i), iend, fgend(i), ks)
        end if
        
        ks%dish_pops(cstat(i), indion+1) = ks%dish_pops(cstat(i), indion+1) + 3.0_q

        ks_list(i) = ksi
      end do
      !!$OMP END PARALLEL DO
    end do

#ifdef ENABLEMPI
    ! MPI_REDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERROR)
    ! Buffer parameters inbuf and inoutbuf must not be aliased
    ! use +0 bypass it.
    if (inp%NPROG > 1) then
      call MPI_REDUCE(ks%dish_pops+0.0_q, ks%dish_pops, ks%ndim*inp%NAMDTIME, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end if
    !call MPI_REDUCE(ks%recom_pops+0.0_q,ks%recom_pops,ks%ndim*inp%NAMDTIME, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (inp%IPROG == 0) then
#endif
    ks%recom_pops = MOD(ks%dish_pops, 1.0_q)
    ks%dish_pops = NINT(ks%dish_pops/3) / REAL(inp%NTRAJ, kind=q)
    ! ks%dish_pops = ks%dish_pops / inp%NTRAJ
    ks%dish_pops(inp%INIBAND, 1) = 1.0_q
    ! ks%recom_pops = ks%recom_pops / inp%NTRAJ
#ifdef ENABLEMPI
    end if
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
  end subroutine


  subroutine printDISH(ks)
    implicit none
    type(TDKS), intent(in) :: ks

    integer :: i, tion, indion, ierr
    character(len=48) :: buf
    character(len=48) :: out_fmt

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f11.6))" )' )  ks%ndim

    open(unit=24, file='SHPROP.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    open(unit=28, file='RECOMB.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: SHPROP file I/O error!"
      stop
    end if

    do indion=1, inp%NAMDTIME
      tion = 2 + MOD(indion+inp%NAMDTINI-1-2,inp%NSW-2)

      write(unit=24,fmt=out_fmt) indion * inp%POTIM, SUM(ks%eigKs(:,tion) * ks%dish_pops(:,indion)), &
                            (ks%dish_pops(i,indion), i=1, ks%ndim)

      write(unit=28,fmt=out_fmt) indion * inp%POTIM, SUM(ks%eigKs(:,tion) * ks%dish_pops(:,indion)), &
                            (ks%recom_pops(i,indion), i=1, ks%ndim)
    end do

    close(24)
    close(28)

  end subroutine

  subroutine printMPDISH(ks)
    implicit none
    type(TDKS), intent(inout) :: ks

    integer :: i, ierr
    integer :: bi, tion, indion
    integer, dimension(inp%NACELE) :: bands

    character(len=48) :: buf
    character(len=48) :: out_fmt

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f11.6))" )' )  inp%NBADNS

    open(unit=52, file='MPSHPROP.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: MPSHPROP file I/O error!"
      stop
    end if

    ks%dish_mppops = 0.0_q
    do i=1, ks%ndim
      bands = inp%BASIS(:,i) - inp%BMIN + 1
      do bi=1, inp%NACELE
        ks%dish_mppops(bands(bi),:) = ks%dish_mppops(bands(bi),:) + &
                                      ks%dish_pops(i,:)
      end do
    end do

    do indion=1, inp%NAMDTIME
      tion = 2 + MOD(indion+inp%NAMDTINI-1-2, inp%NSW-2)
      write(unit=52, fmt=out_fmt) indion * inp%POTIM, & 
                                  SUM(ks%eigKs(:,tion) * ks%dish_pops(:,indion)) / inp%NACELE, &
                                  (ks%dish_mppops(i,indion), i=1, inp%NBADNS)
    end do

    close(52)

  end subroutine

end module
