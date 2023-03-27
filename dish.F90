module dish
  use prec
  use constants
  use fileio
  use utils
  use parallel
  use couplings
  use hamil
  use TimeProp
#ifdef ENABLEMPI
  use mpi
#endif
  implicit none

  private :: calcDect, whichToDec, projector
  private :: printDISH, printMPDISH

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

  subroutine projector(ks, COEFFISQ, which, tion, indion, cstat, iend, fgend, ks0, olap)
    implicit none
    type(TDKS), intent(inout) :: ks, ks0
    type(overlap), intent(in) :: olap
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

    dE = ((olap%Eig(  :  ,tion) + olap%Eig(  :  , tion+1)) - &
          (olap%Eig(cstat,tion) + olap%Eig(cstat, tion+1))) /2.0_q
    if (inp%LHOLE) then
      dE = -dE
    end if
    where ( dE > 0.0_q )
      popBoltz = popBoltz * exp(-dE / kbT)
    end where
    popOnWhich = popBoltz(which)

    if (r <= popOnWhich) then
    !! project in, hop: cstat -> which
      ks%psi_c = con%cero
      ks%psi_c(which) = con%uno

      if ((.NOT. fgend) .AND. which == iend) then
        ks0%dish_pops(cstat, (indion+1):) = ks0%dish_pops(cstat, (indion+1):) + 1.0_q/inp%NTRAJ
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
!      popBoltz(which) = 0.0_q
!      do i=1, ks%ndim
!        if (i == 1) then
!          ! lower = 0.0_q
!          ! upper = COEFFISQ(i)
!          lower = popOnWhich 
!          upper = popOnWhich + popBoltz(i)
!        else
!          ! lower = upper
!          ! upper = upper + COEFFISQ(i)
!          lower = upper
!          upper = upper + popBoltz(i)
!        end if
!
!        if (lower <= r .AND. r < upper) then
!          if ((.NOT. fgend) .AND. i == iend) then
!            ! Here use decimal part of dish_pops to store the recom_pops matrix.
!            ks0%dish_pops(cstat, (indion+1):) = ks0%dish_pops(cstat, (indion+1):) + 1.0_q/inp%NTRAJ
!            fgend = .TRUE.
!          end if
!          cstat = i !! Hop to i
!          exit
!        end if
!      end do
      !end if
    end if

  end subroutine


  subroutine runDISH(ks, olap)
    implicit none
    type(TDKS), intent(inout) :: ks
    type(overlap), intent(inout) :: olap

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
    real(kind=q), dimension(ks%ndim)           :: tmp_pops

    integer :: lbd, ubd, uboundry
    integer :: lower, upper
    integer :: ierr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! A = Diag, B = Hpsi, C = SH, D = write
    ! #PROG  1  A...B..A...B..C....D
    !        2  A...B..A...B..C....
    !        3  A...B..A...B..C....
    !        4  A...B..A...B..C....
    !        5  A...B..A...B..C....
#ifdef ENABLEMPI
    call MPI_BARRIER(inp%COMMUNICATOR, ierr)
#endif

    ! N = ceiling(REAL(inp%NTRAJ)/inp%NPROG)
    call divideTasks(inp%IPROG, inp%NPROG, inp%NTRAJ, lower, upper, N)
    allocate(ks_list(N))
    allocate(fgend(N))
    allocate(cstat(N))

    ! At the first step, current state always equal initial state
    istat = inp%INIBAND
    iend  = inp%LORB
    cstat = istat
    fgend = .FALSE.
    init  = .TRUE.
    lbd   = 2
    ubd   = 1 ! For short verion, ubd <=> indion
    uboundry = MAXSIZE

    deallocate(ks%dish_pops)
    ks_list = ks
    allocate(ks%dish_pops(ks%ndim, MIN(MAXSIZE, inp%NAMDTIME)))
    ks%dish_pops = 0.0_q
    ks%dish_pops(inp%INIBAND, 1) = 1.0_q

    shuffle = [(i, i=1,ks%ndim)]
    edt = inp%POTIM / inp%NELM

    ! long version. First ouput init condition. Then output M-1 at one time.
    ! Beware coup is long enough.
    call setOlapBoundry(olap, inp%NAMDTINI)
    if (inp%IPROG == 0) then
      call printDISH(ks, olap, lbd, ubd, step=0)
      if (inp%LSPACE) call printMPDISH(ks, lbd, ubd, step=0)
    end if

    !! The loop start here. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do indion = 1, inp%NAMDTIME - 1
      if (.NOT. all(fgend)) then
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
            call make_hamil(tion, tele, olap, ks)
            call TrotterMat(ks, edt)
          case (11) ! this method will accumulate error! May renormalize wfc.
            call make_hamil(tion, tele, olap, ks)
            call EulerMat(ks, edt)
          case (10)
            if (init) then
              ks%psi_p = ks%psi_c
              call make_hamil(tion, tele, olap, ks)
              call EulerMat(ks, edt)
            else
              call make_hamil2(tion, tele-1, olap, ks)
              call EulerModMat(ks, edt)
            end if
          case (2)
            call make_hamil(tion, tele, olap, ks)
            call DiagonizeMat(ks, edt)
          end select
          ham = ks%ham_c

          ! TODO init conflict with OMP
          !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ksi) IF (inp%NPARDISH > 1)
          do i = 1, N
            if (fgend(i)) cycle

            select case (inp%ALGO_INT)
            case (10)
              if (.NOT. init) then
                ksi = ks_list(i)
                ksi%psi_n = ksi%psi_p + HamPsi(ham, ksi%psi_p, 'T')
                ksi%psi_p = ksi%psi_c
                ksi%psi_c = ksi%psi_n
                ks_list(i) = ksi
              else
                ks_list(i)%psi_c = ks_list(i)%psi_c + HamPsi(ham, ks_list(i)%psi_c, 'T')
                init = .FALSE.
              end if
            case (11)
              ks_list(i)%psi_c = ks_list(i)%psi_c + HamPsi(ham, ks_list(i)%psi_c, 'T')
            case default
              ks_list(i)%psi_c = HamPsi(ham, ks_list(i)%psi_c, 'N')
            end select
          end do
          !!$OMP END PARALLEL DO
        end do
      
        !!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ksi,norm,which,COEFFISQ,DECOTIME) FIRSTPRIVATE(shuffle) IF (inp%NPARDISH > 1)
        do i = 1, N
          if (fgend(i)) cycle

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
            call projector(ksi, COEFFISQ, which, tion, ubd, cstat(i), iend, fgend(i), ks, olap)
          end if
          
          if ( fgend(i) ) then
            ! recomb completed
            ks%dish_pops(cstat(i), ubd+1:) = ks%dish_pops(cstat(i), ubd+1:) + 3.0_q
          else
            ks%dish_pops(cstat(i), ubd+1) = ks%dish_pops(cstat(i), ubd+1) + 3.0_q
          end if

          ks_list(i) = ksi
        end do
        !!$OMP END PARALLEL DO
      end if

      ubd = ubd + 1
      if (ubd == uboundry) then
#ifdef ENABLEMPI
        if (inp%NPROG > 1) then
          call MPI_REDUCE(ks%dish_pops+0.0_q, ks%dish_pops, ks%ndim*MAXSIZE, MPI_REAL8, MPI_SUM, 0, inp%COMMUNICATOR, ierr)
        end if
        call MPI_BARRIER(inp%COMMUNICATOR, ierr)
#endif
        tmp_pops = MOD(ks%dish_pops(:,ubd), 1.0_q)
        if (inp%IPROG == 0) then
          call printDISH(ks, olap, lbd, ubd, step=1)
          if (inp%LSPACE) call printMPDISH(ks, lbd, ubd, step=1)
        end if
        ! refresh
        lbd = ubd + 1
        ubd = ubd
        uboundry = ubd+MAXSIZE-1
        ! recom: keep accumulating. Copy from the last of the old list.
        ! pops:  keep changing, except those have recombined. 
        !        The first of the new list is wrong, but useless.
        deallocate(ks%dish_pops)
        allocate(ks%dish_pops(ks%ndim, ubd:ubd+MAXSIZE-1))
        ks%dish_pops = 0.0_q
        ks%dish_pops(iend,:) = 3.0_q * COUNT(fgend)
        if (inp%IPROG == 0) then
          do i = ubd, ubd+MAXSIZE-1
            ks%dish_pops(:,i) = ks%dish_pops(:,i) + tmp_pops
          end do
        end if
        call setOlapBoundry(olap, tion+1)
      end if
    end do
    !! The loop end here. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef ENABLEMPI
    ! MPI_REDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERROR)
    ! Buffer parameters inbuf and inoutbuf must not be aliased
    ! use +0 bypass it.
    if (inp%NPROG > 1) then
      call MPI_REDUCE(ks%dish_pops+0.0_q, ks%dish_pops, size(ks%dish_pops), MPI_REAL8, MPI_SUM, 0, inp%COMMUNICATOR, ierr)
    end if
    call MPI_BARRIER(inp%COMMUNICATOR, ierr)
#endif

    if (inp%IPROG == 0) then
      call printDISH(ks, olap, lbd, ubd, step=2)
      if (inp%LSPACE) call printMPDISH(ks, lbd, ubd, step=2)
    end if

    deallocate(ks_list)
    deallocate(fgend)
    deallocate(cstat)
  end subroutine


  subroutine printDISH(ks, olap, lbd, ubd, step)
    implicit none
    type(TDKS), intent(in) :: ks
    type(overlap), intent(in) :: olap
    integer, intent(in) :: lbd, ubd, step !! step = 0|1|2 init|-|stop

    integer :: i, tion, indion, ierr
    character(len=48) :: buf
    character(len=48) :: out_fmt
    real(kind=q), allocatable, dimension(:,:) :: avgene
    real(kind=q), allocatable, dimension(:,:) :: dish_pops, recom_pops

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f11.6))" )' )  ks%ndim

    if (step == 0) then
      indion = 1
      tion = 2 + MOD(indion+inp%NAMDTINI-1-2,inp%NSW-2)
      ! open, in the init step ks%dish_pops is the real pop & ks%recom_pops = 0
      if (inp%LBINOUT) then
        open(unit=27,                                   &
             file='AVGENE.bin.' // trim(adjustl(buf)),  &
             form='unformatted',                        &
             status='unknown',                          &
             access='stream',                           &
             action='write',                            &
             iostat=ierr)
        open(unit=24,                                   &
             file='SHPROP.bin.' // trim(adjustl(buf)),  &
             form='unformatted',                        &
             status='unknown',                          &
             access='stream',                           &
             action='write',                            &
             iostat=ierr)
        open(unit=28,                                   &
             file='RECOMB.bin.' // trim(adjustl(buf)),  &
             form='unformatted',                        &
             status='unknown',                          &
             access='stream',                           &
             action='write',                            &
             iostat=ierr)
        if (ierr /= 0) then
          write(*,*) "[E] IOError: SHPROP.bin file I/O error!"
          stop
        end if

        write(unit=27) [REAL(tion, kind=q), &
                        indion * inp%POTIM, &
                        SUM(olap%Eig(:,tion) * ks%dish_pops(:,indion))]
        write(unit=24) ks%dish_pops(:,indion)
        write(unit=28) (0.0_q, i=1, ks%ndim)
        return
      end if

      open(unit=24, file='SHPROP.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
      open(unit=28, file='RECOMB.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
      if (ierr /= 0) then
        write(*,*) "[E] IOError: SHPROP file I/O error!"
        stop
      end if

      write(unit=24,fmt=out_fmt) indion * inp%POTIM, SUM(olap%Eig(:,tion) * ks%dish_pops(:,indion)), &
                            (ks%dish_pops(i,indion), i=1, ks%ndim)

      write(unit=28,fmt=out_fmt) indion * inp%POTIM, SUM(olap%Eig(:,tion) * ks%dish_pops(:,indion)), &
                            (0.0_q, i=1, ks%ndim)

      return
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(dish_pops(ks%ndim, lbd:ubd))
    allocate(recom_pops(ks%ndim, lbd:ubd))
    dish_pops = NINT(ks%dish_pops(:, lbd:ubd)/3) / REAL(inp%NTRAJ, kind=q)
    recom_pops = MOD(ks%dish_pops(:, lbd:ubd), 1.0_q)

    if (inp%LBINOUT) then
      allocate(avgene(3, lbd:ubd)) ! tion, time, avgene
      do indion=lbd, ubd
        tion = 2 + MOD(indion+inp%NAMDTINI-1-2,inp%NSW-2)
        avgene(:,indion) = [REAL(tion, kind=q), &
                            indion * inp%POTIM, &
                            SUM(olap%Eig(:,tion) * dish_pops(:,indion))]
      end do

      write(unit=27) avgene
      write(unit=24) dish_pops
      write(unit=28) recom_pops
      deallocate(avgene)
    else
      do indion=lbd, ubd
        tion = 2 + MOD(indion+inp%NAMDTINI-1-2,inp%NSW-2)

        write(unit=24,fmt=out_fmt) indion * inp%POTIM, SUM(olap%Eig(:,tion) * dish_pops(:,indion)), &
                              (dish_pops(i,indion), i=1, ks%ndim)

        write(unit=28,fmt=out_fmt) indion * inp%POTIM, SUM(olap%Eig(:,tion) * dish_pops(:,indion)), &
                              (recom_pops(i,indion), i=1, ks%ndim)
      end do
    end if
    deallocate(dish_pops)
    deallocate(recom_pops)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (step == 2) then
      if (inp%LBINOUT) close(27)
      close(24)
      close(28)
    end if

  end subroutine

  subroutine printMPDISH(ks, lbd, ubd, step)
    implicit none
    type(TDKS), intent(inout) :: ks
    integer, intent(in) :: lbd, ubd, step !! step = 0|1|2 init|-|stop

    integer :: i, ierr
    integer :: bi, tion, indion
    integer, dimension(inp%NACELE) :: bands

    character(len=48) :: buf
    character(len=48) :: out_fmt
    real(kind=q), allocatable, dimension(:,:) :: dish_mppops

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f11.6))" )' )  inp%NBADNS

    if (step == 0) then
      indion = 1
      tion = 2 + MOD(indion+inp%NAMDTINI-1-2,inp%NSW-2)

      allocate(dish_mppops(inp%NBADNS, 1))
      dish_mppops = 0.0_q
      do i=1, ks%ndim
        bands = inp%BASIS(:,i) - inp%BMIN + 1
        do bi=1, inp%NACELE
          dish_mppops(bands(bi),:) = dish_mppops(bands(bi),:) + &
                                     ks%dish_pops(i,:)
        end do
      end do
      ! open, in the init step ks%dish_pops is the real pop & ks%recom_pops = 0
      if (inp%LBINOUT) then
        open(unit=52,                                      &
             file='MPSHPROP.bin.' // trim(adjustl(buf)),   &
             form='unformatted',                           &
             status='unknown',                             &
             access='stream',                              &
             action='write',                               &
             iostat=ierr)
        if (ierr /= 0) then
          write(*,*) "[E] IOError: MPSHPROP.bin file I/O error!"
          stop
        end if

        write(unit=52) dish_mppops
        deallocate(dish_mppops)
        return
      end if

      open(unit=52, file='MPSHPROP.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
      if (ierr /= 0) then
        write(*,*) "[E] IOError: MPSHPROP file I/O error!"
        stop
      end if

      write(unit=52, fmt=out_fmt) indion * inp%POTIM, & 
                                  ! SUM(ks%eigKs(:,tion) * ks%dish_pops(:,indion)) / inp%NACELE, &
                                  (dish_mppops(i,indion), i=1, inp%NBADNS)
      deallocate(dish_mppops)
      return
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(dish_mppops(inp%NBADNS, lbd:ubd))
    dish_mppops = 0.0_q
    ks%dish_pops = NINT(ks%dish_pops/3) / REAL(inp%NTRAJ, kind=q)

    do i=1, ks%ndim
      bands = inp%BASIS(:,i) - inp%BMIN + 1
      do bi=1, inp%NACELE
        dish_mppops(bands(bi),:) = dish_mppops(bands(bi),:) + &
                                   ks%dish_pops(i,lbd:ubd)
      end do
    end do

    if (inp%LBINOUT) then
      write(unit=52) dish_mppops
    else
      do indion=lbd, ubd
        tion = 2 + MOD(indion+inp%NAMDTINI-1-2, inp%NSW-2)

        write(unit=52, fmt=out_fmt) indion * inp%POTIM, & 
                                    ! SUM(ks%eigKs(:,tion) * ks%dish_pops(:,indion)) / inp%NACELE, &
                                    (dish_mppops(i,indion), i=1, inp%NBADNS)
      end do
    end if
    deallocate(dish_mppops)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (step == 2) then
      close(52)
    end if

  end subroutine

end module
