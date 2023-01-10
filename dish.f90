module dish
  use prec
  use fileio
  use utils
  use hamil
  use TimeProp
  implicit none

  private :: calcdect, whichtodec, projector
  private :: mpirundish

contains

  subroutine calcdect(ks, DECOTIME, COEFFISQ)
    implicit none
    type(TDKS), intent(in) :: ks

    real(kind=q), dimension(inp%NBASIS), intent(inout) :: DECOTIME
    real(kind=q), dimension(inp%NBASIS), intent(inout) :: COEFFISQ
    integer :: i, j
    real(kind=q) :: tmp

    do j=1, inp%NBASIS
      COEFFISQ(j) = REAL(CONJG(ks%psi_c(j)) * ks%psi_c(j), kind=q)
    end do

    ! Calculating Decoherence Time tau_a(t)
    ! 1/tau_a(t) = SUM_{i!=a}|c_i(t)|^2 * r_ai
    do i=1, inp%NBASIS
      tmp = 0
      do j=1, inp%NBASIS
        tmp = tmp + inp%DEPHMATR(j, i) * COEFFISQ(j)   ! DEPHMATR(j,i)=DEPHMATR(i,j)
      end do
      DECOTIME(i) = 1.0_q / tmp
    end do
    !write(90,*),DECOTIME
  end subroutine

  subroutine whichtodec(DECOTIME, which, decmoment, shuffle)
    implicit none
    real(kind=q), intent(in) :: DECOTIME(inp%NBASIS)
    integer, intent(inout) :: which
    real(kind=q), dimension(inp%NBASIS), intent(inout) :: decmoment
    integer, intent(inout),dimension(inp%NBASIS) :: shuffle
    integer :: tmp
    integer :: i, j
    !! integer :: tm(inp%NBASIS), decstate(inp%NBASIS)
    !! the colomn of decstate is the sum number of decoherence states (n)
    !! the element of decstate is the serial number of decoherence states (which)
    real(kind=q) :: r

    !! n = 0
    !! decstate = 0
    which = 0


    do j=inp%NBASIS,1,-1
      call random_number(r)
      i=INT(inp%NBASIS*r)+1
      tmp=shuffle(i)
      shuffle(i)=shuffle(j)
      shuffle(j)=tmp
    end do


    do j=1,inp%NBASIS
      i=shuffle(j)
      if ( DECOTIME(i) <= decmoment(i) ) then
        which = i
        decmoment(which) = 0.0_q    !! update the decoherence moment
        exit
      end if
    end do
  end subroutine

  subroutine projector(ks, which, tion, indion, cstat, iend, fgend)
    implicit none
    type(TDKS), intent(inout) :: ks
    integer, intent(in) :: which, tion, indion, iend
    integer, intent(inout) :: cstat, fgend

    integer :: i, j
    real(kind=q) :: r, dE, kbT, popmax, lower, upper, norm_c
    real(kind=q) :: popBoltz(inp%NBASIS)

    kbT = inp%TEMP * BOLKEV
    call random_number(r)

    !! hopping probability with Boltzmann factor
    popBoltz(which) = REAL(CONJG(ks%psi_c(which)) * ks%psi_c(which), kind=q)
    !if (Mod(tion,100)==0) then
    !  write(89,*) which,tion,popBoltz(which)
    !end if
    dE = ((ks%eigKs(which,tion) + ks%eigKs(which,tion+1)) - &
          (ks%eigKs(cstat,tion) + ks%eigKs(cstat, tion+1))) /2.0_q
    if (inp%LHOLE) then
      dE = -dE
    end if
    if ( dE > 0.0_q ) then
      popBoltz(which) = popBoltz(which) * exp(-dE / kbT)
    end if
    !write(*,*) which,popBoltz(which)

    !! project in/out
    if (r <= popBoltz(which)) then
      do i=1, inp%NBASIS
        ks%psi_c(i) = con%cero
      end do
      ks%psi_c(which) = con%uno

      if (fgend==0 .AND. which==iend) then
        ks%recom_pops(cstat, indion + 1:inp%NAMDTIME)=ks%recom_pops(cstat, indion + 1:inp%NAMDTIME)+1.0_q
        fgend=-1
      end if

      cstat = which

    else
    !! project out
      ks%psi_c(which) = con%cero
      norm_c=REAL(SUM(CONJG(ks%psi_c)*ks%psi_c), kind=q) !! 1-|c_a|^2
      do i=1, inp%NBASIS
        ks%psi_c(i) = ks%psi_c(i) / DSQRT(norm_c)
      end do
      !ks%psi_c(which) = cero
      !if (cstat==which) then
      !    call random_number(r)
      !    do i=1,inp%NBASIS
      !        if (i==1) then
      !            lower = 0.0_q
      !            upper = DCONJG(ks%psi_c(i)) * ks%psi_c(i)
      !        else
      !            lower = upper
      !            upper = upper+DCONJG(ks%psi_c(i)) * ks%psi_c(i)

      !        end if

      !        if(lower <= r .AND. r< upper) then
      !            if (fgend==0 .AND. i==iend) then

      !                ks%recom_pops(cstat,tion + 1:inp%RTIME)=ks%recom_pops(cstat, tion + 1:inp%RTIME)+1.0_q
      !                fgend=-1
      !            end if
      !            cstat=i
      !            exit
      !        end if

      !    end do

      !end if
    end if

  end subroutine

  subroutine runDISH(ks)
    implicit none

    type(TDKS), intent(inout) :: ks

    integer :: i
    integer :: istat, iend

    ! ks%dish_pops = 0.0_q
    ! ks%recom_pops = 0.0_q
    istat = inp%INIBAND
    if (inp%LHOLE) then
      iend=inp%NBASIS
    else
      iend=1
    end if

    ! initialize the random seed for ramdom number production
    call init_random_seed()
    !$OMP PARALLEL DO SHARED(ks,inp,istat), PRIVATE(i),REDUCTION(+:ks%dish_pops)
    !MPI is under development
    do i=1, inp%NTRAJ
      call mpirundish(ks, istat, iend)
    end do
    !$OMP END PARALLEL DO
    ks%dish_pops = ks%dish_pops / inp%NTRAJ
    ks%dish_pops(istat, 1) = 1.0_q

  end subroutine

  subroutine mpirundish(ks, istat, iend)
    implicit none

    type(TDKS), intent(inout) :: ks
    integer, intent(in)       :: iend
    integer, intent(inout)    :: istat

    integer :: fgend = 0 !! recomb indicator
    integer :: i, j, k, tion, indion
    integer :: which, cstat
    real(kind=q), dimension(inp%NBASIS) :: decmoment !! t_i(t)
    integer, dimension(inp%NBASIS)      :: shuffle
    real(kind=q), dimension(inp%NBASIS) :: COEFFISQ  !! |c_i(t)|^2
    real(kind=q), dimension(inp%NBASIS) :: DECOTIME  !! tau_i(t)
    ! At the first step, current state always equal initial state
    ks%psi_c = con%cero
    ks%psi_c(istat) = con%uno
    cstat = istat
    decmoment = 0.0_q

    do j=1,inp%NBASIS
      shuffle(j)=j
    end do

    indion = 1
    do indion=1, inp%NAMDTIME - 1
      ! NAMDTINI >= 2, tion+inp%NAMDTINI-1 <= NSW(FSSH), dump 1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !          S                  E
      ! E1--D1--E2--D2--E3--D3------EN--DN--XX
      !                         E1--D1--E2--D2--E3--D3------EN--DN--XX
      !                              S                      E
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! for DISH, one cycle has NSW-2 intervals.
      tion = 2 + MOD(indion+inp%NAMDTINI-1-2,inp%NSW-2) ! 2 <= tion <= 2+(NSW-3) = NSW-1

      call PropagationEle(ks, tion)
     
      call calcdect(ks, DECOTIME, COEFFISQ)   

      call whichtodec(DECOTIME, which, decmoment, shuffle)  
      decmoment = decmoment + inp%POTIM
        
      if ( which > 0 ) then
        call projector(ks, which, tion, indion, cstat, iend, fgend)
      end if
      
      ks%dish_pops(cstat, indion + 1) = ks%dish_pops(cstat, indion + 1) + 1.0_q
    end do

  end subroutine

  subroutine printDISH(ks)
    implicit none
    type(TDKS), intent(in) :: ks

    integer :: i, j, tion, indion, ierr
    character(len=48) :: buf
    character(len=48) :: out_fmt

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f9.4))" )' )  ks%ndim

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

end module
