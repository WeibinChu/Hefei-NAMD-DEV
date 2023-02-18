module fssh
  use prec
  use constants
  use fileio
  use utils
  use couplings
  use hamil
  use TimeProp
  implicit none

  private :: whichToHop, calcProb, calcProbwithDBC

  contains

  subroutine whichToHop(indion, ks, cstat, which)
    implicit none

    integer, intent(in) :: indion, cstat
    integer, intent(out) :: which
    type(TDKS), intent(in) :: ks

    integer :: i
    real(kind=q) :: lower, upper, r

    which = 0
    call random_number(r)

    do i=1, ks%ndim
      if (i == 1) then
        lower = 0.0_q
        ! upper = ks%sh_prop(i, indion+1)
        upper = ks%sh_prob(i, cstat, indion)
      else
        lower = upper
        ! upper = upper + ks%sh_prop(i, indion+1)
        upper = upper + ks%sh_prob(i, cstat, indion)
      end if
      if (lower <= r .AND. r < upper) then
        which = i
        exit
      end if
    end do

  end subroutine

  subroutine calcProb(tion, indion, tele, ks, olap)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(overlap), intent(in) :: olap
    integer, intent(in) :: tion, indion, tele

    integer :: i
    real(kind=q) :: Akk_p, Akk_c      !! |c_k|^2
    real(kind=q), dimension(ks%ndim) :: Bkm_p, Bkm_c
    real(kind=q) :: edt

    edt = inp%POTIM / inp%NELM

    ! Akk = 0 => P_ij = NaN
    do i = 1, ks%ndim
      Akk_p = REAL(CONJG(ks%psi_p(i)) * ks%psi_p(i), kind=q)
      Akk_c = REAL(CONJG(ks%psi_c(i)) * ks%psi_c(i), kind=q)

      !Changed the matrix structure for NAC. 
      !Now NAC(i,j) is stored as  olap%Dij(j,i,time)
      Bkm_c = REAL(CONJG(ks%psi_c(i)) * ks%psi_c(:) * &
                   calcDij(olap%Dij(:,i,tion-1),olap%Dij(:,i,tion),olap%Dij(:,i,tion+1),tele), kind=q)
      Bkm_p = REAL(CONJG(ks%psi_p(i)) * ks%psi_p(:) * &
                   calcDij(olap%Dij(:,i,tion-1),olap%Dij(:,i,tion),olap%Dij(:,i,tion+1),tele-1), kind=q)

      ks%sh_prob(:, i, indion) = ks%sh_prob(:, i, indion) + &
                                 (Bkm_p / Akk_p + Bkm_c / Akk_c) * edt
    end do

  end subroutine
    
  subroutine calcProbwithDBC(tion, indion, ks, olap)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(overlap), intent(in) :: olap
    integer, intent(in) :: tion, indion

    integer :: i, j
    real(kind=q) :: dE(ks%ndim), kbT

    kbT = inp%TEMP * BOLKEV

    ! Detail Balance Condition
    ! Beware there is a tricky part: 
    ! P_ij * EXP(-dE/kT) != SUM(P_ij/N * EXP(-dE/N/kT)) = P_ij * EXP(-dE/NkT)
    do i = 1, ks%ndim        
      dE = ((olap%Eig(:,tion) + olap%Eig(:,tion+1)) - &
            (olap%Eig(i,tion) + olap%Eig(i,tion+1))) / 2.0_q
      if (inp%LHOLE) then
        dE = -dE
      end if
      where ( dE > 0.0_q )
        ks%sh_prob(:,i,indion) = ks%sh_prob(:,i,indion) * exp(-dE / kbT)
      end where

      forall (j=1:ks%ndim, ks%sh_prob(j,i,indion) < 0) ks%sh_prob(j,i,indion) = 0.0_q
    end do
  end subroutine


  subroutine runSE(ks, olap)
    implicit none
    type(TDKS), intent(inout) :: ks
    type(overlap), intent(in) :: olap

    integer :: tion, indion, tele
    integer :: istat, cstat
    real(kind=q) :: norm, edt
    logical :: init

    ! init
    istat = inp%INIBAND
    cstat = istat
    edt = inp%POTIM / inp%NELM
    init = .TRUE.

    ! Enter the Loop
    do indion=1, inp%NAMDTIME - 1

      tion = indion + inp%NAMDTINI - 1

      do tele = 1, inp%NELM
        select case (inp%ALGO_INT)
        case (0)
          ks%psi_p = ks%psi_c
          call make_hamil(tion, tele, olap, ks)
          call Trotter(ks, edt)
        case (11) ! this method will accumulate error! May renormalize wfc.
          ks%psi_p = ks%psi_c
          call make_hamil(tion, tele, olap, ks)
          call Euler(ks, edt)
        case (10)
          if (init) then
            ks%psi_p = ks%psi_c
            call make_hamil(tion, tele, olap, ks)
            call Euler(ks, edt)
            init = .FALSE.
          else
            call make_hamil2(tion, tele-1, olap, ks)
            call EulerMod(ks, edt)
          end if
        case (2)
          ks%psi_p = ks%psi_c
          call make_hamil(tion, tele, olap, ks)
          call DiagonizeMat(ks, edt)
          ks%psi_c = HamPsi(ks%ham_c, ks%psi_c, 'N')
        end select
        if (inp%LSHP) call calcProb(tion, indion, tele, ks, olap)
      end do
      norm = REAL(SUM(CONJG(ks%psi_c) * ks%psi_c), kind=q) 
      if ( norm <= 0.99_q .OR. norm >= 1.01_q)  then
        write(*,*) "[E] Error in Electronic Propagation"
        stop
      end if

      if (inp%LSHP) call calcProbwithDBC(tion, indion, ks, olap)
      ks%pop_a(:, indion+1) = REAL(CONJG(ks%psi_c) * ks%psi_c, kind=q)
      ks%psi_a(:, indion+1) = ks%psi_c 

    end do
  end subroutine

  ! calculate surface hopping probabilities
  subroutine runSH(ks, olap)
    implicit none

    type(TDKS), intent(inout) :: ks
    type(overlap), intent(in) :: olap
    integer :: i, tion, indion
    integer :: istat, cstat, which

    istat = inp%INIBAND

    ! initialize the random seed for ramdom number production
    call init_random_seed()

    do i=1, inp%NTRAJ
      ! in the first step, current step always equal initial step
      cstat = istat
      do indion=1, inp%NAMDTIME - 1
        tion = indion + inp%NAMDTINI - 1

        call whichToHop(indion, ks, cstat, which)

        if (which > 0) then
          cstat = which
        end if
        ks%sh_pops(cstat, indion+1) = ks%sh_pops(cstat, indion+1) + 1
      end do
    end do

    ks%sh_pops = ks%sh_pops / inp%NTRAJ
    ks%sh_pops(istat, 1) = 1.0_q

  end subroutine

  subroutine printSE(ks, olap)
    implicit none
    type(TDKS), intent(in) :: ks
    type(overlap), intent(in) :: olap

    integer :: i, ierr
    integer :: tion, indion
    character(len=48) :: buf
    character(len=48) :: out_fmt, out_fmt_cmplx

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f11.6))" )' )  ks%ndim
    write (out_fmt_cmplx, '( "(f13.2,f11.6, ", I5, "(SS,f11.6,SP,f9.6", A3, "))" )' )  ks%ndim, '"i"'
    
    open(unit=25, file='PSICT.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    open(unit=26, file='POPRT.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: PSICT file I/O error!"
      stop
    end if

    do indion=1, inp%NAMDTIME
      tion = indion + inp%NAMDTINI - 1
      write(unit=25, fmt=out_fmt_cmplx) indion * inp%POTIM, &
                                  SUM(olap%Eig(:,tion) * ks%pop_a(:,indion)), &
                                  (ks%psi_a(i,indion), i=1, ks%ndim)
      write(unit=26, fmt=out_fmt) indion * inp%POTIM, &
                                  SUM(olap%Eig(:,tion) * ks%pop_a(:,indion)), &
                                  (ks%pop_a(i,indion), i=1, ks%ndim)
    end do

    close(25)
    close(26)

  end subroutine

  subroutine printSH(ks, olap)
    implicit none
    type(TDKS), intent(in) :: ks
    type(overlap), intent(in) :: olap

    integer :: i, ierr
    integer :: tion, indion
    character(len=48) :: buf
    character(len=48) :: out_fmt

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f11.6))" )' )  ks%ndim
    
    open(unit=24, file='SHPROP.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: SHPROP file I/O error!"
      stop
    end if

    do indion=1, inp%NAMDTIME
      tion = indion + inp%NAMDTINI - 1
      write(unit=24, fmt=out_fmt) indion * inp%POTIM, &
                                  SUM(olap%Eig(:,tion) * ks%sh_pops(:,indion)), &
                                  (ks%sh_pops(i,indion), i=1, ks%ndim)
    end do

    close(24)

  end subroutine

  subroutine printMPFSSH(ks)
    implicit none
    type(TDKS), intent(inout) :: ks

    integer :: i, ierr
    integer :: bi, tion, indion
    integer, dimension(inp%NACELE) :: bands

    character(len=48) :: buf
    character(len=48) :: out_fmt

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f11.6))" )' )  inp%NBADNS
    
    open(unit=51, file='MPPOPRT.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: MPPOPRT file I/O error!"
      stop
    end if

    ks%mppop_a = 0.0_q
    do i=1, ks%ndim
      bands = inp%BASIS(:,i) - inp%BMIN + 1
      do bi=1, inp%NACELE
        ks%mppop_a(bands(bi),:) = ks%mppop_a(bands(bi),:) + &
                                  ks%pop_a(i,:)
      end do
    end do

    do indion=1, inp%NAMDTIME
      tion = indion + inp%NAMDTINI - 1
      write(unit=51, fmt=out_fmt) indion * inp%POTIM, & 
                                  ! SUM(ks%eigKs(:,tion) * ks%pop_a(:,indion)) / inp%NACELE, &
                                  (ks%mppop_a(i,indion), i=1, inp%NBADNS)
    end do

    close(51)

    !! Surfacr Hopping
    if (inp%LSHP) then
      open(unit=52, file='MPSHPROP.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
      if (ierr /= 0) then
        write(*,*) "[E] IOError: MPSHPROP file I/O error!"
        stop
      end if
  
      ks%sh_mppops = 0.0_q
      do i=1, ks%ndim
        bands = inp%BASIS(:,i) - inp%BMIN + 1
        do bi=1, inp%NACELE
          ks%sh_mppops(bands(bi),:) = ks%sh_mppops(bands(bi),:) + &
                                      ks%sh_pops(i,:)
        end do
      end do
  
      do indion=1, inp%NAMDTIME
        tion = indion + inp%NAMDTINI - 1
        write(unit=52, fmt=out_fmt) indion * inp%POTIM, & 
                                    ! SUM(ks%eigKs(:,tion) * ks%sh_pops(:,indion)) / inp%NACELE, &
                                    (ks%sh_mppops(i,indion), i=1, inp%NBADNS)
      end do
  
      close(52)
    end if

  end subroutine

end module
