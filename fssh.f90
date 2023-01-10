module fssh
  use prec
  use fileio
  use utils
  use hamil
  use TimeProp
  implicit none

  private :: whichToHop, calcprop

  contains

  subroutine whichToHop(indion, ks, which)
    implicit none

    integer, intent(in) :: indion
    integer, intent(inout) :: which
    type(TDKS), intent(in) :: ks

    integer :: i
    real(kind=q) :: lower, upper, r

    which = 0
    call random_number(r)


    lower = 0
    upper = ks%sh_prop(i, indion+1)
    do i=1, ks%ndim
      if (lower <= r .AND. r < upper) then
        which = i
        exit
      end if
      lower = upper
      upper = upper + ks%sh_prop(i, indion+1)
    end do

  end subroutine

  subroutine calcprop(tion, indion, cstat, ks)
    implicit none

    type(TDKS), intent(inout) :: ks
    integer, intent(in) :: tion, indion
    integer, intent(in) :: cstat

    integer :: i, j
    real(kind=q) :: Akk      !! |c_k|^2
    real(kind=q) :: dE, kbT

    kbT = inp%TEMP * BOLKEV

    Akk = REAL(CONJG(ks%psi_a(cstat, indion+1)) * ks%psi_a(cstat, indion+1), kind=q)
    ! Bkm = REAL(CONJG(Akm) * Ckm)
    !ks%Bkm = 2. * REAL(CONJG(ks%psi_a(cstat, tion)) * ks%psi_a(:, tion) * &
    !                ks%NAcoup(cstat, :, tion))
    !Changed the matrix structure for NAC. 
    !Now NAC(i,j) is stored as  ks%NAcoup(j,i,time)
    ks%Bkm = 2.0_q * REAL(CONJG(ks%psi_a(cstat, indion+1)) * ks%psi_a(:, indion+1) * &
                    ks%NAcoup(:, cstat, tion), kind=q)


    ks%sh_prop(:, indion+1) = ks%Bkm / Akk * inp%POTIM


   ! if (inp%LHOLE) then
   !   do i = 1, cstat
   !     dE = ks%eigKs(cstat, tion) - ks%eigKs(i,tion)
   !     ks%sh_prop(i,tion) = ks%sh_prop(i,tion) * exp(-dE / kbT)
   !   end do
   ! else
   !   do i=cstat, ks%ndim
   !     dE = ks%eigKs(i,tion) - ks%eigKs(cstat, tion)
   !     ks%sh_prop(i,tion) = ks%sh_prop(i,tion) * exp(-dE / kbT)
   !   end do
   ! end if
       
    do i = 1, ks%ndim
      dE = ((ks%eigKs(i,tion) + ks%eigKs(i,tion+1)) - &
            (ks%eigKs(cstat,tion)+ ks%eigKs(cstat,tion+1)))    &
           / 2.0_q
      if (inp%LHOLE) then
        dE = -dE
      end if
      if (dE>0) then
        ks%sh_prop(i,indion+1) = ks%sh_prop(i,indion+1) * exp(-dE / kbT)
      end if
    end do

    forall (i=1:ks%ndim, ks%sh_prop(i,indion+1) < 0) ks%sh_prop(i,indion+1) = 0.0_q
    ! write(*,*) (ks%Bkm(i), i=1, ks%ndim) 
    ! write(*,*) (ks%sh_prop(i, tion), i=1, ks%ndim) 

  end subroutine

  ! subroutine calcProb(tion, indion, ks)
  !   implicit none

  !   type(TDKS), intent(inout) :: ks
  !   integer, intent(in) :: tion, indion

  !   integer :: i, j
  !   real(kind=q) :: Akk_p, Akk_c      !! |c_k|^2
  !   real(kind=q), dimension(ks%ndim) :: Bkm_p, Bkm_c
  !   real(kind=q) :: dE, kbT, edt

  !   kbT = inp%TEMP * BOLKEV
  !   edt = inp%POTIM / inp%NELM

  !   do i = 1, ks%ndim
  !     Akk_p = REAL(CONJG(ks%psi_p(i)) * ks%psi_p(i), kind=q)
  !     Akk_c = REAL(CONJG(ks%psi_c(i)) * ks%psi_c(i), kind=q)

  !     !Changed the matrix structure for NAC. 
  !     !Now NAC(i,j) is stored as  ks%NAcoup(j,i,time)
  !     Bkm_c = REAL(CONJG(ks%psi_c(i)) * ks%psi_c(:) * &
  !                  ks%NAcoup(:, i, tion), kind=q)
  !     Bkm_p = REAL(CONJG(ks%psi_p(i)) * ks%psi_p(:) * &
  !                  ks%NAcoup(:, i, tion), kind=q)

  !     ks%sh_prob(:, i, indion) = (Bkm_p / Akk_p + Bkm_c / Akk_c) * edt
        
  !     do j = 1, ks%ndim
  !       dE = ((ks%eigKs(j,tion) + ks%eigKs(j,tion+1)) - &
  !             (ks%eigKs(i,tion)+ ks%eigKs(i,tion+1)))    &
  !            interpolate((TELE - 0.5_q) / inp%NELM, &
  !                                 ks%eigKs(i,RTIME), ks%eigKs(i,XTIME))
  !           / 2.0_q
  !       if (inp%LHOLE) then
  !         dE = -dE
  !       end if
  !       if (dE>0) then
  !         ks%sh_prob(j,i,indion) = ks%sh_prob(j,i,indion) * exp(-dE / kbT)
  !       end if
  !     end do

  !     forall (j=1:ks%ndim, ks%sh_prob(j,i,indion) < 0) ks%sh_prob(j,i,indion) = 0.0_q
  !   end do

  ! end subroutine


  subroutine runSE(ks)
    implicit none
    type(TDKS), intent(inout) :: ks
    integer :: tion, indion
    integer :: i, j
    integer :: istat, iend, cstat
    real(kind=q) :: norm

    ! init
    istat = inp%INIBAND
    ! ks%psi_a = con%cero
    ! ks%psi_a(istat, 1) = con%uno
    ! ks%pop_a = 0.0_q 
    ! ks%pop_a(istat, 1) = 1.0_q
    cstat = istat

    ! Enter the Loop
    do indion=1, inp%NAMDTIME - 1

      ks%psi_p = ks%psi_c
      tion = indion + inp%NAMDTINI - 1

      call PropagationEle(ks, tion)

      ks%pop_a(:, indion+1) = REAL(CONJG(ks%psi_c) * ks%psi_c, kind=q)
      ks%psi_a(:, indion+1) = ks%psi_c 

    end do
  end subroutine

  ! calculate surface hopping probabilities
  subroutine runSH(ks)
    implicit none

    type(TDKS), intent(inout) :: ks
    integer :: i, j, tion, indion
    integer :: istat, cstat, which
    real(kind=q), dimension(inp%NBASIS) :: SHPROB

    istat = inp%INIBAND
    ! ks%sh_pops = 0.0_q
    ! ks%sh_prop = 0.0_q

    ! initialize the random seed for ramdom number production
    call init_random_seed()

    do i=1, inp%NTRAJ
      ! in the first step, current step always equal initial step
      cstat = istat
      do indion=1, inp%NAMDTIME - 1
        tion = indion + inp%NAMDTINI - 1
        SHPROB = ks%sh_prob(:, cstat, indion)

        call calcprop(tion, indion, cstat, ks)

        call whichToHop(indion, ks, which)

        if (which > 0) then
          cstat = which
        end if
        ks%sh_pops(cstat, indion+1) = ks%sh_pops(cstat, indion+1) + 1
      end do
    end do

    ks%sh_pops = ks%sh_pops / inp%NTRAJ
    ks%sh_pops(istat, 1) = 1.0_q

    ! do tion=1, inp%NAMDTIME
    !   write(*,*) (ks%sh_pops(i,tion), i=1, ks%ndim)
    ! end do
  end subroutine

  subroutine printSE(ks)
    implicit none
    type(TDKS), intent(in) :: ks

    integer :: i, j, ierr
    integer :: tion, indion
    character(len=48) :: buf
    character(len=48) :: out_fmt

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f9.4))" )' )  ks%ndim
    
    open(unit=25, file='PSICT.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    open(unit=26, file='POPRT.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: PSICT file I/O error!"
      stop
    end if

    do indion=1, inp%NAMDTIME
      tion = indion + inp%NAMDTINI - 1
      write(unit=25, fmt=out_fmt) indion * inp%POTIM, SUM(ks%eigKs(:,tion) * ks%pop_a(:,indion)), &
                                  (ks%psi_a(i,indion), i=1, ks%ndim)
      write(unit=26, fmt=out_fmt) indion * inp%POTIM, SUM(ks%eigKs(:,tion) * ks%pop_a(:,indion)), &
                                  (ks%pop_a(i,indion), i=1, ks%ndim)
    end do

    close(25)
    close(26)

  end subroutine


  subroutine printSH(ks)
    implicit none
    type(TDKS), intent(in) :: ks

    integer :: i, j, ierr
    integer :: tion, indion
    character(len=48) :: buf
    character(len=48) :: out_fmt

    write(buf, *) inp%NAMDTINI
    write (out_fmt, '( "(f13.2,f11.6, ", I5, "(f9.4))" )' )  ks%ndim
    
    open(unit=24, file='SHPROP.' // trim(adjustl(buf)), status='unknown', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: SHPROP file I/O error!"
      stop
    end if

    do indion=1, inp%NAMDTIME
      tion = indion + inp%NAMDTINI - 1
      write(unit=24, fmt=out_fmt) indion * inp%POTIM, &
                                  SUM(ks%eigKs(:,tion) * ks%sh_pops(:,indion)), &
                                  (ks%sh_pops(i,indion), i=1, ks%ndim)
    end do

    close(24)

  end subroutine

end module
