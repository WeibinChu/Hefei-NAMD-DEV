Program main
  use prec
  use fileio
  use couplings
  use hamil
  use TimeProp
  use fssh
  use dish

  implicit none

  type(TDKS) :: ks
  type(overlap) :: olap

  real(kind=q) :: start, fin
  integer :: ns,cr,cm

  CALL system_clock(count_rate=cr)
  CALL system_clock(count_max=cm)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First, get user inputs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call inp%getInstance()
  call inp%getUserInp()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Secondly, get couplings
  ! In the very first run, the following subroutine will calculate the
  ! NA-couplings from WAVECARs and then write it to a binary file called
  ! COUPCAR.  From the second run on, the subroutine will just read the
  ! NA-couplings from the file. However, for a general NAMD run, the file is way
  ! too huge, the solution is to write only the information we need to another
  ! plain text file. If such files exist (set LCPTXT = .TRUE. in the inp), then
  ! we may skip the huge binary file and read the plain text file instead. This
  ! is done in the subroutine 'initTDKS'.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call TDCoupIJ(olap)
  ! write(*,*) "T_coup: ", fin - start
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ns=1, inp%NSAMPLE
    call inp%setIni(ns)
    call inp%printUserInp()

    ! initiate KS matrix
    call cpu_time(start)
    call initTDKS(ks, olap)
    call cpu_time(fin)
    write(*,'(A, F8.2)') "CPU Time in initTDKS [s]:", fin - start

    select case(inp%ALGO)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case ('FSSH')
      ! Time propagation
      start=fin
      call runSE(ks)
      call cpu_time(fin)
      write(*,'(A, F8.2)') "CPU Time in runSE [s]:", fin - start

      start=fin
      call printSE(ks)
      call cpu_time(fin)
      write(*,'(A, F8.2)') "CPU Time in printSE [s]:", fin - start
      ! Run surface hopping
      if (inp%LSHP) then
        start=fin
        call runSH(ks)
        call cpu_time(fin)
        write(*,'(A, F8.2)') "CPU Time in runSH [s]:", fin - start

        start=fin
        call printSH(ks)
        call cpu_time(fin)
        write(*,'(A, F8.2)') "CPU Time in printSH [s]:", fin - start
      end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case ('DISH')
      start=fin
      call runDISH(ks)
      call cpu_time(fin)
      write(*,'(A, F8.2)') "CPU Time in runDISH [s]:", fin - start

      start=fin
      call printDISH(ks)
      call cpu_time(fin)
      write(*,'(A, F8.2)') "CPU Time in printDISH [s]:", fin - start
    end select
  end do

! contains

!   subroutine printTime(func, arg, message)
!     procedure, intent(inout) :: func
!     type(TDKS), intent(inout) :: arg
!     character(len=48), intent(in) :: message
!     real(kind=q) :: start, fin

!     call cpu_time(start)
!     call func(arg)
!     call cpu_time(fin)
!     write(*, 'A, F8.2') message, fin - start
!   end subroutine

end Program
