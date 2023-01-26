Program main
  use prec
  use fileio
  use couplings
  use hamil
  use fssh
  use dish

  implicit none

  type(TDKS) :: ks
  type(overlap) :: olap, olap_sp

  real(kind=q) :: start, fin, startall, finall
  integer :: ns,cr,cm

  call system_clock(count_rate=cr)
  call system_clock(count_max=cm)
  call cpu_time(startall)
  call printWelcome()
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
  call TDCoupIJ(olap, olap_sp)
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
      if (inp%LSPACE) then
        start=fin
        call printMPFSSH(ks, olap_sp)
        write(*,'(A, F8.2)') "CPU Time in printMPFSSH [s]:", fin - start
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
  call cpu_time(finall)
  write(*,'(A)') "------------------------------------------------------------"
  write(*,'(A, F8.2)') "All Time Elapsed [s]:", finall - startall

contains

  subroutine printWelcome()
    ! Varsity
    write(*,'(A)') "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,'(A)') "|  ____      ____  ________  _____       ______    ___   ____    ____  ________    _________    ___         |"
    write(*,'(A)') "| |_  _|    |_  _||_   __  ||_   _|    .' ___  | .'   `.|_   \  /   _||_   __  |  |  _   _  | .'   `.       |"
    write(*,'(A)') "|   \ \  /\  / /    | |_ \_|  | |     / .'   \_|/  .-.  \ |   \/   |    | |_ \_|  |_/ | | \_|/  .-.  \      |"
    write(*,'(A)') "|    \ \/  \/ /     |  _| _   | |   _ | |       | |   | | | |\  /| |    |  _| _       | |    | |   | |      |"
    write(*,'(A)') "|     \  /\  /     _| |__/ | _| |__/ |\ `.___.'\\  `-'  /_| |_\/_| |_  _| |__/ |     _| |_   \  `-'  /      |"
    write(*,'(A)') "|      \/  \/     |________||________| `.____ .' `.___.'|_____||_____||________|    |_____|   `.___.'       |"
    write(*,'(A)') "|                         ____  ____  ________       ____  _____       _       ____    ____  ______     _   |"
    write(*,'(A)') "|                        |_   ||   _||_   __  |     |_   \|_   _|     / \     |_   \  /   _||_   _ `.  | |  |"
    write(*,'(A)') "|                          | |__| |    | |_ \_|______ |   \ | |      / _ \      |   \/   |    | | `. \ | |  |"
    write(*,'(A)') "|                          |  __  |    |  _|  |______|| |\ \| |     / ___ \     | |\  /| |    | |  | | | |  |"
    write(*,'(A)') "|                         _| |  | |_  _| |_          _| |_\   |_  _/ /   \ \_  _| |_\/_| |_  _| |_.' / |_|  |"
    write(*,'(A)') "|                        |____||____||_____|        |_____|\____||____| |____||_____||_____||______.'  (_)  |"
    write(*,'(A)') "| Version: Jan. 2023.                                                                                       |"
    write(*,'(A)') "| Authors:                                                                                                  |"
    write(*,'(A)') "| Big Thanks to:                                                                                            |"
    write(*,'(A)') "|                                                                                                           |"
    write(*,'(A)') "| Supported input paramters:                                                                                |"
    write(*,'(A)') "|     bmin, bmax,                                                                                           |"
    write(*,'(A)') "|     nsample, ntraj, nsw, nelm,                                                                            |"
    write(*,'(A)') "|     temp, namdtime, potim,                                                                                |"
    write(*,'(A)') "|     lhole, lshp, algo, algo_int, lcpext,                                                                  |"
    write(*,'(A)') "|     lspace, nacbasis, nacele,                                                                             |"
    write(*,'(A)') "|     rundir, tbinit, diinit, spinit,                                                                       |"
    write(*,'(A)') "|     debuglevel                                                                                            |"
    write(*,'(A)') "| Supported input files:                                                                                    |"
    write(*,'(A)') "|     inp, INICON, DEPHTIME, ACSPACE                                                                        |"
    write(*,'(A)') "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,'(A)') "                                                                                                             "
  end subroutine

end Program
