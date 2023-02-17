Program main
  use prec
  use fileio
  use couplings
  use hamil
  use fssh
  use dish
#ifdef ENABLEMPI
  use mpi
#endif

  implicit none

  type(TDKS) :: ks
  type(overlap) :: olap, olap_sp

  integer :: ns, cr, cm, t1, t2, ttot1, ttot2
  integer :: nprog = 1, iprog = 0, ierr

#ifdef ENABLEMPI
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprog, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, iprog, ierr)
#endif

  call system_clock(count_rate=cr)
  call system_clock(count_max=cm)
  
  if (iprog == 0) call printWelcome()
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First, get user inputs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call inp%getInstance()
  call inp%getUserInp(nprog, iprog)

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
  call system_clock(ttot1)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ns=1, inp%NSAMPLE
    call inp%setIni(ns)
    if (iprog == 0) call inp%printUserInp()

    ! initiate KS matrix
    call system_clock(t1)
    call initTDKS(ks, olap)
    call system_clock(t2)
    if (iprog == 0) write(*,'(A, T31, F11.3)') "CPU Time in initTDKS [s]:", MOD(t2-t1, cm) / REAL(cr)

    select case(inp%ALGO)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case ('FSSH')
      if (iprog == 0) then
        ! Time propagation
        t1 = t2
        call runSE(ks)
        call system_clock(t2)
        write(*,'(A, T31, F11.3)') "CPU Time in runSE [s]:", MOD(t2-t1, cm) / REAL(cr)

        t1 = t2
        call printSE(ks)
        call system_clock(t2)
        write(*,'(A, T31, F11.3)') "CPU Time in printSE [s]:", MOD(t2-t1, cm) / REAL(cr)
        ! Run surface hopping
        if (inp%LSHP) then
          t1 = t2
          call runSH(ks)
          call system_clock(t2)
          write(*,'(A, T31, F11.3)') "CPU Time in runSH [s]:", MOD(t2-t1, cm) / REAL(cr)

          t1 = t2
          call printSH(ks)
          call system_clock(t2)
          write(*,'(A, T31, F11.3)') "CPU Time in printSH [s]:", MOD(t2-t1, cm) / REAL(cr)
        end if
        if (inp%LSPACE) then
          t1 = t2
          call printMPFSSH(ks)
          write(*,'(A, T31, F11.3)') "CPU Time in printMPFSSH [s]:", MOD(t2-t1, cm) / REAL(cr)
        end if
      end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case ('DISH')
      t1 = t2
      call runDISH(ks, olap)
      call system_clock(t2)
      if (iprog == 0) then
        write(*,'(A, T31, F11.3)') "CPU Time in runDISH [s]:", MOD(t2-t1, cm) / REAL(cr)

        t1 = t2
        call printDISH(ks)
        call system_clock(t2)
        write(*,'(A, T31, F11.3)') "CPU Time in printDISH [s]:", MOD(t2-t1, cm) / REAL(cr)

        if (inp%LSPACE) then
          t1 = t2
          call printMPDISH(ks)
          write(*,'(A, T31, F11.3)') "CPU Time in printMPFSSH [s]:", MOD(t2-t1, cm) / REAL(cr)
        end if
      end if
    end select
  end do
  call system_clock(ttot2)
  if (iprog == 0) write(*,'(A)') "------------------------------------------------------------"
  if (iprog == 0) write(*,'(A, T31, F11.3)') "All Time Elapsed [s]:", MOD(ttot2-ttot1, cm) / REAL(cr)

#ifdef ENABLEMPI
  call MPI_FINALIZE(ierr)
#endif

contains

  subroutine printWelcome()
    implicit none
    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone
    ! Big
    write(*,'(A)') "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,'(A)') "| __          ________ _      _____ ____  __          __ ______   _______ ____   |"
    write(*,'(A)') "| \ \        / /  ____| |    / ____/ __ \|  \        /  |  ____| |__   __/ __ \  |"
    write(*,'(A)') "|  \ \  /\  / /| |__  | |   | |   | |  | | \ \      / / | |__       | | | |  | | |"
    write(*,'(A)') "|   \ \/  \/ / |  __| | |   | |   | |  | | |\ \    / /| |  __|      | | | |  | | |"
    write(*,'(A)') "|    \  /\  /  | |____| |___| |___| |__| | | \ \  / / | | |____     | | | |__| | |"
    write(*,'(A)') "|     \/_ \/ _ |______|______\_____\____/| |  \ \/ /  | |______|  _ |_|  \____/  |"
    write(*,'(A)') "|      | |  | |  ____|    | \ | |   /\   | |   \  /   | |  __ \  | |             |"
    write(*,'(A)') "|      | |__| | |__ ______|  \| |  /  \  | |    \/    | | |  | | | |             |"
    write(*,'(A)') "|      |  __  |  __|______| . ` | / /\ \ | |          | | |  | | | |             |"
    write(*,'(A)') "|      | |  | | |         | |\  |/ ____ \| |          | | |__| | |_|             |"
    write(*,'(A)') "|      |_|  |_|_|         |_| \_/_/    \_\_|          |_|_____/  (_)             |"
    write(*,'(A)') "|                                                                                |"
    write(*,'(A)') "| Version: Jan. 2023.                                                            |"
    write(*,'(A)') "| Authors:                                                                       |"
    write(*,'(A)') "| Big Thanks to:                                                                 |"
    write(*,'(A)') "|                                                                                |"
    write(*,'(A)') "| Supported input paramters:                                                     |"
    write(*,'(A)') "|     BMIN, BMAX,                                                                |"
    write(*,'(A)') "|     NSAMPLE, NTRAJ, NSW, NELM,                                                 |"
    write(*,'(A)') "|     TEMP, NAMDTIME, POTIM,                                                     |"
    write(*,'(A)') "|     LHOLE, LSHP, ALGO, ALGO_INT, LCPTXT,                                       |"
    write(*,'(A)') "|     LSPACE, NACBASIS, NACELE,                                                  |"
    write(*,'(A)') "|     NPARDISH, LBINOUT,                                                         |"
    write(*,'(A)') "|     RUNDIR, TBINIT, DIINIT, SPINIT,                                            |"
    write(*,'(A)') "|     DEBUGLEVEL                                                                 |"
    write(*,'(A)') "| Supported input files:                                                         |"
    write(*,'(A)') "|     inp, COUPCAR, EIGTXT, NATXT, INICON, DEPHTIME, ACSPACE                     |"
    write(*,'(A)') "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,'(A)') "                                                                                  "
    call date_and_time(date, time, zone)
    write(*,'(A21,X,A4,".",A2,".",A2,X,A2,":",A2,":",A6,X,"UTC",A5,/)') 'The program starts at',          &
                                                                        date(1:4), date(5:6), date(7:8),  &
                                                                        time(1:2), time(3:4), time(5:10), &
                                                                        zone
  end subroutine

end Program
