Program main
  use prec
  use fileio
  use utils
  use parallel
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
  integer :: communicator = 0, color = 0
  integer :: lower, upper, ncount
  real(kind=qs) :: tottime
  character(len=256) :: buf

#ifdef ENABLEMPI
  call MPI_INIT(ierr)
  if (ierr /= 0) then
    write(*,*) "[E] MPI initialization failed. Aborting..."
    call MPI_FINALIZE(ierr)
    stop
  end if
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprog, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, iprog, ierr)
#endif

  call system_clock(count_rate=cr)
  call system_clock(count_max=cm)
  
  if (iprog == 0) call printWelcome()
  ! initialize the random seed for ramdom number production
  call init_random_seed(salt=iprog)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First, get user inputs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call inp%getInstance()
  call inp%getUserInp(nprog)

#ifdef ENABLEMPI
  ! MPI_COMM_SPLIT(COMM, COLOR, KEY, NEWCOMM, IERROR)
  color = MODULO(iprog, inp%NPAR)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, iprog, communicator, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  call MPI_COMM_SIZE(communicator, nprog, ierr)
  call MPI_COMM_RANK(communicator, iprog, ierr)
#endif
  call inp%setMPI(iprog, nprog, color, communicator)

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
  ! distribute tasks.
  call divideTasks(color, inp%NPAR, inp%NSAMPLE, lower, upper, ncount)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ns=lower, upper
    call inp%setIni(ns)
    if (iprog == 0) call inp%printUserInp()

    ! initiate KS matrix
    call system_clock(t1)
    call initTDKS(ks)
    call system_clock(t2)
    if (iprog == 0) then
      write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in initTDKS [s]:", &
                                              MODULO(t2-t1, cm) / REAL(cr)
      write(*,*) trim(buf)
    end if

    select case(inp%ALGO)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case ('FSSH')
      if (iprog == 0) then
        ! Time propagation
        t1 = t2
        call runSE(ks, olap)
        call system_clock(t2)
        write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in runSE [s]:", &
                                                MODULO(t2-t1, cm) / REAL(cr)
        write(*,*) trim(buf)

        t1 = t2
        call printSE(ks, olap)
        call system_clock(t2)
        write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in printSE [s]:", &
                                                MODULO(t2-t1, cm) / REAL(cr)
        write(*,*) trim(buf)
        ! Run surface hopping
        if (inp%LSHP) then
          t1 = t2
          call runSH(ks, olap)
          call system_clock(t2)
          write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in runSH [s]:", &
                                                  MODULO(t2-t1, cm) / REAL(cr)
          write(*,*) trim(buf)

          t1 = t2
          call printSH(ks, olap)
          call system_clock(t2)
          write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in printSH [s]:", &
                                                  MODULO(t2-t1, cm) / REAL(cr)
          write(*,*) trim(buf)
        end if
        if (inp%LSPACE) then
          t1 = t2
          call printMPFSSH(ks)
          call system_clock(t2)
          write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in printMPFSSH [s]:", &
                                                  MODULO(t2-t1, cm) / REAL(cr)
          write(*,*) trim(buf)
        end if
      end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    case ('DISH')
      t1 = t2
      call runDISH(ks, olap)
      call system_clock(t2)
      if (iprog == 0) then
        write(buf,'(A3,I3,A5,I4,A,T48,F11.3)') "MPI", color, " TINI", inp%NAMDTINI, ": CPU Time in runDISH [s]:", &
                                                MODULO(t2-t1, cm) / REAL(cr)
        write(*,*) trim(buf)
      end if
    end select
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call system_clock(ttot2)
  tottime = MODULO(ttot2-ttot1, cm) / REAL(cr)
#ifdef ENABLEMPI
  call MPI_REDUCE(tottime+0.0, tottime, 1, MPI_REAL, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
#endif
  
  if (iprog == 0 .AND. color == 0) then
    write(*,'(A)') "------------------------------------------------------------"
    write(buf,'(A,T48,F11.3)') "All Time Elapsed [s]:", tottime
    write(*,*) trim(buf)
  end if

#ifdef ENABLEMPI
  call MPI_COMM_FREE(communicator, ierr)
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
    write(*,'(A)') "|     NPAR, LBINOUT,                                                             |"
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
