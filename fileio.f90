module fileio
  use omp_lib
  use prec

  implicit none

  type, private :: namdInfo
    logical, private :: isCreated = .false.
    integer :: NPROG !! TRAJPAR
    integer :: IPROG !! MPI prog index in one group, start from 0
    integer :: NPAR
    integer :: COLOR
    integer :: COMMUNICATOR

    integer :: BMIN
    integer :: BMAX
    integer :: NBADNS
    integer :: HORB        ! Highest energy orbit, usually N(e) or 1(h)
    integer :: LORB        ! Lowest  energy orbit, usually 1(e) or N(h)
    integer, allocatable, dimension(:,:) :: BASIS !! Active Space Basis Macrostate (nacele, nacbasis)
    integer :: NBASIS      ! No. of adiabatic states as basis
    integer :: INIBAND     ! inititial adiabatic state of excited electron/hole
    integer :: NAMDTINI    ! Initial time step of NAMD
    integer :: NAMDTIME    ! No. of steps of NAMD
    integer, allocatable, dimension(:), private :: NAMDTINI_A    ! No. of steps of NAMD
    integer, allocatable, dimension(:), private :: INIBAND_A     ! No. of steps of NAMD

    integer :: NSW         ! No. of MD steps (ticks)
    integer :: NTRAJ       ! No. of surface hopping trajectories
    integer :: NELM        ! No. of steps of electron wave propagation
    integer :: NSAMPLE     ! No. of steps of electron wave propagation

    integer :: NACBASIS
    integer :: NACELE
    real(kind=q) :: POTIM  ! Time step of MD run
    real(kind=q) :: TEMP   ! MD Temperature

    ! hole or electron surface hopping
    logical :: LHOLE
    ! whether to perform surface hopping, right now the value is .TRUE.
    logical :: LSHP
    character(len=256) :: ALGO !! SH Algorrithm
    integer :: ALGO_INT
    logical :: LBINOUT
    logical :: LCPTXT
    logical :: LSPACE
    ! running directories
    character(len=256) :: RUNDIR
    character(len=256) :: TBINIT
    character(len=256) :: DIINIT   !! input of pure dephasing time matrix
    ! DISH parameters
    real(kind=q), allocatable, dimension(:,:) :: DEPHMATR     !! Pure dephasing rate between any two adiabatic states from the file DEPHTIME, 1/fs unit, symmetric
    integer :: NPARDISH
    integer :: DEBUGLEVEL
  contains
    procedure :: getInstance
    procedure :: getUserInp, printUserInp
    procedure :: setIni, setHLORB, setMPI
    procedure, private :: checkUserInp
  end type

  type(namdInfo), public :: inp

contains

    subroutine getInstance(this)
      class(namdInfo), intent(inout) :: this
      if (this%isCreated) then
        return
      else
        this%isCreated = .TRUE.
      end if
    end subroutine

    subroutine setIni(this, value)
      class(namdInfo), intent(inout) :: this
      integer, intent(in) :: value
      this%INIBAND = this%INIBAND_A(value)
      this%NAMDTINI = this%NAMDTINI_A(value)
    end subroutine

    subroutine setHLORB(this, horb, lorb)
      class(namdInfo), intent(inout) :: this
      integer, intent(in) :: horb, lorb
      this%HORB = horb
      this%LORB = lorb
    end subroutine

    subroutine setMPI(this, iprog, nprog, color, communicator)
      implicit none
      class(namdInfo), intent(inout) :: this
      integer, intent(in) :: iprog, nprog, color, communicator

      this%IPROG = iprog
      this%NPROG = nprog
      this%COLOR = color
      this%COMMUNICATOR = communicator
    end subroutine

    subroutine getUserInp(this, nprog)
      implicit none

      class(namdInfo), intent(inout) :: this
      integer, intent(in) :: nprog

      ! local variables with the same name as those in "inp"
      integer :: npar      = 1

      integer :: bmin
      integer :: bmax

      integer :: nsw
      integer :: ntraj     = 1000
      integer :: nelm      = 10
      integer :: nsample   

      integer :: nacbasis  = 100
      integer :: nacele    = 1

      integer :: namdtime  = 0
      real(kind=q) :: potim= 1.0_q
      real(kind=q) :: temp = 300_q

      ! hole or electron surface hopping
      logical :: lhole          = .FALSE.
      logical :: lshp           = .TRUE.  !! surface hopping?
      character(len=16) :: algo = 'DISH'
      integer :: algo_int       = 0
      logical :: lcptxt         = .TRUE.
      logical :: lspace         = .FALSE.
      logical :: lbinout        = .FALSE.
      ! running directories
      character(len=256) :: rundir = 'run'
      character(len=256) :: tbinit = 'INICON'
      character(len=256) :: diinit = 'DEPHTIME'
      character(len=256) :: spinit = 'ACSPACE'

      character(len=256) :: debuglevel = 'I'

      namelist /NAMDPARA/ bmin, bmax,                          &
                          nsample, ntraj, nsw, nelm,           &
                          temp, namdtime, potim,               &
                          lhole, lshp, algo, algo_int, lcptxt, &
                          lspace, nacbasis, nacele,            &
                          npar, lbinout,                   &
                          rundir, tbinit, diinit, spinit,      &
                          debuglevel

      integer :: ierr, i, j, nthread = 1
      logical :: lext

      ! inp
      open(file="inp", unit=8, status='unknown', action='read', iostat=ierr)
      if ( ierr /= 0 ) then
        write(*,*) "[E] IOError: I/O error with input file: 'inp'"
        stop
      end if

      ! read inp
      read(unit=8, nml=NAMDPARA)
      close(unit=8)

      ! omp
      ! nthread = MIN(npardish, omp_get_max_threads())
      ! call omp_set_num_threads(nthread)

      ! DEBUGLEVEL
      select case (debuglevel)
      case ('D')
        this%DEBUGLEVEL = 0
      case ('I')
        this%DEBUGLEVEL = 1
      case ('W')
        this%DEBUGLEVEL = 2
      case ('E')
        this%DEBUGLEVEL = 3
      end select

      ! INICON
      allocate(this%INIBAND_A(nsample), this%NAMDTINI_A(nsample))
      inquire(file=tbinit, exist=lext)
      if (.NOT. lext) then
        write(*,*) "[E] IOError: File containing initial conditions does NOT exist!"
        stop
      else
        open(unit=9, file=tbinit, action='read')
        do i=1, nsample
          read(unit=9,fmt=*) this%NAMDTINI_A(i), this%INIBAND_A(i)
        end do
        close(9)
        ! assuming BMIN is big enough, deprecated
        if (this%INIBAND_A(1) >= bmin) then
          this%INIBAND_A = this%INIBAND_A - bmin + 1
        end if
      end if

      ! ACSPACE
      if (lspace) then
        allocate(this%BASIS(nacele, nacbasis))
        inquire(file=spinit, exist=lext)
        if (.NOT. lext) then
          write(*,*) "[E] IOError: File containing Active Space does NOT exist!"
          stop
        else
          open(unit=11, file=spinit, action='read')
          do i = 1, nacbasis
            read(unit=11, fmt=*) (this%BASIS(j,i), j=1,nacele)
          end do
          close(11)
          this%NBASIS = nacbasis
        end if
      else
        allocate(this%BASIS(bmax - bmin + 1, 1))
        this%BASIS = reshape([(i, i=bmin, bmax)], shape=[1,bmax-bmin+1])
        this%NBASIS = bmax - bmin + 1
      end if

      ! DEPHTIME
      if (algo == 'DISH') then
        allocate(this%DEPHMATR(this%NBASIS, this%NBASIS))   !! read in the pure dephasing time matrix, the diagonal elements are zero
        inquire(file=diinit, exist=lext)
        if (.NOT. lext) then
          write(*,*) "[E] IOError: File containing initial conditions of DISH does NOT exist!"
          stop
        else
          open(unit=10, file=diinit, action='read')
          read(unit=10, fmt=*) ((this%DEPHMATR(i,j), j=1, this%NBASIS), i=1, this%NBASIS)
          do i = 1, this%NBASIS
            do j = 1, this%NBASIS
              if (i /= j) then
                this%DEPHMATR(j,i) = 1.0_q / this%DEPHMATR(j,i)
              end if
            end do
          end do
          close(10)
        end if
      end if

      ! assign the parameters
      this%NPAR     = MIN(npar, nprog, nsample)

      this%BMIN     = bmin
      this%BMAX     = bmax
      this%NBADNS   = bmax - bmin + 1

      this%NSW      = nsw
      this%NTRAJ    = ntraj
      this%NELM     = nelm
      this%NSAMPLE  = nsample

      this%NAMDTIME = namdtime 
      this%POTIM    = potim
      this%TEMP     = temp

      this%LHOLE    = lhole
      this%LSHP     = lshp
      this%LCPTXT   = lcptxt
      this%ALGO     = algo
      this%ALGO_INT = algo_int
      this%LBINOUT  = lbinout

      this%RUNDIR   = trim(rundir)
      this%TBINIT   = trim(tbinit)
      this%DIINIT   = trim(diinit)

      this%LSPACE   = lspace
      this%NACBASIS = nacbasis
      this%NACELE   = nacele

      call this%checkUserInp()
    end subroutine

    subroutine checkUserInp(this)
      implicit none
      class(namdInfo), intent(in) :: this

      integer :: i 

      ! do some checking...
      ! put the following checks in the future version
      if (this%BMIN <= 0 .OR. this%BMAX <= 0 .OR. this%BMIN >= this%BMAX) then
        write(*,*) "[E] Please specify the correct BMIN/BMAX"
        stop
      end if

      if (this%ALGO == 'FSSH') then
        do i=1, this%NSAMPLE
          if (this%NAMDTINI_A(i) + this%NAMDTIME - 1 > this%NSW) then
            write(*,*) "[E] NAMDTIME too long..."
            stop
          else if (this%NAMDTINI_A(i) == 1) then
            write(*,*) "[E] NAMDTINI should > 1..."
            stop
          end if
        end do
      end if
    end subroutine

    ! Need a subroutine to print out all the input parameters
    subroutine printUserInp(this)
      implicit none
      class(namdInfo), intent(in) :: this
      character(len=4096) :: buf

      write(buf,'(A60,A)') "------------------------------------------------------------", new_line('')
      write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'NPAR',     ' = ', this%NPAR,             new_line('')
      write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'TRAJPAR',  ' = ', this%NPROG,            new_line('')
      write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'MPIGROUP', ' = ', this%COLOR,            new_line('')
      write(buf,'(A,A)')              trim(buf),                                           new_line('')
      write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'BMIN',     ' = ', this%BMIN,             new_line('')
      write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'BMAX',     ' = ', this%BMAX,             new_line('')
      write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'INIBAND',  ' = ', this%INIBAND,          new_line('')
      write(buf,'(A,A)')              trim(buf),                                           new_line('')
      write(buf,'(A,A30,A3,F10.1,A)') trim(buf), 'POTIM',    ' = ', this%POTIM,            new_line('')
      write(buf,'(A,A30,A3,F10.1,A)') trim(buf), 'TEMP',     ' = ', this%TEMP,             new_line('')
      write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'NAMDTINI', ' = ', this%NAMDTINI,         new_line('')
      write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'NAMDTIME', ' = ', this%NAMDTIME,         new_line('')
      write(buf,'(A,A)')              trim(buf) ,                                          new_line('')
      write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'NSW',      ' = ', this%NSW,              new_line('')
      write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'NTRAJ',    ' = ', this%NTRAJ,            new_line('')
      write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'NELM',     ' = ', this%NELM,             new_line('')
      write(buf,'(A,A)')              trim(buf),                                           new_line('')
      write(buf,'(A,A30,A3,L10,A)')   trim(buf), 'LHOLE',    ' = ', this%LHOLE,            new_line('')
      write(buf,'(A,A30,A3,L10,A)')   trim(buf), 'LSHP',     ' = ', this%LSHP,             new_line('')
      write(buf,'(A,A30,A3,A10,A)')   trim(buf), 'ALGO',     ' = ', TRIM(this%ALGO),       new_line('')
      write(buf,'(A,A30,A3,I10,A)')   trim(buf), 'ALGO_INT', ' = ', this%ALGO_INT,         new_line('')
      write(buf,'(A,A30,A3,L10,A)')   trim(buf), 'LBINOUT',  ' = ', this%LBINOUT,          new_line('')
      write(buf,'(A,A30,A3,L10,A)')   trim(buf), 'LCPTXT',   ' = ', this%LCPTXT,           new_line('')
      write(buf,'(A,A)')              trim(buf),                                           new_line('')
      write(buf,'(A,A30,A3,L10,A)')   trim(buf), 'LSPACE',   ' = ', this%LSPACE,           new_line('')
      if (this%LSPACE) then
        write(buf,'(A,A30,A3,I10,A)') trim(buf), 'NACBASIS', ' = ', this%NACBASIS,         new_line('')
        write(buf,'(A,A30,A3,I10,A)') trim(buf), 'NACELE',   ' = ', this%NACELE,           new_line('')
      end if
      write(buf,'(A,A60)') trim(buf), "------------------------------------------------------------"
      write(*,'(A)') trim(buf)
    end subroutine

end module
