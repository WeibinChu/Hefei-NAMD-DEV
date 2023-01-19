module fileio
  use prec

  implicit none

  type, private :: namdInfo
    logical, private :: isCreated = .false.
    integer :: BMIN
    integer :: BMAX
    integer, allocatable, dimension(:,:) :: BASIS !! Active Space Basis Macrostate (nacbasis, nacele)
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
    ! whether to perform DISH, right now the value is .TRUE.
    character(len=256) :: ALGO !! SH Algorrithm
    logical :: LCPTXT
    logical :: LSPACE
    ! running directories
    character(len=256) :: RUNDIR
    character(len=256) :: TBINIT
    character(len=256) :: DIINIT   !! input of pure dephasing time matrix
    real(kind=q), allocatable, dimension(:,:) :: DEPHMATR     !! Pure dephasing rate between any two adiabatic states from the file DEPHTIME, 1/fs unit, symmetric
    integer :: DEBUGLEVEL
  contains
    procedure :: getInstance
    procedure :: getUserInp, printUserInp
    procedure :: setIni
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

    subroutine getUserInp(this)
      implicit none

      class(namdInfo), intent(inout) :: this

      ! local variables with the same name as those in "inp"
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
      logical :: lhole  = .FALSE.
      ! surface hopping?
      logical :: lshp   = .TRUE.
      character(len=256) :: algo = 'DISH'
      logical :: lcpext = .TRUE.
      logical :: lspace = .FALSE.
      ! running directories
      character(len=256) :: rundir = 'run'
      character(len=256) :: tbinit = 'INICON'
      character(len=256) :: diinit = 'DEPHTIME'
      character(len=256) :: spinit = 'ACSPACE'

      character(len=256) :: debuglevel = 'I'

      namelist /NAMDPARA/ bmin, bmax,                &
                          nsample, ntraj, nsw, nelm, &
                          temp, namdtime, potim,     &
                          lhole, algo, lcpext, lshp, &
                          lspace, nacbasis, nacele,  &    !!
                          rundir, tbinit, diinit, spinit, &
                          debuglevel

      integer :: ierr, i, j
      logical :: lext

      ! inp
      open(file="inp", unit=8, status='unknown', action='read', iostat=ierr)
      if ( ierr /= 0 ) then
        write(*,*) "[E] IOError: I/O error with input file: 'inp'"
        stop
      end if

      read(unit=8, nml=NAMDPARA)
      close(unit=8)

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
      
      ! DEPHTIME
      if (algo == 'DISH') then
        allocate(this%DEPHMATR(bmax - bmin + 1, bmax - bmin + 1))   !! read in the pure dephasing time matrix, the diagonal elements are zero
        inquire(file=diinit, exist=lext)
        if (.NOT. lext) then
          write(*,*) "[E] IOError: File containing initial conditions of DISH does NOT exist!"
          stop
        else
          open(unit=10, file=diinit, action='read')
          read(unit=10, fmt=*) ((this%DEPHMATR(i,j), j=1, bmax - bmin + 1), i=1, bmax - bmin + 1)
          do i = 1, bmax - bmin + 1
            do j = 1, bmax - bmin + 1
              if (i /= j) then
                this%DEPHMATR(j,i) = 1.0_q / this%DEPHMATR(j, i)
              end if
            end do
          end do
          ! do i=1, bmax - bmin + 1
          !   do j=1, bmax - bmin + 1
          !     read(unit=10,fmt=*) this%DEPHMATR(i, j)
          !   end do
          ! end do
          close(10)
        end if
      end if

      ! ACSPACE
      if (lspace) then
        allocate(this%BASIS(nacbasis, nacele))
        inquire(file=spinit, exist=lext)
        if (.NOT. lext) then
          write(*,*) "[E] IOError: File containing Active Space does NOT exist!"
          stop
        else
          open(unit=11, file=spinit, action='read')
          read(unit=11, fmt=*) ((this%BASIS(i,j),j=1,nacele),i=1,nacbasis)
          close(11)
          this%NBASIS = nacbasis
        end if
      else
        allocate(this%BASIS(bmax - bmin + 1, 1))
        this%BASIS = reshape([(i, i=bmin, bmax)], shape=[bmax-bmin+1, 1])
        this%NBASIS = bmax - bmin + 1
      end if

      ! if (realtime<=namdtime) then
      !   realtime=namdtime
      ! else
      !   if (algo /= 'DISH') then 
      !     write(*,*) "We do not recommend replicating NAC in FSSH"
      !     write(*,*) "Used NAMDTIME instead of REALTIME"
      !     realtime=namdtime
      !   end if
      ! end if

      ! assign the parameters
      this%BMIN     = bmin
      this%BMAX     = bmax

      this%NSW      = nsw
      this%NTRAJ    = ntraj
      this%NELM     = nelm
      this%NSAMPLE  = nsample

      this%NAMDTIME = namdtime 
      this%POTIM    = potim
      this%TEMP     = temp

      this%LHOLE    = lhole
      this%LSHP     = lshp
      this%LCPTXT   = lcpext
      this%ALGO     = algo

      this%RUNDIR   = trim(rundir)
      this%TBINIT   = trim(tbinit)
      this%DIINIT   = trim(diinit)

      this%LSPACE   = lspace
      this%NACBASIS = nacbasis
      this%NACELE= nacele

      call this%checkUserInp()
    end subroutine

    subroutine checkUserInp(this)
      implicit none
      class(namdInfo), intent(inout) :: this

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

      write(*,'(A)') "------------------------------------------------------------"
      write(*,'(A30,A3,I8)') 'BMIN',     ' = ', this%BMIN
      write(*,'(A30,A3,I8)') 'BMAX',     ' = ', this%BMAX
      write(*,'(A30,A3,I8)') 'INIBAND',  ' = ', this%INIBAND
      write(*,'(A)') ""
      write(*,'(A30,A3,F8.1)') 'POTIM',  ' = ', this%POTIM
      write(*,'(A30,A3,F8.1)') 'TEMP',   ' = ', this%TEMP
      write(*,'(A30,A3,I8)') 'NAMDTINI', ' = ', this%NAMDTINI
      write(*,'(A30,A3,I8)') 'NAMDTIME', ' = ', this%NAMDTIME
      write(*,'(A)') ""
      write(*,'(A30,A3,I8)') 'NSW',      ' = ', this%NSW
      write(*,'(A30,A3,I8)') 'NTRAJ',    ' = ', this%NTRAJ
      write(*,'(A30,A3,I8)') 'NELM',     ' = ', this%NELM
      write(*,'(A)') ""
      write(*,'(A30,A3,L8)') 'LHOLE',    ' = ', this%LHOLE
      write(*,'(A30,A3,L8)') 'LSHP',     ' = ', this%LSHP
      write(*,'(A30,A3,A8)') 'ALGO',     ' = ', TRIM(ADJUSTL(this%ALGO))
      write(*,'(A30,A3,L8)') 'LCPTXT',   ' = ', this%LCPTXT
      write(*,'(A)') ""
      write(*,'(A30,A3,L8)') 'LSPACE',   ' = ', this%LSPACE
      if (this%LSPACE) then
        write(*,'(A30,A3,I8)')'NACBASIS',' = ', this%NACBASIS
        write(*,'(A30,A3,I8)')'NACELE',  ' = ', this%NACELE
      end if
      write(*,'(A)') "------------------------------------------------------------"
    end subroutine

end module
