module couplings
  use prec
  use constants
  use fileio
  implicit none

  private :: readNaEig, initSpace

  type overlap
    logical :: LLONG
    integer :: LBD
    integer :: UBD
    integer :: NBANDS !! NBASIS
    integer :: TSTEPS !! NSW
    real(kind=q) :: dt
    real(kind=q), allocatable, dimension(:,:,:) :: Dij !! anti-hermite, real matrix
    ! complex(kind=q), allocatable, dimension(:,:,:) :: DijR
    ! complex(kind=q), allocatable, dimension(:,:,:) :: DijI
    real(kind=q), allocatable, dimension(:,:) :: Eig
  end type

contains

  subroutine TDCoupIJ(olap, olap_sp)
    implicit none
    type(overlap), intent(out) :: olap
    type(overlap), intent(out) :: olap_sp !! single particle overlap

    ! AIMD, usually NSW < MAXSIZE, here store all coup coeff.
    ! DPMD, usually NSW > MAXSIZE, here store M+2 coeff.
    olap%LLONG    = (inp%NSW > MAXSIZE)
    olap_sp%LLONG = (inp%NSW > MAXSIZE)
    if (inp%ALGO == 'FSSH' .AND. olap%LLONG) then
      write(*,*) "[E] The NSW is too long for a FSSH calculation."
      stop
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialization
    olap%NBANDS = inp%NBASIS
    olap%TSTEPS = inp%NSW
    olap%dt = inp%POTIM

    olap_sp%NBANDS = inp%BMAX - inp%BMIN + 1
    olap_sp%TSTEPS = inp%NSW
    olap_sp%dt = inp%POTIM
    !Trotter factorization integrator is not compatible with complex NAC
    if (.NOT. olap%LLONG) then
      allocate(olap%Dij(olap%NBANDS, olap%NBANDS, olap%TSTEPS))
      allocate(olap%Eig(olap%NBANDS, olap%TSTEPS))
      allocate(olap_sp%Dij(olap_sp%NBANDS, olap_sp%NBANDS, olap_sp%TSTEPS))
      allocate(olap_sp%Eig(olap_sp%NBANDS, olap_sp%TSTEPS))
      ! read coupling from EIGTXT and NATXT
      if (inp%LCPTXT) then
        call readNaEig(olap_sp)
        if (inp%LSPACE) then
          call initSpace(olap, olap_sp)
        else
          olap = olap_sp
        end if
      else
        write(*,*) "[E] This version does not support coupling from COUPCAR, please use NATXT"
        stop
      end if
    else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(olap%Dij(olap%NBANDS, olap%NBANDS, MAXSIZE+2))
      allocate(olap%Eig(olap%NBANDS, MAXSIZE+2))
      allocate(olap_sp%Dij(olap_sp%NBANDS, olap_sp%NBANDS, MAXSIZE+2))
      allocate(olap_sp%Eig(olap_sp%NBANDS, MAXSIZE+2))
      ! read coupling from EIGTXT and NATXT
      if (inp%LCPTXT) then
        call readLongNaEig(olap, olap_sp)
      else
        write(*,*) "[E] This version does not support coupling from COUPCAR, please use NATXT"
        stop
      end if
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! deallocate(olap_sp%Dij, olap_sp%Eig)
  end subroutine

  subroutine readNaEig(olap)
    implicit none

    type(overlap), intent(inout) :: olap
    integer :: i, j, k, N, ierr

    open(unit=22, file='EIGTXT', status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: EIGTXT does NOT exist!"
      stop
    end if
    open(unit=23, file='NATXT', status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: NATXT does NOT exist!"
      stop
    end if

    N = inp%NSW 
    do j=1, N 
      read(unit=22, fmt=*) (olap%Eig(i,j), i=1, olap%NBANDS)
    end do
    do k=1, N
      read(unit=23, fmt=*) ((olap%Dij(j, i, k), j=1, olap%NBANDS), &
                                                   i=1, olap%NBANDS)
    end do

    if (inp%LHOLE) then
      call inp%setHLORB(1, olap%NBANDS)
    else
      call inp%setHLORB(olap%NBANDS, 1)
    end if

    olap%Dij = olap%Dij / (2*inp%POTIM)

    close(unit=22)
    close(unit=23)
  end subroutine

  subroutine readLongNaEig(olap, olap_sp)
    implicit none

    type(overlap), intent(inout) :: olap, olap_sp
    integer :: i, ii, j, k, N, ierr

    !!!!!!!!!!!!!!!!!!!
    ! SAVE:    1 ~  M+2
    !        M+3 ~ 2M+4
    !       2M+5 ~ 3M+6

    open(unit=22, file='EIGTXT', status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: EIGTXT does NOT exist!"
      stop
    end if
    open(unit=23, file='NATXT', status='unknown', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: NATXT does NOT exist!"
      stop
    end if
    open(unit=42, file='EIGTMP', form='unformatted', status='unknown', &
         access='stream', action='write', iostat=ierr)
    open(unit=43, file='NATMP', form='unformatted', status='unknown',  &
         access='stream', action='write', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] IOError: Can not write to EIGTMP file."
      stop
    end if

    ! init
    N = olap_sp%NBANDS

    if (inp%LHOLE) then
      call inp%setHLORB(1, N)
    else
      call inp%setHLORB(N, 1)
    end if

    do ii=1, inp%NSW, MAXSIZE+2

      read(unit=22, fmt=*, iostat=ierr) ((olap_sp%Eig(i,k), i=1, N), k=1, MAXSIZE+2)
      read(unit=23, fmt=*, iostat=ierr) (((olap_sp%Dij(j, i, k), j=1, N), i=1, N), k=1, MAXSIZE+2)

      olap_sp%Dij = olap_sp%Dij / (2*inp%POTIM)

      !! init_space
      if (inp%LSPACE) then
        call initSpace(olap, olap_sp)
      else
        olap = olap_sp
      end if

      ! write
      write(unit=42) olap%Eig
      write(unit=43) olap%Dij
    end do

    close(22)
    close(23)
    close(42)
    close(43)
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! H_ii = SUM_{ik \in ES} Eig_ik - SUM_{ik \in GS} Eig_ik
  ! D_ij = SUM_k d_{ik.jk} PROD_{k'!=k} delta_{ik'.jk'}
  ! i,j : multi-electron state notation
  ! k,k': electron index 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initSpace(olap, olap_sp)
    implicit none

    type(overlap), intent(inout) :: olap
    type(overlap), intent(in)    :: olap_sp

    integer :: i, j, k, l, N, band, ierr
    integer :: horb, lorb
    real(q) :: heig, leig
    integer, dimension(inp%NACELE) :: bi, bj
    logical :: beq = .false.

    open(unit=35, file='ACEIGTXT', status='unknown', action='write', iostat=ierr)
    open(unit=36, file='ACNATXT', status='unknown', action='write', iostat=ierr)

    olap%Eig = 0.0_q
    olap%Dij = 0.0_q

    ! Excite State Energy
    do i=1,inp%NBASIS
      do j=1,inp%NACELE
        band = inp%BASIS(j,i) - inp%BMIN + 1
        olap%Eig(i,:) = olap%Eig(i,:) + olap_sp%Eig(band,:)
      end do
    end do

    heig = -huge(0.0_q)
    leig = huge(0.0_q)
    do i=1,inp%NBASIS
      if (olap%Eig(i,1) > heig) then
        horb = i
        heig = olap%Eig(i,1)
      end if
      if (olap%Eig(i,1) < leig) then
        lorb = i
        leig = olap%Eig(i,1)
      end if
    end do
    if (inp%LHOLE) then
      call inp%setHLORB(lorb, horb)
    else
      call inp%setHLORB(horb, lorb)
    end if

    ! Minus Gound State Energy
    do N=1,inp%NSW
      olap%Eig(:,N) = olap%Eig(:,N) - olap%Eig(inp%LORB,N)
    end do

    ! NAC Matrix Elements, anti-hermite, assuming real Dij = -Dji
    !         N                  N
    ! H_ij = SUM  h_{i_k' j_k'} PROD delta_{i_k'' j_k''}    
    !        k'=1               k''=1,
    !                           k''!=k'
    ! N = NACELE, i,j: states, k',k'': electrons
    ! Dij is stored as Dij(j,i,N)
    do i=1,inp%NBASIS
      do j=i+1,inp%NBASIS
        bi = inp%BASIS(:,i) - inp%BMIN + 1
        bj = inp%BASIS(:,j) - inp%BMIN + 1
        do k=1,inp%NACELE   !! k'
          beq = .true.
          do l=1,inp%NACELE !! k'' != k'
            if (l == k) cycle
            if (bi(l) /= bj(l)) then
              beq = .false.
              exit
            end if
          end do
          if (beq) then
            olap%Dij(j,i,:) = olap%Dij(j,i,:) + olap_sp%Dij(bj(k),bi(k),:)
          end if
        end do 
        olap%Dij(i,j,:) = -olap%Dij(j,i,:)
      end do
    end do

    if (inp%DEBUGLEVEL > 0) return

    do N=1,inp%NSW
      write(unit=35,fmt=*) (olap%Eig(i,N), i=1, olap%NBANDS)
      write(unit=36,fmt=*) ((olap%Dij(j,i,N), j=1, olap%NBANDS), &
                                              i=1, olap%NBANDS)
    end do
  end subroutine

  subroutine setOlapBoundry(olap, start)
    implicit none
    type(overlap), intent(inout) :: olap
    integer, intent(in) :: start
    integer :: ierr

    if (.NOT. olap%LLONG) return

    olap%LBD = start - 1
    olap%UBD = olap%LBD + MAXSIZE + 1

    if (allocated(olap%Eig)) then
      deallocate(olap%Eig)
      deallocate(olap%Dij)
    end if
    allocate(olap%Eig(olap%NBANDS, olap%LBD : olap%UBD))
    allocate(olap%Dij(olap%NBANDS, olap%NBANDS, olap%LBD : olap%UBD))

    ! read 
    open(unit=42, file='EIGTMP', form='unformatted', status='unknown', &
         access='stream', action='read', iostat=ierr)
    open(unit=43, file='NATMP', form='unformatted', status='unknown',  &
         access='stream', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) "[E] Read EIG and NA from TMP file fails."
      stop
    end if

    ! start from 1?
    read(unit=42, iostat=ierr, pos=1+kind(0.0_q)*olap%NBANDS*(olap%LBD-1)) olap%Eig
    read(unit=43, iostat=ierr, pos=1+kind(0.0_q)*olap%NBANDS*olap%NBANDS*(olap%LBD-1)) olap%Dij

    close(42)
    close(43)
  end subroutine

end module
