module parallel
  implicit none
contains

  subroutine divideTasks(iprog, ngroup, ntask, lower, upper, nout)
  !! distributing tasks between processes.
    integer, intent(in)    :: iprog, ntask
    integer, intent(inout) :: ngroup
    integer, intent(out)   :: lower, upper, nout

    integer :: quotient, remainder

    ngroup = MIN(ngroup, ntask)
    quotient  = ntask / ngroup
    remainder = MOD(ntask, ngroup)

    !! q+1 | q+1 | q+1 | q | q ...
    !! For example:
    !! ntask=8, ngroup=5, q=1, r=3
    !!  i lower upper nout
    !!  0    1    2    2
    !!  1    3    4    2
    !!  2    5    6    2
    !!  3    7    7    1
    !!  4    8    8    1
    if (iprog + 1 > remainder) then
      nout = quotient
      lower = remainder + iprog*quotient + 1
      upper = lower + nout - 1
    else
      nout = quotient + 1
      lower = iprog*(quotient+1) + 1
      upper = lower + nout - 1
    end if
  end subroutine

end module parallel