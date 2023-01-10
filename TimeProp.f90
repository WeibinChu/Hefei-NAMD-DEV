module TimeProp
  use prec
  use fileio
  use hamil

  implicit none

contains 
  
  !Use Trotter formula to integrate the Time-dependtent Schrodinger equation
  !The new scheme provides a robust and efficient integration. The electronic
  !time step (POTIM/NELM) may be increased up to the value of the nuclear time step (Not in all)
  !Always check convergence before using a small NELM 
  !This scheme is proposed by Akimov, A. V., & Prezhdo, O. V. J. Chem. Theory Comput. 2014, 10, 2, 789–804
  
  
  !This scheme has been revised by Dr.Li yunhai (liyunhai1016@hotmail.com)
  !To use this scheme, the offdiagonal elements of Hamiltonian should be real numbers (without
  !the the imaginary unit) 
  subroutine Trotter(ks, edt)
    implicit none

    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in)  :: edt
    integer :: jj, kk
    complex(kind=q) :: phi, cos_phi, sin_phi, cjj, ckk


    ! propagate the psi_c according to Liouville-Trotter algorithm
    !   exp[i*(L_{ij}+L_i)*dt/2]
    ! = exp(i*L_{ij}*dt/2) * exp(i*L_i*dt/2) * exp(i*L_i*dt/2) * exp(i*L_{ij}*dt/2)
    ! = exp(i*L_{ij}*dt/2) * exp(i*L_i*dt) * exp(i*L_{ij}*dt/2)

    ! First L_{ij} part
    
    !Changed the matrix structure for NAC.
    !Now Ham_c(i,j) is stored as  ks%ham_c(j,i)
    !Weibin

    do jj = 1, inp%NBASIS
      do kk = jj+1, inp%NBASIS
        phi = 0.5_q * edt * con%miI * ks%ham_c(kk, jj) / con%hbar
        cos_phi = cos(phi)
        sin_phi = sin(phi)
        cjj = ks%psi_c(jj)
        ckk = ks%psi_c(kk)
        ks%psi_c(jj) =  cos_phi * cjj + sin_phi * ckk
        ks%psi_c(kk) = -sin_phi * cjj + cos_phi * ckk
      end do
    end do
    
    ! Li part

    !Changed the matrix structure for NAC.
    !Now Ham_c(i,j) is stored as  ks%ham_c(j,i) ?
    !Weibin


    do jj = 1, inp%NBASIS
      phi = edt * con%miI * ks%ham_c(jj, jj) / con%hbar
      ks%psi_c(jj) = ks%psi_c(jj) * exp(phi)
    end do

    ! Second L_{ij} part
    ! this should be simplified by the symmetry of Dij, I will change it later.
    do jj = inp%NBASIS, 1, -1
      do kk = inp%NBASIS, jj+1, -1
        phi = 0.5_q * edt * con%miI * ks%ham_c(kk, jj) / con%hbar
        cos_phi = cos(phi)
        sin_phi = sin(phi)
        cjj = ks%psi_c(jj)
        ckk = ks%psi_c(kk)
        ks%psi_c(jj) =  cos_phi * cjj + sin_phi * ckk
        ks%psi_c(kk) = -sin_phi * cjj + cos_phi * ckk

      end do
    end do
  end subroutine
  
  subroutine PropagationEle(ks, tion)
    implicit none
    type(TDKS), intent(inout)  :: ks
    integer, intent(in) :: tion

    integer :: tele
    real(kind=q) :: edt, norm

    edt = inp%POTIM / inp%NELM
    
    do tele = 1, inp%NELM
      ! construct hamiltonian matrix
      call make_hamil_rtime(tion, tele, ks)

      call Trotter(ks, edt)
    end do
    norm = REAL(SUM(CONJG(ks%psi_c) * ks%psi_c), kind=q) 
    if ( norm <= 0.99_q .OR. norm >= 1.01_q)  then
        write(*,*) "[E] Error in Electronic Propagation"
        stop
    end if

    if (inp%DEBUGLEVEL <= 0) write(*,'(A, f12.8)') 'NORM', norm

    !ks%psi_a(:,tion) = ks%psi_c
    !ks%pop_a(:,tion) = ks%pop_a(:,tion)/ks%norm(tion) 


    !write(89,*) tion, conjg(ks%psi_c)*ks%psi_c
    !write(89,*) 'sum', sum(conjg(ks%psi_c)*ks%psi_c)
    ! ! end of the INNER loop
    ! ! call cpu_time(fin)
    ! ! write(*,*) "T_ion ", tion, fin - start
  end subroutine

end module
