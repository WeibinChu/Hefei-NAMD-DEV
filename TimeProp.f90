module TimeProp
  use prec
  use constants
  use utils
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

    do jj = 1, ks%ndim
      do kk = jj+1, ks%ndim
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
    !Now Ham_c(i,j) is stored as  ks%ham_c(j,i) 
    !Weibin


    do jj = 1, ks%ndim
      phi = edt * con%miI * ks%ham_c(jj, jj) / con%hbar
      ks%psi_c(jj) = ks%psi_c(jj) * exp(phi)
    end do

    ! Second L_{ij} part
    do jj = ks%ndim, 1, -1
      do kk = ks%ndim, jj+1, -1
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
  
  subroutine TrotterMat(ks, edt)
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in)  :: edt

    integer :: jj, kk
    complex(kind=q) :: cos_phi, sin_phi
    complex(kind=q), dimension(ks%ndim,ks%ndim) :: identity, ham_l, ham_r
    complex(kind=q), dimension(ks%ndim) :: cjj, ckk

    ! Here input Ham_c(j,i), output Ham_c(i,j)
    identity = zDIAG(ks%ndim, [(con%uno, jj=1,ks%ndim)])
    ! phi
    ks%ham_c(:,:) = 0.5_q * edt * con%miI * ks%ham_c(:,:) / con%hbar

    ham_l = identity
    do jj = 1, ks%ndim
      do kk = jj+1, ks%ndim
        cos_phi = cos(ks%ham_c(kk, jj))
        sin_phi = sin(ks%ham_c(kk, jj))
        cjj(:) = cos_phi * ham_l(:,jj) - sin_phi * ham_l(:,kk)
        ckk(:) = sin_phi * ham_l(:,jj) + cos_phi * ham_l(:,kk)
        ham_l(:,jj) = cjj(:)
        ham_l(:,kk) = ckk(:)
        ! ham_r = identity
        ! ham_r(jj,jj) =  cos_phi
        ! ham_r(kk,kk) =  cos_phi
        ! ham_r(jj,kk) =  sin_phi
        ! ham_r(kk,jj) = -sin_phi
        ! ham_l = matmul(ham_l, ham_r)
      end do
    end do

    ! ham_r = identity
    do jj = 1, ks%ndim
      ham_l(:,jj) = ham_l(:,jj) * exp(2 * ks%ham_c(jj,jj))
      ! ham_r(jj,jj) = exp(2 * ks%ham_c(jj,jj))
    end do
    ! ham_l = matmul(ham_l, ham_r)

    do jj = ks%ndim, 1, -1
      do kk = ks%ndim, jj+1, -1
        cos_phi = cos(ks%ham_c(kk, jj))
        sin_phi = sin(ks%ham_c(kk, jj))
        cjj(:) = cos_phi * ham_l(:,jj) - sin_phi * ham_l(:,kk)
        ckk(:) = sin_phi * ham_l(:,jj) + cos_phi * ham_l(:,kk)
        ham_l(:,jj) = cjj(:)
        ham_l(:,kk) = ckk(:)
        ! ham_r = identity
        ! ham_r(jj,jj) =  cos_phi
        ! ham_r(kk,kk) =  cos_phi
        ! ham_r(jj,kk) =  sin_phi
        ! ham_r(kk,jj) = -sin_phi
        ! ham_l = matmul(ham_l, ham_r)
      end do
    end do

    ks%ham_c = ham_l

  end subroutine

  subroutine Euler(ks, edt)  
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in) :: edt

    ks%psi_c = ks%psi_c &
               & - con%I * edt / con%hbar * matmul(transpose(ks%ham_c), ks%psi_c)

  end subroutine

  !! The first step should use Euler method.
  subroutine EulerMod(ks, edt)
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in) :: edt

    !! check, should be deleted.
    if ((inp%DEBUGLEVEL <= 0 .AND. (.NOT. allocated(ks%psi_p)) .OR. (.NOT. allocated(ks%psi_n)))) then
      write(*,*) '[E] Did not allocate memory for psi_p and psi_n'
      stop
    end if

    ! if (.NOT. init) then
      ks%psi_n = ks%psi_p &
                 & - 2 * con%I * edt / con%hbar * matmul(transpose(ks%ham_c), ks%psi_c)
    ! else
    !   ks%psi_n = ks%psi_c &
    !              & - con%I * edt / con%hbar * matmul(transpose(ks%ham_c), ks%psi_c)
    ! end if
    ks%psi_p = ks%psi_c
    ks%psi_c = ks%psi_n

  end subroutine

  subroutine EulerMat(ks, edt)
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in) :: edt
  
    ks%ham_c = con%miI * edt / con%hbar * ks%ham_c
  end subroutine

  subroutine EulerModMat(ks, edt)
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in) :: edt
  
    ks%ham_c = 2 * con%miI * edt / con%hbar * ks%ham_c
  end subroutine

  subroutine DiagonizeMat(ks, edt)
    implicit none
    type(TDKS), intent(inout) :: ks
    real(kind=q), intent(in) :: edt

    integer :: INFO, LWORK
    real(kind=q),    dimension(ks%ndim)         :: eig
    complex(kind=q), dimension(ks%ndim,ks%ndim) :: exp_eig
    complex(kind=q), dimension(2*ks%ndim-1)     :: WORK
    real(kind=q),    dimension(3*ks%ndim-2)     :: RWORK

    complex(kind=q), dimension(ks%ndim,ks%ndim) :: C1, C2
    ! complex(kind=q), dimension(ks%ndim)         :: Y

    !! https://netlib.org/lapack/explore-html
    !! zHEEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
    LWORK = 2*ks%ndim - 1
    call zHEEV('V', 'U',                  &
              ks%ndim, ks%ham_c, ks%ndim, &
              eig,                        &
              WORK, LWORK, RWORK,         &
              INFO)
    if (INFO == 0) then
      exp_eig = zDIAG(ks%ndim, exp(con%miI*edt/con%hbar*eig))
    else
      write(*,*) '[E] Hamiltonian Diagonization Error. Code: ', INFO
      stop
    end if

    !! zGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    !! C2 = matmul(matmul(conjg(ks%ham_c), exp_eig), transpose(ks%ham_c)) !! H_ij
    C1 = con%cero
    C2 = con%cero
    call zGEMM('C', 'N',                                             &
              ks%ndim, ks%ndim, ks%ndim, con%uno, ks%ham_c, ks%ndim, &
              exp_eig, ks%ndim,                                      &
              con%cero, C1, ks%ndim)
    call zGEMM('N', 'T',                                             &
              ks%ndim, ks%ndim, ks%ndim, con%uno, C1, ks%ndim,       &
              ks%ham_c, ks%ndim,                                     &
              con%cero, C2, ks%ndim)

    ks%ham_c = C2
    !! zGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    ! Y = con%cero
    ! call zGEMV('N',                                   &
    !           ks%ndim, ks%ndim, con%uno, C2, ks%ndim, &
    !           ks%psi_c, 1,                            &
    !           con%cero, Y, 1)
    ! ks%psi_c = Y

  end subroutine

  function HamPsi(ham, psi, jobz)
    implicit none
    complex(kind=q), intent(in), dimension(:,:) :: ham
    complex(kind=q), intent(in), dimension(:)   :: psi
    character(len=1), intent(in) :: jobz
    
    complex(kind=q), dimension(size(psi)) :: HamPsi

    select case (jobz)
    case ('N')
      HamPsi = matmul(ham, psi)
    case ('C')
      HamPsi = matmul(transpose(ham), psi)
    end select
  end function

  function HamPsiBlas(ham, psi, jobz)
    implicit none
    complex(kind=q), intent(in), dimension(:,:) :: ham
    complex(kind=q), intent(in), dimension(:)   :: psi
    character(len=1), intent(in) :: jobz
    
    complex(kind=q), dimension(size(psi)) :: HamPsiBlas, Y
    integer :: N

    !! zGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    N = size(psi)
    Y = con%cero
    call zGEMV(jobz,                 &
              N, N, con%uno, ham, N, &
              psi, 1,                &
              con%cero, Y, 1)
    HamPsiBlas = Y
  end function

  subroutine PropagationEle(ks, tion)
    implicit none
    type(TDKS), intent(inout)  :: ks
    integer, intent(in) :: tion

    integer :: tele
    real(kind=q) :: edt, norm

    edt = inp%POTIM / inp%NELM
    
    do tele = 1, inp%NELM
      ! construct hamiltonian matrix
      call make_hamil(tion, tele, ks)

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
