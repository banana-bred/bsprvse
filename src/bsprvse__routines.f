! ================================================================================================================================ !
module bsprvse__routines

  use bsprvse__constants, only: wp

  implicit none

  save

  private

  public :: bound_states_and_resonances

  interface
    !! Explicit interface for the `ZGGEV()` routine from `LAPACK`
    subroutine zggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
      use bsprvse__constants, only: wp
      character(1) :: jobvl
      character(1) :: jobvr
      integer :: n
      integer :: lda
      integer :: ldb
      integer :: ldvl
      integer :: ldvr
      integer :: lwork
      integer :: info
      real(wp) :: rwork(*)
      complex(wp) :: a(lda, *)
      complex(wp) :: b(ldb, *)
      complex(wp) :: alpha(*)
      complex(wp) :: beta(*)
      complex(wp) :: vl(ldvl, *)
      complex(wp) :: vr(ldvl, *)
      complex(wp) :: work(*)
    end subroutine zggev
  end interface


! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine bound_states_and_resonances(R_in, V_in, nWF, nR_wf, R_wf, wf_g, E_for_K     &
                                        , np_in, legpoints_in, order_in, mass_in, B_rot  &
                                        , CAP_exists, CAP_length, CAP_type, CAP_strength &
                                        , left_bc_zero_in, right_bc_zero_in)
    !! The driving routine to calculate wavefunctions given an internuclear potential as a function
    !! of the internuclear distance $R$. A complex absorbing potential (CAP) can be added to calculate
    !! bound continuum states

    use bsprvse__globals,          only: R_mod, V_mod, x_grid, x_weights, x_sectors, B_overlap, H_matrix, basis, bound_energies &
                                       , basis_dimension, legpoints, mass, np, nR_mod, order, V_min, x_begin, x_end, ntotal &
                                       , left_bc_zero, right_bc_zero
    use bsprvse__constants,        only: zero, ci
    use bsprvse__utilities,        only: die
    use bsprvse__quadrature_weights, only: available_points

    real(wp), intent(in) :: R_in(:)
      !! the gid of the potential curve
    real(wp), intent(in) :: V_in(:)
      !! the potential ccurve on the grid
    integer, intent(in) :: nWF
      !! number of bound and resonance energies and wfs to compute
    integer, intent(in) :: nR_wf
      !! the number of grid points for wfs
    real(wp), intent(in) :: R_wf(nR_wf)
      !! the grid points on which wfs should be computed
    complex(wp), intent(out) :: wf_g(nR_wf, nWF)
      !! wave functions
    complex(wp), intent(out) :: E_for_K(nWF)
      !! complex energies E - i Γ / 2
    integer, intent(in) :: np_in
      !! The number of B-spline intervals
    integer, intent(in) :: legpoints_in
      !! The number of Gauss-Legendre quadrature points used to calculate integrals. The available number of points is given
      !! in src/bsprvse__quadrature_weights.f : [ 6, 8, 10, 12, 14, 16, 24, 28, 32, 36, 40, 44, 48, 64, 80 ]
    integer, intent(in) :: order_in
      !! The order of the B-splines
    real(wp), intent(in) :: mass_in
      !! the reduced mass μ of the system in atomic units (not atomic mass unit)
    complex(wp), intent(out) :: B_rot(nWF)
      !! The rotational constant for each vibrational state in atomic units
    logical, intent(in) :: CAP_exists
      !! Add an imaginary potential to the real internuclear potential ?
      !!   .true.  -> yes
      !!   .false. -> no
    real(wp), intent(in) :: CAP_length
      !! The length of the CAP, starting from the larget $R$ value and moving inward. The CAP takes its maximal value at
      !! the end of the $R$ interval. This variable does nothing if CAP_exists is .false.
    integer, intent(in) :: CAP_type
      !! The type of the CAP. See `optpot()` for the implementation.
      !! Allowable CAP types :
      !!   0 : exponential
      !!   1 : linear
      !!   2 : quadratic
      !!   3 : cubic
      !!   4 : quartic
    real(wp), intent(in) :: CAP_strength
      !! The CAP strength (parameters A1, A2, A3, A4, A5 in [1])
      !! This variable does nothing if CAP_exists is .false.
    logical, intent(in) :: left_bc_zero_in
      !! Set the left boundary condition to be zero ?
    logical, intent(in) :: right_bc_zero_in
      !! Set the right boundary condition to be zero ?

    integer :: i, k, l
    integer :: first_sector

    ! -- check array sizes
    nR_mod = size(R_in, 1)
    if(size(V_in, 1) .ne. nR_mod) call die("The supplied R_in and V_in arrays are of different sizes")

    left_bc_zero  = left_bc_zero_in
    right_bc_zero = right_bc_zero_in

    np = np_in
    legpoints = legpoints_in
    order = order_in
    mass = mass_in

    if(.not. any(legpoints .eq. available_points)) call die("Please choose a valid number of points for the Gauss-Legendre &
      &quadrature :  6, 8, 10, 12, 14, 16, 24, 28, 32, 36, 40, 44, 48, 64, 80.")

    R_mod  = R_in
    V_mod  = V_in
    V_min  = minval(V_mod, 1)

    x_begin = R_in(1)
    x_end   = R_in(nR_mod)

    call GetAllGaussFactors

    ! -- assume no boundary conditions (BCs)
    basis_dimension = np + order

    ! ! -- remove the spline for the function and its derivative from each side that has a zero BC
    ! if(left_bc_zero .eqv. .true.)  basis_dimension = basis_dimension - 2
    ! if(right_bc_zero .eqv. .true.) basis_dimension = basis_dimension - 2

    allocate(x_grid(legpoints, np))
    allocate(x_weights(legpoints, np))
    allocate(x_sectors(np + 1))
    allocate(B_overlap(basis_dimension, basis_dimension))
    allocate(H_matrix(basis_dimension, basis_dimension))

    call grid_for_1D_B_splines(np,legpoints,x_begin,x_end,x_grid,x_weights,x_sectors)

    allocate(basis(legpoints, order + 1, basis_dimension))

    call Basis_of_1D_B_splines(np, legpoints, order, x_grid, x_sectors, 0, basis)

    call Matrix_of_Overlap(np, order, legpoints, basis, x_weights, B_overlap)

    call Calc_Matrix_H(CAP_exists, CAP_length, CAP_type, CAP_strength)

    call diag

    E_for_K(1:nWF) = Bound_energies(1:nWF)

    ! -- calculate the wavefunction ψ_ν(R) for each ν
    do i = 1, nWF

      call wf_R_calc(nR_wf, R_wf, i, wf_g(:, i), B_rot(i))

      ! if(i .eq. 1) cycle

      ! ! -- rotational constants should decrease monotonically with increasing vibratinal
      ! !    quantum number ν. However, calculating rotational constants for quasi-continuum
      ! !    states doesn't really make sense. The moment of inertia for such states should
      ! !    be infinite, which implies B = 0, so we set all B_rot = 0 past the first B_rot
      ! !    that increases with increasing ν.
      ! if(B_rot(i) > B_rot(i -1)) B_rot(i) = zero

    enddo

    ! -- deallocate arrays in case this routine is called again during the same program execution
    deallocation: block
      use bsprvse__globals, only: R_mod, V_mod, xwlegall, basis, x_grid, x_weights, x_sectors, B_overlap, psi, H_matrix &
                                , bound_energies
      deallocate(R_mod)
      deallocate(V_mod)
      deallocate(xwLegAll)
      deallocate(basis)
      deallocate(x_grid)
      deallocate(x_weights)
      deallocate(x_sectors)
      deallocate(B_overlap)
      deallocate(Psi)
      deallocate(H_matrix)
      deallocate(Bound_energies)
    end block deallocation

    call interv_reset

  end subroutine bound_states_and_resonances

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine diag
    !! Diagonalize the Hamiltonian, get the wavefunctions

    use iso_fortran_env,    only: stdout => output_unit
    use bsprvse__constants, only: two
    use bsprvse__globals,   only: psi, bound_energies, B_overlap, H_matrix, mass, basis_dimension &
                                , left_bc_zero, right_bc_zero

    integer :: N
    integer :: IL
    integer :: IU
    integer :: M
    integer :: INFO
    integer :: k
    integer :: i1, i2

    real(wp), allocatable :: RWORK(:)

    complex(wp) :: norma
    complex(wp), allocatable :: VL(:,:),work2(:)
    complex(wp), allocatable :: NomE(:)
    complex(wp), allocatable :: DenomE(:)
    complex(wp), allocatable :: S_matrixc(:,:)
    complex(wp), allocatable :: tmpv(:)
    complex(wp), allocatable :: H_matrix2(:,:)

    ! N = np + order
    N = basis_dimension

    ! -- Applying boundary conditions
    i1 = 1
    i2 = N
    if(left_bc_zero) then
      ! -- skip the first 2 splines
      N = N - 2
      i1 = i1 + 2
    endif
    if(right_bc_zero) then
      ! -- skip the last 2 splines
      N = N - 2
      i2 = basis_dimension - 2
    endif

    allocate(psi(N, N))
    allocate(VL(1, N))
    allocate(work2(2 * N))
    allocate(NomE(N))
    allocate(DenomE(N))
    allocate(S_matrixc(N, N))
    allocate(RWORK(8 * N))
    allocate(Bound_energies(N))

    ! -- OP
    ! H_matrix  = H_matrix / (two * mass)
    H_matrix2 = H_matrix(i1 : i2, i1 : i2) / (two * mass)
    S_matrixc = B_overlap(i1 : i2, i1 : i2)
    call zggev('N', 'V', N, H_matrix2, N, S_matrixc, N, NomE, DenomE, VL, 1, psi, N, work2, 2 * N, rwork, info)

    write(stdout, '("Diagonalization done !")')

    ! -- normalization
    do k = 1, N
      tmpv = matmul(B_overlap(i1 : i2,i1 : i2), psi(:,k))
      norma = sum(psi(:,k) * tmpv)
      ! if(k < 10) print *,norma
      psi(:,k) = psi(:,k) / sqrt(norma)
    enddo

    call sort_E2(N, NomE, DenomE, psi, Bound_energies)

  end subroutine diag

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine sort_E2(MatrixDim, NomE, DenomE, psi_R, Bound_energies)
    !! Sort the wavefunctions by the real part of their energies

    use bsprvse__constants, only: zero
    use bsprvse__utilities, only: swapvals, swapcols
    use bsprvse__globals,   only: B_overlap

    integer, intent(in) :: MatrixDim
    complex(wp), intent(in) :: NomE(MatrixDim)
    complex(wp), intent(in) :: DenomE(MatrixDim)
    complex(wp), intent(out) :: Bound_energies(MatrixDim)
    complex(wp), intent(out) :: psi_R(MatrixDim,MatrixDim)
    integer :: i
    integer :: k
    integer :: iloc

    real(wp), parameter :: min_denom = 1e-10_wp

    Bound_energies = zero
    k = 0
    do i = 1, MatrixDim
      if(abs(DenomE(i)) .le. min_denom) cycle
      bound_energies(i) = nomE(i) / denomE(i)
    enddo

    do i = 1, MatrixDim
      do k = i + 1, MatrixDim

        if(Bound_energies(i) % re .le. Bound_energies(k) % re) cycle

        call swapvals(Bound_energies, i, k)

        call swapcols(psi_R, i, k)

      enddo
    enddo

    ! -- setting of bound wavefunctions to be real
    ! do i = 1, MatrixDim
    !   iloc = maxloc(abs(psi_R(:, i)), 1)
    !   psi_R(:, i) = psi_R(:, i) * abs(psi_R(iloc, i)) / psi_R(iloc, i)
    ! enddo

  end subroutine sort_E2

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine wf_R_calc(N_g, R_g, N_wf, wf_g, bv)

    use bsprvse__globals,   only: psi, np, order, x_sectors, basis_dimension, mass &
                                , left_bc_zero, right_bc_zero
    use bsprvse__constants, only: zero, ci

    integer,     intent(in)  :: N_g
    integer,     intent(in)  :: N_wf
    real(wp),    intent(in)  :: R_g(N_g)
    complex(wp), intent(out) :: wf_g(N_g)
    complex(wp),    intent(out) :: bv

    integer :: i, i1, i2, N
    ! real(wp) :: bcoef(basis_dimension)
    real(wp) :: wf_r(N_g)
    real(wp) :: wf_c(N_g)
    real(wp) :: bcoef_r(basis_dimension)
    real(wp) :: bcoef_c(basis_dimension)

    bcoef_r = zero
    bcoef_c = zero
    N = basis_dimension
    i1 = 1
    i2 = N
    if(left_bc_zero .eqv. .true.) then
      i1 = i1 + 2
      N  = N - 2
    endif
    if(right_bc_zero .eqv. .true.) then
      i2 = i2 - 2
      N  = N - 2
    endif
    ! N = basis_dimension - 1
    bcoef_r(i1 : i2) = psi(1 : N, N_wf) % re
    bcoef_c(i1 : i2) = psi(1 : N, N_wf) % im

    do i = 1, N_g
      wf_r(i) = B_Spline_to_R(order, np + 1, R_g(i), x_sectors, bcoef_r)
      wf_c(i) = B_Spline_to_R(order, np + 1, R_g(i), x_sectors, bcoef_c)
      wf_g(i) = cmplx(wf_r(i), wf_c(i), kind = wp)
    enddo

    ! -- Calculate the rotational constant: B = < ψ_ν(R) | 1 / (2μR^2) | ψ_ν(R) >, where the bra is not
    !    conjugated. For a purely real potential, these wavefunctions will be real anyway. With an
    !    added imaginary potential, the normalization requires that the bra not be conjugated.
    !    The rotational constants will be complex. The vibrational energies with a CAP are in general complex,
    !    so complex rotational constants will just add further complex energy corrections if rotation is added
    !    later on with, e.g., the rigid rotor approximation.
    bv = integral_trapezoid(R_g, wf_g * wf_g / ( 2 * mass * R_g * R_g ))
    ! bv = integral_trapezoid(R_g, conjg(wf_g) * wf_g / ( 2 * mass * R_g * R_g ))

  end subroutine wf_R_calc

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine Calc_Matrix_H(Absorbing, AP_length, AP_type, A)

    use bsprvse__constants, only: zero, one, two, ci
    use bsprvse__globals,   only: H_matrix, basis, legpoints, order, np, mass, x_sectors, x_grid, x_weights

    logical,  intent(in) :: Absorbing
    integer,  intent(in) :: AP_type
    real(wp), intent(in) :: AP_length
    real(wp), intent(in) :: A

    integer :: i, i1, i2, l

    real(wp) :: f(legpoints, np)
    real(wp) :: basis1(legpoints, order + 1, np + order)

    ! d^2/dx^2
    call Basis_of_1D_B_splines(np, legpoints, order, x_grid, x_sectors, 1, basis1)

    H_matrix = zero
    f = one

    do i1 = 1, np + order
      do i2 = 1, np + order

        H_matrix(i1, i2) = Overlap_2_bfunction(np, order, legpoints, basis1(:, :, i1), i1, basis1(:, :, i2), i2, x_weights, f)

        ! -- Boundary terms
        if(i1 <= order+1 .and. i2 <= order + 1)  H_matrix(i1, i2) = H_matrix(i1, i2) &
          + basis(1, 1, i1) * x_weights(1, 1) * basis1(1, 1, i2)

        if(i1 >= np .and. i2 >= np) H_matrix(i1,i2) =  H_matrix(i1, i2) &
          - basis(legpoints, np + order + 1 - i1, i1) * x_weights(legpoints, np) * basis1(legpoints, np + order + 1 - i2, i2)

      enddo
    enddo

    ! + V(x)
    do concurrent( i = 1 : np, l = 1 : legpoints )
      f(l, i) = V_potential(x_grid(l, i))
    enddo

    do concurrent ( i1 = 1 : np + order, i2 = 1 : np + order )
      H_matrix(i1, i2) = H_matrix(i1, i2) + two * mass &
        * Overlap_2_bfunction(np, order, legpoints, basis(:, :, i1), i1, basis(:, :, i2), i2, x_weights, f)
    enddo

    if(Absorbing .eqv. .false.) return

    do concurrent( i = 1 : np, l = 1 : legpoints )
     f(l, i) = OptPot(x_grid(l, i), x_grid(legpoints, np), AP_length, AP_type, A)
    enddo

    do concurrent ( i1 = 1 : np + order, i2 = 1 : np + order )
      H_matrix(i1, i2) = H_matrix(i1, i2) + two * mass * ci &
        * Overlap_2_bfunction(np, order, legpoints, basis(:, :, i1), i1, basis(:, :, i2), i2, x_weights, f)
    enddo

  end subroutine Calc_Matrix_H

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine Matrix_of_Overlap(np_local, order_local, legpoints_local, basis_local, x_weights_local, B_overlap_local)

    use bsprvse__constants, only: zero, one

    integer, intent(in) :: order_local
    integer, intent(in) :: np_local
    integer, intent(in) :: legpoints_local
    real(wp), intent(in) :: basis_local(legpoints_local, order_local + 1, np_local + order_local)
    real(wp), intent(out) :: B_overlap_local(np_local + order_local, np_local + order_local)
    real(wp), intent(in) :: x_weights_local(legpoints_local,np_local)
    integer :: i1
    integer :: i2
    real(wp) :: f(legpoints_local,np_local)
    B_overlap_local = zero
    f = one
    do concurrent( i1 = 1 : np_local + order_local, i2 = 1 : np_local + order_local )
        B_overlap_local(i1, i2) = Overlap_2_bfunction(                                                                     &
          np_local, order_local, legpoints_local, basis_local(:, :, i1), i1, basis_local(:, :, i2), i2, x_weights_local, f &
        )
    enddo
  end subroutine Matrix_of_Overlap

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function OptPot(r, r_end, L, potential_type, A) result(res)
    !! Optical potential. See works by
    !! Vibok et al. (J. Phys. Chem. 96, 8712 (1992)) for details and parameters
    use bsprvse__constants, only: zero, two, three, five
    integer, intent(in) :: potential_type
    real(wp), intent(in) :: r
    real(wp), intent(in) :: r_end
    real(wp), intent(in) :: L
    real(wp), intent(in) :: A
    real(wp) :: r0
    real(wp) :: res
    real(wp) :: p
    real(wp) :: x
    p   = 13.22_wp
    res = zero
    r0 = r_end - L
    if(r < r0) return
    if(r > r_end) return
    x = (r - r0) / L
    select case(potential_type)
    case(0)
      res = -A * p * exp(-two / x)
    case(1)
      res = -A * x
    case(2)
      res = -A * three / two * x ** 2
    case(3)
      res = -A * two * x ** 3
    case(4)
      res = -A * five / two * x ** 4
    end select
  end function OptPot

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function Overlap_2_bfunction(np_local, order_local, legpoints_local, u1, i1, u2, i2, x_weights_local, f) result(res)
    use bsprvse__constants, only: zero
    integer, intent(in) :: np_local
    integer, intent(in) :: order_local
    integer, intent(in) :: legpoints_local
    real(wp), intent(in) :: u1(legpoints_local, order_local + 1)
    real(wp), intent(in) :: u2(legpoints_local, order_local + 1)
    real(wp), intent(in) :: x_weights_local(legpoints_local, np_local)
    real(wp), intent(in) :: f(legpoints_local, np_local)
    integer, intent(in) :: i1
    integer, intent(in) :: i2
    real(wp) :: res
    integer :: s1_min
    integer :: s1_max
    integer :: s2_min
    integer :: s2_max
    integer :: s
    res = zero
    s1_min = max(1, i1 - order_local)
    s1_max = min(np_local, i1)
    s2_min = max(1, i2 - order_local)
    s2_max = min(np_local, i2)
    do s = max(s1_min, s2_min), min(s1_max, s2_max)
      res = res + sum( u1(:, s - s1_min + 1) * u2(:, s - s2_min + 1) * x_weights_local(:, s) * f(:, s) )
    enddo
  end function Overlap_2_bfunction

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine Basis_of_1D_B_splines(np_local, legpoints_local, order_local, x_grid_local, x_sectors_local, deriv, u)
    use bsprvse__constants, only: zero
    integer, intent(in) :: np_local
    integer, intent(in) :: legpoints_local
    integer, intent(in) :: order_local
    integer, intent(in) :: deriv
    real(wp), intent(in) :: x_grid_local(legpoints_local, np_local)
    real(wp), intent(in) :: x_sectors_local(np_local + 1)
    real(wp), intent(out) :: u(legpoints_local, order_local + 1, np_local + order_local)
    integer :: i
    integer :: k
    integer :: l
    integer :: first_sector
    u = zero
    do i = 1, np_local + order_local
      do k = 1, order_local + 1
        first_sector = max(1, i - order_local) + k - 1
        if(first_sector > np_local) cycle
        do l = 1, legpoints_local
          u(l, k, i) = One_B_Spline(order_local, deriv, np_local + 1, i, x_grid_local(l, first_sector), &
            x_sectors_local(1 : np_local + 1))
        enddo
      enddo
    enddo
  end subroutine Basis_of_1D_B_splines

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine grid_for_1D_B_splines(np_local, legpoints_local, x_begin_local, x_end_local, x_grid_local, x_weights_local, &
                                   x_sectors_local &
    )
    integer, intent(in)   :: np_local
    integer, intent(in)   :: legpoints_local
    real(wp), intent(in)  :: x_begin_local
    real(wp), intent(in)  :: x_end_local
    real(wp), intent(out) :: x_grid_local(legpoints_local, np_local)
    real(wp), intent(out) :: x_weights_local(legpoints_local, np_local)
    real(wp), intent(out) :: x_sectors_local(np_local + 1)
    integer  :: i
    real(wp) :: dx
    real(wp) :: x_begin2
    dx = (x_end_local - x_begin_local) / np_local
    do i = 1, np_local
      x_begin2 = x_begin_local + dx * (i - 1)
      x_sectors_local(i) = x_begin2
      call make_Legendre_grid(legpoints_local, x_begin2, x_begin2 + dx, x_grid_local(:, i), x_weights_local(:, i))
    enddo
    x_sectors_local(np_local + 1) = x_end_local
  end subroutine grid_for_1D_B_splines

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine make_Legendre_grid(nx, x_begin_local, x_end_local, x, x_weights_local)

    use bsprvse__constants, only: one, two
    use bsprvse__globals, only: xwLegAll

    integer, intent(in)   :: nx
    real(wp), intent(in)  :: x_begin_local
    real(wp), intent(in)  :: x_end_local
    real(wp), intent(out) :: x(nx)
    real(wp), intent(out) :: x_weights_local(nx)
    real(wp) :: coef
    x = xwLegAll(1:nx, 1, nx)
    x_weights_local = xwLegAll(1:nx, 2, nx)
    coef = (x_end_local - x_begin_local) / two
    x = coef * (x + one) + x_begin_local
    x_weights_local = coef * x_weights_local
  end subroutine make_Legendre_grid

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine GetAllGaussFactors
    ! use Siegert_Constants

    use bsprvse__constants,          only: zero
    use bsprvse__globals,            only: xwLegAll
    use bsprvse__quadrature_weights, only: available_points, gaussian_nodes_and_weights

    ! integer :: Points
    integer :: i
    integer :: npoints
    integer :: maxpoints
    real(wp), allocatable :: nodes(:)
    real(wp), allocatable :: weights(:)
    ! open(unit = 7, file = 'Legendre.dat')

    maxpoints = maxval(available_points)

    allocate(xwLegAll(maxpoints, 2, maxpoints))
    xwLegAll = zero

    do i = 1, size(available_points)

      npoints = available_points(i)

      if(npoints .eq. maxpoints) exit

      call gaussian_nodes_and_weights(npoints, nodes, weights)
      xwLegAll(1:npoints, 1, npoints) = nodes(1:npoints)
      xwLegAll(1:npoints, 2, npoints) = weights(1:npoints)

    enddo

  end subroutine GetAllGaussFactors

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function V_potential(x) result(res)
    use bsprvse__auxiliary_routines, only: spline, spl
    use bsprvse__globals, only: nR_mod, R_mod, V_mod
    real(wp), intent(in) :: x
    real(wp) :: cm(nR_mod)
    real(wp) :: res
    call spline(nR_mod, R_mod, V_mod, cm)
    res = spl(nR_mod, R_mod, V_mod, cm, x)
   !   V_potential=1/(x*x)-1./x
  end function V_potential

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  function B_Spline_to_R(order_local, np_local, x, x_grid_local, bcoef) result(res)
    !  np is number of points here (not number of intervals)
    integer, intent(in) :: order_local
    integer, intent(in) :: np_local
    real(wp), intent(in) :: x
    real(wp), intent(in) :: x_grid_local(np_local)
    real(wp), intent(in) :: bcoef(np_local + order_local - 1)
    real(wp) :: res
    real(wp) :: t(np_local + 2 * order_local)
    t(1:order_local)               =  x_grid_local(1)
    t(1 + order_local:np_local + order_local)      =  x_grid_local(1:np_local)
    t(1 + order_local + np_local:np_local + 2*order_local) =  x_grid_local(np_local)
    res = bvalue(t, bcoef, np_local + order_local - 1, order_local + 1, x, 0)
  end function B_Spline_to_R

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  function One_B_Spline(order_local, deriv, np_local, n_spl, x, x_grid_local) result(res)
    use bsprvse__constants, only: zero, one
    integer, intent(in) :: order_local
    integer, intent(in) :: deriv
    integer, intent(in) :: np_local
    integer, intent(in) :: n_spl
    real(wp), intent(in) :: x
    real(wp), intent(in) :: x_grid_local(np_local)
    real(wp) :: bcoef(np_local + order_local - 1)
    real(wp) :: t(np_local + 2 * order_local)
    real(wp) :: res
    bcoef = zero
    bcoef(n_spl) = one
    t(1:order_local)                                           = x_grid_local(1)
    t(1 + order_local : np_local + order_local)                = x_grid_local(1:np_local)
    t(1 + order_local + np_local : np_local + 2 * order_local) = x_grid_local(np_local)
    res = bvalue(t, bcoef, n_spl, order_local + 1, x, deriv)
  end function One_B_Spline

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  function bvalue ( t, bcoef, n, k, x, jderiv ) result(res)
    use bsprvse__constants, only: zero
    ! real(wp), intent(in) :: t(*)
    ! real(wp), intent(in) :: bcoef(*)
    real(wp), intent(in) :: t(:)
    real(wp), intent(in) :: bcoef(:)
    integer,  intent(in) :: n
    integer,  intent(in) :: k
    real(wp), intent(in) :: x
    integer,  intent(in) :: jderiv
    real(wp) :: res
    integer, parameter :: kmax = 20
    integer :: i
    integer :: ilo
    integer :: imk
    integer :: j
    integer :: jc
    integer :: jcmin
    integer :: jcmax
    integer :: jj
    integer :: kmj
    integer :: km1
    integer :: mflag
    integer :: nmi
    integer :: jdrvp1
    real(wp) :: aj(kmax)
    real(wp) :: dl(kmax)
    real(wp) :: dr(kmax)
    real(wp) :: fkmj

    res = zero
    if (jderiv .ge. k) return

    call interv ( t, n + k, x, i, mflag )

    if (mflag .ne. 0)Return

    km1 = k - 1
    if (km1 .gt. 0)                   go to 1
    res = bcoef(i)
    return

    1 jcmin = 1
    imk = i - k
    if (imk .ge. 0)                   go to 8
    jcmin = 1 - imk
    do j = 1, i
      dl(j) = x - t(i + 1 - j)
    enddo
    do j = i, km1
      aj(k - j) = zero
      dl(j) = dl(i)
    enddo
                                     go to 10
    8 do j=1,km1
      dl(j) = x - t(i + 1 - j)
    enddo

    10 jcmax = k
    nmi = n - i
    if (nmi .ge. 0)                   go to 18
    jcmax = k + nmi
    do j = 1, jcmax
      dr(j) = t(i + j) - x
    enddo
    do j = jcmax, km1
      aj(j + 1) = zero
      dr(j) = dr(jcmax)
       enddo
                                         go to 20
    18 do j = 1, km1
          dr(j) = t(i + j) - x
       enddo

    20 do jc = jcmin, jcmax
          aj(jc) = bcoef(imk + jc)
       enddo

       ! -- difference the coefficients  jderiv  times
       if (jderiv .eq. 0)                go to 30
       do j = 1, jderiv
          kmj = k - j
          fkmj = real(kmj, kind = wp)
          ilo = kmj
          do jj = 1, kmj
             aj(jj) = ((aj(jj + 1) - aj(jj)) / (dl(ilo) + dr(jj))) * fkmj
             ilo = ilo - 1
          enddo
       enddo

       ! -- compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
       !    given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
    30 if (jderiv .eq. km1)              go to 39
       jdrvp1 = jderiv + 1
       do j = jdrvp1, km1
          kmj = k-j
          ilo = kmj
          do jj = 1, kmj
             aj(jj) = (aj(jj + 1) * dl(ilo) + aj(jj) * dr(jj)) / (dl(ilo) + dr(jj))
             ilo = ilo - 1
          enddo
       enddo
    39 res = aj(1)

  end function bvalue

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine interv (xt, lxt, x, left, mflag)
    real(wp), intent(in) :: xt(lxt)
    integer,  intent(in) :: lxt
    real(wp), intent(in) :: x
    integer,  intent(out) :: left
    integer,  intent(out) :: mflag
    integer :: ihi
    integer :: istep
    integer :: middle
    integer :: ilo = 1
    ihi = ilo + 1
    if (ihi .lt. lxt)                 go to 20
       if (x .ge. xt(lxt))            go to 110
       if (lxt .le. 1)                go to 90
       ilo = lxt - 1
       ihi = lxt

    20 if (x .ge. xt(ihi))            go to 40
    if (x .ge. xt(ilo))               go to 100

  !              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
    istep = 1
    31    ihi = ilo
       ilo = ihi - istep
       if (ilo .le. 1)                go to 35
       if (x .ge. xt(ilo))            go to 50
       istep = istep*2
                                      go to 31
    35 ilo = 1
    if (x .lt. xt(1))                 go to 90
                                      go to 50
  !              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
    40 istep = 1
    41    ilo = ihi
       ihi = ilo + istep
       if (ihi .ge. lxt)              go to 45
       if (x .lt. xt(ihi))            go to 50
       istep = istep*2
                                      go to 41
    45 if (x .ge. xt(lxt))            go to 110
    ihi = lxt

  !           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
    50 middle = (ilo + ihi) / 2
    if (middle .eq. ilo)              go to 100
  !     note. it is assumed that middle = ilo in case ihi = ilo+1 .
    if (x .lt. xt(middle))            go to 53
       ilo = middle
                                      go to 50
    53    ihi = middle
                                      go to 50
    90 mflag = -1
    left = 1
                                      return
    100 mflag = 0
    left = ilo
                                      return
    110 mflag = 1
    if (x .eq. xt(lxt)) mflag = 0
    left = lxt
    111 if (left .eq. 1)              return
    left = left - 1
    if (xt(left) .lt. xt(lxt))        return
                                      go to 111
                                      return
    ! -- reset ilo to 1 in case the parent routines are called multiple times in one program execution.
    entry interv_reset
    ilo = 1
                                      return

  end subroutine interv

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function integral_trapezoid(x, y) result(r)
    use bsprvse__constants, only: two
    real(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(size(x, 1))
    complex(wp) :: r
    integer :: n
    n = size(x, 1)
    r = sum( (y(2 : n) + y(1 : n - 1)) * (x(2 : n) - x(1 : n - 1)) ) / two
  end function integral_trapezoid

! ================================================================================================================================ !
end module bsprvse__routines
! ================================================================================================================================ !
