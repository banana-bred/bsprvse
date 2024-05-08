! ================================================================================================================================ !
module bsprvse
  !! Solve the (ro)vibrational Schrödinger equation for a diatomic molecule using a DVR and optionally an optical potential :
  !! $$H_{j}(R) = -1/(2  μ)  \frac{d^2}{dR^2} + V(R) + \frac{j(j+1)/}{2  μ  R^2}$$
  !! where $H_{j}(R)$ is the Hamiltonian and $V(R)$ is the internuclear potential. The potential $V(R)$ can be purely real or
  !! a complex absorbing potential (CAP) used to calculate continuum states, as described in [1].
  !!
  !! Original code written by Viatcheslav Kokoouline
  !!
  !! Refactored as an fpm package by Josh Forer
  !!
  !! Reference
  !!
  !! [1]  Vibok et al. (J. Phys. Chem. 96, 8712 (1992))

  implicit none

  private

  public :: solve_RVSE

  interface solve_RVSE
    module procedure solve_RVSE_with_CAP
    module procedure solve_RVSE_without_CAP
  end interface solve_RVSE

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine solve_RVSE_without_CAP(R_vals, V_vals, j, reduced_mass, nwf, nR_wf, R_wf, wf, wf_nrg, B_rot &
                                          , np, legpoints, order)
    !! Wrapper routine to `solve_RVSE_with_CAP` for solving the RVSE without an aborbing potential.
    !! This just sets the `CAP_exists` variable to `.false.` without requiring the CAP variables as input
    !! from the user.

    use bsprvse__constants, only: wp, zero

    real(wp), intent(in)  :: R_vals(:)
      !! Array containing the values of the internuclear distance in atomic units
    real(wp), intent(in)  :: V_vals(:)
      !! Array containing the values of the intermolecular potential in atomic units
    integer, intent(in) :: j
      !! The value $j$ in the Schödinger equation
    real(wp), intent(in)  :: reduced_mass
      !! The system's reduced mass in atomic units (electron mass = 1; not atomic mass units)
    integer, intent(in) :: nwf
      !! The number of wavefunctions to calculate
    integer, intent(in) :: nR_wf
      !! The number of $R$-grid points on which to evaluate the wavefunctions
    real(wp), intent(in) :: R_wf(:)
      !! The $R$-grid on which wavefunctions will be calculated
    complex(wp), intent(out) :: wf(:,:)
      !! Array containing the values of the wavefunctions indexed as (iR, iv), where iR runs over the
      !! internuclear distances and iv runs over the vibrational quantum number ν
    complex(wp), intent(out) :: wf_nrg(:)
      !! The energies of the wavefunctions (in atomic units)
    real(wp), intent(out) :: B_rot(:)
      !! The rotational constant for each vibrational state in atomic units
    integer, intent(in) :: np
      !! The number of B-spline intervals
    integer, intent(in) :: legpoints
      !! The number of Gauss-Legendre quadrature points used to calculate integrals
    integer, intent(in) :: order
      !! The order of the B-splines

    logical, parameter :: CAP_exists    = .false.
    integer, parameter :: CAP_type      = 0
    real(wp), parameter :: CAP_length   = zero
    real(wp), parameter :: CAP_strength = zero

    call solve_RVSE_with_CAP(R_vals, V_vals, j, reduced_mass, nwf, nR_wf, R_wf, wf, wf_nrg, B_rot, np, legpoints, order &
                            , CAP_exists, CAP_length, CAP_type, CAP_strength)

  end subroutine solve_RVSE_without_CAP

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine solve_RVSE_with_CAP(R_vals, V_vals, j, reduced_mass, nwf, nR_wf, R_wf, wf, wf_nrg, B_rot  &
                                       , np, legpoints, order, CAP_exists, CAP_length, CAP_type, CAP_strength )
    !! Solve the (ro)vibrational time-independent Schrödinger equation with an added
    !! imaginary potential on the tail end.

    use bsprvse__routines,  only: bound_states_and_resonances
    use bsprvse__utilities, only: i2char, die
    use bsprvse__constants, only: one, two, wp

    real(wp), intent(in)  :: R_vals(:)
      !! Array containing the values of the internuclear distance in atomic units
    real(wp), intent(in)  :: V_vals(:)
      !! Array containing the values of the intermolecular potential in atomic units
    integer, intent(in) :: j
      !! The value $j$ in the Schödinger equation
    real(wp), intent(in)  :: reduced_mass
      !! The system's reduced mass in atomic units (electron mass = 1; not atomic mass units)
    integer, intent(in) :: nwf
      !! The number of wavefunctions to calculate
    integer, intent(in) :: nR_wf
      !! The number of $R$-grid points on which to evaluate the wavefunctions
    real(wp), intent(in) :: R_wf(:)
      !! The $R$-grid on which wavefunctions will be calculated
    complex(wp), intent(out) :: wf(:,:)
      !! Array containing the values of the wavefunctions indexed as (iR, iv), where iR runs over the
      !! internuclear distances and iv runs over the vibrational quantum number ν
    complex(wp), intent(out) :: wf_nrg(:)
      !! The energies of the wavefunctions (in atomic units)
    real(wp), intent(out) :: B_rot(:)
      !! The rotational constant for each vibrational state in atomic units
    integer, intent(in) :: np
      !! The number of B-spline intervals
    integer, intent(in) :: legpoints
      !! The number of Gauss-Legendre quadrature points used to calculate integrals
    integer, intent(in) :: order
      !! The order of the B-splines
    logical, intent(in) :: CAP_exists
      !! Add an imaginary potential to the real internuclear potential ?
      !!   .true.  -> yes
      !!   .false. -> no
    real(wp), intent(in) :: CAP_length
      !! The length of the CAP, starting from the largest R-value and moving inward. The CAP takes its maximal value at the final
      !! R-value. This variable does nothing if CAP_exists is .false.
    real(wp), intent(in) :: CAP_strength
      !! The CAP strength (parameters A1, A2, A3, A4, A5 in [1])
      !! This variable does nothing if CAP_exists is .false.
    integer, intent(in) :: CAP_type
      !! The type of the CAP. See `optpot()` for the implementation.
      !! Allowable CAP types :
      !!   0 : exponential
      !!   1 : linear
      !!   2 : quadratic
      !!   3 : cubic
      !!   4 : quartic

    integer :: i
    integer :: io
    integer :: nR
      !! The number of input $R$-values
    integer :: funit

    real(wp) :: massa
    real(wp) :: dR
    real(wp), allocatable :: V_vals2(:)

    complex(wp), allocatable :: E_for_K(:)

    ! -- check dimensions
    nR = size(R_vals, 1)
    if(size(wf, 1) .ne. nR_wf) call die("nR_wf " // i2char(nR_wf) //" and size(wf, 1) " // i2char(size(wf, 1)) // " don't match")
    if(size(wf, 2) .ne. nwf) call die("nwf " // i2char(nwf) // " and size(wf, 2) " // i2char(size(wf, 2)) // " don't match")
    if(size(wf_nrg) .ne. nwf) call die("nwf " // i2char(nwf) // " and size(wf_nrg) " // i2char(size(wf_nrg)) // " don't match")
    if(nR .ne. size(V_vals, 1)) call die("size(R_vals) " // i2char(size(R_vals, 1)) // " and size(V_vals) " &
      // i2char(size(V_vals, 1)) // " mismatch")

    select case(CAP_type)
      case(0:4)    ; continue
      case default ; call die("The CAP type " // i2char(CAP_type) // " is not supported. Please choose between 0 – 5")
    end select

    massa = reduced_mass

    ! -- includes the centrifugral potential, if needed
    V_vals2 = V_vals + (j * (j + one)) / (two * massa * R_vals ** 2)

    ! -- calculating bound states and resonances
    call bound_states_and_resonances(R_vals, V_vals2, nWF, nR_wf, R_wf, wf, wf_nrg, np, legpoints, order, massa, B_rot &
                                    , CAP_exists, CAP_length, CAP_type, CAP_strength )

  end subroutine solve_RVSE_with_CAP

! ================================================================================================================================ !
end module bsprvse
! ================================================================================================================================ !
