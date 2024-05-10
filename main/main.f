! ================================================================================================================================ !
program main
  !! Solve the (ro)vibrational Schrödinger equation (RVSE) for a diatomic molecule using a DVR and optionally an optical potential :
  !! $$H_{j}(R) = -1/(2  μ)  \frac{d^2}{dR^2} + V(R) + \frac{j(j+1)}{2  μ  R^2}$$
  !! where $H_{j}(R)$ is the Hamiltonian and $V(R)$ is the internuclear potential. The potential $V(R)$ can be purely real or
  !! a complex absorbing potential (CAP) used to calculate continuum states, as described in [1].
  !!
  !! Reference
  !!
  !! [1]  Vibok et al. (J. Phys. Chem. 96, 8712 (1992))

  use iso_fortran_env, only: wp => real64, stdin => input_unit, iostat_end, stdout => output_unit

  use bsprvse, only: solve_RVSE
  use bsprvse__constants, only: au2ev, au2amu, one, ryd
  use bsprvse__utilities, only: die, append, i2char, ndigits

  implicit none

  integer, parameter :: iostat_ok = 0

  integer :: iR
  integer :: iv
  integer :: io
  integer :: funit
  real(wp) :: dR
  real(wp) :: R
  real(wp) :: V
  real(wp), allocatable :: R_vals(:)
  real(wp), allocatable :: V_vals(:)
  real(wp), allocatable :: R_wf(:)
    !! The $R$-grid on which wavefunctions will be calculated
  complex(wp), allocatable :: wf(:,:)
  complex(wp), allocatable :: wf_nrg(:)
  character(:), allocatable :: filename

  ! -- namelist variables
  integer :: j
    !! The quantum number $j$ in the Hamiltonian
  logical :: CAP_exists
    !! Use a comlex absorbing potential (CAP) ?
    !! yes -> CAP_exists = .true.
    !! no  -> CAP_exists = .false.
  integer :: CAP_type
    !! The CAP type:
    !!   0: exponential
    !!   1: linear
    !!   2: quadratic
    !!   3: cubic
    !!   4: quartic
  integer :: np
    !! The number of B-spline intervals
  integer :: nwf
    !! The number of wavefunctions to calculate, from ν = 0 to ν = nwf - 1
  integer :: legpoints
    !! The number of Gauss-Legendre quadrature points used to calculate integrals
  integer :: order
    !! The order of the B-splines
  integer :: nR_wf
    !! The number of $R$-values for the wavefunctions
  real(wp) :: reduced_mass
    !! The reduced mass of the molecule, in atomic mass units, i.e., the mass of
    !! carbon-12 is 12 atomic mass units.
  real(wp) :: R_max
    !! The largest value to consider for the internuclear potential
  real(wp) :: CAP_length
    !! The length of the CAP, in atomic units. The CAP starts at R_max - CAP_length and
    !! continues to R_max, where it takes its maximal value. R_max is the
    !! largest value at which the internuclear potential was calculated
  real(wp) :: CAP_strength
    !! The CAP_strength $A$, in atomic units.
  real(wp), allocatable :: B_rot(:)
    !! The rotational constant for each vibrational state in atomic units
  character(:), allocatable :: potential_file
    !! The location for the input internuclear potential in the format
    !!   R1 V1
    !!   R2 V2
    !!   R3 V3
    !!    .  .
    !!    .  .
    !!    .  .
    !!
    !! where R and V are both in atomic units
  character(:), allocatable :: output_directory
    !! The directory in which to write the wavefunctions and energies. This directory
    !! should already exist

  namelist / input_parameters /                 &
                                j,              &
                                np,             &
                                nwf,            &
                                legpoints,      &
                                order,          &
                                reduced_mass,   &
                                R_max,          &
                                nR_wf,          &
                                CAP_exists,     &
                                CAP_length,     &
                                CAP_strength,   &
                                CAP_type,       &
                                potential_file, &
                                output_directory

  ! allocate(character(1000) :: potential_file)
  ! allocate(character(1000) :: output_directory)

  call read_input

  ! -- convert reduced mass to atomic units
  reduced_mass = reduced_mass / au2amu

  ! -- read the inout R-grid and potential
  open(newunit = funit, file = potential_file)
  do

    read(funit, *, iostat = io) R, V

    select case(io)
      case(iostat_end) ; exit
      case(iostat_ok)  ; continue
      case default     ; call die("Unknown error reading the file " // potential_file)
    end select

    ! -- consider only R-values smaller than R_max
    if(R .gt. R_max) exit

    call append(R_vals, R)
    call append(V_vals, V)

  enddo

  allocate(wf_nrg(nwf))
  allocate(wf(nR_wf, nwf))
  allocate(B_rot(nwf))

  ! -- the grid points on which wavefunctions should be computed
  dR = (R_vals(size(R_vals, 1)) - R_vals(1)) / (nR_wf - one)
  R_wf = [( R_vals(1) + (iR - one) * dR, iR = 1, nR_wf )]

  ! -- call the procedure to solve the 1D RVSE with or without a CAP
  call solve_RVSE(R_vals, V_vals, j, reduced_mass, nwf, nR_wf, R_wf, wf, wf_nrg, B_rot, np, legpoints, order, &
    CAP_exists, CAP_length, CAP_type, CAP_strength)

  ! -- example of calling the procedure without supplying the CAP variables for a bound-state calculation
  ! call solve_RVSE(R_vals, V_vals, j, reduced_mass, R_max, nwf, nR_wf, wf, wf_nrg, np, legpoints, order)

  ! -- write the output wavefunctions
  do iv = 0, nwf - 1
    filename =  output_directory // "/wf_v" // i2char(iv) // ".dat"
    open(newunit = funit, file = filename)
    do iR = 1, nR_wf

      write(funit, '(3e30.20)') R_wf(iR), wf(iR, iv + 1)

    enddo
    close(funit)
  enddo

  ! -- write the output energies and rotational constants
  filename =  output_directory // "/energies.dat"
  open(newunit = funit, file = filename)
  do iv = 0, nwf - 1
    write(funit, '(A, I' // i2char(ndigits(nwf - 1)) // ', A, 2(e30.20,X), "i")') "Energy for ν = ", iv, " (eV) : " &
      , wf_nrg(iv + 1) * au2ev
  enddo
  close(funit)
  filename =  output_directory // "/rotational.constants.dat"
  open(newunit = funit, file = filename)
  do iv = 0, nwf - 1
    write(funit, '(A, I' // i2char(ndigits(nwf - 1)) // ', A, e30.20)') "Rotational constant B for ν = ", iv, " (invcm) : "  &
      , B_rot(iv + 1) * ryd
  enddo
  close(funit)

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  subroutine read_input

    implicit none

    allocate(character(1000) :: potential_file)
    allocate(character(1000) :: output_directory)

    ! -- read the namelist variables
    read(stdin, input_parameters, iostat = io)

    potential_file = trim(potential_file)
    output_directory= trim(output_directory)

    write(stdout,'(A)') "The input parameters for this run are as follows :"
    write(stdout,*)
    write(stdout, input_parameters)
    write(stdout,*)

  end subroutine read_input


! ================================================================================================================================ !
end program main
! ================================================================================================================================ !
