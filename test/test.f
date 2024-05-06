program test
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

  use iso_fortran_env, only: wp => real64

  implicit none

  integer :: i
  integer :: id
  integer :: j
  integer :: nwf
  integer :: nR_wf
  integer :: nR_in
  integer :: ncolumn_R
  integer :: ncolumn_V
  integer :: funit

  real(wp), parameter :: zero = 0.0_wp
  real(wp), parameter :: one  = 1.0_wp
  real(wp), parameter :: two  = 2.0_wp
  real(wp), parameter :: Ryd  = 219474.6313710_wp
  real(wp) :: R_max
  real(wp) :: massa
  real(wp) :: R_in(10000)
  real(wp) :: V_in(10000)
  real(wp), allocatable :: R_wf(:)

  complex(wp), allocatable :: wf_g(:,:)
  complex(wp), allocatable :: E_for_K(:)
  character(:), allocatable :: input_filename
  character(:), allocatable :: file_with_curve

  R_max = 10.0_wp

  ! input_filename = 'bound_spl.in'
  ! open(newunit = funit, file = input_filename)
  ! call FindStr(input_filename, funit, '#system')
  ! read(funit,*) massa
  ! close(funit)

  ! ! rewind funit

  ! open(newunit = funit, file = input_filename)
  ! call FindStr(input_filename, funit, '#Pot_di')
  ! R_in = zero
  ! V_in = zero
  ! allocate(character(500) :: file_with_curve)
  ! read(funit,*) file_with_curve, ncolumn_R, ncolumn_V
  ! file_with_curve = trim(file_with_curve)
  ! close(funit)

  massa = 13413.54_wp
  file_with_curve = "CFp_pot_30.dat"
  ncolumn_R = 1
  ncolumn_V = 2

  open(newunit = funit, file = file_with_curve)
  id = 0
  do
    id = id + 1
    read(funit, *, end=99) (R_in(id), i=1, ncolumn_R), (V_in(id), i=ncolumn_R+1, ncolumn_V)
    if(R_in(id) > R_max) exit
  enddo
  99 nR_in = id-1
  close(funit)
  j = 0
  ! -- includes the centrifugral potential, if needed
  V_in(1:nR_in) = V_in(1:nR_in) + (j * (j + one)) / (two * massa * R_in(1:nR_in) * R_in(1:nR_in))

  ! -- the grid points on which  wfs should be computed
  nwf = 64; nR_wf = 1000
  allocate(R_wf(nR_wf),wf_g(nR_wf,nwf),E_for_K(nwf))

  ! -- calculating bound states and resonances
  call bound_states_and_resonances(nR_in,R_in(1:nR_in),V_in(1:nR_in),nWF,nR_wf,R_wf,wf_g,E_for_K)
  ! call Siegert_states(nR_in, R_in(1:nR_in), V_in(1:nR_in), nwf, nR_wf, R_wf, wf_g, E_for_K)

  open(newunit = funit, file='energies.dat')
  do i = 1, nwf
    write(funit, '(10e20.12)') E_for_K(i) % re, E_for_K(i) % im
  enddo
  close(funit)

  open(newunit = funit, file='wf.dat')
  do i = 1, nR_wf
    write(funit, '(1000e20.12)') R_wf(i), wf_g(i, 1:nwf)
  enddo
  close(funit)

  deallocate(R_wf,wf_g,E_for_K)

end program test
