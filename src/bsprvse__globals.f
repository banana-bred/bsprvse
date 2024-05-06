! ================================================================================================================================ !
module bsprvse__globals

  use bsprvse__constants, only: wp

  implicit none

  private

  save

  integer, public :: basis_dimension
  integer, public :: nR_mod
  integer, public :: np
  integer, public :: legpoints
    !!The number of Gauss-Legendre quadrature points to calculate integrals
  integer, public :: order
  integer, public :: n_s_energies

  real(wp), public :: mass
  real(wp), public :: x_begin
  real(wp), public :: x_end
  real(wp), public :: V_min

  real(wp), allocatable, public :: R_mod(:)
  real(wp), allocatable, public :: V_mod(:)
  real(wp), allocatable, public :: xwLegAll(:,:,:)
  real(wp), allocatable, public :: basis(:,:,:)
  real(wp), allocatable, public :: x_grid(:,:)
  real(wp), allocatable, public :: x_weights(:,:)
  real(wp), allocatable, public :: x_sectors(:)
  real(wp), allocatable, public :: B_overlap(:,:)
  complex(wp), allocatable, public :: Psi(:,:)
  complex(wp), allocatable, public :: H_matrix(:,:)
  complex(wp), allocatable, public :: Bound_energies(:)

! ================================================================================================================================ !
end module bsprvse__globals
! ================================================================================================================================ !
